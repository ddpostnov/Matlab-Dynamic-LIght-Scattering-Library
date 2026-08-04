%runInternalCycle  Collapse a high-frame-rate *.rls recording to one mean period
%
%   runInternalCycle(s,fNames) detects an endogenous periodic signal
%   (cardiac pulse) in each raw LSCI acquisition file
%   (*.rls), rejects outlier cycles by multiple feature criteria, and
%   averages the accepted cycles into a single X×Y×T contrast cube.  The
%   resulting cube represents one “canonical” period sampled at
%   s.interpFactor points per frame.  Preview JPEGs and three MAT-files are
%   written per recording:
%
%       *_c_K_d.mat   SOURCE – mean-period contrast cube (data) and
%                                 its time vector in seconds, source.time
%                                 [T x 1] - the library-wide orientation,
%                                 shared with runContrastFromRLS/runBolus/
%                                 runIntensity and with results.sData's
%                                 [nT x nSeg] frames-down-the-rows layout
%       *_c_K_r.mat   RESULTS – masks, cycle metrics, timestamps
%       *_c_K_s.mat   SETTINGS – copy of the parameter structure *s*, including
%                                 decimationSpaceUsed – the decimation actually
%                                 applied, which equals decimationSpace unless the
%                                 pre-processing cube did not fit in memory
%       *_rep_cycle-detect.jpg   – cycle-rejection overview
%       *_rep_cycle-average.jpg  – mean-period contrast image, the propagated
%                                 mask, and the BFI time-course
%
%   INPUTS
%     s        parameter structure (fields include sizeT, framesToAverage,
%              decimationSpace, contrastKernel*, fps correction, spectral
%              search bounds, rejection coefficients,
%              method = {'sLSCIMM','tLSCIMM','ltLSCIMM','sLSCIMMM'}, etc.).
%              TWO MASKS are parametrised separately, because they answer
%              different questions:
%                maskLimitsK, maskLimitsI  – mask A, the contrast and
%                     intensity range of the pixels averaged into the
%                     reference time course the heartbeat is DETECTED from.
%                     It does not leave this function.
%                trustLimitsK, trustLimitsI, minTrust – mask B, written to
%                     results.mask and read by every downstream step.  Same
%                     names and the same meaning as in runContrastFromRLS,
%                     but measured on the averaged cycle and on the frames
%                     that went into it.
%     fNames   cell array of full paths to *.rls files produced by the
%              in-house acquisition software.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUTS
%     None – function acts via side-effects (files listed above).
%
%   EXAMPLE
%     p = defaultInternalCycleParams();
%     D = dir(fullfile(dataRoot,'*.rls'));
%     runInternalCycle(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     readRLS, getK, getFFT.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
%
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
% s.contrastKernelS=5; %contrast kernel for spatial (sLSCI) processing method
% s.maxFrqIni=20; % initial max frequency of the activity of interest, Hz
% s.minFrqIni=1; % initial min frequency of the activity of interest, Hz
%
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - PULSE DETECTION
% s.maskLimitsK=[0.01,0.3]; %contrast range of the pixels averaged to detect the pulse
% s.maskLimitsI=[5,250]; %intensity range of the pixels averaged to detect the pulse
%
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - OUTPUT MASKING
% s.trustLimitsK=[0.001,0.99]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast
% s.trustLimitsI=[1,254]; %minimum (first value) and maximum (second value) of expected intensity.
% s.minTrust=[0.99,0.99]; %per-pixel trust limits in relation to the portion of frames with minimum (0) or maximum (usually 255) intensity.
%
% %ADJUSTED IF NECESSARY - EXCLUSION CRITERIA
% s.excludeFirstNCycles=0; %reject given number of cycles
% s.coeffsSTD=[3,2,2,2,2,3,3,2,2]; %pulses rejection coefficients relative to the feature standard deviation
% s.coeffsRel=[0.5,0.1]; %pulses rejection coefficients relative to the feature value
% s.coeffsAbs=2; %pulses rejection coefficients relative to the absolute feature value
%
% %ADJUSTED IF NECESSARY - CYCLE CALCULATION
% s.method='sLSCIMM';%,'tLSCIMM','ltLSCIMM' %Typically 'sLSCIMM' is recommended. For high quality data 'ltLSCIMM' will produce better results. Other options are 'tLSCIMM' and 'sLSCIMMM'.
% % method refers to spatial, temporal or lossless contrast calculation,
% % while the MM or MMM refers to minimum to minimum stretching or minimum to
% % maximum + maximum to minimum stretching.
% s.decimationSpace=4; %spatial decimation used to conserve memory in the pre-processing steps
% s.framesToAverage=1; %allows averaging multiple raw frames to artificially increase expsoure time
% s.contrastKernelT=25; %contrast kernel for temporal (tLSCI) and lossless (ltLSCI) processing methods
% s.contrastKernelPreproc=s.contrastKernelS; %contrast kernel used in preprocessing (spatial)
% s.rangeFrq=1;%1/2; % relative frq range around the central frequency, Hz
% s.interpFactor=10; %Sets the number of points that will replace two consequitive values during the interpolation sequence.
% s.smoothCoef1=1/3; %in respect to minimum points per cycle value
% s.minPromCoef=1/4;%1/2; % in respect to the std of the signal

function runInternalCycle(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'.rls'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".rls".');
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Internal cycle',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        % The blanket close-every-figure call that used to start each iteration is
        % gone.  It existed only to mop up the report figures this wrapper leaked, and
        % it closed UIFIGURES too - it could take the Processing Workbench window down
        % mid-run.  reportSave deletes its own figure on every path.
        clearvars results source settings
        s.fName=char(fNames{fidx});
        reportFile(rep,fidx,s.fName);

        %read the raw file meta data and open a stream on the first frame.  readRLS
        %owns the header: the parse that used to live here compared a uint8 vector to
        %'Ver.', so it never found the version block and read every uint16 recording as
        %uint8 - while pass 2, which goes through readRLS, read the same file correctly.
        %The stream also owns the handle, which the fopen this replaces never closed.
        [rawTail,tsTail,stream]=readRLS(s.fName,'KeepOpen',true,'FramesToRead',1);
        s.sizeX=stream.sizeX;
        s.sizeY=stream.sizeY;
        s.sizeT=stream.sizeT-1;
        s.exposureTime=single(stream.sampling);
        s.fps=1000./s.exposureTime; %converting time between frames to fps. DOES NOT PROVIDE AN ACCURATE FPS VALUES, RE-EVALUATED BELOW.
        s.dataSize=stream.dataSize;
        s.dataType=stream.dataType;

        s.sizeT=s.sizeT-s.framesToAverage+1;

        %The decimated cube is by far the largest array this function holds -
        %4*fileSize/decimationSpace^2 bytes, 25 GB for a 100 GB recording at the
        %default 4 - and it is only fully written once the whole recording has been
        %streamed, so a machine that cannot hold it finds out an hour into the read.
        %decimationSpace stays the user's choice and is only ever raised, never
        %lowered: the cube feeds tsK, a mask-weighted spatial average, and a coarser
        %sampling of an average is still that average.  The loop needs no bound: numel(1:d:n)
        %falls to 1, so the cube shrinks monotonically to a floor of s.sizeT*4 bytes - a few
        %hundred kB even for a 100 GB recording - which any machine running MATLAB clears.
        [~,memAvail]=memory;
        memAvail=memAvail.PhysicalMemory.Available;
        decimSpace=s.decimationSpace;
        while numel(1:decimSpace:s.sizeY)*numel(1:decimSpace:s.sizeX)*s.sizeT*4 > memAvail/2
            decimSpace=decimSpace+1;
        end
        if decimSpace>s.decimationSpace
            warning(['s.decimationSpace=%d needs %.1f GB for the pre-processing cube ' ...
                'and half the free memory is %.1f GB; raised to %d (%.1f GB).'], ...
                s.decimationSpace, ...
                numel(1:s.decimationSpace:s.sizeY)*numel(1:s.decimationSpace:s.sizeX)*s.sizeT*4/2^30, ...
                memAvail/2/2^30,decimSpace, ...
                numel(1:decimSpace:s.sizeY)*numel(1:decimSpace:s.sizeX)*s.sizeT*4/2^30);
        end
        s.decimationSpaceUsed=decimSpace;

        data=zeros(numel(1:decimSpace:s.sizeY),numel(1:decimSpace:s.sizeX),s.sizeT,'single');
        meanI=zeros(1,s.sizeT,'single');
        imgI=zeros(s.sizeY,s.sizeX,'single','gpuArray');
        imgK=zeros(s.sizeY,s.sizeX,'single','gpuArray');
        timeStamps=zeros(s.sizeT,1,'int64');

        %perform data pre-processing, store spacially decimated information.  The
        %recording is streamed in batches and the contrast never leaves the GPU until it
        %has been decimated - the full-resolution contrast of a 100 GB recording is
        %400 GB, and all this pass wants from it is one spatial average per frame.
        kernelDecim=gpuArray(ones(decimSpace,decimSpace,'single'));
        gpuDev=gpuDevice;
        %256 frames per batch measured fastest on a 24 GB card (64-256 is flat, past
        %that the allocator rather than the arithmetic decides); the memory bound is
        %what actually applies on a small card or a large frame
        batchSize=min([256, ...
            max(1,floor(gpuDev.AvailableMemory.*0.2./(s.sizeY.*s.sizeX.*4))),s.sizeT]);
        iOut=0;
        while iOut<s.sizeT
            nOut=min(batchSize,s.sizeT-iOut);
            %a sliding frame average needs framesToAverage-1 frames of overlap, so the
            %tail of the previous batch is carried into this one.  The loop this
            %replaces read framesToAverage FRAMES after a SINGLE timestamp, but the
            %layout is [ts][frame][ts][frame]... - so it consumed timestamps as pixels,
            %drifted 8 bytes per iteration and fed garbage into the framerate correction.
            [raw,ts,stream]=readRLS(stream,'FramesToRead',nOut+s.framesToAverage-1-size(rawTail,3));
            raw=cat(3,rawTail,raw);
            ts=cat(1,tsTail,ts);                 %cat, not [;] - it is a fresh bounded vector each round, not a growing one
            timeStamps(iOut+(1:nOut))=int64(ts(1:nOut));
            tmp=single(gpuArray(raw));
            if s.framesToAverage>1
                tmp=movmean(tmp,[0,s.framesToAverage-1],3,'Endpoints','discard');
            end
            meanI(iOut+(1:nOut))=gather(mean(tmp,[1,2]));
            imgI=imgI+sum(tmp,3);
            tmp=getK(tmp,'spatial','kernelSize',s.contrastKernelPreproc, ...
                'procType','gpu','outputType','gpu');
            imgK=imgK+sum(tmp,3);
            tmp=convn(tmp,kernelDecim,'same')./(decimSpace.*decimSpace);
            data(:,:,iOut+(1:nOut))=gather(tmp(1:decimSpace:end,1:decimSpace:end,:));
            rawTail=raw(:,:,end-s.framesToAverage+2:end);
            tsTail=ts(end-s.framesToAverage+2:end);
            iOut=iOut+nOut;
        end
        fclose(stream.fId);
        imgI=gather(imgI)./s.sizeT;
        imgK=gather(imgK)./s.sizeT;
        %results.imgI stays the WHOLE-RECORDING mean intensity, deliberately unlike
        %mask B below, which measures only the accepted-cycle frames.  This image is a
        %background for setRegions and runSegmentation to draw and threshold on, and
        %what it should show is the field of view as it was recorded - not the subset
        %of it the cycle rejection happened to keep.
        results.imgI=imgI;
        %mask A - the pixels averaged into the reference time course below.  It decides
        %which pixels the heartbeat is DETECTED from and goes no further; the mask that
        %travels with the result is mask B, built from the averaged cycle at the end.
        mask=imgK<s.maskLimitsK(2) & imgK>s.maskLimitsK(1) & isfinite(imgK) & imgI<s.maskLimitsI(2) & imgI>s.maskLimitsI(1);
        mask=imerode(mask,[0,1,0;1,1,1;0,1,0]);
        mask=(convn(mask,ones(decimSpace),'same')./(decimSpace.*decimSpace))>0;

        mask=mask(1:decimSpace:end,1:decimSpace:end);
        data(~isfinite(data))=0;                 %isfinite is false for NaN too - the isnan line that stood here swept the whole cube for nothing
        data=reshape(data,size(data,1)*size(data,2),s.sizeT);
        tsK=(single(mask(:)).'*data)./nnz(mask); %a matrix-vector product; sum(data.*mask(:),1) allocated a second copy of the largest array in the function
        tsBFI=1./(tsK.*tsK);

        %correct the framerate
        s.fps=1000./(median(double(timeStamps(2:end)-timeStamps(1:end-1))));
        time=((1:1:length(tsK))-1)./s.fps;

        %identify central frequency
        [fftPow,~,f]=getFFT(tsBFI,s.fps,2.^(nextpow2(s.fps/s.maxFrqIni*50)),'cpu');
        idxLimits=[find(f(:)>s.minFrqIni,1,'first'),find(f(:)>s.maxFrqIni,1,'first')-1];
        [~,idx,~,prm]=findpeaks(fftPow(idxLimits(1):idxLimits(2)));
        [~,idxTMP]=max(prm);
        idx=idx(idxTMP);
        idx=idx+idxLimits(1)-1;
        s.centralFrq=f(idx);
        s.minFrq=max(s.centralFrq*(1-s.rangeFrq),s.minFrqIni);
        s.maxFrq=min(s.centralFrq*(1+s.rangeFrq),s.maxFrqIni);
        s.meanPointsPerCycle=floor(s.fps/s.centralFrq);

        tsBFIF=smooth(tsBFI,max(floor(s.meanPointsPerCycle.*s.smoothCoef1),1),'loess');
        s.minProm=s.minPromCoef.*std(tsBFIF);
        idxPeak=idx;                            % the spectral peak; idx is a loop temp below

        [~,locsMin]=findpeaks(-tsBFIF,'MinPeakDistance',floor(s.fps/s.maxFrq),'MinPeakProminence',s.minProm);
        locsMax=zeros(length(locsMin)-1,1);
        for i=1:1:length(locsMin)-1
            [~,idx]=max(tsBFIF(locsMin(i):locsMin(i+1)));
            locsMax(i)=locsMin(i)+idx-1;

        end
        pulsesList=zeros(length(locsMax),3);
        for i=1:1:length(locsMin)-1
            pulsesList(i,:)=[locsMin(i),locsMax(i),locsMin(i+1)];
        end

        pulsesFeatures=zeros(size(pulsesList,1),9);
        for i=1:1:size(pulsesFeatures,1)
            pulsesFeatures(i,1)=mean(tsBFIF(squeeze(pulsesList(i,:))));
            pulsesFeatures(i,2)=std(tsBFIF(squeeze(pulsesList(i,:))));
            pulsesFeatures(i,3)=max(tsBFIF(squeeze(pulsesList(i,:))))-min(tsBFIF(squeeze(pulsesList(i,:))));
            pulsesFeatures(i,4)=max(squeeze(pulsesList(i,:)))-min(squeeze(pulsesList(i,:)));
            pulsesFeatures(i,5)=abs(tsBFIF(squeeze(pulsesList(i,1)))-tsBFIF(squeeze(pulsesList(i,3))));
            pulsesFeatures(i,6)=sum((tsBFIF(pulsesList(i,1)+1:pulsesList(i,2))-tsBFIF(pulsesList(i,1):pulsesList(i,2)-1))<0);
            pulsesFeatures(i,7)=sum((tsBFIF(pulsesList(i,2)+1:pulsesList(i,3))-tsBFIF(pulsesList(i,2):pulsesList(i,3)-1))>0);
            pulsesFeatures(i,8)=mean((tsBFIF(pulsesList(i,2)+1:pulsesList(i,3))-tsBFIF(pulsesList(i,2):pulsesList(i,3)-1)));
            pulsesFeatures(i,9)=mean(meanI(pulsesList(i,1)+1:pulsesList(i,3)));
        end

        %One row per rejection rule, one column per cycle.  Rules 1-9 compare a cycle's
        %feature to the median and spread of the OTHER cycles and are taken in closed
        %form by featureOutliers below; the other four never depended on the other
        %cycles at all, so they read straight off the feature columns.  Rule 12 asks
        %about the METHOD, not about the cycle, so its strcmp is outside everything.
        pulsesToReject=zeros(14,size(pulsesList,1));
        pulsesToReject(1:size(pulsesFeatures,2),:)=featureOutliers(pulsesFeatures,s.coeffsSTD);
        pulsesToReject(10,(pulsesFeatures(:,5)./pulsesFeatures(:,3)>s.coeffsRel(1)).')=1;
        cycleFrq=s.fps./pulsesFeatures(:,4);
        pulsesToReject(11,(cycleFrq<s.minFrq | cycleFrq>s.maxFrq).')=1;
        if any(strcmp(s.method,'tLSCIMM'))
            %Rule 12 guards BOTH ends of the recording.  tLSCIMM's pass-2 read is padded
            %by half a temporal kernel on either side so getK's centred window is full
            %for every frame of an accepted cycle; the trailing pad has always been
            %protected here, the LEADING one was not.  A cycle starting within kernelPad
            %of frame 1 then asks readRLS for a negative FramesToSkip, which it answers
            %with a misaligned block rather than an error - so that cycle was averaged
            %from the wrong frames.  Rejecting it costs one cycle and keeps every window
            %in the average full, which clamping the read would not: getK's edge
            %normalisation would then engage on the leading window alone, and one cycle
            %of the average would be taken over a shorter kernel than all the others.
            pulsesToReject(12,(pulsesList(:,3)>(s.sizeT-s.contrastKernelT) ...
                | pulsesList(:,1)<=floor(s.contrastKernelT*s.framesToAverage/2)).')=1;
        end
        pulsesToReject(13,(abs(1-pulsesFeatures(:,1)./median(pulsesFeatures(:,1)))>s.coeffsRel(2)).')=1;

        if s.excludeFirstNCycles>0
            pulsesToReject(1:end,1:s.excludeFirstNCycles)=1;
        end

        %reject/accept based on number of consequtive accepts/rejects
        a=max(pulsesToReject,[],1);
        a=movmean(a,[s.coeffsAbs(1),s.coeffsAbs(1)]);
        pulsesToReject(14,:)=round(a);

        % hold on
        % for i=1:1:size(pulsesList,1)
        %     plot([time(pulsesList(i,1)),time(pulsesList(i,3))],[tsBFIF(pulsesList(i,1)),tsBFIF(pulsesList(i,3))],'ok')
        %     if sum(pulsesToReject(:,i))==0
        %         plot(time(pulsesList(i,1):pulsesList(i,3)),tsBFIF(pulsesList(i,1):pulsesList(i,3)),'g')
        %     else
        %         plot(time(pulsesList(i,1):pulsesList(i,3)),tsBFIF(pulsesList(i,1):pulsesList(i,3)),'r')
        %     end
        %     for ii=1:1:size(pulsesToReject,1)
        %         if pulsesToReject(ii,i)==0
        %             plot([time(pulsesList(i,1)),time(pulsesList(i,3))],[min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(pulsesToReject,1)/2,min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(pulsesToReject,1)/2],'g');
        %         else
        %             plot([time(pulsesList(i,1)),time(pulsesList(i,3))],[min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(pulsesToReject,1)/2,min(tsBFIF)-(max(tsBFIF)-min(tsBFIF)).*ii/size(pulsesToReject,1)/2],'r');
        %         end
        %
        %     end
        % end
        % hold off

        pulsesListFinal=pulsesList(~any(pulsesToReject,1),:);

        % --- report page: what the cycle detection saw and what it kept ---------
        % Drawn here rather than while the pulses are being sorted, so nothing in
        % the loop above has to hold a figure open across a hundred lines of
        % arithmetic.  Same three panels as before.
        fh=reportFigure(rep,'cycle-detect');
        tl=tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl,1);
        plot(ax,time(1:s.meanPointsPerCycle*2),tsBFIF(1:s.meanPointsPerCycle*2))
        hold(ax,'on')
        plot(ax,time(1:s.meanPointsPerCycle*2),1./(tsK(1:s.meanPointsPerCycle*2).*tsK(1:s.meanPointsPerCycle*2)))
        hold(ax,'off')
        xlabel(ax,'Time, s')
        ylabel(ax,'BFI')
        ax=nexttile(tl,3);
        plot(ax,f,fftPow)
        hold(ax,'on')
        plot(ax,f(idxPeak),fftPow(idxPeak),'or')
        hold(ax,'off')
        xlabel(ax,'Frq, Hz')
        ylabel(ax,'Amplitude')
        title(ax,['Central frq=',num2str(s.centralFrq)]);
        ax=nexttile(tl,2,[2,1]);
        plot(ax,time,1./(tsK.*tsK));
        xlabel(ax,'Time,s')
        ylabel(ax,'BFI')
        title(ax,['acceptance rate=',num2str(size(pulsesListFinal,1)./size(pulsesList,1))]);
        reportSave(rep,fh,'cycle-detect');
        clearvars data;

        s.descendTimePts=round(median(pulsesListFinal(:,3)-pulsesListFinal(:,2)));
        s.ascendTimePts=round(median(pulsesListFinal(:,2)-pulsesListFinal(:,1)));
        s.cycleTimePts=s.ascendTimePts+s.descendTimePts-1;

        %Pass 2 re-reads the accepted cycles in CONTIGUOUS RUNS.  Cycles are
        %min-to-min, so two neighbours share a frame whenever nothing between them was
        %rejected: one read per run keeps the pass sequential over a 100 GB file and
        %gives getK a batch worth calling it on - called once per cycle it costs more
        %than the loop it replaces.  Each run's contrast then stays on the GPU as
        %[pixels x frames] and is split back into its cycles by a column slice.
        %Cycles of EQUAL LENGTH are summed on their own frame grid and interpolated
        %onto the common phase grid ONCE PER DISTINCT LENGTH.  Interpolating every
        %cycle separately built a temporary the size of the whole interpolated cube
        %per cycle and was two thirds of this step's runtime; linear interpolation is
        %linear, so the sum of the interpolations IS the interpolation of the sum and
        %nothing about the result changes.  The ltLSCIMM branch already grouped by
        %duration - this is its own idea applied to the other three.
        nPx=s.sizeY.*s.sizeX;
        cyclesDurList=pulsesListFinal(:,3)-pulsesListFinal(:,1)+1;

        %Mask B's intensity statistics, accumulated by addTrustStats over the frames
        %the accepted cycles actually read - NOT over pass 1's whole-recording imgI,
        %because the mask has to describe what went into the average.  Declared here
        %so all four branches write into one accumulator, resolved once after the
        %chain, and so the two constants that decide WHICH frames of a read are
        %counted can be set once: kernelPad is 0 for the three branches whose read is
        %not padded, and iFirst is 1 unless a run opens on a frame the previous one
        %already counted.  The accumulate line is then the same in all four.
        kernelPad=0;
        iFirst=1;
        t=struct('n',0,'sumI',zeros(s.sizeY,s.sizeX,'single','gpuArray'), ...
            'nZero',zeros(s.sizeY,s.sizeX,'single','gpuArray'), ...
            'nSat', zeros(s.sizeY,s.sizeX,'single','gpuArray'));

        %SLSCI method Min Max Min
        if any(strcmp(s.method,'sLSCIMMM'))
            %Ascend and descend land on different phase grids, so each is binned by
            %ITS OWN length.  They share the peak frame - deliberately - and the
            %1:end-1 in the concatenation below is what removes the duplicate.
            ascendDurList=pulsesListFinal(:,2)-pulsesListFinal(:,1)+1;
            descendDurList=pulsesListFinal(:,3)-pulsesListFinal(:,2)+1;
            [ascendDur,ascendCol,ascendBin]=makeCycleBins(ascendDurList,nPx,gpuDev);
            [descendDur,descendCol,descendBin]=makeCycleBins(descendDurList,nPx,gpuDev);
            runs=cycleRuns(pulsesListFinal,max(1,batchSize-s.framesToAverage+1));
            for i=1:1:size(runs,1)
                iOff=pulsesListFinal(runs(i,1),1)-1;
                span=pulsesListFinal(runs(i,2),3)-iOff;
                [data,raw]=readAveraged(s,iOff+1,span);
                %Every frame of a run belongs to an accepted cycle, and within a run the
                %min-to-min boundary two neighbours share is read once - that is what the
                %run-based read bought.  ACROSS runs it is not: a run cut by the frame cap
                %rather than by a rejected cycle opens on the frames the previous one
                %ended on, so iFirst steps past whatever has already been counted.  Left
                %alone, the statistics - and therefore mask B - would depend on batchSize,
                %i.e. on the size of the card.  framesToAverage>1 makes the overlap
                %framesToAverage frames rather than one, hence the arithmetic rather than
                %a test for equality.
                if i>1, iFirst=max(1,pulsesListFinal(runs(i-1,2),3)+s.framesToAverage-iOff); end
                t=addTrustStats(t,raw(:,:,kernelPad+(iFirst:span+s.framesToAverage-1)));
                data=getK(data,'spatial','kernelSize',s.contrastKernelS, ...
                    'procType','gpu','outputType','gpu');
                data=reshape(data,nPx,[]);
                for ii=runs(i,1):1:runs(i,2)
                    iCol=ascendCol(ii)+(0:ascendDurList(ii)-1);
                    ascendBin(:,iCol)=ascendBin(:,iCol) ...
                        +data(:,(pulsesListFinal(ii,1):pulsesListFinal(ii,2))-iOff);
                    iCol=descendCol(ii)+(0:descendDurList(ii)-1);
                    descendBin(:,iCol)=descendBin(:,iCol) ...
                        +data(:,(pulsesListFinal(ii,2):pulsesListFinal(ii,3))-iOff);
                end
            end
            dataASLSCI=interpolateBins(ascendBin,ascendDur,s.ascendTimePts*s.interpFactor);
            dataDSLSCI=interpolateBins(descendBin,descendDur,s.descendTimePts*s.interpFactor);
            dataASLSCI=reshape(gather(dataASLSCI),s.sizeY,s.sizeX,[]);
            dataDSLSCI=reshape(gather(dataDSLSCI),s.sizeY,s.sizeX,[]);
            source.data=cat(3,dataASLSCI(:,:,1:end-1),dataDSLSCI)./size(pulsesListFinal,1);
            source.time=((0:1:size(source.data,3)-1)./s.fps./s.interpFactor)';

        elseif any(strcmp(s.method,'sLSCIMM'))
            %SLSCI method Min Min
            [cyclesDur,cycleCol,cyclesBin]=makeCycleBins(cyclesDurList,nPx,gpuDev);
            runs=cycleRuns(pulsesListFinal,max(1,batchSize-s.framesToAverage+1));
            for i=1:1:size(runs,1)
                iOff=pulsesListFinal(runs(i,1),1)-1;
                span=pulsesListFinal(runs(i,2),3)-iOff;
                [data,raw]=readAveraged(s,iOff+1,span);
                if i>1, iFirst=max(1,pulsesListFinal(runs(i-1,2),3)+s.framesToAverage-iOff); end
                t=addTrustStats(t,raw(:,:,kernelPad+(iFirst:span+s.framesToAverage-1)));
                data=getK(data,'spatial','kernelSize',s.contrastKernelS, ...
                    'procType','gpu','outputType','gpu');
                data=reshape(data,nPx,[]);
                for ii=runs(i,1):1:runs(i,2)
                    iCol=cycleCol(ii)+(0:cyclesDurList(ii)-1);
                    cyclesBin(:,iCol)=cyclesBin(:,iCol) ...
                        +data(:,(pulsesListFinal(ii,1):pulsesListFinal(ii,3))-iOff);
                end
            end
            dataCycle=interpolateBins(cyclesBin,cyclesDur,s.cycleTimePts*s.interpFactor);
            source.data=reshape(gather(dataCycle),s.sizeY,s.sizeX,[])./size(pulsesListFinal,1);
            source.time=((0:1:size(source.data,3)-1)./s.fps./s.interpFactor)';

        elseif any(strcmp(s.method,'tLSCIMM'))
            %TLSCI method Min Min
            %use with care TLSCI might be contaminated by the rejected pulses
            %The read stays padded by half a temporal kernel, but the TRIM moves: the
            %movstd this replaces applied a FORWARD window to the padded stretch, so
            %the cycle sat at 1:cycleLength, while getK applies a CENTRED one, which
            %puts it at kernelPad+(1:cycleLength).  With an odd kernel every one of
            %those windows is full, so getK's edge normalisation never engages.
            %The LEADING pad is guaranteed to exist: rejection rule 12 drops any cycle
            %starting within kernelPad of frame 1, which is what keeps the FramesToSkip
            %below non-negative.
            kernelT=s.contrastKernelT*s.framesToAverage;
            kernelPad=floor(kernelT/2);
            [cyclesDur,cycleCol,cyclesBin]=makeCycleBins(cyclesDurList,nPx,gpuDev);
            runs=cycleRuns(pulsesListFinal,max(1,batchSize-kernelT-s.framesToAverage+2));
            for i=1:1:size(runs,1)
                iOff=pulsesListFinal(runs(i,1),1)-1;
                span=pulsesListFinal(runs(i,2),3)-iOff;
                [data,raw]=readAveraged(s,iOff+1-kernelPad,span+kernelT-1);
                %the kernelPad frames on either side exist so getK's centred window has
                %data; they belong to no accepted cycle, so the trim keeps them out of
                %the statistics - the same expression the unpadded branches use
                if i>1, iFirst=max(1,pulsesListFinal(runs(i-1,2),3)+s.framesToAverage-iOff); end
                t=addTrustStats(t,raw(:,:,kernelPad+(iFirst:span+s.framesToAverage-1)));
                data=getK(data,'temporal','kernelSize',kernelT, ...
                    'procType','gpu','outputType','gpu');
                data=reshape(data,nPx,[]);
                for ii=runs(i,1):1:runs(i,2)
                    iCol=cycleCol(ii)+(0:cyclesDurList(ii)-1);
                    cyclesBin(:,iCol)=cyclesBin(:,iCol) ...
                        +data(:,(pulsesListFinal(ii,1):pulsesListFinal(ii,3))-iOff+kernelPad);
                end
            end
            dataCycle=interpolateBins(cyclesBin,cyclesDur,s.cycleTimePts*s.interpFactor);
            source.data=reshape(gather(dataCycle),s.sizeY,s.sizeX,[])./size(pulsesListFinal,1);
            source.time=((0:1:size(source.data,3)-1)./s.fps./s.interpFactor)';

        elseif any(strcmp(s.method,'ltLSCIMM'))
            %LTLSCI method Min Min
            %This branch takes the contrast ACROSS cycles rather than across a kernel,
            %so getK does not apply and the estimator is unchanged.  What changes is
            %that the cycles are accumulated as running sums instead of being stacked
            %into one sizeY*sizeX*duration*nCycles array - 53 GB for a thousand cycles
            %at 800x800.  The sums are taken about the group's first cycle so that the
            %difference below subtracts two quantities of the same order.
            nOut=s.cycleTimePts*s.interpFactor;
            dataCycle=zeros(nPx,nOut,'single','gpuArray');
            cyclesDur=unique(cyclesDurList);
            cyclesN=0;
            for i=1:1:length(cyclesDur)
                pulsesIdxs=find(cyclesDurList==cyclesDur(i));
                if length(pulsesIdxs)>=s.contrastKernelT
                    sumD =zeros(s.sizeY,s.sizeX,cyclesDur(i),'single','gpuArray');
                    sumD2=zeros(s.sizeY,s.sizeX,cyclesDur(i),'single','gpuArray');
                    shiftD=sumD;
                    span=cyclesDur(i);
                    for ii=1:1:length(pulsesIdxs)
                        [data,raw]=readAveraged(s,pulsesListFinal(pulsesIdxs(ii),1),span);
                        %accumulated INSIDE the if, so an under-populated duration group
                        %- whose cycles never enter the average - contributes nothing to
                        %the statistics either.  iFirst stays 1 here, unlike the three
                        %run-reading branches: this one reads and averages per CYCLE, so
                        %the min frame two neighbouring accepted cycles share genuinely
                        %entered the average twice and is counted twice.  Nothing about
                        %that depends on the machine, which is why it is kept.
                        t=addTrustStats(t,raw(:,:,kernelPad+(iFirst:span+s.framesToAverage-1)));
                        if ii==1, shiftD=data; end
                        data=data-shiftD;
                        sumD=sumD+data;
                        sumD2=sumD2+data.*data;
                    end
                    nCycle=length(pulsesIdxs);
                    data=sumD./nCycle;              % mean of the shifted values
                    data=sqrt(max((sumD2-sumD.*data)./(nCycle-1),0))./(data+shiftD);
                    dataCycle=dataCycle+nCycle.*(reshape(data,nPx,[]) ...
                        *gpuArray(interpWeights(cyclesDur(i),nOut)).');
                    cyclesN=cyclesN+nCycle;
                end
            end
            source.data=reshape(gather(dataCycle),s.sizeY,s.sizeX,[])./cyclesN;
            source.time=((0:1:size(source.data,3)-1)./s.fps./s.interpFactor)';

        else
            error('Unknown processing method requested');
        end

        %The accepted-cycle statistics, resolved once.  These are
        %intensityMetrics(:,:,3), (:,:,1) and (:,:,2) of getContrastFromRLS - same
        %quantities, same sign convention - restricted to the frames that went into
        %the average.  Gathered here so results.mask below is a host logical array and
        %not a gpuArray that would have to be unpacked by every step that loads it.
        imgIcycle   =gather(t.sumI ./t.n);
        trustNonZero=gather(1-t.nZero./t.n);     % fraction of non-zero frames
        trustNonSat =gather(1-t.nSat ./t.n);     % fraction of non-saturated frames

        %mask B - what travels with the result.  Same criteria and the same parameter
        %names as runContrastFromRLS, so a downstream step reads one vocabulary whether
        %it was handed a contrast product or a cardiac one.  The difference is the
        %input: every frame here is a PHASE of the averaged cycle, and the intensity
        %statistics come only from the frames that were accepted into that average.
        %min(...,[],3) means EVERY phase must be in range - systole is the extreme, so
        %a pixel that leaves the range only at peak flow is excluded, which is the same
        %reading runContrastFromRLS takes of its own time axis.
        results.mask=min(~isnan(source.data),[],3)...
            & min(source.data>=s.trustLimitsK(1),[],3)...
            & min(source.data<=s.trustLimitsK(2),[],3)...
            & imgIcycle>=s.trustLimitsI(1)...
            & imgIcycle<=s.trustLimitsI(2)...
            & trustNonZero>s.minTrust(1)...
            & trustNonSat >s.minTrust(2);

        %A mask nobody sees is a mask nobody checks, so the page carries it beside the
        %image it was computed from - the same reason runContrastFromRLS shows its own.
        fh=reportFigure(rep,'cycle-average');
        tl=tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl,1);
        img=squeeze(mean(source.data,3));
        imagesc(ax,img);
        clim(ax,[prctile(img(:),5),prctile(img(:),99)]);
        axis(ax,'image')
        %the method belongs to the panel it describes, and it is read from s - the
        %page used to carry a hardcoded 'sLSCIMMM' above it, which was simply wrong
        %on the other three methods.  reportSave titles the page with the recording.
        title(ax,['Contrast, ',char(s.method)])
        ax=nexttile(tl,2);
        imagesc(ax,results.mask);
        clim(ax,[0,1])                          %without it a mask that kept EVERYTHING and
        axis(ax,'image')                        %one that kept NOTHING are the same flat tile
        title(ax,['Mask, kept ',num2str(round(100.*nnz(results.mask)./numel(results.mask))),'%'])
        ax=nexttile(tl,3,[1,2]);
        plot(ax,source.time,1./(squeeze(mean(source.data,[1,2]))).^2);
        xlabel(ax,'Time,s')
        ylabel(ax,'BFI')
        reportSave(rep,fh,'cycle-average');

        reportWriting(rep);
        settings.runInternalCycle=reportSettings(s);
        results.imgK=squeeze(mean(source.data,3,'omitmissing'));
        results.time=source.time;
        %-nocompression is the library-wide save rule (README, 'File format'): deflate
        %is single-threaded and cost 13.2 s of the 14.2 s spent writing this step's
        %0.6 GB cube, for 19% of the disk.  Measured here first, now applied to every
        %save in the library - so do not "restore" the default on one file.
        %This is the ENTRY step, so this is where a raw path becomes a product path, and
        %the only place that has to know a project may keep its results apart from its
        %recordings.  s.fName is NOT reassigned - it stays the raw recording.  With no
        %results folder set the name comes back verbatim.
        outBase=getResultsPath(s.fName,s);
        save(strrep(outBase,'.rls','_c_K_d.mat'),'source','-v7.3','-nocompression');
        save(strrep(outBase,'.rls','_c_K_r.mat'),'results','-v7.3','-nocompression');
        save(strrep(outBase,'.rls','_c_K_s.mat'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);

end

% =====================================================================
function reject=featureOutliers(pulsesFeatures,coeffsSTD)
%featureOutliers  Rejection rules 1-9: a cycle whose feature stands out from the rest.
%   reject(ii,i) is true when |f - median(f without i)| > coeffsSTD(ii)*std(f without
%   i) for feature ii of cycle i.  One row per feature, one column per cycle.
%
%   The leave-one-out median and spread are taken in CLOSED FORM rather than
%   recomputed for every (cycle,feature) pair: doing that built an N-element index
%   vector 9N times and cost 5.5 s at 7000 cycles, growing as N^2.
%
%   With mu=mean(x) and Q=sum((x-mu).^2) over all N values of one feature,
%       mu_i  = (N*mu-x_i)/(N-1)
%       Q_i   = Q-(x_i-mu).^2*N/(N-1)
%       std_i = sqrt(max(Q_i,0)/(N-2))        std of N-1 samples normalises by N-2
%   because sum_{j~=i}(x_j-mu_i)^2 = sum_{j~=i}(x_j-mu)^2-(N-1)*(mu_i-mu)^2 and
%   mu_i-mu = (mu-x_i)/(N-1), so the correction is (x_i-mu)^2*(1+1/(N-1)).  Q is the
%   CENTRED sum: sum(x.^2)-N*mu^2 would cancel two digits on feature 9, whose values
%   are intensities near 100 with a small spread.
%
%   The median comes off the sorted column - dropping the value of rank r leaves a
%   sorted array whose j-th order statistic is xs(j+(j>=r)) - and the even case is
%   combined the way median itself combines it, a+(b-a)/2 except across a sign change.
nF=size(pulsesFeatures,2);
nCycles=size(pulsesFeatures,1);
reject=false(nF,nCycles);
for ii=1:1:nF
    x=pulsesFeatures(:,ii);
    dev=x-mean(x);
    devQ=dev.*dev;
    devSum=sum(devQ);
    if nCycles<=2 || ~isfinite(devSum)
        %Two exceptions, both rare, and neither has a closed form.  With two cycles
        %the leave-one-out set is a single sample, whose std MATLAB defines as 0 while
        %the identity above reads 0/0.  And Q not being a usable number covers the
        %rest: any NaN or Inf among the values reaches it through the mean, and so
        %does the overflow of a column of ~1e155.  A non-finite value is the case that
        %matters - median and std then return NaN for every cycle EXCEPT the one
        %holding it, whose own leave-one-out set has just dropped it, and no sum or
        %order statistic reproduces that asymmetry.  Feature 8 reaches it: it is the
        %mean of an empty difference whenever a cycle's maximum lands on its last frame.
        for i=1:1:nCycles
            xLoo=x([1:i-1,i+1:end]);
            reject(ii,i)=abs(x(i)-median(xLoo))>coeffsSTD(ii)*std(xLoo);
        end
        continue
    end

    sdLoo=sqrt(max(devSum-devQ.*(nCycles/(nCycles-1)),0)./(nCycles-2));

    [xs,ord]=sort(x);
    rnk=zeros(nCycles,1);
    rnk(ord)=1:nCycles;
    if rem(nCycles,2)==0                        % the leave-one-out set is odd
        medLoo=xs(nCycles/2+(rnk<=nCycles/2));
    else                                        % even - two order statistics to combine
        jLo=(nCycles-1)/2;
        a=xs(jLo+(rnk<=jLo));
        b=xs(jLo+1+(rnk<=jLo+1));
        medLoo=a+(b-a)./2;
        acrossZero=sign(a)~=sign(b);            % median's own guard against overflow
        medLoo(acrossZero)=(a(acrossZero)+b(acrossZero))./2;
    end

    reject(ii,:)=(abs(x-medLoo)>coeffsSTD(ii).*sdLoo).';
end
end

% =====================================================================
function [data,raw]=readAveraged(s,firstFrame,nFrames)
%readAveraged  nFrames averaged frames from firstFrame, as single on the GPU.
%   framesToAverage>1 is a SLIDING window that discards the partial one at the end -
%   the same semantics pass 1 uses, so both passes see the same frames.  The read is
%   nFrames+framesToAverage-1 raw frames because of it.  Reading straight onto the
%   card also avoids movmean's 8x blow-up: on an INTEGER input it returns double.
%   The second output is that raw block, still in its file class and already on the
%   card - mask B counts zero and saturated SAMPLES, so it needs the frames before the
%   cast and before the frame average.  Callers that do not ask for it pay nothing.
raw=gpuArray(readRLS(s.fName,'FramesToSkip',firstFrame-1, ...
    'FramesToRead',nFrames+s.framesToAverage-1));
data=single(raw);
if s.framesToAverage>1
    data=movmean(data,[0,s.framesToAverage-1],3,'Endpoints','discard');
end
end

% =====================================================================
function t=addTrustStats(t,raw)
%addTrustStats  Mask B's per-pixel intensity statistics, over accepted-cycle frames.
%   Counts of zero and saturated SAMPLES plus a running intensity sum, taken on the
%   raw block so they describe the exposure rather than the frame average.  The three
%   quantities are intensityMetrics(:,:,3), (:,:,1) and (:,:,2) of getContrastFromRLS
%   restricted to the accepted cycles - same sign convention, because mask B's
%   thresholds are the ones runContrastFromRLS uses.
t.n    =t.n    +size(raw,3);
t.sumI =t.sumI +sum(single(raw),3);
t.nZero=t.nZero+sum(single(raw==0),3);
%underlyingType, not classUnderlying (discouraged since R2023b) and not class, which
%answers 'gpuArray' and would make intmax throw.
t.nSat =t.nSat +sum(single(raw==intmax(underlyingType(raw))),3);
end

% =====================================================================
function [cyclesDur,cycleCol,cyclesBin]=makeCycleBins(cyclesDurList,nPx,gpuDev)
%makeCycleBins  One accumulator block per DISTINCT cycle length.
%   Cycles of equal length share a frame grid, so they can be summed BEFORE they are
%   interpolated - which is what turns the interpolation from a per-cycle cost into a
%   per-length one.  cyclesBin is [nPx x sum(cyclesDur)]: cycle i writes into the
%   columns cycleCol(i)+(0:cyclesDurList(i)-1).
[cyclesDur,~,binIdx]=unique(cyclesDurList);
binStart=cumsum([1;cyclesDur(1:end-1)]);
cycleCol=binStart(binIdx);
binBytes=nPx*sum(cyclesDur)*4;
%A stable heart rate leaves a handful of distinct lengths.  Many of them means the
%cycle rejection has not converged, and averaging those cycles would be wrong however
%it is computed - hence an error rather than a fallback to a per-cycle path.
if binBytes>gpuDev.AvailableMemory*0.5
    error('runInternalCycle:TooManyCycleLengths', ...
        ['%d distinct accepted cycle lengths need %.1f GB of GPU memory. ' ...
        'Rejection is not converging - check coeffsSTD(4) and rangeFrq.'], ...
        numel(cyclesDur),binBytes/2^30);
end
cyclesBin=zeros(nPx,sum(cyclesDur),'single','gpuArray');
end

% =====================================================================
function dataCycle=interpolateBins(cyclesBin,cyclesDur,nOut)
%interpolateBins  Each length's accumulated block onto the common phase grid, summed.
%   Once per distinct length, not once per cycle.  Linear interpolation is linear, so
%   this is exactly the sum of the per-cycle interpolations it replaces.
dataCycle=zeros(size(cyclesBin,1),nOut,'single','gpuArray');
binStart=cumsum([1;cyclesDur(1:end-1)]);
for i=1:1:length(cyclesDur)
    dataCycle=dataCycle+cyclesBin(:,binStart(i)+(0:cyclesDur(i)-1)) ...
        *gpuArray(interpWeights(cyclesDur(i),nOut)).';
end
end

% =====================================================================
function W=interpWeights(nIn,nOut)
%interpWeights  The linear-interpolation operator interp1 would apply.
%   W is [nOut x nIn] and D*W.' equals interp1(D.',linspace(1,nIn,nOut)).' for a
%   row-wise D.  Writing it out is what lets equal-length cycles be summed BEFORE
%   they are interpolated - the sum of interpolations is the interpolation of the
%   sum - and it makes the interpolation a matrix product, so no permute is needed
%   to reach it and none is needed to come back.
if nIn<2
    error('runInternalCycle:ShortCycle','A cycle needs at least two frames to interpolate.');
end
x=linspace(1,nIn,nOut);
i0=min(floor(x),nIn-1);
w=x-i0;
W=zeros(nOut,nIn,'single');
W(sub2ind([nOut,nIn],(1:nOut).',i0.'))  =1-w.';
W(sub2ind([nOut,nIn],(1:nOut).',i0.'+1))=w.';
end

% =====================================================================
function runs=cycleRuns(pulsesListFinal,maxFrames)
%cycleRuns  Accepted cycles grouped into contiguous, length-capped runs.
%   Cycles are min-to-min, so two neighbours in the list share a frame whenever
%   nothing between them was rejected.  Reading such a stretch in one go keeps the
%   pass sequential over a 100 GB file and gives getK a batch worth calling it on -
%   called once per cycle it costs 109 ms/cycle, which is slower than the loop it
%   replaces.  Each row is [firstCycle,lastCycle].
runs=zeros(0,2);
iFirst=1;
for i=1:1:size(pulsesListFinal,1)
    if i==size(pulsesListFinal,1) ...
            || pulsesListFinal(i+1,1)~=pulsesListFinal(i,3) ...
            || (pulsesListFinal(i+1,3)-pulsesListFinal(iFirst,1)+1)>maxFrames
        runs(end+1,:)=[iFirst,i]; %#ok<AGROW>
        iFirst=i+1;
    end
end
end