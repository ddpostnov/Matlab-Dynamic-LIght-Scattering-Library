%runExternalCycle  Extract and average stimulus-locked LSCI epochs
%
%   runExternalCycle(s,fNames) processes each *_K_r.mat file in fNames that
%   originates from the basic contrast pipeline.  For a periodic external
%   stimulus (e.g. whisker stimulation in neuro-vascular coupling studies)
%   the function
%       • locates every epoch relative to the stimulus onset
%       • builds baseline / epoch / finale BFI images
%       • rejects artefactual epochs by multiple criteria (baseline drift,
%         peak outliers, image-similarity, time-loss, etc.) with optional
%         manual override
%       • averages accepted epochs into an X×Y×T contrast cube
%       • derives the corresponding BFI time-series
%       • saves preview JPEGs and three MAT-files per recording:
%
%            *_e_K_d.mat   SOURCE   – averaged contrast cube + source.time [T x 1],
%                                     the library-wide orientation (= results.time)
%            *_e_K_r.mat   RESULTS  – masks, epoch metrics, timestamps
%            *_e_K_s.mat   SETTINGS – copy of parameter struct *s*
%            *_rep_epochs.jpg        – epoch-rejection overview
%            *_rep_epoch-average.jpg – averaged BFI time-series
%
%   INPUTS
%     s        parameter structure (fields: stimStartType, stimOffset,
%              epochsN, epochDurationSec, epochBaselineSec, epochFinaleSec,
%              reject*Coefs, maskType, enablelRejectionModification, etc.)
%              • stimDurationSec - HOW LONG THE STIMULUS LASTS, seconds.  This step
%                does not use it: the epoching is driven by the epoch duration, not
%                by the stimulus, so the field is OPTIONAL here.  It is recorded
%                beside the rest of the stimulus geometry because that geometry is
%                protocol, and the steps that DO need the duration - runFitNVC above
%                all, where it is required - then inherit it from this recording's
%                own settings instead of asking for it to be retyped.  Given, it is
%                checked against the epoch.
%     fNames   cell array of full paths to CONTRAST products - *_t_K_r.mat or
%              *_s_K_r.mat, as written by runContrastFromRLS.  A cardiac product
%              (*_c_K_r.mat) is not an input: it is one averaged period and carries
%              no results.timeStamp, so there is no recording clock to place the
%              stimulus on.  Rejected with a named error rather than a missing-field
%              one; the workbench never offers the combination (requires 'contrast').
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUTS
%     None – all results are written to disk (see above).
%
%   EXAMPLE
%     p = defaultExternalCycleParams();
%     files = dir(fullfile(dataRoot,'*_K_r.mat'));
%     runExternalCycle(p, fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     Functions from the LSCI processing library: getEdgeSizeSLSCI,
%     image-mask helpers, etc.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
%
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - STIM PARAMETERS
% s.enablelRejectionModification=1;
% %Type of stimulation start information. Use either 'manual' for a list of
% %starting times to be used in the 'HH:mm:ss.SSS' format,or 'offset' for
% %stimulation that starts at a fixed time from the recording start
% %Note: stim start time corresponds to start of the first epoch, not the
% %stimulation itself
% s.stimStartType='offset';
% %Set the stim offset (for offset stim start mode)
% s.stimOffset=125; %seconds
% %Set the list of stimulation start timestamps (for manual stim start mode)
% s.stimStartAll{1}='09:23:31.346'; %'HH:mm:ss.SSS' %set for multiple files
% as{1} {2} etc
%
% %define epochs (repeated stimulations) parameters
% s.epochsN=20;
% s.epochDurationSec=30;
% %duration of single epoch, seconds
% %time from start of the epoch considered to be baseline. Example [0,5]
% %means that baseline starts with the epoch start (thus 0) and ends in 5
% %seconds.
% s.epochBaselineSec=[0,5];
% s.epochStimStartSec=5; %time when stimulation actually starts
% s.stimDurationSec=5; %how long the stimulation lasts, seconds. Not used here - it
% %is recorded so the response fit can read the protocol instead of being told it again
%
% %time from the end of the epoch when flow is expected to return to baseline
% %Example: [-5,0] means that finale starts 5 seconds before the end of the
% %epoch and ends when the epoch ends.
% s.epochFinaleSec=[-5,0];
%
% s.maskType='selection'; %'basic','cMask','selection';
%
% %ADJUSTED IF NECESSARY - QUALITY CHECK
% s.rejectBlCoef=1; %use Inf to disable rejection by this parameter
% s.rejectEpochCoef=1; %use Inf to disable rejection by this parameter
% s.rejectFinCoef=1; %use Inf to disable rejection by this parameter
% s.rejectPeakCoef=1; %use Inf to disable rejection by this parameter
% s.rejectBlSimCoef=1; %use Inf to disable rejection by this parameter
% s.rejectSimCoef=1; %use Inf to disable rejection by this parameter
% s.rejectTimeLoss=0.5; %allowed time loss due to grabbing faluere in seconds per epoch
% s.rejectFirstEpoch=1; %always reject the first epoch

function runExternalCycle(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_K_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_K_r.mat'').']);
end

%THE STIMULUS DURATION IS RECORDED HERE AND CHECKED HERE, and used by nobody until
%the response fit reads it back (see the header).  The check is worth its four lines:
%a stimulus that outlasts its epoch means the epochs overlap the stimulus and every
%average below is wrong, and hearing that now costs a second instead of sixty files.
if isfield(s,'stimDurationSec') && ~isempty(s.stimDurationSec)
    d=double(s.stimDurationSec);
    if ~(isscalar(d) && isfinite(d) && d>0)
        error('runExternalCycle:stimDuration', ...
            's.stimDurationSec must be a positive number of seconds.');
    end
    if s.epochStimStartSec+d > s.epochDurationSec
        error('runExternalCycle:stimOutlastsEpoch', ...
            ['The stimulus starts %g s into the epoch and lasts %g s, which runs ' ...
             'past the end of the %g s epoch.  Either the stimulus duration or the ' ...
             'epoch duration is wrong.'],s.epochStimStartSec,d,s.epochDurationSec);
    end
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  This step can
% block on an epoch-rejection GUI, so cancel is only checked between files.
rep=reportOpen(s,'External cycle',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        load(getProductPath(s.fName,'d'),'source')
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');
        data=source.data; source=rmfield(source,"data");
        %source.time is [T x 1] library-wide; normalised here anyway so no expression
        %below can be silently orientation-dependent (it used to need a ROW, which no
        %producer writes - see the implicit expansion a few lines down).
        time=source.time(:);

        %This step needs the ABSOLUTE start of the recording to place the stimulus, and
        %only a step that read the raw recording clock can supply it.  runInternalCycle
        %writes no timeStamp - its product is one AVERAGED cardiac period, a cube with no
        %position on the recording timeline - so a *_c_K_r.mat cannot be epoched, by
        %construction rather than by omission (wbStepRegistry: requires 'contrast').
        if ~isfield(results,'timeStamp')
            error('runExternalCycle:noTimeStamp', ...
                ['%s carries no results.timeStamp, so the stimulus cannot be located ' ...
                 'on the recording clock.  Epoch a contrast product (*_t_K_r.mat or ' ...
                 '*_s_K_r.mat); an averaged cardiac cycle (*_c_K_r.mat) has no timeline.'], ...
                s.fName);
        end
        s.rlsStartTime=datetime(results.timeStamp,'ConvertFrom','epochtime','Epoch',datetime(1970,1,1),'TicksPerSecond',1e3,'Format', 'HH:mm:ss.SSS');
        if strcmp(s.stimStartType,'manual')
            s.stimStart=datetime(s.stimStartAll{fidx},'InputFormat','HH:mm:ss.SSS');
        elseif strcmp(s.stimStartType,'offset')
            s.stimStart=s.rlsStartTime+seconds(s.stimOffset);
        end
        %% Accessing timestamps and epoch timing
        %allocate memory for epochs
        epochStartSec=zeros(1,s.epochsN);
        %find time of the stimulation start relative to the recording start
        epochStartSec(1)=((hour(s.stimStart)-hour(s.rlsStartTime))*60+(minute(s.stimStart)-minute(s.rlsStartTime)))*60+(second(s.stimStart)-second(s.rlsStartTime));
        if epochStartSec(1)<0
            error('Recording started later than stimulation?')
        end
        epochStartSec(2:end)=epochStartSec(1) + (1:1:(s.epochsN-1))*s.epochDurationSec;
        epochEndSec=epochStartSec+s.epochDurationSec;
        %[T x 1] against [1 x epochsN] expands to [T x epochsN]; the nearest frame to
        %each epoch edge is then the min DOWN the columns.
        [~,epochStartFrame]=min(abs(time-epochStartSec),[],1);
        [~,epochEndFrame]=min(abs(time-epochEndSec),[],1);

        %% Epoch rejection and calculation of average epoch
        %calculate average epoch time-step and time loss
        timeStep=median(time(2:end)-time(1:end-1));
        timeLoss=[timeStep;time(2:end)-time(1:end-1)]-timeStep;

        %calculate epochs baseline, finale and timeloss values
        switch s.maskType
            case 'basic'
                mask=results.mask;
            case 'cMask'
                mask=results.cMask>1;
            case 'selection'
                f=figure(1);
                f.WindowState='maximized';
                t=tiledlayout(1,2,"TileSpacing",'compact','Padding','compact');
                t1=nexttile(t);
                img=1./(mean(data,3).^2);
                imagesc(img,'Parent', t1)
                clim(prctile(img(:),[5,99]))
                axis image
                t2=nexttile(t);
                idx=find(time>s.stimOffset,1,'first');
                img=mean(1./data(:,:,idx:end).^2-mean(1./data(:,:,1:idx).^2,3,'omitnan'),3,'omitnan');
                imagesc(img,'Parent', t2)
                clim(prctile(img(:),[5,99]))
                axis image
                prctile(img,[5,99])
                mask=roipoly;
                mask=mask.*results.mask;
                delete(f);      % the picker is done with; it used to be left open
            otherwise
                error('Unrecognized maskType')
        end
        data=reshape(data,[],size(data,3));
        tsBFI=mean(data(mask(:)==1,:),1,'omitnan');
        tsBFI=1./(tsBFI.^2);
        data=reshape(data,size(mask,1),size(mask,2),numel(tsBFI));


        baseBFI=zeros(size(data,1),size(data,2),s.epochsN);
        epochBFI=baseBFI;
        epochFinBFI=zeros(size(data,1),size(data,2),s.epochsN);
        timeLossSum=zeros(1,s.epochsN);
        peakBFImean=zeros(1,s.epochsN);
        for i=1:1:s.epochsN
            blStartSec=epochStartSec(i)+s.epochBaselineSec(1);
            blEndSec=epochStartSec(i)+s.epochBaselineSec(2);
            [~,blStartFrame]=min(abs(time-blStartSec));
            [~,blEndFrame]=min(abs(time-blEndSec));
            baseBFI(:,:,i)=1./(squeeze(mean(data(:,:,blStartFrame:blEndFrame).*mask,3)).^2);


            finStartSec=epochEndSec(i)+s.epochFinaleSec(1);
            finEndSec=epochEndSec(i)+s.epochFinaleSec(2);
            [~,finStartFrame]=min(abs(time-finStartSec));
            [~,finEndFrame]=min(abs(time-finEndSec));
            epochFinBFI(:,:,i)=1./(squeeze(mean(data(:,:,finStartFrame:finEndFrame).*mask,3)).^2);

            epochBFI(:,:,i)=1./(mean(data(:,:,epochStartFrame(i):epochEndFrame(i)),3).^2);
            peakBFImean(i)=max(tsBFI(epochStartFrame(i):epochEndFrame(i)));
            timeLossSum(i)=sum(timeLoss(epochStartFrame(i):epochEndFrame(i)));
        end
        epochsBlBFImean=squeeze(sum(baseBFI.*mask,[1,2],'omitnan')./sum(mask(:)));
        epochsFinBFImean=squeeze(sum(epochFinBFI.*mask,[1,2],'omitnan')./sum(mask(:)));
        epochsBFImean=squeeze(sum(epochBFI.*mask,[1,2],'omitnan')./sum(mask(:)));

        %EPOCH REJECTION RULES
        %Change the array size depending on number of rejection rules
        epochsToReject=zeros(7,s.epochsN);
        %Each rule is one ROW, and the rows are the order they are drawn in on the
        %epochs page.  The comments used to name rows 2 and 3 the other way round -
        %the code always matched its own coefficient (rejectFinCoef with the finale
        %mean, rejectEpochCoef with the epoch mean), so only the reading was wrong.
        %based on baseline value deviation
        epochsToReject(1,abs(epochsBlBFImean-median(epochsBlBFImean))>s.rejectBlCoef*std(epochsBlBFImean))=1;
        %based on finale value deviation
        epochsToReject(2,abs(epochsFinBFImean-median(epochsFinBFImean))>s.rejectFinCoef*std(epochsFinBFImean))=1;
        %based on epoch value deviation
        epochsToReject(3,abs(epochsBFImean-median(epochsBFImean))>s.rejectEpochCoef*std(epochsBFImean))=1;
        %based on peak value deviation
        epochsToReject(4,abs(peakBFImean-median(peakBFImean))>s.rejectPeakCoef*std(peakBFImean))=1;
        %based on baseline images similarity between epochs (motion artifacts sensitive)
        blSimilarity=zeros(s.epochsN,s.epochsN);
        for i=1:1:s.epochsN
            for j=1:1:s.epochsN
                blSimilarity(i,j)=ssim(baseBFI(:,:,i),baseBFI(:,:,j));
            end
        end
        blSimilarity=(squeeze(sum(blSimilarity,2))-1)./(s.epochsN-1);
        epochsToReject(5,(median(blSimilarity)-blSimilarity)>s.rejectBlSimCoef*std(blSimilarity))=1;
        %based on epoch images similarity between epochs (motion artifacts sensitive)
        %The twin of the rule above on the epoch image instead of the baseline one, and
        %it had never once fired: the matrix was allocated [1 x epochsN] and the inner
        %loop OVERWROTE it, so only j=epochsN survived, and the reduction then summed
        %that row to a single number.  A scalar makes median(x)-x and std(x) both 0, so
        %the test was 0>0 for every coefficient and s.rejectSimCoef did nothing.
        epochSimilarity=zeros(s.epochsN,s.epochsN);
        for i=1:1:s.epochsN
            for j=1:1:s.epochsN
                epochSimilarity(i,j)=ssim(epochBFI(:,:,i),epochBFI(:,:,j));
            end
        end
        %sum the row, drop the epoch's similarity to itself (which is 1), and the rest
        %is its mean similarity to the other epochs
        epochSimilarity=(squeeze(sum(epochSimilarity,2))-1)./(s.epochsN-1);
        epochsToReject(6,(median(epochSimilarity)-epochSimilarity)>s.rejectSimCoef*std(epochSimilarity))=1;
        %based on the amount of time lost during the epoch
        epochsToReject(7,timeLossSum>=s.rejectTimeLoss)=1;
        %You can add more rejection rules here or remove existing ones

        if (s.rejectFirstEpoch==1)
            epochsToReject(:,1)=1;
        end

        %visualize accepted and rejected epochs
        % The picker is a VISIBLE figure the operator clicks on; the report page
        % below is drawn fresh from the FINAL decisions on the canonical canvas.
        % They used to be one figure, which is why the page's size depended on the
        % monitor and why the window was still open when the next file started.
        if s.enablelRejectionModification==1
            hInt=figure('Name','runExternalCycle - accept or reject epochs','NumberTitle','off');
            axInt=axes('Parent',hInt);
            drawEpochs(axInt,time,tsBFI,timeLoss,epochStartSec,epochEndSec,epochsToReject,s);
            title(axInt,'Epochs (green - accepted, red - rejected). Click on the epoch to change decision. Enter to accept')
            %plot required timeseries, without timevector
            [x,~]=ginput();
            for i=1:1:length(x)
                epochsToReject(:,x(i)>=epochStartSec & x(i)<epochEndSec)=1-max(epochsToReject(:,x(i)>=epochStartSec & x(i)<epochEndSec));
            end
            delete(hInt);
        end

        fh=reportFigure(rep,'epochs','single');
        tl=tiledlayout(fh,1,1,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl);
        drawEpochs(ax,time,tsBFI,timeLoss,epochStartSec,epochEndSec,epochsToReject,s);
        title(ax,'Epochs (green - accepted, red - rejected)')
        reportSave(rep,fh,'epochs');

        %calculate average epoch images
        source.data=zeros(size(data,1),size(data,2),max(epochEndFrame-epochStartFrame));
        epochCounts=zeros(1,1,max(epochEndFrame-epochStartFrame));
        for i=1:1:s.epochsN
            %epochsToReject holds ONE ROW PER RULE, so epoch i's verdict is its whole
            %COLUMN and an epoch is kept only when no rule rejected it - the same
            %reading the report page above draws.  The single subscript this used to
            %carry linear-indexed down column 1 instead, which averaged epochs the page
            %had marked red and, whenever column 1 was fully rejected and there were no
            %more than 7 epochs (rejectFirstEpoch on a short protocol), accumulated
            %nothing at all and left the 0/0 all-NaN cube below.
            if ~any(epochsToReject(:,i))
                for ii=1:1:(epochEndFrame(i)-epochStartFrame(i))
                    source.data(:,:,ii)=source.data(:,:,ii)+data(:,:,epochStartFrame(i)+ii-1);
                    epochCounts(ii)=epochCounts(ii)+1;
                end
            end
        end
        %0/0 is NaN and says nothing; an averaged epoch with no epochs in it is a
        %protocol or rejection-setting problem the operator has to see.
        if ~any(epochCounts)
            error('runExternalCycle:noEpochsAccepted', ...
                ['Every epoch of %s was rejected, so there is nothing to average.  ' ...
                 'Relax the reject* coefficients, or keep an epoch in the picker.'], s.fName);
        end
        source.data=source.data./epochCounts;
        results.time=((0:1:(max(epochEndFrame-epochStartFrame)-1)).*timeStep)';   %[T x 1]
        source.time=results.time;

        % Average contrast time series
        epochTsK=squeeze(sum(source.data.*mask,[1,2])./sum(mask(:)));
        % Simplified conversion to blood flow index
        epochTsBFI=1./(epochTsK.^2);

        fh=reportFigure(rep,'epoch-average','single');
        tl=tiledlayout(fh,1,1,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl);
        plot(ax,results.time,epochTsBFI)
        xlabel(ax,'Time, s')
        ylabel(ax,'BFI')
        axis(ax,'tight')
        grid(ax,'on')
        reportSave(rep,fh,'epoch-average');

        % Save the settings and results
        reportWriting(rep);
        settings.externalCycle=reportSettings(s);
        results.time=source.time;
        % the epoch triplet INSERTS its own flag into the input's name; the other two
        % members are named from the RESULTS one so the three cannot drift apart
        epochName=strrep(s.fName,'_K_r.mat','_e_K_r.mat');
        save(getProductPath(epochName,'d'),'source','-v7.3','-nocompression');
        save(epochName,'results','-v7.3','-nocompression');
        save(getProductPath(epochName,'s'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);

end

% =====================================================================
function drawEpochs(ax,time,tsBFI,timeLoss,epochStartSec,epochEndSec,epochsToReject,s)
%drawEpochs  The epoch overview: BFI trace with one accept/reject bar per rejection
%   rule under each epoch, and the time loss on the right axis.  One function so the
%   interactive picker and the report page cannot drift apart - the old code drew it
%   twice, differing only in the colour of the epoch-start lines.
yyaxis(ax,'left')
plot(ax,time,tsBFI)
hold(ax,'on')
for i=1:1:s.epochsN
    plot(ax,[epochStartSec(i),epochStartSec(i)],[min(tsBFI),max(tsBFI)],'--k')
    for ii=1:1:size(epochsToReject,1)
        y=(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI);
        if epochsToReject(ii,i)==0
            plot(ax,[epochStartSec(i),epochEndSec(i)],[y,y],'-g','LineWidth',1.5)
        else
            plot(ax,[epochStartSec(i),epochEndSec(i)],[y,y],'-r','LineWidth',1.5)
        end
    end
end
plot(ax,[epochEndSec(end),epochEndSec(end)],[min(tsBFI),max(tsBFI)],'--r')
ylabel(ax,'BFI')
hold(ax,'off')
yyaxis(ax,'right')
plot(ax,time,timeLoss)
xlabel(ax,'Time, s')
ylabel(ax,'Time loss, s')
end
