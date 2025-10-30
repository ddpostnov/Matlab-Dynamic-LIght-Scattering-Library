%getExternalCycle  Extract and average stimulus-locked LSCI epochs
%
%   getExternalCycle(s,fNames) processes each *_K_d.mat file in fNames that
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
%            *_e_K_d.mat   SOURCE   – averaged contrast cube (s,results.time)
%            *_e_K_r.mat   RESULTS  – masks, epoch metrics, timestamps
%            *_e_K_s.mat   SETTINGS – copy of parameter struct *s*
%            *_ec.jpg      – epoch-rejection overview
%            *_ec2.jpg     – averaged BFI time-series
%
%   INPUTS
%     s        parameter structure (fields: stimStartType, stimOffset,
%              epochsN, epochDurationSec, epochBaselineSec, epochFinaleSec,
%              reject*Coefs, maskType, enablelRejectionModification, etc.)
%     fNames   cell array of full paths to *_K_d.mat files (same naming
%              convention as produced by getContrast).
%
%   OUTPUTS
%     None – all results are written to disk (see above).
%
%   EXAMPLE
%     p = defaultExternalCycleParams();
%     files = dir(fullfile(dataRoot,'*_K_d.mat'));
%     getExternalCycle(p, fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     Functions from the LSCI processing library: getEdgeSizeSLSCI,
%     image-mask helpers, etc.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: Header generated with ChatGPT and may contain minor inconsistencies.
%   ----------------------------------------------------------------------

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
% stimStart{1}='09:23:31.346'; %'HH:mm:ss.SSS'
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

function getExternalCycle(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

for fidx=1:1:length(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
     load(s.fName,'source')
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');
    data=source.data; source=rmfield(source,"data");
    time=source.time;

    s.rlsStartTime=datetime(results.timeStamp,'ConvertFrom','epochtime','Epoch',datetime(1970,1,1),'TicksPerSecond',1e3,'Format', 'HH:mm:ss.SSS');
    if strcmp(s.stimStartType,'manual')
        s.stimStart=datetime(stimStart{fidx},'InputFormat','HH:mm:ss.SSS');
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
    [~,epochStartFrame]=min(abs(time'-epochStartSec),[],1);
    [~,epochEndFrame]=min(abs(time'-epochEndSec),[],1);

    %% Epoch rejection and calculation of average epoch
    %calculate average epoch time-step and time loss
    timeStep=median(time(2:end)-time(1:end-1));
    timeLoss=[timeStep,time(2:end)-time(1:end-1)]-timeStep;

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
    %based on baseline value deviation (
    epochsToReject(1,abs(epochsBlBFImean-median(epochsBlBFImean))>s.rejectBlCoef*std(epochsBlBFImean))=1;
    %based on epoch value deviation
    epochsToReject(2,abs(epochsFinBFImean-median(epochsFinBFImean))>s.rejectFinCoef*std(epochsFinBFImean))=1;
    %based on finale value deviation
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
    %based on images similarity to baseline within epochs (motion artifacts sensitive)
    epochSimilarity=zeros(1,s.epochsN);
    for i=1:1:s.epochsN
        for j=1:1:s.epochsN
            epochSimilarity(i)=ssim(epochBFI(:,:,i),epochBFI(:,:,j));
        end
    end
    epochSimilarity=(squeeze(sum(epochSimilarity,2))-1)./(s.epochsN-1);
    epochsToReject(6,(median(epochSimilarity)-epochSimilarity)>s.rejectSimCoef*std(epochSimilarity))=1;
    %based on the amount of time lost during the epoch
    epochsToReject(7,timeLossSum>=s.rejectTimeLoss)=1;
    %You can add more rejection rules here or remove existing ones

    if (s.rejectFirstEpoch==1)
        epochsToReject(:,1)=1;
    end

    %visualize accepted and rejected epochs
    h=figure;
    yyaxis left
    plot(time,tsBFI)
    hold on
    for i=1:1:s.epochsN
        plot([epochStartSec(i),epochStartSec(i)],[min(tsBFI),max(tsBFI)],'--k')

        for ii=1:1:size(epochsToReject,1)
            if epochsToReject(ii,i)==0
                plot([epochStartSec(i),epochEndSec(i)],[(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI),(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI)],'-g','LineWidth',1.5)
            else
                plot([epochStartSec(i),epochEndSec(i)],[(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI),(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI)],'-r','LineWidth',1.5)
            end
        end
    end
    plot([epochEndSec(end),epochEndSec(end)],[min(tsBFI),max(tsBFI)],'--r')
    ylabel('BFI')
    hold off
    yyaxis right
    plot(time,timeLoss)
    xlabel('Time, s')
    ylabel('Time loss, s')
    if s.enablelRejectionModification==1
        title('Epochs (green - accepted, red - rejected). Click on the epoch to change decision. Enter to accept')
        %plot required timeseries, without timevector
        [x,~]=ginput();
        for i=1:1:length(x)
            epochsToReject(:,x(i)>=epochStartSec & x(i)<epochEndSec)=1-max(epochsToReject(:,x(i)>=epochStartSec & x(i)<epochEndSec));
        end
        yyaxis left
        plot(time,tsBFI)
        hold on
        for i=1:1:s.epochsN
            plot([epochStartSec(i),epochStartSec(i)],[min(tsBFI),max(tsBFI)],'--r')
            for ii=1:1:size(epochsToReject,1)
                if epochsToReject(ii,i)==0
                    plot([epochStartSec(i),epochEndSec(i)],[(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI),(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI)],'-g','LineWidth',1.5)
                else
                    plot([epochStartSec(i),epochEndSec(i)],[(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI),(1+ii*.015)*min(tsBFI)-ii*0.015*max(tsBFI)],'-r','LineWidth',1.5)
                end
            end
        end
        plot([epochEndSec(end),epochEndSec(end)],[min(tsBFI),max(tsBFI)],'--r')
    end
    title('Epochs (green - accepted, red - rejected)')
    print(h,strrep(s.fName,'.mat','_ec.jpg'), '-djpeg', '-r300');

    %calculate average epoch images
    source.data=zeros(size(data,1),size(data,2),max(epochEndFrame-epochStartFrame));
    epochCounts=zeros(1,1,max(epochEndFrame-epochStartFrame));
    for i=1:1:s.epochsN
        if epochsToReject(i)==0
            for ii=1:1:(epochEndFrame(i)-epochStartFrame(i))
                source.data(:,:,ii)=source.data(:,:,ii)+data(:,:,epochStartFrame(i)+ii-1);
                epochCounts(ii)=epochCounts(ii)+1;
            end
        end
    end
    source.data=source.data./epochCounts;
    results.time=(0:1:(max(epochEndFrame-epochStartFrame)-1)).*timeStep;
    source.time=results.time;

    % Average contrast time series
    epochTsK=squeeze(sum(source.data.*mask,[1,2])./sum(mask(:)));
    % Simplified conversion to blood flow index
    epochTsBFI=1./(epochTsK.^2);

    h=figure;
    plot(results.time,epochTsBFI)
    xlabel('Time, s')
    ylabel('BFI')
    axis tight
    grid on
    set(gcf,'Color','w')
    drawnow
    print(h,strrep(s.fName,'.mat','_ec2.jpg'), '-djpeg', '-r300');

    % Save the settings and results
    disp('Saving the results');
    settings.externalCycle=s;
    results.time=source.time;
    save(strrep(s.fName,'_K_d.mat','_e_K_d.mat'),'source','-v7.3');
    save(strrep(s.fName,'_K_d.mat','_e_K_r.mat'),'results','-v7.3');
    save(strrep(s.fName,'_K_d.mat','_e_K_s.mat'),'settings','-v7.3');
end
end


