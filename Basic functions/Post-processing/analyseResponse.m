% analyseResponse - calculates variety of features of signal increase.
% Syntax:  output1 = function_name(requiredInput1,requiredInput2,requiredInput3,requiredInput4,requiredInput5,requiredInput6,requiredInput7)
%
% Required inputs:
%   data          - data to analyze. Can be either 1d,2d or 3d, but time
%                   has to be the last dimension
%   timeIni       - time vector
%   baselineSec   - vector of size 2, which corresponds to the first and
%                   the last seconds of the baseline (e.g. [0,4])
%   method        - method for the response start and end identification:
%                   'crossing' finds the last value before the peak and the
%                   first after the peak that cross the baseline, 'minima'
%                   finds the last and the first local minima connected to
%                   the peak, which are below the designeted baseline
%                   threshod
%   coef          - coefficient that offsets the baseline threshold based
%                   on the baseline standard deviation. Use 0 for default
%                   mean baseline value.
%   filtSec       - filter size in seconds, used to smooth the numerical
%                   derivative in the maximum acceleration and
%                   decelleration calculations.
%   thresholds    - vector of additional thresholds for area and onset
%                   calculations. E.g. [0.1,0.3,0.9], but exact value
%                   depends on the thresholdType.
%   thresholdsType- type of the additional thresholds. 'magnitude'
%                   identifies thresholds relative to the maximum increase
%                   (as baseline+(max-baseline)*threshold), 'rBFI' uses
%                   relative blood flow index value (as baseline*threshold)
%                   , 'totalRBFI' uses total relative blood flow index
%                   change as a reference for thresholding
%
%
% Outputs:
%  r                     - structure containing all results:
%   r.baseBFI            - baseline BFI
%   r.baseSTD            - standard deviation of the baseline BFI
%   r.maxRBFI            - maximum relative BFI increase
%   r.maxDelay           - delay (from the baseline) until BFI increase reaches max
%   r.totalChangeBFI     - total BFI change during the post baseline period
%   r.totalChangeRBFI    - total rBFI change during the post baseline period
%   r.responseChangeBFI  - total BFI change during the identified response period
%   r.responseChangeRBFI - total rBFI change during the identified response period
%   r.responseDuration   - duration of the identified response period
%   r.responseDelay      - delay (from the baseline) of the identified response period
%   r.accelMaxBFI        - BFI change per second at maximum accelaration point
%   r.accelMaxRBFI       - rBFI change per second at maximum accelaration point
%   r.accelMaxDelay      - delay until BFI increase reaches maximum accelearation
%   r.deccelMaxBFI       - BFI change per second at maximum deccelaration point
%   r.deccelMaxRBFI      - rBFI change per second at maximum deccelaration point
%   r.deccelMaxDelay     - delay until BFI increase reaches maximum deccelearation
%   r.thresholdValue     - identified threshold values for thresholds vector
%   r.thresholdChangeBFI - total BFI change within the threshold limits
%   r.thresholdChangeRBFI- total rBFI change within the threshold limits
%   r.thresholdDuration  - duration of BFI being within the threshold limits
%   r.thresholdDelay     - delay until BFI reaches the threshold value
%
% Example:
%    [results]=analyseResponse(data,time,[0,4],'crossing',1,1,[0.1,0.9],'magnitude');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 10-Feb-2023

function [r]=analyseResponse(data,timeIni,baselineSec,method,coef,filtSec,thresholds,thresholdsType)
fs=1./median(diff(timeIni));
filtPoints=max(1,floor(filtSec.*fs));
[~,baselineFrames(1)]=min(squeeze(abs(squeeze(timeIni')-baselineSec(1))),[],1);
[~,baselineFrames(2)]=min(squeeze(abs(squeeze(timeIni')-baselineSec(2))),[],1);

switch sum((size(data))>1)
    case 1
        data=reshape(data,1,1,length(data(:)));
    case 2
        data=reshape(data,1,size(data,1),size(data,2));
end

r.baseBFI=mean(data(:,:,baselineFrames(1):baselineFrames(2)),3);
r.baseSTD=std(data(:,:,baselineFrames(1):baselineFrames(2)),0,3);
data=data(:,:,baselineFrames(2)+1:end);
time=timeIni(baselineFrames(2)+1:end);
[maxVal,maxIdx]=max(data,[],3);
r.maxRBFI=maxVal./r.baseBFI;
r.maxDelay=time(maxIdx)-timeIni(baselineFrames(2));


r.totalChangeBFI=trapz(time,data-r.baseBFI,3);
r.totalChangeRBFI=trapz(time,data./r.baseBFI-1,3);


responseFrames=zeros(size(maxVal,1),size(maxVal,2),2);
r.responseChangeBFI=zeros(size(maxVal));
r.responseChangeRBFI=zeros(size(maxVal));
r.responseDuration=zeros(size(maxVal));
r.responseDelay=zeros(size(maxVal));
r.accelMaxBFI=zeros(size(maxVal));
r.accelMaxRBFI=zeros(size(maxVal));
r.deccelMaxBFI=zeros(size(maxVal));
r.deccelMaxRBFI=zeros(size(maxVal));
r.accelMaxDelay=zeros(size(maxVal));
r.deccelMaxDelay=zeros(size(maxVal));

% identify start and end points of the response curve using 3 different
% methods: 1 - find the first and the last points crossing the r.baseBFI
% value connected to the maximum response, 2 - find the first rise and the
% last drop connected to the maximum value and taking into the account the
% r.baseBFI variability
thresh=min(r.baseBFI+r.baseSTD*coef,maxVal);
for y=1:1:size(data,1)
    for x=1:1:size(data,2)
        ts=squeeze(data(y,x,:));
        if ~isnan(sum(ts)) && sum(ts)<inf
            switch method
                case 'crossing'
                    idxL=find(ts(1:maxIdx(y,x))<r.baseBFI(y,x),1,'last');
                    idxL(isempty(idxL))=1;
                    responseFrames(y,x,1)=idxL;
                    idxR=find(ts(maxIdx(y,x):end)<r.baseBFI(y,x),1,'first');
                    if ~isempty(idxR)
                        responseFrames(y,x,2)=idxR+maxIdx(y,x)-1;
                    else
                        responseFrames(y,x,2)=length(ts);
                    end
                case 'minima'
                    [~,locs]=findpeaks(-ts);
                    idxL=find(ts(1:maxIdx(y,x))<thresh(y,x),1,'last');
                    if ~isempty(idxL) && ~isempty(locs) && ~isempty(find(locs<=idxL,1,'last'))
                        responseFrames(y,x,1)=locs(find(locs<=idxL,1,'last'));
                    else
                        responseFrames(y,x,1)=1;
                    end

                    idxR=find(ts(maxIdx(y,x):end)<thresh(y,x),1,'first');
                    if ~isempty(idxR) && ~isempty(locs) && ~isempty(find(locs>=idxR+maxIdx(y,x)-1,1,'first'))
                        responseFrames(y,x,2)=locs(find(locs>=idxR+maxIdx(y,x)-1,1,'first'));
                    else
                        responseFrames(y,x,2)=length(ts);
                    end
            end

            r.responseChangeBFI(y,x)=trapz(time(responseFrames(y,x,1):responseFrames(y,x,2)),ts(responseFrames(y,x,1):responseFrames(y,x,2))-r.baseBFI(y,x));
            r.responseChangeRBFI(y,x)=trapz(time(responseFrames(y,x,1):responseFrames(y,x,2)),ts(responseFrames(y,x,1):responseFrames(y,x,2))./r.baseBFI(y,x)-1);
            r.responseDuration(y,x)=time(responseFrames(y,x,2))-time(responseFrames(y,x,1));
            r.responseDelay(y,x)=time(responseFrames(y,x,1))-timeIni(baselineFrames(2));


            tsDiff=[0;diff(ts)];
            tsDiff=smooth(tsDiff,filtPoints);
            if length(responseFrames(y,x,1):maxIdx(y,x))>3
            [pksP,locsP]=findpeaks(tsDiff(responseFrames(y,x,1):maxIdx(y,x)),'SortStr', 'descend', 'NPeaks', 1);
            if ~isempty(pksP)
                locsP=locsP+responseFrames(y,x,1)-1;
                r.accelMaxBFI(y,x)=pksP*fs;
                r.accelMaxRBFI(y,x)=pksP./r.baseBFI(y,x)*fs;
                r.accelMaxDelay(y,x)=time(locsP)-timeIni(baselineFrames(2));
            end
            end

            if length(maxIdx(y,x):responseFrames(y,x,2))>3
            [pksN,locsN]=findpeaks(-tsDiff((maxIdx(y,x):responseFrames(y,x,2))),'SortStr', 'descend', 'NPeaks', 1);
            if ~isempty(pksN)
                locsN=locsN+maxIdx(y,x)-1;
                r.deccelMaxBFI(y,x)=-pksN*fs;
                r.deccelMaxRBFI(y,x)=-pksN./r.baseBFI(y,x)*fs;
                r.deccelMaxDelay(y,x)=time(locsN)-timeIni(baselineFrames(2));
            end
            end
        end
    end
end


thresholdFrames=zeros(size(maxVal,1),size(maxVal,2),length(thresholds),2);
r.thresholdChangeBFI=zeros(size(maxVal,1),size(maxVal,2),length(thresholds));
r.thresholdChangeRBFI=zeros(size(maxVal,1),size(maxVal,2),length(thresholds));
r.thresholdDuration=zeros(size(maxVal,1),size(maxVal,2),length(thresholds));
r.thresholdDelay=zeros(size(maxVal,1),size(maxVal,2),length(thresholds));
r.thresholdValue=zeros(size(maxVal,1),size(maxVal,2),length(thresholds));

% identify crossing points with set thresholds (relative to the max or
% r.baseBFI)
for i=1:1:length(thresholds)
if strcmp(thresholdsType,'totalRBFI')
    for y=1:1:size(data,1)
        for x=1:1:size(data,2)
            ts=squeeze(data(y,x,:));
            csTs=cumsum((ts./r.baseBFI(y,x)-1));
            r.thresholdValue(y,x,i)=max(csTs).*thresholds(i);
            idxL=find(csTs./max(csTs)<=thresholds(i),1,'last');
            idxL(isempty(idxL))=length(ts);
            idxL(idxL==length(ts))=idxL-1;
            thresholdFrames(y,x,i,1)=idxL;
            thresholdFrames(y,x,i,2)=length(ts);
            r.thresholdChangeBFI(y,x,i)=trapz(time(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2)),ts(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2))-r.baseBFI(y,x));
            r.thresholdChangeRBFI(y,x,i)=trapz(time(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2)),ts(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2))./r.baseBFI(y,x)-1);
            r.thresholdDuration(y,x,i)=time(thresholdFrames(y,x,i,2))-time(thresholdFrames(y,x,i,1));
            r.thresholdDelay(y,x,i)=time(thresholdFrames(y,x,i,1))-timeIni(baselineFrames(2));

        end            
    end

else
    switch thresholdsType
        case 'magnitude'
            thresh=min(r.baseBFI+(maxVal-r.baseBFI).*thresholds(i),maxVal);
        case 'rBFI'
            thresh=min(r.baseBFI+r.baseBFI.*thresholds(i),maxVal);
    end
    r.thresholdValue(:,:,i)=tresh;
    for y=1:1:size(data,1)
        for x=1:1:size(data,2)
            ts=squeeze(data(y,x,:));
            idxL=find(ts(1:maxIdx(y,x))<=thresh(y,x),1,'last');
            idxL(isempty(idxL))=1;
            thresholdFrames(y,x,i,1)=idxL;
            idxR=find(ts(maxIdx(y,x):end)<=thresh(y,x),1,'first');
            if ~isempty(idxR)
                thresholdFrames(y,x,i,2)=idxR+maxIdx(y,x)-1;
            else
                thresholdFrames(y,x,i,2)=length(ts);
            end
            r.thresholdChangeBFI(y,x,i)=trapz(time(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2)),ts(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2))-r.baseBFI(y,x));
            r.thresholdChangeRBFI(y,x,i)=trapz(time(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2)),ts(thresholdFrames(y,x,i,1):thresholdFrames(y,x,i,2))./r.baseBFI(y,x)-1);
            r.thresholdDuration(y,x,i)=time(thresholdFrames(y,x,i,2))-time(thresholdFrames(y,x,i,1));
            r.thresholdDelay(y,x,i)=time(thresholdFrames(y,x,i,1))-timeIni(baselineFrames(2));

        end            
    end
end
end


r.baseBFI=squeeze(r.baseBFI);
r.baseSTD=squeeze(r.baseSTD);
r.maxRBFI=squeeze(r.maxRBFI);
r.maxDelay=squeeze(r.maxDelay);
r.totalChangeBFI=squeeze(r.totalChangeBFI);
r.totalChangeRBFI=squeeze(r.totalChangeRBFI);
r.responseChangeBFI=squeeze(r.responseChangeBFI);
r.responseChangeRBFI=squeeze(r.responseChangeRBFI);
r.responseDuration=squeeze(r.responseDuration);
r.responseDelay=squeeze(r.responseDelay);
r.accelMaxBFI=squeeze(r.accelMaxBFI );
r.accelMaxRBFI=squeeze(r.accelMaxRBFI );
r.accelMaxDelay=squeeze(r.accelMaxDelay );
r.deccelMaxBFI=squeeze(r.deccelMaxBFI);
r.deccelMaxRBFI=squeeze(r.deccelMaxRBFI);
r.deccelMaxDelay=squeeze(r.deccelMaxDelay);
r.thresholdValue=squeeze(r.thresholdValue);
r.thresholdChangeBFI=squeeze(r.thresholdChangeBFI);
r.thresholdChangeRBFI=squeeze(r.thresholdChangeRBFI);
r.thresholdDuration=squeeze(r.thresholdDuration);
r.thresholdDelay=squeeze(r.thresholdDelay);

end
