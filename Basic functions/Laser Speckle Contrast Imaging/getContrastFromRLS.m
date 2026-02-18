%getContrastFromRLS - reads .rls data and calculates spatial or temporal
% contrast in batches. This function is used when the amount of available
% RAM memory is insufficient to store the whole LSCI sequence.
%
% Syntax:  [output1,output2,output3] =
% function_name(input1,input2,input3,input4,input5,input6,input7,input8);
%
% Required Inputs:
%   fileName        - path to the .rls file
%   contrastType    - contrast calculation method, use 'spatial' or 'temporal'
%   contrastKernel  - size of the kernel to use for the contrast calculation
%   batchSize       - size of the batch of data to be loaded for analysis
%   procType        - choose the processor type, use:
%                       'cpu', 'gpu', 'fastcpu', 'fastgpu'
%   decimation      - decimation factor. Definines number of frames to be
%                       averaged and kept during decimation. Default 1 (no
%                       decimation). Should be adjusted depending on the
%                       output frequency requirements.
%   saveContrast    - flag to save the processed data in .mat file
%   selectROI       - flag to choose an ROI to process. Only recommended
%                       when parts of field of view have to be discarded.
%
% Outputs:
%   dataLSCI        - processed data as [y,x,t] 3d matrix
%   time            - relative time vector in seconds
%   timeStamp       - absolute timeStamp of first frame
%
% Example:
%
%
% Other m-files required: readRLS, getSLSCI, getTLSCI
% Subfunctions: readRLS, getSLSCI, getTLSCI
% MAT-files required: none
%
% See also:

% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 03-Oct-2025

%------------- BEGIN CODE --------------
function [data,time,timeStamp,intensityMetrics]=getContrastFromRLS(rlsFileName,contrastType, varargin)
p = inputParser;
p.KeepUnmatched = false;
addRequired(p, 'rlsFileName');
addRequired(p, 'contrastType', @(x) any(validatestring(x, {'spatial', 'temporal'})));
addParameter(p, 'kernelSize', [], @isnumeric);
addParameter(p, 'procType', 'gpu', @(x) any(validatestring(x, {'cpu', 'gpu'})));
addParameter(p, 'decimFactor', 1, @isnumeric);
addParameter(p, 'decimMethod', 'leaking', @(x) any(validatestring(x, {'leaking', 'sharp'})));
addParameter(p, 'memoryCoef', 0.8, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
parse(p, rlsFileName, contrastType, varargin{:});

kernelSize  = p.Results.kernelSize;
procType    = p.Results.procType;
decimFactor = p.Results.decimFactor;
decimMethod = p.Results.decimMethod;
memoryCoef  = p.Results.memoryCoef;

if isempty(kernelSize)
    if strcmpi(contrastType, 'temporal')
        kernelSize = 25;
    elseif strcmpi(contrastType, 'spatial')
        kernelSize = 5;
    end
end

if strcmpi(decimMethod, 'sharp')
    if strcmpi(contrastType, 'spatial')
        error('getContrastFromRLS:InvalidMethod', ...
            "The 'sharp' decimation method is only valid for 'temporal' contrast.");
    end
    if strcmpi(contrastType, 'temporal') && decimFactor > 1
        if mod(decimFactor, kernelSize) ~= 0
            error('getContrastFromRLS:InvalidDecimation', ...
                "For 'sharp' temporal contrast, decimFactor (%d) must be an integer multiple of kernelSize (%d).", ...
                decimFactor, kernelSize);
        end
    end
end


[tmp1,timeStamp,stream]=readRLS(rlsFileName,'KeepOpen',true,'FramesToRead',1);
sz = [stream.sizeY,stream.sizeX,stream.sizeT];
outT = floor(sz(3) / decimFactor);
outSz = [sz(1), sz(2), outT];

[~, memAvailRAM] = memory;
memAvailRAM = memAvailRAM.PhysicalMemory.Available;
if memAvailRAM*memoryCoef < (prod(outSz)*4+sz(1)*sz(2)*256*4)
    error('Insufficient memory available for keeping the processed file.');
end

data = zeros(outSz,'single');
time = zeros(outT,1,'double');
intensityMetrics=zeros(sz(1),sz(2),4);
intensityMetrics(:,:,1:2)=1;



batchSize=(memAvailRAM*memoryCoef-prod(outSz)*4)./(sz(1)*sz(2)*(5+6*4/decimFactor)); % logic: dataOUT_single+batch_dataIN_uint8+batch_dataOUT_intermediate_x4_single=memAvailRAM*memoryCoef
batchSize = min(floor(batchSize / decimFactor) * decimFactor,sz(3));
if batchSize == 0, batchSize = decimFactor; end

% histI=zeros(sz(1),sz(2),256,'single');
%     %"Fast" histogram - still takes as long as the entire calculation
%     [Y, X, T] = size(chunkIntensity);
% numPixels = Y * X;
% 
% % 1. Create a spatial offset for every pixel: [1, 257, 513, ...]
% % Reshape to match the spatial dimensions of the batch
% offsets = reshape(uint32(0:numPixels-1) * 256, Y, X);
% 
% % 2. Apply offsets to the temporal stack
% % This shifts each pixel's 0-255 range into a unique global bin range
% shiftedData = uint32(chunkIntensity) + 1 + offsets;
% 
% % 3. Use histcounts on the entire shifted 3D volume at once
% % This is the "fast" global call you were looking for
% globalHist = histcounts(shiftedData, 1:(numPixels * 256 + 1));
% 
% % 4. Reshape the flat global histogram back into per-pixel bins
% batchHist = reshape(globalHist, 256, Y, X); 
% batchHist = permute(batchHist, [2, 3, 1]); % Result: [Y, X, 256] 

batchNum = ceil(sz(3) / batchSize);
curSize = 1;
nextStartIdx = 1;

useOverlap = strcmpi(contrastType, 'temporal') && strcmpi(decimMethod, 'leaking');
tailIntensity = [];
tailTime = [];
H = floor(kernelSize/2);

tic
for i=1:batchNum
    curBatchSize=min(batchSize, sz(3)-(i-1)*batchSize);
    if i==1
        chunkIntensity=zeros(sz(1),sz(2),curBatchSize,class(tmp1));
        chunkTime=zeros(curBatchSize,1,'double');
        chunkIntensity(:,:,1)=tmp1;
        chunkTime(1)=timeStamp;
        [chunkIntensity(:,:,2:end),chunkTime(2:end),stream]=readRLS(stream,'FramesToRead',curBatchSize-1);
    else
        [chunkIntensity,chunkTime,stream]=readRLS(stream,'FramesToRead',curBatchSize);
    end
    intensityMetrics(:,:,1)=intensityMetrics(:,:,1)-double(sum(chunkIntensity==0,3))./sz(3);
    intensityMetrics(:,:,2)=intensityMetrics(:,:,2)-double(sum(chunkIntensity==intmax(class(chunkIntensity)),3))./sz(3);
    intensityMetrics(:,:,3)=intensityMetrics(:,:,3)+sum(single(chunkIntensity),3)./sz(3);
    intensityMetrics(:,:,4)=intensityMetrics(:,:,4)+sum(single(chunkIntensity).*single(chunkIntensity),3)./sz(3);


    if ~isempty(tailIntensity)
        chunkIntensity = cat(3, tailIntensity, chunkIntensity);
        chunkTime = cat(1, tailTime, chunkTime);
    end
    [tmp1,tmp2]=getK(chunkIntensity,contrastType,'time',chunkTime,'kernelSize',kernelSize,'procType', procType,'decimFactor',decimFactor,'decimMethod',decimMethod,'memoryCoef',memoryCoef);

    if useOverlap
        startIdx = nextStartIdx;
        if i == batchNum
            endIdx = size(tmp1, 3);
        else            
            validInputLimit = size(chunkIntensity, 3) - H;
            endIdx = floor(validInputLimit / decimFactor);
        end
        nextFrameStart = (endIdx * decimFactor) + 1;
        rawStartNeeded = nextFrameStart - H;
        alignedStart = floor((rawStartNeeded - 1) / decimFactor) * decimFactor + 1;
        alignedStart = max(1, alignedStart);
        framesSavedInPrevBatch = endIdx - (alignedStart - 1) / decimFactor;
        nextStartIdx = framesSavedInPrevBatch + 1;
        tailIntensity = chunkIntensity(:, :, alignedStart:end);
        tailTime = chunkTime(alignedStart:end);
    else
        startIdx = 1;
        endIdx = size(tmp1, 3);
    end

    numToSave = endIdx - startIdx + 1;
    if numToSave > 0
        saveEnd = curSize + numToSave - 1;
        if saveEnd > outT
            numToSave = outT - curSize + 1;
            endIdx = startIdx + numToSave - 1;
        end
        data(:, :, curSize : curSize+numToSave-1) = tmp1(:, :, startIdx:endIdx);
        time(curSize : curSize+numToSave-1) = tmp2(startIdx:endIdx);
        curSize = curSize + numToSave;
    end
    fprintf('Batch %d/%d processed. Elapsed: %.2fs\n', i, batchNum, toc);
end
time=(time-time(1))./1000; %conversion to seconds
fclose(stream.fId);
end

%     %"Fast" histogram - still takes as long as the entire calculation
%     [Y, X, T] = size(chunkIntensity);
% numPixels = Y * X;
% 
% % 1. Create a spatial offset for every pixel: [1, 257, 513, ...]
% % Reshape to match the spatial dimensions of the batch
% offsets = reshape(uint32(0:numPixels-1) * 256, Y, X);
% 
% % 2. Apply offsets to the temporal stack
% % This shifts each pixel's 0-255 range into a unique global bin range
% shiftedData = uint32(chunkIntensity) + 1 + offsets;
% 
% % 3. Use histcounts on the entire shifted 3D volume at once
% % This is the "fast" global call you were looking for
% globalHist = histcounts(shiftedData, 1:(numPixels * 256 + 1));
% 
% % 4. Reshape the flat global histogram back into per-pixel bins
% batchHist = reshape(globalHist, 256, Y, X); 
% batchHist = permute(batchHist, [2, 3, 1]); % Result: [Y, X, 256] 


%
%
%
%
% tic
% [tmp1,tmp2,s]=readRLS("C:\Dropbox\Work\Data\fake.rls",'KeepOpen',true,'FramesToRead',1);
% data2=zeros(s.sizeY,s.sizeX,s.sizeT,'uint8');
% data2(:,:,1)=tmp1;
% timeStamps=zeros(s.sizeT,1,'uint64');
% timeStamps(1)=tmp2;
% batchSize=1000;
% for i=2:batchSize:s.sizeT
% framesToRead=min(batchSize,s.sizeT-i+1);
% [data2(:,:,i:i+framesToRead-1),timeStamps(i:i+framesToRead-1),s]=readRLS('Stream',s,'FramesToRead',framesToRead);
% end
% fclose(s.fId);
% toc
%
%
%
%
%
% %read metadata
% [dataLSCI,~,timeStamp,rawFramesN]=readRLS(rlsFileName,0,1);
% disp(['Number of frames to read: ',num2str(rawFramesN)]);
% ROI=[1,size(dataLSCI,1);1,size(dataLSCI,2)];
%
% %Time conversion (to seconds) factor:
% timeScale=1/1000; %for timestamps written in ms
%
%
% %batch size test
% if rawBatchSize<decimation
%     error('Batch size smaller than decimation factor, aborting');
% elseif mod(rawBatchSize,decimation)~=0
%     rawBatchSize=single(decimation*floor(rawBatchSize/decimation));
%     warning(['Changing the batch size to be a multiple of decimation. New batch size is ',num2str(rawBatchSize)])
% end
% outBatchSize=floor(rawBatchSize/decimation);
%
%
% %ROI selection if required
% if selectROI==1
%     disp('ROI selection was opened in a separate figure');
%     subdata=getSLSCI(dataLSCI,5,'fastcpu');
%     figure
%     sgtitle(["Select ROI to process","(double-click to confirm)"])
%     subplot(1,2,1)
%     img=dataLSCI;
%     imagesc(img)
%     clim([prctile(img(:),1),prctile(img(:),99)])
%     title('Raw data')
%     axis image
%     subplot(1,2,2)
%     img=subdata;
%     imagesc(img)
%     clim([prctile(img(img(:)>0 & img(:)<1),1),prctile(img(img(:)>0 & img(:)<1),99)])
%     title('Spatial contrast')
%     axis image
%     roi=images.roi.Rectangle(gca,'Position',[size(img,1)/4,size(img,2)/4,size(img,1)/2,size(img,2)/2]);
%     wait(roi);
%     ROI=round([roi.Position(2),roi.Position(2)+roi.Position(4)-1;roi.Position(1),roi.Position(1)+roi.Position(3)-1]);
% end
%
% %contrast calculation
% tic
% switch contrastType
%     case 'spatial'
%         procFramesN=floor(rawFramesN/decimation);
%         dataLSCI=zeros(ROI(1,2)-ROI(1,1)+1,ROI(2,2)-ROI(2,1)+1,procFramesN,'single');
%         %Uncertainty mask
%         trustMatrix=zeros(size(dataLSCI,1),size(dataLSCI,2),4);
%         trustMatrix(:,:,1:3)=rawFramesN;
%         time=zeros(1,procFramesN);
%         batchesN=ceil(rawFramesN/rawBatchSize);
%         for i=1:1:batchesN
%             framesToRead=min(rawFramesN-(i-1)*rawBatchSize,rawBatchSize);
%             [subdata,~,timeStamps,~]=readRLS(rlsFileName,(i-1)*rawBatchSize,framesToRead,ROI);
%             trustMatrix(:,:,1)=trustMatrix(:,:,1)-sum(subdata==0,3);
%             trustMatrix(:,:,2)=trustMatrix(:,:,2)-sum(subdata==intmax(class(subdata)),3);
%             trustMatrix(:,:,4)=trustMatrix(:,:,4)+mean(single(subdata),3)./batchesN; % saving intensity
%             subdata=getSLSCI(subdata,contrastKernel,procType);
%             trustMatrix(:,:,3)=trustMatrix(:,:,3)-sum(subdata==0,3);
%             subdata=movmean(subdata,[0,decimation-1],3);
%             timeStamps=double(timeStamps-timeStamp);
%             timeStamps=movmean(timeStamps,[0,decimation-1]);
%             time(((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(framesToRead/decimation)))=timeStamps(1:decimation:floor(end./decimation)*decimation);
%             dataLSCI(:,:,(i-1)*outBatchSize+1:1:(i-1)*outBatchSize+floor(framesToRead/decimation))= subdata(:,:,1:decimation:floor(end./decimation)*decimation);
%             disp(['Processed batch ',num2str(i),' out of ',num2str(batchesN),'. Time elapsed ',num2str(toc)]);
%         end
%     case 'temporal'
%         procFramesN=floor((rawFramesN-contrastKernel+1)/decimation);
%         dataLSCI=zeros(ROI(1,2)-ROI(1,1)+1,ROI(2,2)-ROI(2,1)+1,procFramesN,'single');
%         %Uncertainty mask
%         trustMatrix=zeros(size(dataLSCI,1),size(dataLSCI,2),4);
%         trustMatrix(:,:,1:3)=rawFramesN;
%         time=zeros(1,procFramesN);
%         batchesN=ceil((rawFramesN-contrastKernel+1)/rawBatchSize);
%         for i=1:1:batchesN
%             framesToRead=min(rawFramesN-(i-1)*rawBatchSize,rawBatchSize+contrastKernel-1);
%             [subdata,~,timeStamps,~]=readRLS(rlsFileName,(i-1)*rawBatchSize,framesToRead,ROI);
%             trustMatrix(:,:,1)=trustMatrix(:,:,1)-sum(subdata==0,3);
%             trustMatrix(:,:,2)=trustMatrix(:,:,2)-sum(subdata==intmax(class(subdata)),3);
%             trustMatrix(:,:,4)=trustMatrix(:,:,4)+mean(single(subdata),3)./batchesN; % saving intensity
%             subdata=getTLSCI(subdata,contrastKernel,procType);
%             trustMatrix(:,:,3)=trustMatrix(:,:,3)-sum(subdata==0,3);
%             timeStamps=double(timeStamps-timeStamp);
%             timeStamps=movmean(timeStamps,[0,decimation-1]);
%
%             if mod(decimation,contrastKernel)==0
%                 subdata=subdata(:,:,1:contrastKernel:end);
%                 timeStamps=timeStamps(1:contrastKernel:end);
%                 subdata=movmean(subdata,[0,decimation./contrastKernel-1],3,'Endpoints','shrink');
%                 dataLSCI(:,:,((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(size(subdata,3)/(decimation./contrastKernel))))= subdata(:,:,1:(decimation./contrastKernel):floor(end./(decimation./contrastKernel))*(decimation./contrastKernel));
%                 time(((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(size(subdata,3)/(decimation./contrastKernel))))=timeStamps(1:(decimation./contrastKernel):floor(end./(decimation./contrastKernel))*(decimation./contrastKernel));
%             else
%                 timeStamps=timeStamps(1:end-contrastKernel+1);
%                 subdata=movmean(subdata,[0,decimation-1],3,'Endpoints','shrink');
%                 subdata=subdata(:,:,1:end-contrastKernel+1);
%                 dataLSCI(:,:,((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(size(subdata,3)/decimation)))= subdata(:,:,1:decimation:floor(end./decimation)*decimation);
%                 time(((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(size(subdata,3)/decimation)))=timeStamps(1:decimation:floor(end./decimation)*decimation);
%             end
%
%             disp(['Processed batch ',num2str(i),' out of ',num2str(batchesN),'. Time elapsed ',num2str(toc)]);
%         end
%         if size(dataLSCI,3)>procFramesN
%             dataLSCI(:,:, procFramesN+1:end)=[];
%             time(procFramesN+1:end)=[];
%         end
% end
% time=time*timeScale; %convert time to seconds
%
% trustMatrix(:,:,[1,2,3])=trustMatrix(:,:,[1,2,3])./rawFramesN;
%
% settings.contrastCalculation.contrastType=contrastType;
% settings.contrastCalculation.contrastKernel=contrastKernel;
% settings.contrastCalculation.rlsFileName=rlsFileName;
% settings.contrastCalculation.rawBatchSize=rawBatchSize;
% settings.contrastCalculation.procType=procType;
% settings.contrastCalculation.decimation=decimation;
% settings.contrastCalculation.selectROI=selectROI;
% settings.contrastCalculation.ROI=ROI;
%
%
% disp('Processing has been compelted.');
% end
% %------------- END OF CODE --------------
% %Comments: For the high-speed cameras and other cases with abnormal
% %timestamps - pay attention to the time scaling,