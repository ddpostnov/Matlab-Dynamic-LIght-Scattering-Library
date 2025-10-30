%getContrastFromMRAW - reads .MRAW data and calculates spatial or temporal
% contrast in batches. This function is used when the amount of available
% RAM memory is insufficient to store the whole LSCI sequence. Requires
% cihx file to be present in the same folder as the MRAW
%
% Syntax:  [output1,output2,output3] =
% function_name(input1,input2,input3,input4,input5,input6,input7,input8);
%
% Required Inputs:
%   fileName        - path to the .mraw file
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
%   timeStamp       - dummy timeStamp (set 0 for MRAW files)
%
% Example:
%
%
% Other m-files required: readMRAW, getSLSCI, getTLSCI
% Subfunctions: readMRAW, getSLSCI, getTLSCI
% MAT-files required: none
%
% See also:

% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 25-July-2022

%------------- BEGIN CODE --------------
function [dataLSCI,time,timeStamp,trustMatrix,settings]=getContrastFromMRAW(mrawFileName,contrastType,contrastKernel,rawBatchSize,procType,decimation,saveContrast,selectROI,rawAveraging)
timeStamp=0;
%read metadata
[dataLSCI,fps,~,rawFramesN]=readMRAW(mrawFileName,0,1);
disp(['Number of frames to read: ',num2str(rawFramesN)]);
ROI=[1,size(dataLSCI,1);1,size(dataLSCI,2)];

%batch size test
if rawBatchSize<decimation
    error('Batch size smaller than decimation factor, aborting');
elseif mod(rawBatchSize,decimation)~=0
    rawBatchSize=single(decimation*floor(rawBatchSize/decimation));
    warning(['Changing the batch size to be a multiple of decimation. New batch size is ',num2str(rawBatchSize)])
end
outBatchSize=floor(rawBatchSize/decimation);


%ROI selection if required
if selectROI
    disp('ROI selection was opened in a separate figure');
    subdata=getSLSCI(dataLSCI,5,'fastcpu');
    figure
    sgtitle(["Select ROI to process","(double-click to confirm)"])
    subplot(1,2,1)
    img=dataLSCI;
    imagesc(img)
    caxis([prctile(img(:),1),prctile(img(:),99)])
    title('Raw data')
    axis image
    subplot(1,2,2)
    img=subdata;
    imagesc(img)
    caxis([prctile(img(img(:)>0 & img(:)<1),1),prctile(img(img(:)>0 & img(:)<1),99)])
    title('Spatial contrast')
    axis image
    roi=images.roi.Rectangle(gca,'Position',[size(img,1)/4,size(img,2)/4,size(img,1)/2,size(img,2)/2]);
    wait(roi);
    ROI=round([roi.Position(2),roi.Position(2)+roi.Position(4)-1;roi.Position(1),roi.Position(1)+roi.Position(3)-1]);
end

%contrast calculation
tic
switch contrastType
    case 'spatial'
        procFramesN=floor(rawFramesN/decimation);
        dataLSCI=zeros(ROI(1,2)-ROI(1,1)+1,ROI(2,2)-ROI(2,1)+1,procFramesN,'single');
        %Uncertainty mask
        trustMatrix=ones(size(dataLSCI,1),size(dataLSCI,2),3).*rawFramesN;
        batchesN=ceil(rawFramesN/rawBatchSize);
        for i=1:1:batchesN
            framesToRead=min(rawFramesN-(i-1)*rawBatchSize,rawBatchSize+rawAveraging-1);
            subdata=readMRAW(mrawFileName,(i-1)*rawBatchSize,framesToRead,ROI);
            subdata=movmean(single(subdata),[0,rawAveraging-1],3,'omitnan','Endpoints','discard');
            trustMatrix(:,:,1)=trustMatrix(:,:,1)-sum(subdata==0,3);
            trustMatrix(:,:,2)=trustMatrix(:,:,2)-sum(subdata==255,3);
            subdata=getSLSCI(subdata,contrastKernel,procType);
            trustMatrix(:,:,3)=trustMatrix(:,:,3)-sum(subdata==0,3);
            subdata=movmean(subdata,[0,decimation-1],3);
            dataLSCI(:,:,(i-1)*outBatchSize+1:1:(i-1)*outBatchSize+floor(framesToRead/decimation))= subdata(:,:,1:decimation:floor(end./decimation)*decimation);
            disp(['Processed batch ',num2str(i),' out of ',num2str(batchesN),'. Time elapsed ',num2str(toc)]);
        end
    case 'temporal'
        procFramesN=floor((rawFramesN-contrastKernel+1)/decimation);
        dataLSCI=zeros(ROI(1,2)-ROI(1,1)+1,ROI(2,2)-ROI(2,1)+1,procFramesN,'single');
        %Uncertainty mask
        trustMatrix=ones(size(dataLSCI,1),size(dataLSCI,2),3).*rawFramesN;
        batchesN=ceil((rawFramesN-contrastKernel+1-rawAveraging+1)/rawBatchSize);
        for i=1:1:batchesN
            framesToRead=min(rawFramesN-(i-1)*rawBatchSize,rawBatchSize+contrastKernel-1+rawAveraging-1);
            subdata=readMRAW(mrawFileName,(i-1)*rawBatchSize,framesToRead,ROI);
            subdata=movmean(single(subdata),[0,rawAveraging-1],3,'omitnan','Endpoints','discard');
            trustMatrix(:,:,1)=trustMatrix(:,:,1)-sum(subdata==0,3);
            trustMatrix(:,:,2)=trustMatrix(:,:,2)-sum(subdata==255,3);
            subdata=getTLSCI(subdata,contrastKernel,procType);
            trustMatrix(:,:,3)=trustMatrix(:,:,3)-sum(subdata==0,3);
            subdata=movmean(subdata,[0,decimation-1],3,'Endpoints','shrink');
            subdata=subdata(:,:,1:end-contrastKernel+1);
            dataLSCI(:,:,((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(size(subdata,3)/decimation)))= subdata(:,:,1:decimation:floor(end./decimation)*decimation);
            
            disp(['Processed batch ',num2str(i),' out of ',num2str(batchesN),'. Time elapsed ',num2str(toc)]);
        end
end
time=((1:1:size(dataLSCI,3))-1).*decimation./fps; %convert time to seconds

trustMatrix=trustMatrix./rawFramesN;

settings.contrastCalculation.contrastType=contrastType;
settings.contrastCalculation.contrastKernel=contrastKernel;
settings.contrastCalculation.mrawFileName=mrawFileName;
settings.contrastCalculation.rawBatchSize=rawBatchSize;
settings.contrastCalculation.procType=procType;
settings.contrastCalculation.decimation=decimation;
settings.contrastCalculation.selectROI=selectROI;
settings.contrastCalculation.ROI=ROI;


if saveContrast
    disp('Saving the results');
    save(strrep(mrawFileName,'.mraw',['_',contrastType,'.mat']),'dataLSCI','time','timeStamp','trustMatrix','settings','-v7.3');
end

disp('Processing has been compelted.');
end
%------------- END OF CODE --------------
%Comments: For the high-speed cameras and other cases with abnormal
%timestamps - pay attention to the time scaling,