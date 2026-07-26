%getContrastFromMRAW - Batch spatial/temporal contrast from a Photron .mraw file.
%
%   Reads a raw high-speed recording (.mraw plus its companion .cihx) in batches
%   and computes spatial or temporal laser speckle contrast with getK, so
%   recordings larger than RAM can be processed.  Supports optional pre-averaging
%   of raw frames, temporal decimation, ROI selection and a per-pixel trust map.
%
% Syntax:
%    [dataLSCI,time,timeStamp,trustMatrix,settings] = getContrastFromMRAW( ...
%        mrawFileName, contrastType, contrastKernel, rawBatchSize, procType, ...
%        decimation, saveContrast, selectROI, rawAveraging)
%
% Inputs:
%    mrawFileName   - path to the .mraw file (.cihx must sit next to it).
%    contrastType   - 'spatial' or 'temporal'.
%    contrastKernel - contrast kernel size (e.g. 5 spatial, 25 temporal).
%    rawBatchSize   - raw frames loaded per batch (must be >= decimation).
%    procType       - 'cpu' or 'gpu'.
%    decimation     - frames averaged and kept during decimation (1 = none).
%    saveContrast   - logical; save the result next to the .mraw as a .mat.
%    selectROI      - logical; interactively pick an ROI to process.
%    rawAveraging   - raw frames averaged before contrast (1 = none).
%
% Outputs:
%    dataLSCI    - contrast data, 3-D [y x t] single.
%    time        - relative time vector, seconds.
%    timeStamp   - absolute start timestamp (0 for .mraw).
%    trustMatrix - [y x 3] fraction of valid frames per pixel (non-zero,
%                  non-saturated, non-zero-contrast).
%    settings    - struct recording the processing parameters used.
%
% Example:
%    K = getContrastFromMRAW('rec.mraw','spatial',5,2000,'gpu',10,true,false,1);
%
% Dependencies: readMRAW, getK; Parallel Computing Toolbox (procType 'gpu').
% See also: getContrastFromRLS, getK, readMRAW
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function [dataLSCI,time,timeStamp,trustMatrix,settings]=getContrastFromMRAW(mrawFileName,contrastType,contrastKernel,rawBatchSize,procType,decimation,saveContrast,selectROI,rawAveraging)
timeStamp=0;
procType=erase(lower(procType),'fast');       % accept legacy 'fastcpu'/'fastgpu' -> 'cpu'/'gpu'
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
    subdata=getK(dataLSCI,'spatial','kernelSize',5,'procType','cpu');
    figure
    sgtitle(["Select ROI to process","(double-click to confirm)"])
    subplot(1,2,1)
    img=dataLSCI;
    imagesc(img)
    clim([prctile(img(:),1),prctile(img(:),99)])
    title('Raw data')
    axis image
    subplot(1,2,2)
    img=subdata;
    imagesc(img)
    clim([prctile(img(img(:)>0 & img(:)<1),1),prctile(img(img(:)>0 & img(:)<1),99)])
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
            subdata=getK(subdata,'spatial','kernelSize',contrastKernel,'procType',procType);
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
            subdata=getK(subdata,'temporal','kernelSize',contrastKernel,'procType',procType);
            trustMatrix(:,:,3)=trustMatrix(:,:,3)-sum(subdata==0,3);
            subdata=movmean(subdata,[0,decimation-1],3,'Endpoints','shrink');
            subdata=subdata(:,:,1:end-contrastKernel+1);
            dataLSCI(:,:,((i-1)*outBatchSize+1):1:((i-1)*outBatchSize+floor(size(subdata,3)/decimation)))= subdata(:,:,1:decimation:floor(end./decimation)*decimation);
            
            disp(['Processed batch ',num2str(i),' out of ',num2str(batchesN),'. Time elapsed ',num2str(toc)]);
        end
end
time=((1:1:size(dataLSCI,3))-1).*decimation./fps; %convert time to seconds

trustMatrix=trustMatrix./rawFramesN;

settings.getContrast.contrastType=contrastType;
settings.getContrast.contrastKernel=contrastKernel;
settings.getContrast.mrawFileName=mrawFileName;
settings.getContrast.rawBatchSize=rawBatchSize;
settings.getContrast.procType=procType;
settings.getContrast.decimation=decimation;
settings.getContrast.selectROI=selectROI;
settings.getContrast.ROI=ROI;


if saveContrast
    disp('Saving the results');
    save(strrep(mrawFileName,'.mraw',['_',contrastType,'.mat']),'dataLSCI','time','timeStamp','trustMatrix','settings','-v7.3');
end

disp('Processing has been completed.');
end
%------------- END OF CODE --------------
%Comments: For the high-speed cameras and other cases with abnormal
%timestamps - pay attention to the time scaling,