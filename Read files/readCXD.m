%readCXD - Read a Hamamatsu .cxd recording via the Bio-Formats library.
%
%   Reads a Hamamatsu .cxd image sequence using Bio-Formats, with optional frame
%   skipping, a frame count and an ROI.
%
% Syntax:
%    [data,sampling,timeStamps,sizeT,startT] = readCXD(fileName)
%    [data,sampling,timeStamps,sizeT,startT] = readCXD(fileName, framesToSkip, framesToRead, ROI)
%
% Inputs:
%    fileName     - path to the .cxd file.
%    framesToSkip - (optional) frames skipped at the start (default 0).
%    framesToRead - (optional) frames to read ([] = all remaining, default).
%    ROI          - (optional) [firstRow lastRow; firstCol lastCol]
%                   (default whole frame).
%
% Outputs:
%    data       - raw image data, 3-D [y x t].
%    sampling   - frame time increment, ms.
%    timeStamps - per-frame timestamps (offset by framesToSkip).
%    sizeT      - total number of frames in the file.
%    startT     - acquisition start date/time from the metadata.
%
% Example:
%    [data,dt,ts,sizeT] = readCXD('rec.cxd',10,100,[10 100;20 200]);
%
% Dependencies: the Bio-Formats MATLAB library (bfGetReader, bfGetPlane).
% See also: readDV, readMRAW, getK
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function [data,sampling,timeStamps,sizeT,startT]=readCXD(fileName,varargin)

%check if the fileName ends with .cxd, add .cxd otherwise
C = strsplit(fileName,'.');
if ~strcmp(C{end},'cxd')
    fileName=[fileName,'.cxd'];
end

%initializing default values for optional and hidden parameters
framesToSkip=0;
framesToRead=[];
ROI=[];

%checking that inputs are in the correct format:
for iVar = 1:length(varargin)
    if iVar==1
        framesToSkip = varargin{1}; 
        if ~(round(framesToSkip)==framesToSkip)
            framesToSkip=round(framesToSkip);
            warning('Rounding number of frames to nearest integer')
        end
    end
    if iVar==2 && ~isempty(varargin{2})
        framesToRead = varargin{2};
        if ~(round(framesToRead)==framesToRead)
            framesToRead=round(framesToRead);
            
            warning('Rounding number of frames to nearest integer')
        end
    end
    if iVar==3 
        ROI = varargin{3};
        if ~all(size(ROI)==[2,2])
            error(['ROI format is incorrect, please check that the matrix dimensions are: ',char(13),...
                    '[firstRow    --> , <--   lastRow   --> ; <--',char(13),...
                    ' firstColumn --> , <--   lastColumn      ]'])
        end
    end
end


%read the meta data
reader = bfGetReader(fileName);
omeMeta = reader.getMetadataStore();
sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
sizeT=double(omeMeta.getPlaneCount(0)); %
startT=omeMeta.getImageAcquisitionDate(0);
dataType=char(omeMeta.getPixelsType(0).toString());
sampling=double(omeMeta.getPixelsTimeIncrement(0).value()).*1000; %ms

%correct default values based on meta data
if isempty(ROI), ROI = [1,sizeY;1,sizeX]; end
if isempty(framesToRead), framesToRead = sizeT-framesToSkip; end

%pre-allocate memory for arrays
timeStamps=zeros(framesToRead,1,'int64');
data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),framesToRead,dataType);


%read data and close the file
for t=1:1:framesToRead
    timeStamps(t)=(t+framesToSkip-1).*sampling;
    frame = bfGetPlane(reader, t+framesToSkip);
    data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
end

reader.close();
end
%------------- END OF CODE --------------
% Comments: this code can be used only with files that are several times
% smaller than the available RAM memory. Use the batch version of readRLS
% (batchReadRLS) to process the data using less RAM memory.