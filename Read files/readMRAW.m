%readMRAW - Read a Photron .mraw high-speed recording.
%
%   Reads a raw Photron recording (.mraw) using the companion .cihx metadata
%   file, which must sit next to it.  Handles 8- and 16-bit stored data and
%   supports frame skipping, a frame count and an ROI.
%
% Syntax:
%    [data,fps,time,sizeT] = readMRAW(fileName)
%    [data,fps,time,sizeT] = readMRAW(fileName, framesToSkip, framesToRead, ROI)
%
% Inputs:
%    fileName     - path to the .mraw file (.cihx must sit next to it).
%    framesToSkip - (optional) frames skipped at the start (default 0).
%    framesToRead - (optional) frames to read ([] = all remaining, default).
%    ROI          - (optional) [firstRow lastRow; firstCol lastCol]
%                   (default whole frame).
%
% Outputs:
%    data  - raw data, 3-D [y x t].
%    fps   - recording frame rate, frames/second.
%    time  - relative time vector, seconds.
%    sizeT - total number of frames in the file.
%
% Example:
%    [data,fps,time,sizeT] = readMRAW('rec.mraw',10,100,[10 100;20 200]);
%
% Dependencies: a companion .cihx metadata file.
% See also: getContrastFromMRAW, getK, readRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function [data,fps,time,sizeT]=readMRAW(fileName,varargin)

%check if the fileName ends with .mraw, add .mraw otherwise
C = strsplit(fileName,'.');
if ~strcmp(C{end},'mraw')
    fileName=[fileName,'.mraw'];
end

%initializing default values for optional and hidden parameters
framesToSkip=0;
framesToRead=[];
ROI=[];

%checking that inputs are in the correct format:
for iVar = 1:length(varargin)
    if iVar==1 && ~isempty(varargin{1})
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
    if iVar==3 && ~isempty(varargin{3})
        ROI = varargin{3};
        if ~all(size(ROI)==[2,2])
            error(['ROI format is incorrect, please check that the matrix dimensions are: ',char(13),...
                '[firstRow    --> , <--   lastRow   --> ; <--',char(13),...
                ' firstColumn --> , <--   lastColumn      ]'])
        end
    end
end


%read the meta data


cihx=readlines(strrep(fileName,'mraw','cihx'));
for i=1:1:length(cihx)
    tmp=cihx{i};
    if contains(tmp,'<resolution>')
        tmp=cihx{i+1};
        sizeY=str2double(strrep(strrep(tmp,'<width>',''),'</width>',''));
        tmp=cihx{i+2};
        sizeX=str2double(strrep(strrep(tmp,'<height>',''),'</height>',''));
    end
    if contains(tmp,'<totalFrame>')
        sizeT=str2double(strrep(strrep(tmp,'<totalFrame>',''),'</totalFrame>',''));
    end
    if contains(tmp,'<recordRate>')
        fps=str2double(strrep(strrep(tmp,'<recordRate>',''),'</recordRate>',''));
    end
    if contains(tmp,'<bit>')
        dataSize=str2double(strrep(strrep(tmp,'<bit>',''),'</bit>',''))/8;
    end
end

%correct default values based on meta data
if isempty(ROI), ROI = [1,sizeY;1,sizeX]; end
if isempty(framesToRead), framesToRead = sizeT; end
switch dataSize
    case 1
        dataType='uint8';
    case 2
        dataType='uint16';
    otherwise
        error('Unindentified data type')
end

%pre-allocate memory for arrays
time=((1:1:sizeT)-1)/fps;

%move to the first timeStamp/frame location
fileReadId=fopen(fileName,'r');
fseek(fileReadId,sizeX*sizeY*uint64(framesToSkip)*dataSize,'bof');

if isequal(ROI,[1,sizeY;1,sizeX])
    data=fread(fileReadId,sizeX*sizeY*framesToRead,['*',dataType]);
    data=reshape(data,sizeY,sizeX,framesToRead);
else
    data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),framesToRead,dataType);
    %read data and close the file
    for t=1:1:framesToRead
        frame=fread(fileReadId,[sizeY,sizeX],['*',dataType]);
        data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
    end
end
fclose(fileReadId);

end
%------------- END OF CODE --------------
% Comments: this code can be used only with files that are several times
% smaller than the available RAM memory. Use the batch version of readMRAW
% (batchReadMRAW) to process the data using less RAM memory.