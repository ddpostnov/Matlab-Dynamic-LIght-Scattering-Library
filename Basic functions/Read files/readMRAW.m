%readMRAW - reads (.mraw) data files, requires cihx file to be
% present in the same location. Works with 16 or 8 bit data, not the real
% bit
%
% Syntax:  output1 = function_name(requiredInput1,optionalInput1,optionalInput2,optionalInput3)
%
% Required inputs:
%   fileName      - path to the .mraw file
% Optional Inputs:
%   framesToSkip  - number of frames to skip at the begining of the recording.
%                   Default is 0. Use 0 when you don't want to skip frames
%                   but need to enter other optional arguments
%   framesToRead  - number of frames to read. Use [] when you want to read
%                   all frames and use ROI. By default reads all frames
%                   except the skipped ones.
%   ROI           - Region of Interest:
%                   [firstRow,      lastRow ;
%                    firstColumn,   lastColumn]
%                   Default is the whole frame.
% Outputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    sampling   - sampling time in ms
%    timeStamps - time stamps in microseconds
%    sizeT      - total number of frames in the file
%
% Example 1:
%    [data,sampling,timeStamps,sizeT]=readMRAW('D:\stroke_20201001.mraw',10,100,[10,100;20,200]);
%
% Example 2:
%    [data,sampling,timeStamps,sizeT]=readMRAW('D:\stroke_20201001.mraw');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getTLSCI.m, getSLSCI.m

% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 25-July-2022

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
        strrep(strrep(tmp,'<width>',''),'</width>','');
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

if ROI==[1,sizeY;1,sizeX]
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