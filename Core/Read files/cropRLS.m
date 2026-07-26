%cropRLS - Crop a raw laser speckle (.rls) recording in frames and space.
%
%   Reads a .rls file and writes a new .rls keeping only the frames in fLim and
%   the pixels in ROI, without loading the whole recording into memory.
%
% Syntax:
%    cropRLS(fileName, fLim, ROI, suffix)
%
% Inputs:
%    fileName - path to the source .rls file ('.rls' appended if missing).
%    fLim     - [firstFrame lastFrame] to keep ([] = all frames).
%    ROI      - [firstRow lastRow; firstCol lastCol] to keep ([] = whole frame).
%    suffix   - appended to the output name ([] -> '_crop').
%
% Outputs:
%    (none)   - writes <fileName><suffix>.rls next to the source.
%
% Example:
%    cropRLS('rec.rls',[100 2000],[1 512;1 512],'_roi');
%
% See also: readRLS, fixMetaRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function cropRLS(fileName,fLim,ROI,suffix)

%check if the fileName ends with .rls, add .rls otherwise
C = strsplit(fileName,'.');
if ~strcmp(C{end},'rls')
    fileName=[fileName,'.rls'];
end

if isempty(suffix)
fileNameSave=strrep(fileName,'.rls','_crop.rls');
else
fileNameSave=strrep(fileName,'.rls',[suffix,'.rls']);
end
fileReadId = fopen(fileName, 'r');
fileWriteId = fopen(fileNameSave, 'w');



%read the meta data
fseek(fileReadId,0*1024,-1 );
dataSize=1;
sizeX=fread(fileReadId,1,'*uint64');
sizeY=fread(fileReadId,1,'*uint64');
sizeT=single(fread(fileReadId,1,'*uint64'));
if isempty(fLim)
    fLim=[1,sizeT];
end

sampling=single(fread(fileReadId,1,'*uint64'));
version=fread(fileReadId,4,'*ubit8')';
if strcmp(char(version),'Ver.')
    nVer = fread(fileReadId,1,'*uint64');
    if nVer>1
        dataSize=fread(fileReadId,1,'*uint64');
    end
end

%correct default values based on meta data
if isempty(ROI), ROI = [1,sizeY;1,sizeX]; end
switch dataSize
    case 1
        dataType='uint8';
    case 2
        dataType='uint16';
    otherwise
        error('Unindentified data type')
end

%write empty meta data
fseek(fileWriteId,0*1024,-1 );
for k=1:(30*1024)
    fwrite(fileWriteId,0 , 'uint64');
end
fseek(fileWriteId,0*1024,-1 );
fwrite(fileWriteId,uint64(ROI(2,2)-ROI(2,1)+1), 'uint64');
fwrite(fileWriteId,uint64(ROI(1,2)-ROI(1,1)+1), 'uint64');
fwrite(fileWriteId,uint64(fLim(2)-fLim(1)+1), 'uint64');
fwrite(fileWriteId,uint64(sampling), 'uint64');
fwrite(fileWriteId,'Ver.','ubit8');             % version tag, so the data type is preserved
fwrite(fileWriteId,uint64(2),'uint64');         % nVer (>1 so dataSize is stored)
fwrite(fileWriteId,uint64(dataSize),'uint64');  % dataSize (1=uint8, 2=uint16)

fseek(fileWriteId,30*1024,-1 );
fseek(fileReadId,30*1024,-1 );
%read data and close the file
for i=1:1:fLim(2)
timeStamp=fread(fileReadId,1,'*uint64');
frame=fread(fileReadId,[sizeY,sizeX],['*',dataType]);
frame=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
if i>=fLim(1)
fwrite(fileWriteId,timeStamp, 'uint64'); 
fwrite(fileWriteId,frame,dataType);
end
end
fclose(fileReadId);
fclose(fileWriteId);
end
%------------- END OF CODE --------------
