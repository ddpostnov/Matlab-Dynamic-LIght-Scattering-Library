%fixMetaRLS - Rewrite a .rls recording with corrected header metadata.
%
%   Reads a .rls file whose header metadata is wrong or missing and writes a new
%   .rls with the supplied sizeY, sizeX, sizeT and sampling.  The frame data is
%   copied verbatim; only the header is rebuilt.
%
% Syntax:
%    fixMetaRLS(fileName, sizeY, sizeX, iSizeT, sampling, dataType)
%
% Inputs:
%    fileName - path to the source .rls file ('.rls' appended if missing).
%    sizeY    - frame height, pixels.
%    sizeX    - frame width, pixels.
%    iSizeT   - number of frames to write ([] = count them while copying).
%    sampling - sampling / frame rate stored in the header.
%    dataType - pixel type of the source data, 'uint8' or 'uint16'.
%
% Outputs:
%    (none)   - writes <fileName>_fixed.rls next to the source.
%
% Example:
%    fixMetaRLS('rec.rls',1000,1000,[],100,'uint8');
%
% See also: readRLS, cropRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function fixMetaRLS(fileName,sizeY,sizeX,iSizeT,sampling,dataType)

%check if the fileName ends with .rls, add .rls otherwise
C = strsplit(fileName,'.');
if ~strcmp(C{end},'rls')
    fileName=[fileName,'.rls'];
end

fileNameSave=strrep(fileName,'.rls','_fixed.rls');
fileReadId = fopen(fileName, 'r');
fileWriteId = fopen(fileNameSave, 'w');
%write empty meta data
fseek(fileWriteId,0*1024,-1 );
for k=1:(30*1024)
    fwrite(fileWriteId,0 , 'uint64');
end
fseek(fileWriteId,30*1024,-1 );
fseek(fileReadId,30*1024,-1 );
%read data and close the file
sizeT=1;
timeStamp=fread(fileReadId,1,'*uint64');
frame=fread(fileReadId,[sizeY,sizeX],['*',dataType]);
while ~isempty(timeStamp)
fwrite(fileWriteId,timeStamp, 'uint64'); %dummy timestamps
fwrite(fileWriteId,frame,dataType);
timeStamp=fread(fileReadId,1,'*uint64');
frame=fread(fileReadId,[sizeY,sizeX],['*',dataType]);
sizeT=sizeT+1;
end
fclose(fileReadId);

fseek(fileWriteId,0*1024,-1 );
fwrite(fileWriteId,uint64(sizeX), 'uint64');
fwrite(fileWriteId,uint64(sizeY), 'uint64');
if isempty(iSizeT)
fwrite(fileWriteId,uint64(sizeT-1), 'uint64');
else 
    fwrite(fileWriteId,uint64(iSizeT), 'uint64');
end
fwrite(fileWriteId,uint64(sampling), 'uint64');
dataSize=1; if strcmpi(dataType,'uint16'), dataSize=2; end
fwrite(fileWriteId,'Ver.','ubit8');             % version tag, so the data type is preserved
fwrite(fileWriteId,uint64(2),'uint64');         % nVer (>1 so dataSize is stored)
fwrite(fileWriteId,uint64(dataSize),'uint64');  % dataSize (1=uint8, 2=uint16)
fclose(fileWriteId);
end
%------------- END OF CODE --------------
