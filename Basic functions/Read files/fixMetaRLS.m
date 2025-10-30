%fixMetaRLS - reads raw laser speckle (.rls) data file and re-saves it with
%a new meta
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 07-October-2022

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
fclose(fileWriteId);
end
%------------- END OF CODE --------------
