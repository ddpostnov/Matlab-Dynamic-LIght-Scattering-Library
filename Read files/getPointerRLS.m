%getPointerRLS - Open a .rls file and return a file pointer at the first frame.
%
%   Opens a raw laser speckle (.rls) file, reads its metadata into a struct and
%   leaves the file open with the read position at the first frame (after any
%   skipped frames), so frames can be streamed manually with fread.  The caller
%   must close the returned pointer with fclose.
%
% Syntax:
%    [fPointer,s] = getPointerRLS(fName)
%    [fPointer,s] = getPointerRLS(fName, framesToSkip)
%
% Inputs:
%    fName        - path to the .rls file ('.rls' is appended if missing).
%    framesToSkip - (optional) frames skipped before the read position (default 0).
%
% Outputs:
%    fPointer - open file identifier, positioned at the first frame to read.
%    s        - metadata struct (sizeX, sizeY, sizeT, sampling, version,
%               dataSize, dataType, fName, framesToSkip).
%
% Example:
%    [fp,s] = getPointerRLS('rec.rls');
%    ts     = fread(fp,1,'*uint64');                         % frame timestamp
%    frame  = fread(fp,[s.sizeY s.sizeX],['*',s.dataType]);  % frame data
%    fclose(fp);
%
% See also: readRLS, getContrastFromRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function [fPointer,s]=getPointerRLS(fName,varargin)

%check if the s.fName ends with .rls, add .rls otherwise
C = strsplit(fName,'.');
if ~strcmp(C{end},'rls')
    fName=[fName,'.rls'];
end
s.fName=fName;

%initializing default values for optional and hidden parameters
s.framesToSkip=0;
s.dataSize=1; %initiate data size as default for 'uint8'

%checking that inputs are in the correct format:
for iVar = 1:length(varargin)
    if iVar==1
        s.framesToSkip = varargin{1}; 
        if ~(round(s.framesToSkip)==s.framesToSkip)
            s.framesToSkip=round(s.framesToSkip);
            warning('Rounding number of frames to nearest integer')
        end
    end
end


%read the meta data
fPointer = fopen(s.fName, 'r');
fseek(fPointer,0*1024,-1 );
s.sizeX=fread(fPointer,1,'*uint64');
s.sizeY=fread(fPointer,1,'*uint64');
s.sizeT=single(fread(fPointer,1,'*uint64'));
s.sampling=single(fread(fPointer,1,'*uint64'));
s.version=fread(fPointer,4,'*ubit8')';
if strcmp(char(s.version),'Ver.')
    nVer = fread(fPointer,1,'*uint64');
    if nVer>1
        s.dataSize=fread(fPointer,1,'*uint64');
    end
end

%correct default values based on meta data
switch s.dataSize
    case 1
        s.dataType='uint8';
    case 2
        s.dataType='uint16';
    otherwise
        error('Unindentified data type')
end

%move to the first timeStamp/frame location
firstByte=30*1024+s.sizeX*s.sizeY*uint64(s.framesToSkip)*s.dataSize+8*uint64(s.framesToSkip);
fseek(fPointer,firstByte,-1 );
end
%------------- END OF CODE --------------
% Comments: this code can be used only with files that are several times
% smaller than the available RAM memory. Use the batch s.version of readRLS
% (batchReadRLS) to process the data using less RAM memory.