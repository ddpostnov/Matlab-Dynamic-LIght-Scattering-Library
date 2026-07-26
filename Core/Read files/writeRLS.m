%writeRLS - Write a 3-D stack to a raw laser speckle (.rls) file.
%
%   Writes a [y x t] stack in the .rls layout read by readRLS/getPointerRLS:
%   a 30 KB header (sizeX, sizeY, sizeT, sampling as uint64, the 'Ver.' tag,
%   nVer=2 and dataSize) followed, from byte 30*1024, by one uint64 time
%   stamp and one column-major frame per time point.
%
% Syntax:
%    writeRLS(fileName, data)
%    writeRLS(fileName, data, sampling, timeStamps)
%
% Inputs:
%    fileName   - output path ('.rls' appended if missing).
%    data       - stack, 3-D [y x t], class uint8 (8-bit) or uint16 (16-bit).
%    sampling   - (optional) acquisition sampling stored in the header.
%                 Default 0.
%    timeStamps - (optional) per-frame time stamps, length t. Default 0:t-1.
%
% Outputs:
%    (none)     - writes fileName.
%
% Example:
%    writeRLS('sim.rls', uint8(stack), 100);
%
% See also: readRLS, cropRLS, getPointerRLS, fixMetaRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 20-July-2026

%------------- BEGIN CODE --------------
function writeRLS(fileName, data, sampling, timeStamps)

C = strsplit(fileName,'.');
if ~strcmp(C{end},'rls')
    fileName = [fileName,'.rls'];
end

[sizeY, sizeX, sizeT] = size(data);
switch class(data)
    case 'uint8',  dataSize = 1;
    case 'uint16', dataSize = 2;
    otherwise
        error('writeRLS:InvalidType','data must be uint8 or uint16 (got %s).', class(data));
end

if nargin < 3 || isempty(sampling),   sampling   = 0;            end
if nargin < 4 || isempty(timeStamps), timeStamps = 0:sizeT-1;    end
if numel(timeStamps) ~= sizeT
    error('writeRLS:BadTimeStamps','timeStamps must have %d elements.', sizeT);
end

fId = fopen(fileName,'w');
if fId == -1, error('writeRLS:CannotOpen','Cannot open file for writing: %s', fileName); end

% Header (matches cropRLS): fields then zero-padding up to byte 30*1024.
fseek(fId,0,-1);
fwrite(fId, uint64(sizeX),    'uint64');
fwrite(fId, uint64(sizeY),    'uint64');
fwrite(fId, uint64(sizeT),    'uint64');
fwrite(fId, uint64(sampling), 'uint64');
fwrite(fId, 'Ver.',           'ubit8');          % version tag so dataSize is stored
fwrite(fId, uint64(2),        'uint64');         % nVer (>1 so dataSize is read back)
fwrite(fId, uint64(dataSize), 'uint64');         % 1 = uint8, 2 = uint16
headerBytes = 8*4 + 4 + 8 + 8;                   % = 52
fwrite(fId, zeros(30*1024-headerBytes,1,'uint8'), 'uint8');

% Data: one time stamp + one column-major frame per time point.
dataType = class(data);
for t = 1:sizeT
    fwrite(fId, uint64(timeStamps(t)), 'uint64');
    fwrite(fId, data(:,:,t), dataType);
end
fclose(fId);
end
%------------- END OF CODE --------------
