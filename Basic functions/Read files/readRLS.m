%readRLS - reads raw laser speckle (.rls) data files.
%
% Syntax:  [data,timeStamps,s] = readRLS(Input1, Name, Value)
%
% Required inputs:
%   Input1        - Path to the .rls file (String).
%                   Required for New File mode.
%                   Omit entirely if 'Stream' structure is provided.
%
% Optional Name-Value Pair Inputs:
%   'Stream'      - Structure containing metadata and fileID from a previous
%                   call. Use this to continue reading from an open stream
%                   without re-parsing the header.
%   'KeepOpen'    - Boolean (true/false). If true, the file remains open after
%                   reading. Required for batch/stream processing. Default is false.
%   'FramesToSkip'- Number of frames to skip at the beginning of the recording.
%                   Default is 0. Ignored (triggers warning) if 'Stream' is provided.
%   'FramesToRead'- Number of frames to read.
%                   - In "New File" mode: Default is all remaining frames.
%                   - In "Stream" mode: Must be specified explicitly.
%
% Outputs:
%    data       - Raw laser speckle data as 3d [y,x,t] matrix.
%    timeStamps - Vector of time stamps (int64).
%    s          - Stream structure containing metadata and the active fileID.
%
% Example 1 (Standard Read):
%    [data, timeStamps, s] = readRLS('data.rls', 'FramesToSkip', 10, 'FramesToRead', 100);
%
% Example 2 (Stream/Batch Processing):
%    % 1. Open and read first batch (KeepOpen = true)
%    [d1, t1, s] = readRLS('data.rls', 'FramesToRead', 100, 'KeepOpen', true);
%
%    % 2. Read next batch using 's' (Notice: fileName is omitted entirely)
%    [d2, t2, s] = readRLS('Stream', s, 'FramesToRead', 100);
%
%    % 3. Close manually when done
%    fclose(s.fId);
%
% See also: getTLSCI.m, getSLSCI.m
% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 25-January-2026

%------------- BEGIN CODE --------------
function [data,timeStamps,s]=readRLS(varargin)
p = inputParser;
p.KeepUnmatched = false;

if nargin == 3 && isstruct(varargin{1}) && strcmpi(varargin{2}, 'FramesToRead')
    %%STREAM MODE
    s= varargin{1};
    s.framesToRead = varargin{3};
    % Validation
    if ~isfield(s, 'fId') || isempty(s.fId)
        error('Invalid stream structure: fId is missing.');
    end
    s.keepOpen = true;
elseif nargin>=1 && (ischar(varargin{1}) || isstring(varargin{1}))
    %%NEW FILE MODE
    parserInput = varargin;
    addRequired(p, 'fileName', @(x) ischar(x) || isstring(x));

    % Optional Arguments
    addParameter(p, 'KeepOpen', false, @(x) islogical(x) || x==1 || x==0);
    addParameter(p, 'FramesToSkip', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'FramesToRead', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    parse(p, parserInput{:});
    fileName = p.Results.fileName;
    [~, ~, ext] = fileparts(fileName);
    if ~strcmpi(ext, '.rls')
        error('Input file name should have .rls extension');
    end
    s = struct();
    s.framesToSkip = p.Results.FramesToSkip;
    s.framesToRead = p.Results.FramesToRead;
    s.keepOpen = p.Results.KeepOpen;


    s.fileName=fileName;
    s.fId = fopen(s.fileName, 'r');
    if s.fId == -1, error('Cannot open file: %s', fileName); end
    fseek(s.fId,0*1024,-1 );
    s.sizeX=double(fread(s.fId,1,'*uint64'));
    s.sizeY=double(fread(s.fId,1,'*uint64'));
    s.sizeT=double(fread(s.fId,1,'*uint64'));
    s.sampling=double(fread(s.fId,1,'*uint64'));
    s.version=fread(s.fId,4,'*ubit8')';
    s.dataSize=1; %default uint8
    if strcmp(char(s.version), 'Ver.') %strcmp(s.version,'Ver.')
        nVer = fread(s.fId,1,'*uint64');
        if nVer>1
            s.dataSize=double(fread(s.fId,1,'*uint64'));
        end
    end

    %correct default values based on meta data
    if isempty(s.framesToRead), s.framesToRead = s.sizeT-s.framesToSkip; end
    switch s.dataSize
        case 1
            s.dataType='uint8';
        case 2
            s.dataType='uint16';
        otherwise
            error('Unindentified data type')
    end

    %move to the first timeStamp/frame location
    firstByte=30*1024+s.sizeX*s.sizeY*s.framesToSkip*s.dataSize+8*s.framesToSkip;
    fseek(s.fId,firstByte,-1 );
else
    % Fallback for invalid inputs
    error('readRLS:InvalidInput', ...
        'Inputs must be either a filename string or a Stream struct with FramesToRead specified.');
end

timeStamps=zeros(s.framesToRead,1,'double');
data=zeros(s.sizeY,s.sizeX,s.framesToRead,s.dataType);

sY=s.sizeY;
sX=s.sizeX;
dT=s.dataType;
fId=s.fId;
for t=1:1:s.framesToRead
    if feof(fId), break; end
    timeStamps(t)=double(fread(fId,1,'*uint64'));
    data(:,:,t)=fread(fId,[sY,sX],['*',dT]);
end


s.fId=fId;
if ~s.keepOpen
    fclose(s.fId);
end
end
%------------- END OF CODE --------------
