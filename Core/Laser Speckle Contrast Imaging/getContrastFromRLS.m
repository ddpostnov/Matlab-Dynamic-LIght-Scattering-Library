%getContrastFromRLS - Batch spatial/temporal contrast from a raw .rls recording.
%
%   Reads a raw laser speckle (.rls) file in batches and computes spatial or
%   temporal contrast with getK, so recordings larger than RAM can be processed.
%   Temporal contrast uses overlap-aware batching ('leaking') so the result is
%   identical to processing the whole recording at once.
%
% Syntax:
%    [data,time,timeStamp,intensityMetrics] = getContrastFromRLS(rlsFileName,contrastType)
%    [...] = getContrastFromRLS(rlsFileName,contrastType, Name,Value,...)
%
% Inputs:
%    rlsFileName  - path to the .rls file.
%    contrastType - 'spatial' or 'temporal'.
%
% Name-Value pairs (optional):
%    'kernelSize'  - contrast kernel size, ODD (default 5 spatial, 25 temporal).
%    'procType'    - 'cpu' or 'gpu' (default 'gpu').
%    'decimFactor' - frames averaged and kept during decimation (default 1).
%                    A recording too short to fill even one decimated frame is
%                    refused here (getContrastFromRLS:TooShort) rather than left
%                    to fail later on an empty time vector.
%    'decimMethod' - 'leaking' (overlapping, default) or 'sharp' (temporal only).
%    'memoryCoef'  - fraction of free RAM used for batching (default 0.8).
%
%   IT SAYS NOTHING WHILE IT WORKS.  Reporting belongs to the wrapper: it names
%   the recording before this is called and closes it with an elapsed time after,
%   and a per-batch line underneath that is noise in a sixty-file run.  So there
%   is no progress sink and no command-window print here - a core neither reports
%   nor takes a reporting argument.
%
% Outputs:
%    data             - contrast data, 3-D [y x t] single.
%    time             - relative time vector, seconds.
%    timeStamp        - absolute timestamp of the first frame.
%    intensityMetrics - [y x 3] per-pixel quality map: fraction of non-zero
%                       frames, fraction of non-saturated frames, and mean
%                       intensity.
%
% Example:
%    [K,t] = getContrastFromRLS('rec.rls','spatial','kernelSize',5,'procType','gpu');
%
% Dependencies: readRLS, getK; Parallel Computing Toolbox (procType 'gpu').
% See also: getContrastFromMRAW, getK, readRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

%------------- BEGIN CODE --------------
function [data,time,timeStamp,intensityMetrics]=getContrastFromRLS(rlsFileName,contrastType, varargin)
p = inputParser;
p.KeepUnmatched = false;
addRequired(p, 'rlsFileName');
addRequired(p, 'contrastType', @(x) any(validatestring(x, {'spatial', 'temporal'})));
addParameter(p, 'kernelSize', [], @isnumeric);
addParameter(p, 'procType', 'gpu', @(x) any(validatestring(x, {'cpu', 'gpu'})));
addParameter(p, 'decimFactor', 1, @isnumeric);
addParameter(p, 'decimMethod', 'leaking', @(x) any(validatestring(x, {'leaking', 'sharp'})));
addParameter(p, 'memoryCoef', 0.8, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
parse(p, rlsFileName, contrastType, varargin{:});

kernelSize  = p.Results.kernelSize;
procType    = p.Results.procType;
decimFactor = p.Results.decimFactor;
decimMethod = p.Results.decimMethod;
memoryCoef  = p.Results.memoryCoef;

if isempty(kernelSize)
    if strcmpi(contrastType, 'temporal')
        kernelSize = 25;
    elseif strcmpi(contrastType, 'spatial')
        kernelSize = 5;
    end
end

%getK owns this rule and states it in full; it is asked here as well so a bad kernel
%stops before the recording is opened rather than one batch into the loop.
if ~isscalar(kernelSize) || ~isfinite(kernelSize) || kernelSize < 1 || mod(kernelSize,2) ~= 1
    error('getContrastFromRLS:EvenKernel', ...
        'kernelSize must be an odd whole number of pixels or frames, and %s is not: a contrast window has to be centred on the pixel it describes.', ...
        mat2str(kernelSize));
end

if strcmpi(decimMethod, 'sharp')
    if strcmpi(contrastType, 'spatial')
        error('getContrastFromRLS:InvalidMethod', ...
            "The 'sharp' decimation method is only valid for 'temporal' contrast.");
    end
    if strcmpi(contrastType, 'temporal') && decimFactor > 1
        if mod(decimFactor, kernelSize) ~= 0
            error('getContrastFromRLS:InvalidDecimation', ...
                "For 'sharp' temporal contrast, decimFactor (%d) must be an integer multiple of kernelSize (%d).", ...
                decimFactor, kernelSize);
        end
    end
end


%'FramesToRead',0 opens the stream sitting on the first frame without consuming it, so
%every batch below is read the same way and the first one is not a special case.
[~,~,stream]=readRLS(rlsFileName,'KeepOpen',true,'FramesToRead',0);
%The handle is this function's for the rest of its life, including the lives it does not
%finish: a batch that throws left the fclose below unreached and the recording locked
%until MATLAB was restarted.  readRLS carries the same fId through every read (KeepOpen),
%so the scalar captured here is the one to close - and capturing the scalar rather than
%the stream keeps the whole struct out of the closure's workspace.
fIdToClose = stream.fId;
closeStream = onCleanup(@() fclose(fIdToClose));

sz = [stream.sizeY,stream.sizeX,stream.sizeT];
outT = floor(sz(3) / decimFactor);
outSz = [sz(1), sz(2), outT];

%A recording shorter than one decimated frame produces nothing, and everything below
%assumes at least one: the time base is anchored on time(1) and the batch loop writes
%into a zero-length cube.  Said here, in the recording's own terms, rather than as
%MATLAB's index message from the last line of the function.
if outT < 1
    error('getContrastFromRLS:TooShort', ...
        ['The recording is %d frames long and %d frames are averaged into each output ' ...
        'frame, so there is nothing to write. Use a longer recording or average fewer ' ...
        'frames.'],sz(3),decimFactor);
end

[~, memAvailRAM] = memory;
memAvailRAM = memAvailRAM.PhysicalMemory.Available;
if memAvailRAM*memoryCoef < (prod(outSz)*4+sz(1)*sz(2)*256*4)
    error('Insufficient memory available for keeping the processed file.');
end

useGPU     = strcmpi(procType, 'gpu');
useOverlap = strcmpi(contrastType, 'temporal') && strcmpi(decimMethod, 'leaking');
H = floor(kernelSize/2);        % the exact half-width of getK's window - the kernel is odd
satValue = intmax(stream.dataType);

data = zeros(outSz,'single');
time = zeros(outT,1,'double');
intensityMetrics=zeros(sz(1),sz(2),3);
intensityMetrics(:,:,1:2)=1;
if useGPU
    %reduced on the card alongside the chunk and brought home once, after the loop
    intensityMetrics=gpuArray(intensityMetrics);
end

%The chunk getK is handed is the batch PLUS the overlap tail carried over from the
%batch before it, so the tail has to come out of the batch budget rather than be
%discovered afterwards.  With C the chunk length the tail is
%C-decimFactor*floor((C-H)/decimFactor)+decimFactor*ceil(H/decimFactor), and that is
%bounded by 2*H+2*decimFactor.  The bound is arithmetic on H and decimFactor alone, so
%the odd-kernel rule above does not move it; what the rule buys is that H is the window's
%true half-width, which is what makes this tail the right thing to carry.
if useOverlap
    maxTail = 2*H + 2*decimFactor;
else
    maxTail = 0;
end

batchSize=(memAvailRAM*memoryCoef-prod(outSz)*4)./(sz(1)*sz(2)*(5+6*4/decimFactor)); % logic: dataOUT_single+batch_dataIN_uint8+batch_dataOUT_intermediate_x4_single=memAvailRAM*memoryCoef
if useGPU
    %The chunk is uploaded whole and stays on the card while getK works on it, so the
    %host formula above is only the upper bound.  Two device limits sit under it.  A
    %gpuArray cannot hold more than 2^31-1 elements - on a 1496x896 camera that alone
    %caps the chunk at 1602 frames, and it is what makes "upload the whole batch"
    %impossible without a cap.  And the chunk plus the one same-sized temporary the
    %metrics below make of it may take only a share of the card: half the coefficient,
    %because getK's own splitter then budgets against what is left.  Wrappers/
    %runInternalCycle.m sizes its batch from AvailableMemory the same way.
    gpuDev = gpuDevice;
    bytesPerFrame = sz(1)*sz(2)*(2*stream.dataSize+1);
    batchSize = min([batchSize, ...
        floor((2^31-1)./(sz(1)*sz(2))) - maxTail, ...
        floor(gpuDev.AvailableMemory*memoryCoef/2./bytesPerFrame) - maxTail]);
end
batchSize = min(floor(batchSize / decimFactor) * decimFactor,sz(3));
if batchSize <= 0, batchSize = decimFactor; end

batchNum = ceil(sz(3) / batchSize);
curSize = 1;
nextStartIdx = 1;

%The chunk buffer is exactly as long as the chunk getK is handed, so the upload inside
%the loop takes the whole of it and never a slice of it.  Empty rather than [] so its
%third size is 0 and the reallocation test below cannot match by accident.
chunkIntensity = zeros(sz(1),sz(2),0,stream.dataType);
chunkTime      = zeros(0,1,'double');
tailIntensity  = zeros(sz(1),sz(2),0,stream.dataType);
tailTime       = zeros(0,1,'double');
tailLen        = 0;
timeStamp      = 0;

for i=1:batchNum
    curBatchSize=min(batchSize, sz(3)-(i-1)*batchSize);
    nChunk = tailLen + curBatchSize;
    %The chunk length settles after the first batch, so this reallocates twice in a
    %run rather than once per batch - and nothing here duplicates a whole chunk while
    %the original is still live, which the cat() it replaced did every round.
    if size(chunkIntensity,3) ~= nChunk
        chunkIntensity = zeros(sz(1),sz(2),nChunk,stream.dataType);
        chunkTime      = zeros(nChunk,1,'double');
    end
    if tailLen > 0
        chunkIntensity(:,:,1:tailLen) = tailIntensity;
        chunkTime(1:tailLen)          = tailTime;
    end
    [chunkIntensity(:,:,tailLen+1:nChunk),chunkTime(tailLen+1:nChunk),stream] = ...
        readRLS(stream,'FramesToRead',curBatchSize);
    if i==1, timeStamp = chunkTime(1); end

    %ONE upload per batch: the metrics and getK below read the same copy, which is what
    %keeps getK's own batch loop off the host - it slices the input per sub-batch, and
    %on a resident input that slice is a device copy rather than a host one.
    if useGPU
        chunkProc = gpuArray(chunkIntensity);
    else
        chunkProc = chunkIntensity;
    end

    %Only the frames read THIS round are counted - the tail in front of them was
    %counted when it was read.  sum() over an integer type accumulates in double and
    %returns double; the single() copy this replaced cost four bytes per element,
    %was made twice, and was the less accurate of the two.
    if tailLen > 0
        newFrames = chunkProc(:,:,tailLen+1:nChunk);
    else
        newFrames = chunkProc;
    end
    intensityMetrics(:,:,1)=intensityMetrics(:,:,1)-sum(newFrames==0,3)./sz(3);
    intensityMetrics(:,:,2)=intensityMetrics(:,:,2)-sum(newFrames==satValue,3)./sz(3);
    intensityMetrics(:,:,3)=intensityMetrics(:,:,3)+sum(newFrames,3)./sz(3);
    clear newFrames

    [tmp1,tmp2]=getK(chunkProc,contrastType,'time',chunkTime,'kernelSize',kernelSize,'procType', procType,'decimFactor',decimFactor,'decimMethod',decimMethod,'memoryCoef',memoryCoef);
    clear chunkProc

    if useOverlap
        startIdx = nextStartIdx;
        if i == batchNum
            endIdx = size(tmp1, 3);
        else
            validInputLimit = nChunk - H;
            endIdx = floor(validInputLimit / decimFactor);
        end
        nextFrameStart = (endIdx * decimFactor) + 1;
        rawStartNeeded = nextFrameStart - H;
        alignedStart = floor((rawStartNeeded - 1) / decimFactor) * decimFactor + 1;
        alignedStart = max(1, alignedStart);
        framesSavedInPrevBatch = endIdx - (alignedStart - 1) / decimFactor;
        nextStartIdx = framesSavedInPrevBatch + 1;
        tailIntensity = chunkIntensity(:, :, alignedStart:nChunk);
        tailTime = chunkTime(alignedStart:nChunk);
        tailLen = nChunk - alignedStart + 1;
    else
        startIdx = 1;
        endIdx = size(tmp1, 3);
    end

    numToSave = endIdx - startIdx + 1;
    if numToSave > 0
        saveEnd = curSize + numToSave - 1;
        if saveEnd > outT
            numToSave = outT - curSize + 1;
            endIdx = startIdx + numToSave - 1;
        end
        data(:, :, curSize : curSize+numToSave-1) = tmp1(:, :, startIdx:endIdx);
        time(curSize : curSize+numToSave-1) = tmp2(startIdx:endIdx);
        curSize = curSize + numToSave;
    end
end
time=(time-time(1))./1000; %conversion to seconds
if useGPU, intensityMetrics = gather(intensityMetrics); end
end
