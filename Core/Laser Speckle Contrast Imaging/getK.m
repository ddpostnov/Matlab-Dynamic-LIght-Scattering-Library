%getK - Calculates Temporal or Spatial Contrast 
%
% Syntax:
%    [data, time] = getK(dataIn, contrastType, 'Name', Value, ...)
%
% Description:
%    Computes contrast images from a 3D intensity matrix. It supports both
%    Temporal (contrast over time per pixel) and Spatial (contrast over
%    space per frame) methods. It includes options for temporal decimation,
%    specific decimation methods (leaking vs sharp), and CPU/GPU computation.
%
% Inputs:
%    dataIn       - Raw  data as a 3D numeric matrix [Y, X, Time].
%    contrastType - String specifying the analysis type: 'spatial' or 'temporal'.
%
% Optional Name-Value Pair Arguments:
%    'kernelSize'  - (Numeric) Size of the window used for contrast calculation.
%                    Must be ODD: the window is centred on the pixel it describes,
%                    and an even size has no centre - the 'cpu' and 'gpu' branches
%                    would place it on opposite sides and return different contrast.
%                    Default: 25 for 'temporal', 5 for 'spatial'.
%    'procType'    - (String) Processor to use: 'cpu' or 'gpu'. Default: 'gpu'.
%    'decimFactor' - (Numeric) Factor for temporal decimation/averaging.
%                    Default: 1 (no decimation).%                    
%    'decimMethod' - (String) Method for temporal decimation:
%                    'leaking' (Default): Uses a sliding window (convolution)
%                    followed by averaging. Preserves temporal resolution but
%                    introduces some temporal cross-talk for temporal contrast.
%                    'sharp': strictly non-overlapping block averaging.
%                    Only valid for 'temporal' contrast. Requires 'decimFactor'
%                    to be an integer multiple of 'kernelSize'.
%    'memoryCoef'  - (Numeric) Fraction of available RAM/GPU memory to use
%                    for batch processing (0 to 1). Default: 0.8.
%    'time'        - (Vector) Input time vector matching size(dataIn, 3).
%                    If provided, the function returns a decimated time vector.
%    'outputType'  - 'cpu' (default) returns a normal array; 'gpu' returns a gpuArray
%                    and skips the gather.  Use it when the caller reduces or
%                    decimates the contrast immediately and never needs the
%                    full-resolution cube in RAM - a 100 GB recording cannot hold it.
%                    Requires procType 'gpu'.
%
% Outputs:
%    data         - Processed contrast data as a 3D matrix [Y, X, OutputTime].
%                   OutputTime = floor(InputTime / decimFactor).
%    time         - (Optional) Decimated time vector corresponding to the
%                   output frames, [OutputTime x 1] - a COLUMN whatever the
%                   orientation of the 'time' input and whatever decimFactor.
%
% Examples:
%    % 1. Standard Temporal Contrast (GPU, default kernel=25):
%    K = getK(rawStack, 'temporal');
%
%    % 2. Spatial Contrast on CPU with 7x7 kernel:
%    K = getK(rawStack, 'spatial', 'procType', 'cpu', 'kernelSize', 7);
%
%    % 3. Temporal Contrast with "Sharp" Decimation (Kernel=25, Decim=50):
%    [K, tOut] = getK(rawStack, 'temporal', 'decimFactor', 50, ...
%                     'decimMethod', 'sharp', 'time', tIn);
%
% See also: getContrastFromRLS, getContrastFromMRAW, readRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026
%------------- BEGIN CODE --------------
function [data, time] = getK(dataIn,contrastType, varargin)
p = inputParser;
p.KeepUnmatched = false;

addRequired(p, 'dataIn', @isnumeric);
addRequired(p, 'contrastType', @(x) any(validatestring(x, {'spatial', 'temporal'})));

addParameter(p, 'kernelSize', [], @isnumeric);
addParameter(p, 'procType', 'gpu', @(x) any(validatestring(x, {'cpu', 'gpu'})));
addParameter(p, 'decimFactor', 1, @isnumeric);
addParameter(p, 'decimMethod', 'leaking', @(x) any(validatestring(x, {'leaking', 'sharp'})));
addParameter(p, 'memoryCoef', 0.8, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'time', [], @isnumeric);
addParameter(p, 'outputType', 'cpu', @(x) any(validatestring(x, {'cpu', 'gpu'})));

parse(p, dataIn, contrastType, varargin{:});

kernelSize  = p.Results.kernelSize;
procType    = p.Results.procType;
decimFactor  = p.Results.decimFactor;
decimMethod = p.Results.decimMethod;
memoryCoef  = p.Results.memoryCoef;
timeIn      = p.Results.time;
outputType  = p.Results.outputType;
time        = [];

if strcmpi(outputType,'gpu') && strcmpi(procType,'cpu')
    error('getK:InvalidOutput', ...
        "outputType 'gpu' requires procType 'gpu'.");
end

if isempty(kernelSize)
    if strcmpi(contrastType, 'temporal')
        kernelSize = 25;
    elseif strcmpi(contrastType, 'spatial')
        kernelSize = 5;
    end
end

%A contrast window is CENTRED on the pixel it describes, and an even size has no centre.
%The two branches below then stop agreeing: conv2/convn 'same' takes an even window as K
%taps placed to one side, movstd's [floor(K/2) floor(K/2)] takes it as K+1 taps placed to
%the other.  Measured on a real recording, 'gpu' and 'cpu' return the same contrast to
%rounding at K=25 and K=49 (1.6e-5, 1.3e-5) and differ by 21% at K=50 and 35% at K=24.
%The callers assume the centre too - getContrastFromRLS sizes its overlap tail from
%floor(kernelSize/2), runContrastInternalCycle pads its read with it - so the size is refused
%here rather than centred one way and left to disagree with everything downstream.
if ~isscalar(kernelSize) || ~isfinite(kernelSize) || kernelSize < 1 || mod(kernelSize,2) ~= 1
    error('getK:EvenKernel', ...
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

sz = size(dataIn);
outT = floor(sz(3) / decimFactor);
outSz = [sz(1), sz(2), outT];
framesToProc = outT * decimFactor;

%The output is the thing being guarded, so the guard follows it: a 'cpu' output has
%to fit in RAM, a 'gpu' output has to fit on the card.  memory() costs ~50 ms, so it
%is called once here and reused by the cpu branches below.
memAvailRAM = [];
if strcmpi(outputType,'gpu')
    memAvailGPU = gpuDevice;
    if memAvailGPU.AvailableMemory*memoryCoef<prod(outSz)*4
        error('Insufficient GPU memory available for keeping the processed file.');
    end
else
    [~, memAvailRAM] = memory;
    memAvailRAM = memAvailRAM.PhysicalMemory.Available;
    if memAvailRAM*memoryCoef<prod(outSz)*4
        error('Insufficient memory available for keeping the processed file.');
    end
end

%The time base is a COLUMN, [T x 1], everywhere in this library - frames run down the
%rows, matching source.time and results.sData's [nT x nSeg].  Both branches say so
%explicitly: the decimated one returned a ROW from mean(...,1) while the undecimated one
%passed whatever the caller handed in, so the SAME function's second output changed
%orientation with decimFactor.
if ~isempty(timeIn)
    if length(timeIn) ~= sz(3)
        error('Time vector length must match the number of frames (size(dataIn,3)).');
    end
    if decimFactor > 1
        time = mean(reshape(timeIn(1:framesToProc), decimFactor, outT),1).';
    else
        time = timeIn(:);
    end
end


switch contrastType
    case 'temporal'
        dataIn=reshape(dataIn,sz(1)*sz(2),sz(3));
        data=allocOut([sz(1)*sz(2),outT],outputType);
        switch procType
            case 'gpu'
                memAvailGPU = gpuDeviceBudget(dataIn);
                batchNum=ceil(max((7*numel(dataIn)*4)./(memAvailGPU*memoryCoef),1)); % The actual memory use is approx 3*numel, we use 7*numel because for high-end consumer GPU performance is better with small batches.
                rowEdges=rowBlockEdges(size(dataIn,1),batchNum);
                kernel=gpuArray(ones(1,kernelSize,'single'))./kernelSize;
                normMap = conv2(gpuArray.ones(1, sz(3), 'single'), kernel, 'same');
                for b=1:numel(rowEdges)-1
                    i =rowEdges(b)+1;
                    i2=rowEdges(b+1);
                    if decimFactor > 1 && strcmpi(decimMethod, 'sharp')
                        dataGPU=single(gpuArray(dataIn(i:i2,1:framesToProc)));
                        numSubBlocks = decimFactor / kernelSize;
                        currentBatchSize = size(dataGPU, 1);
                        dataGPU = reshape(dataGPU, currentBatchSize, kernelSize, numSubBlocks, outT);
                        m = mean(dataGPU, 2);
                        dataGPU = mean(dataGPU.^2, 2);
                        dataGPU = sqrt(max(dataGPU - m.^2, 0));
                        dataGPU = dataGPU ./ m;
                        dataGPU = mean(dataGPU, 3);
                        dataGPU = reshape(dataGPU, currentBatchSize, outT);
                    else
                        dataGPU=single(gpuArray(dataIn(i:i2,:)));
                        m=conv2(dataGPU, kernel, 'same')./normMap;
                        dataGPU=conv2(dataGPU.^2, kernel, 'same')./normMap;
                        dataGPU = sqrt(max(dataGPU - m.^2, 0));
                        dataGPU=dataGPU./m;
                        if decimFactor > 1
                            dataGPU=dataGPU(:,1:framesToProc);                            
                            currentBatchSize = size(dataGPU,1);
                            dataGPU = reshape(dataGPU, currentBatchSize, decimFactor, outT);
                            dataGPU = mean(dataGPU, 2);
                            dataGPU = reshape(dataGPU, currentBatchSize, outT);
                        end
                    end
                    if strcmpi(outputType,'gpu')
                        data(i:i2,:)=dataGPU;
                    else
                        data(i:i2,:)=gather(dataGPU);
                    end
                end
            case 'cpu'
                if isempty(memAvailRAM)                 % the 'gpu' outputType branch skipped it above
                    [~, memAvailRAM] = memory;
                    memAvailRAM = memAvailRAM.PhysicalMemory.Available;
                end
                batchNum=ceil(max((4*numel(dataIn)*4)./(memAvailRAM*memoryCoef),1)); % operations below should take aprox 2x batch array size, but, we allocate 3x just for safety
                rowEdges=rowBlockEdges(size(dataIn,1),batchNum);
                for b=1:numel(rowEdges)-1
                    i =rowEdges(b)+1;
                    i2=rowEdges(b+1);
                    chunk = single(dataIn(i:i2, :));
                    if decimFactor > 1 && strcmpi(decimMethod, 'sharp')
                        chunk = chunk(:, 1:framesToProc);
                        numSubBlocks = decimFactor / kernelSize;
                        currentBatchSize = size(chunk, 1);
                        chunk = reshape(chunk, currentBatchSize, kernelSize, numSubBlocks, outT);
                        chunk=std(chunk, 1, 2)./mean(chunk, 2);
                        chunk = mean(chunk, 3);
                        chunk = reshape(chunk, currentBatchSize, outT);
                    else
                        chunk = movstd(chunk, [floor(kernelSize/2),floor(kernelSize/2)], 1, 2, 'Endpoints','shrink')./movmean(chunk, [floor(kernelSize/2),floor(kernelSize/2)], 2, 'Endpoints','shrink');
                        if decimFactor > 1
                            chunk = chunk(:, 1:framesToProc);
                            currentBatchSize = size(chunk,1);
                            chunk = reshape(chunk, currentBatchSize, decimFactor, outT);
                            chunk = mean(chunk, 2);
                            chunk = reshape(chunk, currentBatchSize, outT);
                        end
                    end
                    data(i:i2,:) = chunk;
                end

        end
        data=reshape(data,sz(1),sz(2),outT);
    case 'spatial'
        kernel=single(ones(kernelSize,kernelSize,1)./(kernelSize^2));
        normMap=single(conv2(ones(size(dataIn,1),size(dataIn,2)), gather(kernel), 'same'));

        switch procType
            case 'gpu'
                data=allocOut(outSz,outputType);
                memAvailGPU = gpuDeviceBudget(dataIn);
                batchNum=ceil(max((7*numel(dataIn)*4)./(memAvailGPU*memoryCoef),1)); % The actual memory use is approx 3*numel, we use 7*numel because for high-end consumer GPU performance is better with small batches.
                batchSize=floor(size(dataIn,3)./batchNum);

                batchSize = floor(batchSize/decimFactor) * decimFactor;
                if batchSize == 0, batchSize = decimFactor; end

                kernel=gpuArray(kernel);
                normMap=gpuArray(normMap);
                for i=1:batchSize:framesToProc
                    i2=min(i+batchSize-1,framesToProc);
                    dataGPU=single(gpuArray(dataIn(:,:,i:i2)));
                    m=convn(dataGPU, kernel, 'same')./normMap;
                    dataGPU=convn(dataGPU.^2, kernel, 'same')./normMap;
                    dataGPU = sqrt(max(dataGPU - m.^2, 0));
                    dataGPU=dataGPU./m;

                    if decimFactor > 1
                        [h, w, t] = size(dataGPU);
                        dataGPU = reshape(dataGPU, h, w, decimFactor, t/decimFactor);
                        dataGPU = mean(dataGPU, 3);
                        dataGPU = reshape(dataGPU, h, w, t/decimFactor);
                    end

                    j = (i-1)/decimFactor + 1;
                    j2 = j + size(dataGPU,3) - 1;
                    if strcmpi(outputType,'gpu')
                        data(:,:,j:j2)=dataGPU;
                    else
                        data(:,:,j:j2)=gather(dataGPU);
                    end
                end
            case 'cpu' %convn can be replaced with imboxfilt for a minor (less than 5% boost in performance)
                data=allocOut(outSz,outputType);
                if isempty(memAvailRAM)                 % the 'gpu' outputType branch skipped it above
                    [~, memAvailRAM] = memory;
                    memAvailRAM = memAvailRAM.PhysicalMemory.Available;
                end
                batchNum=ceil(max((4*numel(dataIn)*4)./(memAvailRAM*memoryCoef),1)); % operations below should take aprox 2x batch array size, but, we allocate 3x just for safety
                batchSize=floor(size(dataIn,3)./batchNum);

                batchSize = floor(batchSize/decimFactor) * decimFactor;
                if batchSize == 0, batchSize = decimFactor; end

                for i=1:batchSize:framesToProc
                    i2=min(i+batchSize-1,framesToProc);
                    chunk = single(dataIn(:,:,i:i2));
                    m=convn(chunk, kernel, 'same')./normMap;
                    chunk=convn(chunk.^2, kernel, 'same')./normMap;
                    chunk = sqrt(max(chunk - m.^2, 0));
                    chunk=chunk./m;

                    if decimFactor > 1
                        [h, w, t] = size(chunk);
                        chunk = reshape(chunk, h, w, decimFactor, t/decimFactor);
                        chunk = mean(chunk, 3);
                        chunk = reshape(chunk, h, w, t/decimFactor);
                    end

                    j = (i-1)/decimFactor + 1;
                    j2 = j + size(chunk,3) - 1;
                    data(:,:,j:j2) = chunk;
                end
        end
end
end

% =====================================================================
function e = rowBlockEdges(nRows,batchNum)
%rowBlockEdges  Pixel-row block boundaries, all of the same height to within one row.
%   Block k is the rows e(k)+1 : e(k+1).
%
%   The two temporal branches used to step 1:floor(nRows/batchNum):nRows, which leaves a
%   RAGGED final block whenever batchNum does not divide the row count - on a 1496x896
%   camera at batchNum 3 that block is exactly ONE row.  Height 1 is the one height at
%   which the decimation mean(reshape(x,rows,decimFactor,outT),2) rounds differently:
%   measured on an RTX 4090 on 2026-08-04, heights 2, 3, 4, 7, 8, 16, 17, 1000, 4999,
%   5000, 6667 and 10000 are all BIT-IDENTICAL to the whole-height answer while height 1
%   is not (max|d| 1.5e-05, 1.6e-07 relative).  So the last pixel of the frame was
%   computed by different arithmetic from every other pixel - and batchNum comes from
%   gpuDevice.AvailableMemory, so which pixels that hit depended on how much of the card
%   happened to be free.  Two runs of the same recording could disagree there.
%   An even split never makes that block, which is what makes the result independent of
%   the split rather than merely less often wrong.  (conv2, the decimFactor 1 path, is
%   invariant at every height including 1, so nothing that does not decimate moves.)
%
%   Capping batchNum at half the row count is the same guarantee stated for the other
%   end: it also removes the missing batchSize==0 guard - the two SPATIAL branches have
%   had one for as long as they have had a split - which on a batchNum above the row
%   count gave a step of 0 and a loop that never advanced.
nB = max(1,min(batchNum,floor(nRows/2)));
e  = round(linspace(0,nRows,nB+1));
end

% =====================================================================
function m = gpuDeviceBudget(dataIn)
%gpuDeviceBudget  What the batch loops above may spend, in bytes.
%   AvailableMemory answers what is free RIGHT NOW, so a dataIn that is already on the
%   card has been subtracted from it - while the 7*numel*4 estimate the callers divide
%   by it is a whole-footprint figure that still counts dataIn.  Left alone the two
%   disagree by the input's own size, the loop splits finer than the card needs, and
%   every extra split is one more real device copy of a sub-batch.  So a resident input
%   is added back, and a host input is not.
m = gpuDevice;
m = m.AvailableMemory;
if isa(dataIn,'gpuArray')
    m = m + numel(dataIn)*numel(typecast(cast(0,underlyingType(dataIn)),'uint8'));
end
end

% =====================================================================
function d = allocOut(sz,outputType)
%allocOut  The result buffer, on the device the caller asked for.  'gpu' skips the
%   gather in the batch loops below - which is the whole point of the option for
%   callers that reduce or decimate the contrast immediately (runContrastInternalCycle).
if strcmpi(outputType,'gpu')
    d = zeros(sz,'single','gpuArray');
else
    d = zeros(sz,'single');
end
end
