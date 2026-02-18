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
%
% Outputs:
%    data         - Processed contrast data as a 3D matrix [Y, X, OutputTime].
%                   OutputTime = floor(InputTime / decimFactor).
%    time         - (Optional) Decimated time vector corresponding to the
%                   output frames.
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
% See also: readRLS.m
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 30-January-2026
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

parse(p, dataIn, contrastType, varargin{:});

kernelSize  = p.Results.kernelSize;
procType    = p.Results.procType;
decimFactor  = p.Results.decimFactor;
decimMethod = p.Results.decimMethod;
memoryCoef  = p.Results.memoryCoef;
timeIn      = p.Results.time;
time        = [];


if isempty(kernelSize)
    if strcmpi(contrastType, 'temporal')
        kernelSize = 25;
    elseif strcmpi(contrastType, 'spatial')
        kernelSize = 5;
    end
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

[~, memAvailRAM] = memory;
memAvailRAM = memAvailRAM.PhysicalMemory.Available;
if memAvailRAM*memoryCoef<prod(outSz)*4
    error('Insufficient memory available for keeping the processed file.');
end

if ~isempty(timeIn)
    if length(timeIn) ~= sz(3)
        error('Time vector length must match the number of frames (size(dataIn,3)).');
    end
    if decimFactor > 1
        time = mean(reshape(timeIn(1:framesToProc), decimFactor, outT),1);
    else
        time = timeIn;
    end
end


switch contrastType
    case 'temporal'
        dataIn=reshape(dataIn,sz(1)*sz(2),sz(3));
        data=zeros(sz(1)*sz(2),outT,'single');
        switch procType
            case 'gpu'
                memAvailGPU = gpuDevice;
                memAvailGPU = memAvailGPU.AvailableMemory;
                batchNum=ceil(max((7*numel(dataIn)*4)./(memAvailGPU*memoryCoef),1)); % The actual memory use is approx 3*numel, we use 7*numel because for high-end consumer GPU performance is better with small batches. 
                batchSize=floor(size(dataIn,1)./batchNum);
                kernel=gpuArray(ones(1,kernelSize,'single'))./kernelSize;
                normMap = conv2(gpuArray.ones(1, sz(3), 'single'), kernel, 'same');
                for i=1:batchSize:size(dataIn,1)
                    i2=min(i+batchSize-1,size(dataIn,1));
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
                    data(i:i2,:)=gather(dataGPU);
                end
            case 'cpu'
                [~, memAvailRAM] = memory;
                memAvailRAM = memAvailRAM.PhysicalMemory.Available;
                batchNum=ceil(max((4*numel(dataIn)*4)./(memAvailRAM*memoryCoef),1)); % operations below should take aprox 2x batch array size, but, we allocate 3x just for safety
                batchSize=floor(size(dataIn,1)./batchNum);
                for i=1:batchSize:size(dataIn,1)
                    i2=min(i+batchSize-1,size(dataIn,1));
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
                data=zeros(outSz,'single');
                memAvailGPU = gpuDevice;
                memAvailGPU = memAvailGPU.AvailableMemory;
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
                    data(:,:,j:j2)=gather(dataGPU);
                end
            case 'cpu' %convn can be replaced with imboxfilt for a minor (less than 5% boost in performance)
                data=zeros(outSz,'single');
                [~, memAvailRAM] = memory;
                memAvailRAM = memAvailRAM.PhysicalMemory.Available;
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
