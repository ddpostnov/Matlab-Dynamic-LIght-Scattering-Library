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
%    'kernelSize'  - contrast kernel size (default 5 spatial, 25 temporal).
%    'procType'    - 'cpu' or 'gpu' (default 'gpu').
%    'decimFactor' - frames averaged and kept during decimation (default 1).
%    'decimMethod' - 'leaking' (overlapping, default) or 'sharp' (temporal only).
%    'memoryCoef'  - fraction of free RAM used for batching (default 0.8).
%    'progressFcn' - progress sink, progressFcn(frac,label).  Absent (the default)
%                    the per-batch line is printed to the command window exactly as
%                    it always has been, so every existing caller is unaffected; a
%                    wrapper passes reportProgress's handle to route it instead.
%
% Outputs:
%    data             - contrast data, 3-D [y x t] single.
%    time             - relative time vector, seconds.
%    timeStamp        - absolute timestamp of the first frame.
%    intensityMetrics - [y x 4] per-pixel quality map: fraction of non-zero and
%                       of non-saturated frames, mean intensity, and mean squared
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
% Last revision: 07-July-2026

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
addParameter(p, 'progressFcn', [], @(x) isempty(x) || isa(x,'function_handle'));
parse(p, rlsFileName, contrastType, varargin{:});

kernelSize  = p.Results.kernelSize;
procType    = p.Results.procType;
decimFactor = p.Results.decimFactor;
decimMethod = p.Results.decimMethod;
memoryCoef  = p.Results.memoryCoef;
progressFcn = p.Results.progressFcn;

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


[tmp1,timeStamp,stream]=readRLS(rlsFileName,'KeepOpen',true,'FramesToRead',1);
sz = [stream.sizeY,stream.sizeX,stream.sizeT];
outT = floor(sz(3) / decimFactor);
outSz = [sz(1), sz(2), outT];

[~, memAvailRAM] = memory;
memAvailRAM = memAvailRAM.PhysicalMemory.Available;
if memAvailRAM*memoryCoef < (prod(outSz)*4+sz(1)*sz(2)*256*4)
    error('Insufficient memory available for keeping the processed file.');
end

data = zeros(outSz,'single');
time = zeros(outT,1,'double');
intensityMetrics=zeros(sz(1),sz(2),4);
intensityMetrics(:,:,1:2)=1;



batchSize=(memAvailRAM*memoryCoef-prod(outSz)*4)./(sz(1)*sz(2)*(5+6*4/decimFactor)); % logic: dataOUT_single+batch_dataIN_uint8+batch_dataOUT_intermediate_x4_single=memAvailRAM*memoryCoef
batchSize = min(floor(batchSize / decimFactor) * decimFactor,sz(3));
if batchSize == 0, batchSize = decimFactor; end

batchNum = ceil(sz(3) / batchSize);
curSize = 1;
nextStartIdx = 1;

useOverlap = strcmpi(contrastType, 'temporal') && strcmpi(decimMethod, 'leaking');
tailIntensity = [];
tailTime = [];
H = floor(kernelSize/2);

tic
for i=1:batchNum
    curBatchSize=min(batchSize, sz(3)-(i-1)*batchSize);
    if i==1
        chunkIntensity=zeros(sz(1),sz(2),curBatchSize,class(tmp1));
        chunkTime=zeros(curBatchSize,1,'double');
        chunkIntensity(:,:,1)=tmp1;
        chunkTime(1)=timeStamp;
        [chunkIntensity(:,:,2:end),chunkTime(2:end),stream]=readRLS(stream,'FramesToRead',curBatchSize-1);
    else
        [chunkIntensity,chunkTime,stream]=readRLS(stream,'FramesToRead',curBatchSize);
    end
    intensityMetrics(:,:,1)=intensityMetrics(:,:,1)-double(sum(chunkIntensity==0,3))./sz(3);
    intensityMetrics(:,:,2)=intensityMetrics(:,:,2)-double(sum(chunkIntensity==intmax(class(chunkIntensity)),3))./sz(3);
    intensityMetrics(:,:,3)=intensityMetrics(:,:,3)+sum(single(chunkIntensity),3)./sz(3);
    intensityMetrics(:,:,4)=intensityMetrics(:,:,4)+sum(single(chunkIntensity).*single(chunkIntensity),3)./sz(3);


    if ~isempty(tailIntensity)
        chunkIntensity = cat(3, tailIntensity, chunkIntensity);
        chunkTime = cat(1, tailTime, chunkTime);
    end
    [tmp1,tmp2]=getK(chunkIntensity,contrastType,'time',chunkTime,'kernelSize',kernelSize,'procType', procType,'decimFactor',decimFactor,'decimMethod',decimMethod,'memoryCoef',memoryCoef);

    if useOverlap
        startIdx = nextStartIdx;
        if i == batchNum
            endIdx = size(tmp1, 3);
        else            
            validInputLimit = size(chunkIntensity, 3) - H;
            endIdx = floor(validInputLimit / decimFactor);
        end
        nextFrameStart = (endIdx * decimFactor) + 1;
        rawStartNeeded = nextFrameStart - H;
        alignedStart = floor((rawStartNeeded - 1) / decimFactor) * decimFactor + 1;
        alignedStart = max(1, alignedStart);
        framesSavedInPrevBatch = endIdx - (alignedStart - 1) / decimFactor;
        nextStartIdx = framesSavedInPrevBatch + 1;
        tailIntensity = chunkIntensity(:, :, alignedStart:end);
        tailTime = chunkTime(alignedStart:end);
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
    if isempty(progressFcn)
        fprintf('Batch %d/%d processed. Elapsed: %.2fs\n', i, batchNum, toc);
    else
        progressFcn(i/batchNum, 'Speckle contrast');
    end
end
time=(time-time(1))./1000; %conversion to seconds
fclose(stream.fId);
end
