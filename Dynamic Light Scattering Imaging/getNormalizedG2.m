%getNormalizedG2 - Normalized intensity autocorrelation g2 from an image stack.
%
% Syntax:
%    g2 = getNormalizedG2(dataIn)
%    g2 = getNormalizedG2(dataIn, 'Name', Value, ...)
%
% Description:
%    Computes the temporal, per-pixel normalized intensity autocorrelation
%    function
%
%        g2(tau) = <I(t) I(t+tau)> / <I(t)>^2,
%
%    where the averages <.> run over a fixed base window of the recording.
%    g2 is the quantity fitted by the DLSI models (see getFitDLSI) and
%    the input to the simplified decorrelation-time estimate getTauC.
%
%    The stack is processed in contiguous spatial tiles so that recordings
%    larger than the available GPU/CPU memory can be handled; the tile size
%    is derived automatically from 'memoryCoef' unless 'blocksN' is given.
%    Tiling does not affect the result (each pixel is independent).
%
% Inputs:
%    dataIn - Raw intensity stack as a 3-D numeric matrix [Y, X, Time].
%
% Optional Name-Value Pair Arguments:
%    'lagMax'       - (Numeric) Largest lag, in frames. g2 has lagMax+1
%                     points spanning lags 0..lagMax. Default:
%                     min(100, size(dataIn,3)-2).
%    'sampleLength' - (Numeric) Number of leading frames used to build the
%                     averaging window. The base window has
%                     sampleLength-lagMax frame pairs. Default: [] (use all
%                     frames, i.e. size(dataIn,3)).
%    'procType'     - (String) 'gpu' (default) or 'cpu'. 'gpu' falls back to
%                     'cpu' automatically when no GPU is available.
%    'blocksN'      - (Numeric) [nRowTiles, nColTiles] spatial tiling. Default:
%                     [] (row tiles chosen automatically to fit memory).
%    'memoryCoef'   - (Numeric) Fraction of free device memory used per tile
%                     (0..1). Default: 0.5.
%
% Outputs:
%    g2 - Normalized intensity autocorrelation, single [Y, X, lagMax+1].
%         g2(:,:,1) is the lag-0 term (~1+beta) and g2(:,:,end) the longest
%         lag. Convert lag index k to time with (k-1)/fps seconds.
%
% Examples:
%    % 1. Default (all frames, 100 lags, GPU with CPU fallback):
%    g2 = getNormalizedG2(rawStack);
%
%    % 2. 151 lags on the CPU:
%    g2 = getNormalizedG2(rawStack, 'lagMax', 151, 'procType', 'cpu');
%
%    % 3. Lag times for plotting:
%    g2   = getNormalizedG2(rawStack, 'lagMax', 151);
%    lags = (0:size(g2,3)-1)/fps;   % seconds
%
% Dependencies: Parallel Computing Toolbox (only for 'procType' 'gpu').
% See also: getTauC, getFitDLSI, getDynamicSpeckles, getK
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 22-July-2026

%------------- BEGIN CODE --------------
function g2 = getNormalizedG2(dataIn, varargin)

p = inputParser;
p.KeepUnmatched = false;
addRequired(p, 'dataIn', @(x) isnumeric(x) && ndims(x)==3);
addParameter(p, 'lagMax', [], @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'sampleLength', [], @(x) isempty(x) || (isscalar(x) && x>=2));
addParameter(p, 'procType', 'gpu', @(x) any(validatestring(x, {'cpu','gpu'})));
addParameter(p, 'blocksN', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'memoryCoef', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
parse(p, dataIn, varargin{:});

sizeY = size(dataIn,1);
sizeX = size(dataIn,2);
sizeT = size(dataIn,3);

lagMax = p.Results.lagMax;
if isempty(lagMax)
    lagMax = min(100, sizeT-2);
end
lagMax = round(lagMax);

sampleLength = p.Results.sampleLength;
if isempty(sampleLength)
    sampleLength = sizeT;
end
sampleLength = round(min(sampleLength, sizeT));

% Number of frame pairs in the base averaging window.
windowL = sampleLength - lagMax;
if windowL < 2
    error('getNormalizedG2:tooFewFrames', ...
        'sampleLength (%d) must exceed lagMax (%d) by at least 2 frames.', ...
        sampleLength, lagMax);
end

% Device selection (gpu requested but unavailable -> cpu).
useGPU = strcmpi(p.Results.procType, 'gpu');
if useGPU
    try
        useGPU = gpuDeviceCount > 0;
    catch
        useGPU = false;
    end
end
if useGPU
    toDev   = @(x) gpuArray(single(x));
    gd      = gpuDevice;
    memAvail = gd.AvailableMemory;
else
    toDev   = @(x) single(x);
    [~, mem] = memory;
    memAvail = mem.PhysicalMemory.Available;
end

% Spatial tiling: automatic row tiles unless blocksN is given.
blocksN = p.Results.blocksN;
if isempty(blocksN)
    % Peak device use per row tile ~ rows*sizeX*4 * (sampleLength + 2*windowL).
    bytesPerRow = sizeX * 4 * (sampleLength + 2*windowL);
    rowsPerTile = max(1, floor(memAvail*p.Results.memoryCoef / bytesPerRow));
    nRowTiles = max(1, ceil(sizeY / rowsPerTile));
    nColTiles = 1;
else
    nRowTiles = max(1, round(blocksN(1)));
    nColTiles = max(1, round(blocksN(2)));
end
yEdges = round(linspace(0, sizeY, nRowTiles+1));
xEdges = round(linspace(0, sizeX, nColTiles+1));

g2 = zeros(sizeY, sizeX, lagMax+1, 'single');
for jy = 1:nRowTiles
    rows = yEdges(jy)+1 : yEdges(jy+1);
    for jx = 1:nColTiles
        cols = xEdges(jx)+1 : xEdges(jx+1);
        subData = toDev(dataIn(rows, cols, :));
        base    = subData(:,:,1:windowL);          % fixed base window
        denom   = mean(base, 3).^2;                % <I>^2 over the base window
        for lag = 0:lagMax
            g2(rows, cols, lag+1) = gather( ...
                mean(base .* subData(:,:,1+lag:windowL+lag), 3) ./ denom);
        end
    end
end
end
