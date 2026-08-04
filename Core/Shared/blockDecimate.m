%blockDecimate - Non-overlapping block-mean decimation of a vector.
%
%   blockDecimate(v,avgN) low-passes v with a boxcar of length avgN and then keeps
%   every avgN-th sample, so the result is the mean of consecutive non-overlapping
%   blocks.  It is the library's decimation rule for a TIME SERIES that is about to be
%   analysed at a lower rate, and it exists as one function so a decimated time axis
%   and the decimated trace on it can never be built by two slightly different rules.
%
%   IT LOW-PASSES BEFORE IT SUB-SAMPLES, and that is the whole point.  A plain stride
%   (v(1:avgN:end)) is not decimation, it is aliasing: the cardiac and respiratory
%   content of an LSCI trace folds straight down into the slow band a drug response or
%   a vasomotion analysis lives in, and it arrives looking exactly like signal.
%
%   THE WINDOW SHOULD BE ODD, which is what blockDecimateWindow returns, so the block
%   mean is centred on the sample it replaces and the decimated time axis is the
%   decimated clock rather than a clock shifted by half a block.
%
%   'Endpoints','discard' means a partial block at the end is dropped rather than
%   averaged over fewer samples - a final point computed from three samples where every
%   other point came from twenty-five is a point with a different noise level, and
%   nothing downstream would know.
%
% Syntax:
%    avgN = blockDecimate('window', fsRaw, tgtFS)
%    v    = blockDecimate(v, avgN)
%
% Inputs:
%    fsRaw - the series' own sampling rate, Hz.
%    tgtFS - the wanted rate, Hz.  Empty, non-finite, non-positive or at/above fsRaw
%            all mean "keep every sample" and give avgN = 1.
%    v     - the vector to decimate.
%    avgN  - the block length (1 returns v unchanged).
%
% Outputs:
%    avgN  - ('window') the odd block length that takes fsRaw down to about tgtFS.
%    v     - the decimated vector, floor((numel(v)-avgN)/avgN)+1 samples long.
%
% Example:
%    avgN = blockDecimate('window', 25, 1);      % 25 Hz -> about 1 Hz
%    tD   = blockDecimate(results.time, avgN);
%    yD   = blockDecimate(results.sData(:,7), avgN);
%
% See also: fitVasoreactivity, runFitVasoreactivity, movmean
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

%------------- BEGIN CODE --------------
function out = blockDecimate(arg1,arg2,arg3)

if (ischar(arg1)||isstring(arg1)) && strcmp(char(arg1),'window')
    out=windowFor(arg2,arg3);
    return
end
out=decimate(arg1,arg2);
end

% =====================================================================
function avgN=windowFor(fsRaw,tgtFS)
%windowFor  The odd block length that takes fsRaw down to about tgtFS (1 = no change).
avgN=1;
if isempty(tgtFS) || ~isscalar(tgtFS) || ~isfinite(tgtFS) || tgtFS<=0, return, end
if ~isfinite(fsRaw) || fsRaw<=0 || tgtFS>=fsRaw, return, end
avgN=max(1,floor(fsRaw./tgtFS/2)*2+1);
end

% =====================================================================
function v=decimate(v,avgN)
%decimate  The block mean itself (a window of 1 is the identity, not a movmean).
if isempty(avgN) || avgN<=1, return, end
v=movmean(v,avgN,'Endpoints','discard');
v=v(1:avgN:end);
end
