%cutMyographIntervals  Slice a whole-recording diameter result into sub-intervals
%
%   intervals = cutMyographIntervals(big,ivT,names) takes BIG, the single
%   whole-recording interval produced by getMyographDiameter (a struct with
%   fields name/measures/time/idxL/idxR/diameter/mask/valid), and cuts it into the
%   time windows ivT = [start end; ...] (seconds) with the matching NAMES.  The
%   result is a struct array in exactly the same format as getMyographDiameter
%   would return for predefined intervals - no re-detection is performed, the
%   already measured diameter is simply sliced by time.
%
%   The cut is along FRAMES only.  idxL/idxR/diameter are [frames x nY x 3] (the
%   outer / wall-centre / luminal measures, see getMyographDiameter) and all three
%   are carried through untouched; the 3rd dimension is never sliced.
%
%   If ivT is empty BIG is returned unchanged (whole recording as one interval).
%   If NAMES is omitted or empty, intervals are auto-named interval1, interval2...
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

function intervals = cutMyographIntervals(big,ivT,names)

if isempty(ivT)
    intervals=big;                       % nothing selected -> keep whole recording
    return;
end
if nargin<3 || isempty(names)
    names=arrayfun(@(k)sprintf('interval%d',k),(1:size(ivT,1))','UniformOutput',false);
end

t=big.time(:);
hasValid = isfield(big,'valid') && ~isempty(big.valid);
if isfield(big,'measures') && ~isempty(big.measures), meas=big.measures; else, meas={'outer','mid','inner'}; end
intervals=struct('name',{},'measures',{},'time',{},'idxL',{},'idxR',{},'diameter',{},'mask',{},'valid',{});
for j=1:size(ivT,1)
    sel = t>=ivT(j,1) & t<=ivT(j,2);
    intervals(j).name     = char(names{j});
    intervals(j).measures = meas;
    intervals(j).time     = big.time(sel);
    intervals(j).idxL     = big.idxL(sel,:,:);        % keep all three measures
    intervals(j).idxR     = big.idxR(sel,:,:);
    intervals(j).diameter = big.diameter(sel,:,:);
    intervals(j).mask     = big.mask(sel,:);
    if hasValid, intervals(j).valid = big.valid(sel); else, intervals(j).valid = true(nnz(sel),1); end
end
end
