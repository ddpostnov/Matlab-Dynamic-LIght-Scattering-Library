%myographChannelSamples  One wire-myograph channel's samples, however the product holds them
%
%   [data,time] = myographChannelSamples(c) returns a channel's recorded signal and
%   its absolute times as two column vectors, whether the channel still carries the
%   WHOLE recording or only the windows that were kept of it.
%
%   A WIRE RECORDING IS STORED TWO WAYS, IN THIS ORDER.  runLabChart writes every
%   channel whole - .data and .time over the entire file - because nothing has been
%   chosen yet.  The intervals step then cuts each channel to the windows the
%   operator defined, moves the samples into .intervals(k).samples and empties .data
%   and .time: what falls outside the analysed windows is discarded, and keeping a
%   whole-recording copy beside the windows would be exactly the copy a myograph
%   product exists not to keep.  Both shapes are legitimate; which one a file is in
%   says how far through the pipeline it has come, not whether it is right.
%
%   SO EVERY READER ASKS HERE INSTEAD OF TESTING .data ITSELF.  A reader that looks
%   at .data alone sees an empty channel after the cut and concludes the recording
%   holds nothing; a reader that looks at the windows alone cannot open a file that
%   has not been through the intervals step yet.  One answer, in one place.
%
%   THE WINDOWS ARE LAID END TO END, IN STORED ORDER, AND THE GAPS BETWEEN THEM ARE
%   NOT FILLED.  The result therefore has a time base that JUMPS wherever a stretch
%   was discarded, which is what every consumer of it already expects: the editor
%   breaks its drawn line on a gap (contiguousRuns), and each window is cut back out
%   of it by TIME.  Nothing here interpolates and nothing here re-bases: a sample
%   keeps the absolute second it was recorded at.
%
%   Syntax:
%      [data,time] = myographChannelSamples(c)
%
%   INPUTS
%     c        one element of results.channel: .data .time and .intervals(k).samples.
%
%   OUTPUTS
%     data     [n x 1] double - the channel's samples, in time order.
%     time     [n x 1] double - their absolute seconds from the start of the file.
%              Both are empty for a channel that carries neither shape, which is a
%              channel of a recording that was never read rather than an error.
%
%   EXAMPLE
%     [y,t] = myographChannelSamples(results.channel(2));
%     inWindow = y(t>=iv.tStart & t<=iv.tEnd);
%
% See also: cutMyographIntervals, getMyographTrace, runLabChart, myographProduct,
%           runMyographVasomotion
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function [data,time] = myographChannelSamples(c)

data=[]; time=[];
if ~isstruct(c) || ~isscalar(c), return; end

% THE WHOLE-RECORDING COPY WINS WHILE IT IS THERE.  It is what the reader wrote and
% it has no gaps; only the intervals step empties it, and only after moving the
% samples it kept into the windows.
d=double(fieldOr(c,'data',[]));
t=double(fieldOr(c,'time',[]));
if ~isempty(d)
    % SAMPLES WITHOUT THEIR TIMES ARE NOT A SIGNAL.  Every consumer cuts by time, so
    % a channel whose two vectors disagree comes back empty rather than half-usable:
    % the caller then says the channel carries nothing, which is true of it.
    if numel(d)~=numel(t), return; end
    data=d(:); time=t(:);
    return
end

ivs=fieldOr(c,'intervals',[]);
if isempty(ivs) || ~isstruct(ivs), return; end
dp=cell(1,numel(ivs)); tp=cell(1,numel(ivs));
for k=1:1:numel(ivs)
    sm=fieldOr(ivs(k),'samples',[]);
    if ~isstruct(sm) || ~isscalar(sm), continue; end
    dk=double(fieldOr(sm,'data',[])); tk=double(fieldOr(sm,'time',[]));
    if isempty(dk) || numel(dk)~=numel(tk), continue; end
    dp{k}=dk(:); tp{k}=tk(:);
end
data=cat(1,dp{:}); time=cat(1,tp{:});
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
