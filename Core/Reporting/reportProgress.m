%reportProgress - One rate-limited progress tick, or a parfor-safe queue that emits them.
%
%   A progress indicator only has to be monotone and roughly right, so it is
%   throttled at the SEND and by TWO independent bounds: a tick costs a line only
%   when it is at least ~1 % further on than the last one AND at least one second
%   later.  A 30 s loop then prints about thirty times and a two-hour loop about a
%   hundred - both bounded, neither noisy.  The 100 % tick is always emitted, so a
%   bar that has been throttled all the way still finishes full.
%
%   ONE TRANSPORT, TWO RENDERINGS.  The command window rewrites a single line in
%   place (the \r idiom registerRetinaLSCI already uses); the host hook gets the
%   fraction and the label and updates its own label and cell percentage.  A
%   uitextarea cannot do \r, and appending a tick per percent would put a hundred
%   lines per file into a four-hundred-line log - so progress NEVER reaches the
%   host log, only its progress indicator.
%
%   THE QUEUE FORM exists because a parfor body cannot print, and it OWNS THE SEND
%   THROTTLE so no caller has to invent one.  It returns a parallel.pool.DataQueue
%   plus the batch size to send at: workers report how many items they have
%   finished, the client accumulates and ticks.  Throttling only the PRINTING is
%   not enough - every send costs a client-side afterEach callback whether or not
%   it produces a line, and one send per pixel in a 10^6-iteration parfor is the
%   one place in this library where reporting measurably impedes processing.
%   sendEvery is total/200, so a loop of any size costs about two hundred
%   callbacks.  The count is then approximate and monotone, which is all a
%   progress indicator needs; emit a final reportProgress(rep,1,label) after the
%   loop so the bar always reaches 100 %.
%
% Syntax:
%    reportProgress(rep, frac)
%    reportProgress(rep, frac, label)
%    [dq,sendEvery] = reportProgress(rep, 'queue', total)
%    [dq,sendEvery] = reportProgress(rep, 'queue', total, label)
%
% Inputs:
%    rep   - the context from reportOpen.
%    frac  - progress in [0 1]; values outside are clamped, non-finite ignored.
%    label - (optional) what is being counted ('per-pixel vasomotion').  Defaults
%            to the call's procedure label.
%    total - queue form: how many items make 100 %.
%
% Outputs:
%    dq       - queue form only: the DataQueue to send counts to.  Empty otherwise.
%    sendEvery- queue form only: send once per this many items (>=1).
%
% Example:
%    reportProgress(rep, fidx/numel(fNames));
%
%    [dq,sendEvery] = reportProgress(rep,'queue',nPix,'per-pixel vasomotion');
%    parfor p = 1:nPix
%        ...
%        if mod(p,sendEvery)==0, send(dq,sendEvery); end
%    end
%    reportProgress(rep,1,'per-pixel vasomotion');   % forced final tick
%
% Dependencies: Parallel Computing Toolbox (queue form only).
% See also: reportOpen, reportStage, reportFile
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function [dq, sendEvery] = reportProgress(rep, frac, varargin)

dq        = [];
sendEvery = 1;

% ---- queue form -----------------------------------------------------------
if (ischar(frac)||isstring(frac)) && strcmpi(char(frac),'queue')
    total = 1; label = '';
    if ~isempty(varargin) && ~isempty(varargin{1}), total = double(varargin{1}); end
    if numel(varargin)>1, label = char(string(varargin{2})); end
    if ~isscalar(total) || ~isfinite(total) || total<=0, total = 1; end
    sendEvery = max(1, floor(total/200));   % ~200 client callbacks, whatever the size
    cnt = 0;
    dq  = parallel.pool.DataQueue;
    afterEach(dq, @accumulate);
    return
end

% ---- single tick ----------------------------------------------------------
label = '';
if ~isempty(varargin), label = char(string(varargin{1})); end
tick(rep, double(frac), label);

    function accumulate(n)
    %accumulate  One worker report: n more items done.  The count is approximate
    %   and monotone, which is all a progress indicator needs.
        if nargin<1 || isempty(n), n = 1; end
        cnt = cnt + double(n);
        tick(rep, min(1,cnt/total), label);
    end
end

% =====================================================================
function tick(rep, frac, label)
%tick  Spend a line only if BOTH bounds are cleared - or if this is the 100 % tick
%   and 100 % has not been reported yet, which is always worth a line.
if ~isscalar(frac) || ~isfinite(frac), return; end
frac = max(0, min(1, frac));

st   = rep.state;
tNow = toc(st('t0'));
force = frac>=1 && st('lastTickF')<1;
if ~force && (tNow - st('lastTickT') < 1 || frac - st('lastTickF') < 0.01)
    return
end
st('lastTickT') = tNow;
st('lastTickF') = frac;
st('nTicks')    = st('nTicks') + 1;              % st is a handle Map, not a copy

rep.progressFcn(frac, label);                    % host sink: label + cell percentage

% A hosted run needs the event queue flushed or the label and the cell percentage
% the hook just set are not painted until the wrapper returns.  It is inside the
% throttle, so it costs at most one flush per tick, and 'limitrate' lets MATLAB
% skip it when it is already keeping up.
if rep.hosted, drawnow limitrate; end

if strcmp(rep.level,'quiet'), return; end        % command window: ONE line, in place
if isempty(label), label = rep.procLabel; end

% Both fields only ever grow, so \r alone always covers the previous line and no
% padding is needed.  The line is left unterminated until 100 % so the next tick
% can overwrite it; reportOpen's emit closes it if anything else prints first.
fprintf('\r  %s: %d%%  [%.0fs]', label, floor(100*frac), tNow);
if frac>=1
    fprintf('\n'); st('openTick') = false; %#ok<NASGU>  Map is a handle
else
    st('openTick') = true;                 %#ok<NASGU>  Map is a handle
end
end
