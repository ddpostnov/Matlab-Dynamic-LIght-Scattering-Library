%myographDiameterBranch  One window's diameter block, built the one way it is built
%
%   d = myographDiameterBranch(t,nY,meas,lines,wallL,wallR,measured,valid,keepArrays)
%   assembles results.intervals(k).diameter: the measurement of ONE analysis window,
%   over the rows that were really measured.
%
%   THERE ARE TWO PRODUCERS AND ONE BUILDER.  runMyographDiameter writes this block
%   when a recording is measured, and cutMyographIntervals writes it again when a
%   boundary moves - and they are the same struct or they are a bug waiting for the
%   one code path that compares them.  Two copies kept in step by hand is how the
%   field ORDER drifts, so the copy lives here and both callers hand over arrays
%   they have already sliced.
%
%   WHAT IS STORED, and why the diameter is stored rather than derived.  The last
%   two lines of getMyographDiameter>postProcessWalls are 'R=max(R,L); D=R-L;', so
%   the per-line diameter IS the difference of the two post-processed wall
%   positions, exactly, in single.  It is written anyway, for three reasons in this
%   order: the explorer reads leaves and does not compute; a quantity every consumer
%   re-derives is the "two copies that can disagree" risk with the arithmetic moved
%   to the reader; and the third of the file it costs is not the size problem this
%   product was reshaped to solve.
%
%   THE ARRAYS ARE THE FILLED ONES - the numbers the trace is the mean of,
%   interpolated fill included.  NaN-ing the fill would make the pooled form of
%   .lines disagree with .<measure> (measured on the reference recording: up to
%   1.19 px, 0.3%), and the two claiming to be one quantity is the whole point.  How
%   much of the window was really measured is not lost: it is .measured, and its
%   summary is .stats.<measure>.measuredFraction.
%
%   THE TRACE IS POOLED HERE AND NOWHERE ELSE, which is what makes slicing a window
%   and re-measuring it agree: trace(i) is the mean along the vessel of frame i, so
%   cutting rows out of the trace and pooling the cut rows of the array are the same
%   arithmetic on the same numbers.  A window that was written WITHOUT its arrays
%   (see keepArrays) is narrowed by handing its stored trace in as its own single
%   'line' - the pooling of one column is that column, exactly, with no cast on the
%   way - so there is still only one builder.
%
%   Syntax:
%      d = myographDiameterBranch(t,nY,meas,lines,wallL,wallR,measured,valid)
%      d = myographDiameterBranch(t,nY,meas,lines,wallL,wallR,measured,valid,keepArrays)
%
%   INPUTS
%     t         [frames x 1] the window's own times, ABSOLUTE seconds.
%     nY        the number of measured lines.  It is size(lines{m},2) whenever the
%               arrays are there, and the count the window recorded when they are
%               not - which is why it is passed rather than measured off the data.
%     meas      1 x m cellstr of measure names, in the order the cores return them
%               ({'outer','mid','inner'}).
%     lines     1 x m cell, each [frames x lines] - every line's own diameter.  With
%               keepArrays false it may instead be the [frames x 1] trace itself.
%     wallL     1 x m cell, each [frames x lines] - the left wall position, px.
%     wallR     1 x m cell, each [frames x lines] - the right wall position, px.
%               Both are read only when keepArrays is true.
%     measured  [frames x lines] logical - 1 = detected, 0 = interpolated fill.  ONE
%               array, shared by the three measures: it is the wall-centre measure's
%               (getMyographDiameter), and the three measuredFractions on the
%               reference set are identical.
%     valid     [frames x 1] logical - false = a wall had dilated out of the field
%               of view, so that frame's diameter is a lower bound.
%     keepArrays  true (default) stores the per-point arrays; false stores the trace
%               and its statistics alone.  IT IS NOT FREE: with no source member the
%               arrays are the only copy, so a window written without them can be
%               compared as a trace and can never be analysed for propagation or per
%               line.
%
%   OUTPUT
%     d        the diameter block, fields in this order:
%                .time     [frames x 1] absolute seconds
%                .nY       number of measured lines
%                .wallL.<measure>   [frames x lines] single, px
%                .wallR.<measure>   [frames x lines] single, px
%                .lines.<measure>   [frames x lines] single - the measurement
%                .measured [frames x lines] logical
%                .valid    [frames x 1]     logical
%                .stats.<measure>   mean std min max amplitude measuredFraction
%                                   validFraction
%                .<measure>         [frames x 1] the line-averaged trace
%              With keepArrays false, .wallL .wallR .lines and .measured are absent.
%
%   EXAMPLE
%     d = myographDiameterBranch(iv.time,numel(r),iv.measures, ...
%             {L1,L2,L3},{WL1,WL2,WL3},{WR1,WR2,WR3},iv.mask(:,r),iv.valid,true);
%
% See also: runMyographDiameter, cutMyographIntervals, getMyographDiameter,
%           myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function d = myographDiameterBranch(t,nY,meas,lines,wallL,wallR,measured,valid,keepArrays)

if nargin<9 || isempty(keepArrays), keepArrays=true; end
keepArrays=logical(keepArrays(1));
meas=reshape(cellstr(meas),1,[]);

% THE FIELD ORDER IS DECLARED HERE, ONCE, and the per-measure traces are appended
% after it by the loop - so a block built by either producer is the same struct and
% not merely an equivalent one.
d=struct('time',[],'nY',0,'wallL',struct(),'wallR',struct(),'lines',struct(), ...
    'measured',[],'valid',[],'stats',struct());

d.time=double(t(:));
d.nY=double(nY);
d.measured=logical(measured);
d.valid=logical(valid(:));

for m=1:1:numel(meas)
    X=lines{m};
    % POOL BEFORE CASTING.  The cores measure in single, so single(X) is a copy and
    % not a decision; a narrowed window hands its stored DOUBLE trace in here, and
    % casting first would round it for no reason.
    trace=mean(double(X),2,'omitnan');
    if keepArrays
        d.lines.(meas{m})=single(X);
        d.wallL.(meas{m})=single(wallL{m});
        d.wallR.(meas{m})=single(wallR{m});
    end
    d.(meas{m})=trace;
    d.stats.(meas{m})=traceStats(trace,d.valid,d.measured);
end

% A run that did not keep them leaves no empty container behind: an ABSENT .lines is
% what tells a reader - and exScan - that the trace is the only form of the quantity
% in this file, which an empty struct sitting there would not.
if ~keepArrays
    d=rmfield(d,{'wallL','wallR','lines','measured'});
end
end

% =====================================================================
function st=traceStats(trace,valid,measured)
%traceStats  What a protocol compares: the level, the swing and how much of it was
%   actually seen.  OFF-FOV FRAMES ARE EXCLUDED from the level and the swing - a
%   wall that dilated out of view leaves an edge-clamped lower bound behind, and
%   averaging it in would report a constriction that never happened.  How many such
%   frames there were is not hidden: it is validFraction.
if isempty(valid) || numel(valid)~=numel(trace), valid=true(size(trace)); end
v=trace(valid); v=v(isfinite(v));
st=struct('mean',NaN,'std',NaN,'min',NaN,'max',NaN,'amplitude',NaN, ...
    'measuredFraction',NaN,'validFraction',NaN);
if ~isempty(v)
    st.mean=mean(v);  st.std=std(v);
    st.min=min(v);    st.max=max(v);
    st.amplitude=st.max-st.min;
end
if ~isempty(measured), st.measuredFraction=mean(double(measured(:))); end
if ~isempty(trace),    st.validFraction=mean(double(valid)); end
end
