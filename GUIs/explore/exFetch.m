%exFetch - The numbers behind one plot, in a shape no renderer has to interpret.
%
%   payload = exFetch(results, d, plotId, opts) resolves one descriptor into the
%   arrays a renderer draws.  Four jobs, in this order: find the leaf, decide which
%   UNITS of it the user asked for, reduce over those units, and slice away the
%   dimensions the chosen plot does not draw.  What comes back is the same small
%   struct whatever container the numbers came from - which is the point: a renderer
%   never has to know whether it is looking at a per-segment time series, a
%   per-pixel map or a myograph window.
%
%   items = exFetch('where', results, d) is the selection list for that descriptor,
%   BUILT FROM WHAT THE FILE ACTUALLY CONTAINS.  On the reference contrast set
%   sMetrics.category is [1 2 3 5] - parenchyma, unsegmented, outer wall, lumen -
%   with no category 4, so no inner-wall entry is offered.  Offering one would be
%   offering an empty selection, and that is the difference between a list built
%   from a constant and a list built from the data.
%
%   selectRows, weightCol AND tableDomain CAME OUT OF guiExplore and are now the
%   only copies.  They encode category rules that were arrived at against real
%   recordings and cannot be re-derived from first principles: lumen is category 5
%   and every vessel selection is restricted to it; parenchyma is category 1; a
%   label that exists only on parenchyma segments falls back to category 1 rather
%   than returning nothing; a row whose idx is not finite is an empty table slot.
%   THE ONE RULE TO KEEP IN MIND WHEN EDITING selectRows: its otherwise-arm falls
%   back to "all lumen", so a selection this module offers and forgets to name there
%   becomes a plot of every vessel - the outer wall silently drawn as the whole
%   vasculature.  Every entry the contents-driven list can produce is named in that
%   switch, the four category ones included.
%
%   A SELECTION MAY BE SEVERAL (H).  opts.select takes a list, and the row masks are
%   ORed.  For a curve that would be wrong - several selections compete for one axes
%   and each is drawn as its own series - but a PICTURE combines them spatially:
%   arteries, veins and parenchyma together are the whole map, and painting only the
%   first draws a fragment and calls it the map.  The union makes the otherwise-arm
%   more dangerous, not less, because an unnamed selection ORed into the mask would
%   quietly add every vessel to a picture that would still look plausible - so
%   exFetch('named',sel) exposes the question "does the switch name this one?" and
%   testExplorePlan asks it of every entry the Where list can produce.
%
%   getNested IS LIFTED AND GENERALISED, in two ways a descriptor path needs and a
%   guiExplore path never did: a component of the form name(<digits>) takes that
%   element of a struct array, so intervals(2).diameter.outer resolves; and a final
%   component that names a table VARIABLE returns the column, so sMetrics.CLR does.
%   A component is an index only when it is name(<digits>), which is why the metrics
%   column 'std(diameter)' can never be mistaken for one.
%
%   THE AREA-WEIGHTED MEAN IS THE DEFAULT AND ITS WEIGHTS ARE REAL.  A pooled
%   reduction over segments weights by area (or by length when there is no area
%   column), never by ones - a 4000-pixel artery and a 40-pixel capillary do not
%   contribute equally to "the mean over arteries".  The one place ones is correct
%   is pixels, which are equal by construction, and the one place it is unavoidable
%   is a recording that carries no metrics table at all - and there the payload says
%   so in .note rather than pretending.
%
%   WHICH IS WHY A PIXEL AND A SEGMENT DO NOT GET THE SAME UNIT RECORD.  A segment is
%   a named thing - it has a label, a type, an area and a row in a metrics table, and
%   a reader has to see which one a point belongs to - while a pixel is a position in
%   an image that no renderer ever names and that weighs exactly what every other
%   pixel weighs; so the record for a pixel leaf carries the selection and nothing
%   else, and its names, weights and ids are built only if some plot actually asks
%   for them.  On the reference maps that is the difference between 1 340 416
%   sprintf calls per render and none: the whole image is .all rather than a subscript
%   vector of every index, and both the leaf copy and the NaN canvas go with it.
%
%   image.fromSegments IS THE INVERSE OF A ROW SELECTION and belongs here rather
%   than in a renderer.  Each selected unit's own number is painted into results.sMap
%   (dvsMap for tracked vessels) through the segment's idx, and every pixel that
%   belongs to no selected segment - or to one whose value is not finite - is left
%   NaN so it renders transparent.  It is the only way a per-segment scalar reaches
%   an image, and the only reason a segment index never has to become a coordinate.
%
%   AND video.fromSegments IS THAT SAME OPERATION ONCE PER SAMPLE (D9) - the video an
%   LSCI recording actually has, against a mask stack, which is what Kind = Video used
%   to offer.  Two things decide how it is built:
%
%     THE CUBE IS NEVER MATERIALISED.  2448 x 1496 x 896 doubles is 26 TB, so the
%     payload carries .paint - the [nUnit x nFrame] matrix and the look-up table that
%     is the SAME for every frame - and exFetch('frame',payload,k) paints frame k on
%     demand.  The renderer paints the frame the slider is on; the export paints one
%     frame at a time as it writes.  Neither ever holds more than two pictures.
%
%     A LONG RECORDING IS WALKED AT A STRIDE.  gsData is 61 208 samples, and a slider
%     with 61 208 positions is a defect however correct each frame is, so the slider
%     is capped at scrubCap() evenly spaced positions.  Each is a REAL sample - not an
%     average of a window - labelled in seconds off the leaf's own time base, and the
%     stride is applied on the way IN, as a subscript into the leaf, so the 725 MB
%     that gsData would otherwise become is never built.  scrubSpecFor and the payload
%     call the same scrubPositions, so the slider and the picture cannot disagree.
%
%   THE PER-PIXEL PULSE CYCLE IS NOT RECONSTRUCTED AND IS NOT INVENTED.  pulsatility.
%   ppx.timeVectors is a runPulsatility output that the run behind the reference set
%   did not write; the ppx group there holds scalars only.  If a per-pixel pulse cycle
%   is wanted that is a wrapper decision, not an explorer one.
%
%   THE TREE IS ALWAYS LOADED (D7).  A file past a size limit used to arrive as a
%   handle whose leaves were read in HDF5 slabs, and the price was that everything
%   BESIDE the leaf - the metrics table above all - could not be read at all: a
%   table is an opaque compound to HDF5, so a 3.63 GB recording offered no vessel
%   types, no labels and no weights, and said every unit counted equally.  There is
%   no partial route to a table, so the file is loaded instead and the slab reader
%   is gone.  Everything here now asks a struct, once.
%
%   A LEAF'S COMPANIONS COME BACK WITH IT.  exSchema declares that some leaves ride
%   with others in a named role - speedCI is an INTERVAL on speed, slope is the FIT
%   of lagByRow, confidenceText is its CAPTION - and those are resolved here, into
%   .companions, rather than by a renderer reaching into a results tree for them.
%   The payload already has .note for a caveat, and a companion is not a caveat: it
%   is part of the answer, so it gets a place of its own.  A companion the file does
%   not carry simply does not come back - a window whose propagation failed carries
%   no confidence sentence, and borrowing another window's would be worse than
%   drawing none.  NOTHING IS RECOMPUTED: a caption is a finished sentence the
%   producer wrote and is passed through as it stands.
%
%   THE SLIDER'S AXIS IS ANSWERED HERE TOO.  exFetch('scrub',...) reports which
%   dimension a scrubbed plot leaves for a slider, its label, its coordinate values
%   and its length.  The window could not work that last one out for itself: the
%   length of a dimension is a fact about the leaf and the axis registry together,
%   and size(...,3) alone does not say which coordinate that dimension is.
%
%   Syntax:
%      payload = exFetch(results, d, plotId)
%      payload = exFetch(results, d, plotId, opts)
%      items   = exFetch('where', results, d)
%      spec    = exFetch('scrub', results, d, plotId, opts)
%      key     = exFetch('slicekey', dimName)
%      img     = exFetch('frame', payload, k)
%      tf      = exFetch('named', selection)
%
%   INPUTS
%     results  the RESULTS member of a triplet, loaded.
%     d        one descriptor from exScan.
%     plotId   one of the ids exPlotRules('ids') lists.
%     payload  a payload this function returned, for the 'frame' action.
%     k        which position along the scrub axis, for the 'frame' action.
%     opts     struct, every field optional:
%                .select    the Where entry, e.g. 'Arteries (lumen)', 'Label: MCA',
%                           '(whole image)', 'Outer wall' - or a CELLSTR of several,
%                           whose row masks are ORed.  Default: everything.
%                .units     'pooled' | 'points' | 'curves'.  Default: whatever the
%                           plot asks for - points for a box, pooled for the rest.
%                .rows      [lo hi] range of positions along a vessel, for line units
%                .interval  which analysed window is wanted.  The descriptor's path
%                           already pins one, so this is a GUARD: a mismatch returns
%                           an empty payload rather than quietly drawing the window
%                           the cascade has moved away from.
%                .slice     struct with any of .f .pct .harmonic .frame - the index
%                           along a dimension the plot does not draw.  Default 1.
%                .path      the _r.mat this tree came from, so the percentile LEVELS
%                           can be read from the sibling _s.mat.  That one sibling
%                           read is the only file this tool opens beside the _r -
%                           see the exAxes header for why it is not slicing.
%                .scrub     which dimension a scrubbed plot walks, when the leaf
%                           offers more than one and the user has chosen.  Ignored
%                           when the leaf does not have it, so a stale choice cannot
%                           empty a picture.
%
%   OUTPUTS
%     payload  struct:
%                .ok      false when there is nothing to draw; .note says why
%                .note    a caveat worth showing in the status bar, '' when none
%                .plot .kind
%                .x       [nx x 1] the x coordinates ([] for a box or an image)
%                .Y       [nx x nSeries] curve values, or [nUnit x 1] point values
%                .w       the weights behind them
%                .names   {1 x nSeries} what each series or point is
%                .img     [Y x X] a picture, NaN where transparent
%                .cdata   [ny x nx] a heatmap - rows are y, columns are x
%                .yvals   [ny x 1] the y coordinates of a heatmap
%                .frames  [Y x X x nFrame] a video cube
%                .paint   the kit a RECONSTRUCTED video is painted from, in place of
%                         a cube: .kit (the map's look-up, the same every frame) and
%                         .S [nUnit x nFrame].  Empty for every other plot.  Ask
%                         exFetch('frame',payload,k) rather than reading it
%                .fvals   [nFrame x 1] where each frame sits on the scrub axis
%                .xlab .ylab .clab .flab   the axis titles
%                .xlog    true when x is logarithmic
%                .xscale  'linear' | 'log' | 'ordinal' | 'categorical'
%                .nUnits  how many units survived the selection
%                .companions  struct array .leaf .role .value .label .unit - the
%                         leaves that ride WITH this one rather than being drawn as
%                         variables of their own.  .role is what a renderer switches
%                         on: 'interval' (a [lo hi] pair -> an error bar), 'fit' (a
%                         slope -> the line through the scatter), 'caption' (a
%                         finished sentence -> shown beside the plot).  Empty when
%                         the leaf has none or the file does not carry them.  The
%                         fifth role, 'pooled', never appears here: a pooled form is
%                         what opts.units already produces from the host itself
%     spec     struct: .name .label .values .n .choices - the dimension a slider
%              walks for this plot, empty .name when the plot has none.  For a
%              reconstructed video .n is the number of SLIDER positions, which is
%              fewer than the samples on a long recording, and .values are those
%              samples' own coordinates
%     key      the opts.slice field name for one dimension: 'f' | 'pct' |
%              'harmonic' | 'frame'
%     img      [Y x X] one frame of a reconstructed video, NaN where transparent
%     tf       true when selectRows NAMES that selection rather than falling back
%
%   EXAMPLE
%      S = load('LSCI_20240809_1ADCF08BP_c_BFI_r.mat');
%      D = exScan(S.results);
%      d = D(strcmp({D.path},'pulsatility.sData.scalars.psPI'));
%      p = exFetch(S.results, d, 'box', struct('select','Arteries (lumen)'));
%      numel(p.Y)                      % one pulsatility index per artery
%
%   DEPENDS ON
%     exAxes, exPlotRules (GUIs/explore).
%
% See also: exPlotRules, exScan, exAxes, exSchema, exModality, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function out = exFetch(varargin)

a1 = varargin{1};
if ischar(a1) || isstring(a1)
    switch lower(char(a1))
        case 'where',    out = whereList(varargin{2:end});
        case 'scrub',    out = scrubSpecFor(varargin{2:end});
        case 'slicekey', out = sliceKey(varargin{2});
        case 'frame',    out = frameOf(varargin{2:end});
        case 'named',    out = isNamedSelection(varargin{2});
        otherwise
            error('exFetch:action','Unknown action ''%s''.', char(a1));
    end
    return
end
if nargin<3, error('exFetch:args','results, a descriptor and a plot id are required.'); end
if nargin<4, varargin{4} = struct(); end
out = fetchOne(varargin{1:4});
end

% ================================================================== THE FETCH

function P = fetchOne(R, d, plotId, opts)
plotId = char(plotId);
P = blankPayload(plotId);
opts = fillOpts(opts);

% ---- 0. the guards -------------------------------------------------------------
if ~isempty(d.pairedWith)
    P.note = sprintf('%s is the %s of %s and is drawn with it, not on its own.', ...
        riderName(d), roleNoun(d), d.pairedWith); return
end
if ~any(strcmp(exPlotRules(d), plotId))
    P.note = sprintf('%s cannot be drawn as %s.', d.leaf, plotId); return
end
if ~isempty(opts.interval) && ~intervalMatches(d.path, opts.interval)
    P.note = 'This variable belongs to a different analysed window.'; return
end

A  = exAxes(R, scopeOf(d.path), opts.path);
ax = exPlotRules('axes', d, plotId, A, opts.scrub);

P.kind   = ax.kind;
P.xlab   = ax.x.label;
P.ylab   = ax.y.label;
P.clab   = ax.c.label;
P.xlog   = ax.x.log;
P.xscale = ax.x.scale;

% ---- 0b. what rides with this leaf ----------------------------------------------
% Resolved before the numbers, because a companion is a fact about the DESCRIPTOR
% and the tree, not about the plot: the same interval belongs to a speed whether it
% is drawn as a box or as a bar.
P.companions = companionsFor(R, d);

% ---- 1. the leaf, sliced on the way in when the file is lazy --------------------
sl = sliceIndices(d, ax, opts, A);
[V, why] = readLeaf(R, d, sl);
if isempty(V), P.note = why; return; end

if strcmp(d.layout,'xyPairs')
    P = xyPayload(P, V, ax);
    P.ok = hasContent(P);
    return
end
M = canonical(V, d);                       % [nUnit x n1 x n2], sliced dims singleton

% ---- 2. the units ---------------------------------------------------------------
U = unitsFor(R, d, opts, size(M,1));
P.note = U.note;
if U.n==0
    P.note = trim([U.note ' Nothing is selected here.']); return
end
% A selection of everything has no subscript vector, so there is nothing to index
% by: M(1:end,:,:) is a copy of M made to reproduce M, and on a 1496x896 map that
% copy is the second largest thing this function does.
if ~U.all, M = M(U.idx,:,:); end
P.nUnits = U.n;

% ---- 3. reduce over units, and 4. over the dims this plot does not draw ---------
red = reduceMode(opts.units, ax);
switch plotId
    case {'box','bar'},                P = scalarPayload(P, M, U, red);
    case {'curve.time','curve.f','curve.pct','curve.harmonic','curve.position'}
                                       P = curvePayload(P, M, U, red, d, ax);
    case 'curves.family',              P = familyPayload(P, M, U, d, ax);
    case 'image',                      P = imagePayload(P, M, U, d);
    case 'image.fromSegments',         P = paintPayload(P, M, U, R, d);
    case {'image.frame','image.timeAverage'}
                                       P = imagePayload(P, reduceDims(M,d,ax), U, d);
    case {'heat.ft','heat.fpct'},      P = heatPayload(P, M, U, d, ax);
    case 'heat.positionTime',          P = kymoPayload(P, M, U, d, ax);
    case 'video',                      P = videoPayload(P, M, U, d, ax);
    case 'video.fromSegments',         P = paintVideoPayload(P, M, U, R, d, ax, A);
    otherwise
        P.note = sprintf('No fetch is defined for %s.', plotId); return
end
% A payload that could not build its array said why in .note; it must not come back
% claiming to be drawable just because it reached the end of the switch.
P.ok = hasContent(P);
end

function tf = hasContent(P)
% .paint counts: a reconstructed video has no cube by design, and a payload that can
% draw every frame of one is not empty just because it is not holding 26 TB.
tf = ~isempty(P.Y) || ~isempty(P.img) || ~isempty(P.cdata) || ~isempty(P.frames) ...
     || ~isempty(P.paint);
end

% ================================================================== THE SCRUB AXIS

function s = scrubSpecFor(R, d, plotId, opts)
%scrubSpecFor  WHAT A SLIDER UNDER THE AXES HAS TO WALK: the dimension, its label,
%   its coordinate values in the recording's own units, and how many positions it
%   has.  It is answered here rather than in the window because the LENGTH of a
%   dimension is a fact about the leaf and the axis registry together - a window that
%   worked it out from size(...,3) would be guessing, and would guess wrong on a lazy
%   file it has never read.  Nothing is loaded: only the descriptor and the registry
%   are consulted.
s = struct('name','','label','','values',[],'n',0,'choices',{{}},'choiceLabels',{{}});
if nargin<4, opts = struct(); end
opts = fillOpts(opts);
A  = exAxes(R, scopeOf(d.path), opts.path);
ax = exPlotRules('axes', d, plotId, A, opts.scrub);
if isempty(ax.scrub.name), return; end
n = axisLength(A, ax.scrub.name, d);
if ~isfinite(n) || n<1, n = 0; end
% A PLOT MAY WALK ITS AXIS AT A STRIDE.  The positions are computed by the same
% function the payload uses, so a slider position and the frame it paints cannot
% describe different samples - and the VALUES are the recording's own coordinates at
% those positions, never a renumbering of them.
pos = scrubPositions(n, plotId);
v   = axisValues(ax.scrub, n);
s.name    = ax.scrub.name;
s.label   = ax.scrub.label;
s.n       = numel(pos);
s.values  = v(pos);
s.choices = ax.scrubChoices;
% The candidates are named by asking the rules for each of them in turn, so the words
% beside a choice are the words on the axis if it is taken.
s.choiceLabels = cell(1,numel(s.choices));
for i = 1:numel(s.choices)
    q = exPlotRules('axes', d, plotId, A, s.choices{i});
    s.choiceLabels{i} = q.scrub.label;
end
end

% ================================================================== PER PLOT

function P = scalarPayload(P, M, U, red)
%scalarPayload  One number per unit, or one number for the lot.  A box wants the
%   distribution, a bar wants the pooled value - which is why the two ids exist.
v = double(M(:,1,1));
w = unitWeights(U);
if strcmp(red,'pooled')
    P.Y = wmean(v, w, 1);
    P.w = sum(w);
    P.names = {''};
else
    P.Y = v(:);
    P.w = w;
    P.names = unitNames(U);
end
end

% =====================================================================
function P = curvePayload(P, M, U, red, d, ax)
%curvePayload  One curve, pooled over the units, or one curve per unit.  The x is
%   the recording's own coordinate vector whenever the registry has it; a bare index
%   is the fallback and never the first choice.
M = reduceDims(M, d, ax);
k = dimPos(d, ax.x.name);
if k==0, k = 1; end                                   % the unit is the x (a line)
w = unitWeights(U);
if isempty(ax.x.name) || strcmp(ax.unitAxis,'line')
    Y = double(M(:,1,1));                             % one value per position
    P.Y = Y(:);
    P.x = axisValues(ax.x, numel(Y));
    P.w = w;
    P.names = {''};
    return
end
Y = squeezeTo2D(M, k);                                % [nUnit x nx]
n = size(Y,2);
if strcmp(red,'pooled')
    P.Y = reshape(wmean(Y, w, 1), n, 1);
    P.w = sum(w);
    P.names = {''};
else
    P.Y = Y.';                                        % [nx x nUnit]
    P.w = w;
    P.names = unitNames(U);
end
P.x = axisValues(ax.x, n);
end

% =====================================================================
function P = familyPayload(P, M, U, d, ax)
%familyPayload  Several curves on one axes.  The family is the percentile bins when
%   the leaf has them - one curve per bin of ampMeanPct - and otherwise the unit,
%   which is only meaningful because a position along a vessel is an ordering.
M = reduceDims(M, d, ax);
kx = dimPos(d, ax.x.name);
w  = unitWeights(U);
if isempty(ax.family.name) || strcmp(ax.family.name,'line') || dimPos(d,ax.family.name)==0
    Y = squeezeTo2D(M, kx);                           % [nUnit x nx]
    P.Y = Y.';
    P.x = axisValues(ax.x, size(Y,2));
    P.w = w;
    P.names = unitNames(U);
    return
end
kf = dimPos(d, ax.family.name);
Yu = wmean(M, w, 1);                                  % [1 x n1 x n2] pooled
Y  = permute(Yu, [kx+1, kf+1, 1]);                    % [nx x nFam]
Y  = reshape(Y, size(Y,1), size(Y,2));
P.Y = double(Y);
P.x = axisValues(ax.x, size(Y,1));
P.w = sum(w);
fv  = axisValues(ax.family, size(Y,2));
P.names = arrayfun(@(z) sprintf('%g%%', z), fv(:).', 'UniformOutput', false);
end

% =====================================================================
function P = imagePayload(P, M, U, d)
%imagePayload  A picture.  Every pixel outside the selection is NaN, which is what
%   the renderer turns into transparency rather than into a hole of zeros.  When the
%   selection IS the whole image there is no outside, so the values are the picture
%   and neither the canvas nor the scatter has to be built.
sz = imageSize(d);
if isempty(sz), P.note = trim([P.note ' This recording declares no image size.']); return; end
v = double(M(:,1,1));
if U.all && numel(v)==prod(sz)
    P.img = reshape(v, sz);
    return
end
img = nan(sz);
img(unitIdx(U)) = v;
P.img = img;
end

% =====================================================================
function P = paintPayload(P, M, U, R, d)
%paintPayload  image.fromSegments: each selected unit's own number, painted back
%   into the map by its idx.
[K, why] = paintKit(R, U, d);
if isempty(K), P.note = trim([P.note ' ' why]); return; end
P.img = paintFrame(K, M(:,1,1));
end

% =====================================================================
function P = paintVideoPayload(P, M, U, R, d, ax, A)
%paintVideoPayload  video.fromSegments: the same painting, once per sample - and NOT
%   into a cube.  [nUnit x nFrame] and one look-up table is all a renderer needs to
%   make any frame, which is the difference between 5 MB and 26 TB; the frame itself
%   comes from exFetch('frame',...) when something actually asks to see it.
k = dimPos(d, ax.scrub.name);
if k==0, P.note = trim([P.note ' There is nothing to scrub along.']); return; end
[K, why] = paintKit(R, U, d);
if isempty(K), P.note = trim([P.note ' ' why]); return; end
S = squeezeTo2D(M, k);                                % [nUnit x nFrame], decimated
if isempty(S) || size(S,2)<1
    P.note = trim([P.note ' This leaf has no samples to walk.']); return
end
% The leaf arrived already cut to the slider's positions - sliceIndices puts the
% stride in as a subscript - so the coordinates are read at those same positions,
% off the FULL axis.  A length nothing could measure falls back to the frames the
% leaf turned out to have, rather than asking the registry for 1:Inf.
nF  = size(S,2);
n   = axisLength(A, ax.scrub.name, d);
P.fvals = axisValues(ax.scrub, nF);
if isfinite(n) && n>=nF
    pos = scrubPositions(n, ax.plot);
    v   = axisValues(ax.scrub, n);        % exactly n long, so pos indexes it safely
    if numel(pos)>=nF, P.fvals = v(pos(1:nF)); end
end
P.paint = struct('kit',K, 'S',double(S));
P.flab  = ax.scrub.label;
end

% =====================================================================
function [K, why] = paintKit(R, U, d)
%paintKit  EVERYTHING ABOUT PAINTING THAT DOES NOT CHANGE FROM FRAME TO FRAME: which
%   map, which pixels of it address a selected unit, and where in the look-up table
%   each unit writes.  Built once and reused, which is the whole reason a
%   reconstructed video costs a frame rather than a cube.  The look-up is by index
%   rather than by a comparison per segment - 1480 segments against 300000 pixels is
%   not a loop worth writing.
K = []; why = '';
mapName = 'sMap';
if strcmp(d.unit,'dvs'), mapName = 'dvsMap'; end
map = getNested(R, mapName);
if isempty(map) || ~isnumeric(map)
    why = ['This recording has no ' mapName ' to paint into.']; return
end
ids = unitIds(U);
if isempty(ids) || ~all(isfinite(ids))
    why = 'The units of this recording carry no index to paint by.'; return
end
m   = double(map);
top = max(ids);
ok  = isfinite(m) & m==round(m) & m>=0 & m<=top;
K = struct('sz',size(m), 'ok',ok, 'code',m(ok)+1, 'slot',round(ids(:))+1, 'top',top);
end

function img = paintFrame(K, v)
%paintFrame  One number per unit into one picture.  Every pixel that belongs to no
%   selected unit - or to one whose value is not finite - is left NaN, so it renders
%   transparent rather than as a hole of zeros.
lut = nan(K.top+1,1);
lut(K.slot) = double(v(:));
img = nan(K.sz);
img(K.ok) = lut(K.code);
end

function img = frameOf(P, k)
%frameOf  ONE FRAME OF A RECONSTRUCTED VIDEO, painted when it is asked for.  The
%   renderer asks for the frame under the slider and the export asks for each in
%   turn, so neither ever holds more than two pictures at a time.
img = [];
if ~isstruct(P) || ~isfield(P,'paint') || isempty(P.paint), return; end
k = round(double(k));
n = size(P.paint.S,2);
if ~isfinite(k) || k<1, k = 1; end
if k>n, k = n; end
img = paintFrame(P.paint.kit, P.paint.S(:,k));
end

% =====================================================================
function n = scrubCap(plotId)
%scrubCap  HOW MANY SLIDER POSITIONS A PLOT MAY OFFER.  A cube that is already in
%   memory is walked whole - its length is bounded by the fact that it fits.  A
%   RECONSTRUCTED video is not: gsData is 61 208 samples, and dragging a slider
%   through 61 208 positions is a defect however correct each frame is, so it is
%   walked at a stride.  400 is a slider a hand can land on and a video of about
%   fifty seconds at the default frame rate.
switch char(plotId)
    case 'video.fromSegments', n = 400;
    otherwise,                 n = Inf;
end
end

function pos = scrubPositions(nTotal, plotId)
%scrubPositions  Which samples the slider stops on.  Evenly spaced REAL samples -
%   nothing is averaged into a frame and no frame is between two samples - so the
%   coordinate under the slider is a time the recording actually has.
nTotal = floor(double(nTotal));
if ~isfinite(nTotal) || nTotal<1, pos = zeros(0,1); return; end
cap = scrubCap(plotId);
if ~isfinite(cap) || nTotal<=cap, pos = (1:nTotal).'; return; end
pos = (1:ceil(nTotal/cap):nTotal).';
end

% =====================================================================
function P = heatPayload(P, M, U, d, ax)
%heatPayload  Rows are the y axis, columns the x - frequency up the side and time or
%   percentile along the bottom, which is the orientation guiExplore already draws.
kx = dimPos(d, ax.x.name); ky = dimPos(d, ax.y.name);
if kx==0 || ky==0, P.note = trim([P.note ' This leaf has no two axes to spread out.']); return; end
w  = unitWeights(U);
Yu = wmean(M, w, 1);
C = permute(Yu, [ky+1, kx+1, 1]);
P.cdata = double(reshape(C, size(C,1), size(C,2)));
P.x     = axisValues(ax.x, size(P.cdata,2));
P.yvals = axisValues(ax.y, size(P.cdata,1));
P.w     = sum(w);
end

% =====================================================================
function P = kymoPayload(P, M, U, d, ax)
%kymoPayload  The position-time map.  The UNIT is the y axis here, so nothing is
%   reduced over units - this is the one plot where a line's identity survives.
M  = reduceDims(M, d, ax);
kx = dimPos(d, ax.x.name);
C  = squeezeTo2D(M, kx);                              % [nUnit x nx]
P.cdata = double(C);
P.x     = axisValues(ax.x, size(C,2));
P.yvals = (1:size(C,1)).';
if ~isempty(ax.y.values), P.yvals = clipTo(ax.y.values, size(C,1)); end
P.w     = unitWeights(U);
end

% =====================================================================
function P = videoPayload(P, M, U, d, ax)
%videoPayload  A cube whose third dimension is the one the slider walks.  A whole
%   image needs no canvas and no scatter per frame: the values already ARE the cube,
%   in the order the frames want them.
sz = imageSize(d);
if isempty(sz), P.note = trim([P.note ' This recording declares no image size.']); return; end
k = dimPos(d, ax.scrub.name);
if k==0, P.note = trim([P.note ' There is nothing to scrub along.']); return; end
S = squeezeTo2D(M, k);                                % [nUnit x nFrame]
nF = size(S,2);
if U.all && size(S,1)==prod(sz)
    V = reshape(S, [sz nF]);
else
    idx = unitIdx(U);
    V = nan([sz nF]);
    for j = 1:nF
        fr = nan(sz);
        fr(idx) = S(:,j);
        V(:,:,j) = fr;
    end
end
P.frames = V;
P.fvals  = axisValues(ax.scrub, nF);
P.flab   = ax.scrub.label;
end

% =====================================================================
function P = xyPayload(P, V, ax)
%xyPayload  lagByRow: column 1 IS the coordinate and column 2 the measurement, so
%   there is nothing to reduce and nothing to look up.
V = double(V);
if size(V,2) < 2, P.note = 'This leaf does not carry an x and a y.'; return; end
P.x = V(:,1);
P.Y = V(:,2);
P.w = ones(size(V,1),1);
P.names = {''};
P.nUnits = size(V,1);
P.xlab = ax.x.label;
P.ok = true;
end

% ============================================================== THE COMPANIONS

function K = companionsFor(R, d)
%companionsFor  THE LEAVES THAT RIDE WITH THIS ONE, resolved out of the tree so no
%   renderer has to.  exSchema names them - which leaf, on which host, in what role
%   - and the host's own path supplies the container they are siblings in, so the
%   lookup is the host path with its last component swapped.
%
%   A COMPANION THE FILE DOES NOT CARRY IS SIMPLY ABSENT.  Nothing is padded and
%   nothing is invented: a window whose propagation never fitted has no slope and no
%   sentence, and the plot then shows what it always showed.
K = emptyCompanions();
if ~isstruct(d) || ~isfield(d,'leaf') || isempty(d.leaf) || isempty(d.family), return; end
spec = exSchema('companions', d.family, d.leaf);
if isempty(spec), return; end
base = regexprep(char(d.path), '[^.]+$', '');      % '...propagation.outer.'
for i = 1:numel(spec)
    pth = [base spec(i).leaf];
    v = getNested(R, pth);
    if isempty(v), continue; end
    lab = spec(i).leaf; un = '';
    r = exSchema('match', pth);
    if ~isempty(r), lab = r.label; un = r.yUnit; end
    K(end+1) = struct('leaf',spec(i).leaf, 'role',spec(i).role, 'value',v, ...
                      'label',lab, 'unit',un); %#ok<AGROW>
end
end

function K = emptyCompanions()
K = struct('leaf',{},'role',{},'value',{},'label',{},'unit',{});
end

function s = roleNoun(d)
%roleNoun  A rider's role, in the words the message it appears in needs.
r = ''; if isstruct(d) && isfield(d,'pairRole'), r = char(d.pairRole); end
switch r
    case 'interval', s = 'confidence interval';
    case 'fit',      s = 'fitted line';
    case 'caption',  s = 'caption';
    case 'pooled',   s = 'pooled form';
    otherwise,       s = 'companion';
end
end

function s = riderName(d)
%riderName  What to call a rider in that message.  A WHOLE CONTAINER MAY RIDE, and
%   then there is no leaf name to use - the myograph diameter's measure is the signal
%   token of its path - so the phrase the reader already saw stands in for it.
s = '';
if isstruct(d) && isfield(d,'leaf'), s = char(d.leaf); end
if isempty(s) && isstruct(d) && isfield(d,'label'), s = char(d.label); end
if isempty(s), s = 'This leaf'; end
end

% ================================================================== THE UNITS

function U = unitsFor(R, d, opts, nTotal)
%unitsFor  Which items of the leaf the user asked for, what they weigh, and what
%   they are called.  Four kinds of unit and one of them - the pixel - is selected
%   by geometry rather than by a table row.
%
%   ONLY .n AND .all ARE ALWAYS TRUE OF THE RECORD; everything else is an override.
%   .all says the selection is every unit of the leaf, which is the common case and
%   the one that used to be the expensive one: on a 1496x896 map, "all of them"
%   spelled out as a subscript vector, a weight of one per pixel and a name per pixel
%   is three arrays of 1.34 million elements built to say nothing.  So an arm that
%   selects a SUBSET writes .idx and turns .all off, an arm that has real weights or
%   real names writes those, and what is not written is derived on demand by
%   unitIdx / unitWeights / unitNames / unitIds - which is why a pixel leaf, whose
%   renderers are pictures and never show a unit name, pays for none of it.
%
% Built field by field, not through struct(): a cell passed to struct() makes a
% struct ARRAY of that cell's size, which is how one line here would quietly turn
% one unit list into nUnits of them.
U = struct('unit','', 'idx',[], 'w',[], 'names',{{}}, 'ids',[], 'note','', ...
           'n',0, 'all',true);
U.unit = char(d.unit);
U.n    = double(nTotal);

switch d.unit
    case {'seg','dvs'}
        T = metricsFor(R, d.unit);
        if isempty(T)
            % The only remaining reason to be here is that the recording carries no
            % metrics table.  It used to be reached far more often, by any file over
            % a size limit, and the sentence blamed the reading rather than the file.
            U.note = 'This recording carries no metrics table, so every unit counts equally.';
            return
        end
        if height(T) ~= nTotal
            U.note = sprintf(['The metrics table has %d rows against %d in the data, ' ...
                'so the selection was not applied.'], height(T), nTotal);
            return
        end
        rows = selectRows(T, opts.select, tableDomain(T));
        w    = weightCol(T);
        U.idx   = find(rows);
        U.n     = numel(U.idx);
        U.all   = false;
        U.w     = w(rows);
        U.ids   = idColumn(T, rows);
        U.names = rowNames(T, rows);

    case 'pixel'
        % A PIXEL IS NOT A TABLE ROW.  Nothing here but the selection: no name,
        % because no picture shows one; no weight, because pixels are equal by
        % construction; and no second construction of either.
        [keep, U.note] = pixelMask(R, opts, nTotal);
        if ~isempty(keep)
            U.idx = find(keep);
            U.n   = numel(U.idx);
            U.all = false;
        end

    case 'line'
        if numel(opts.rows)==2
            lo = max(1, floor(opts.rows(1))); hi = min(nTotal, ceil(opts.rows(2)));
            if hi>=lo, U.idx = (lo:hi)'; else, U.idx = zeros(0,1); end
            U.n     = numel(U.idx);
            U.all   = false;
            U.names = arrayfun(@(i) sprintf('position %d',i), U.idx(:).', 'UniformOutput', false);
        end
end
% A names list an arm built for a different count is worse than none: dropping it
% lets unitNames rebuild the defaults, and only if a plot ever asks for them.
if ~isempty(U.names) && numel(U.names)~=U.n, U.names = {}; end
end

% =====================================================================
function idx = unitIdx(U)
%unitIdx  The subscripts of the selected units, spelled out.  Only a caller that
%   genuinely has to scatter values into a canvas needs this; everything else asks
%   .all first and skips the indexing altogether.
if U.all, idx = (1:U.n)'; else, idx = U.idx(:); end
end

function w = unitWeights(U)
%unitWeights  What each selected unit weighs.  Ones unless a metrics table said
%   otherwise - and materialised HERE rather than in the record, so a picture, which
%   weights nothing, never builds a million of them.
if isempty(U.w), w = ones(U.n,1); else, w = U.w(:); end
end

function nm = unitNames(U)
%unitNames  What each selected unit is called, built on demand.  A box and a set of
%   curves need these; an image and a video have nowhere to put them.
if isempty(U.names), nm = defaultNames(U.unit, U.n); else, nm = U.names; end
end

function ids = unitIds(U)
%unitIds  The index each unit is known by in the map it was segmented from, which is
%   its row position when the recording carries no id column.
if isempty(U.ids), ids = unitIdx(U); else, ids = U.ids(:); end
end

% =====================================================================
function [keep, note] = pixelMask(R, opts, nTotal)
%pixelMask  The whole image, or the pixels of the segments a vessel selection picks.
%   For pixels the selection IS the reduction - there is no separate pooling step -
%   so this is where "Arteries (lumen)" turns into a region of the map.
%
%   AN EMPTY MASK MEANS EVERY PIXEL, and it is returned instead of true(nTotal,1)
%   because the whole image is the case this tool spends its time in: a mask of all
%   trues only ever becomes find(...), which is 1.34 million subscripts saying "all
%   of them" and then a copy of the leaf to reproduce the leaf.  An empty mask is not
%   the same answer as an all-FALSE one - that second one means the selection matched
%   nothing, and it still comes back as a mask so the caller can say so.
%
%   The per-pixel test below is the one that IS worth its cost: comparing every pixel
%   of the map against the selected segment ids is what "Arteries (lumen)" means, and
%   there is no cheaper way to ask it.
%
%   SEVERAL SELECTIONS ARE ONE REGION.  The masks combine spatially, so the union is
%   what "arteries and veins" means on a picture - and it is had for free, because
%   selectRows is where the OR lives.  One of them naming the whole image swallows
%   the rest, which is the same answer by a shorter road.
note = '';
keep = [];
sel  = asSelections(opts.select);
if all(cellfun(@isBlankSel, sel)) || any(cellfun(@isEverythingSel, sel)), return; end

map = getNested(R,'sMap');
T   = metricsFor(R,'seg');
if isempty(map) || isempty(T) || numel(map)~=nTotal
    note = 'This recording has no segment map, so the whole image is shown.';
    return
end
rows = selectRows(T, sel, tableDomain(T));
ids  = idColumn(T, rows);
if isempty(ids)
    note = 'No segment of this recording matches that selection.';
    keep = false(nTotal,1);
    return
end
keep = ismember(double(map(:)), double(ids(:)));
if ~any(keep)
    note = 'No pixel of this recording belongs to that selection.';
end
end

% =====================================================================
function items = whereList(R, d)
%whereList  What this file can actually be asked for.  Everything here comes from
%   the contents: a category that is not in the table is not offered, a vessel type
%   nothing carries is not offered, and a label no segment has is not offered.
items = {};
switch d.unit
    case 'whole', items = {'(the whole vessel)'}; return
    case 'line',  items = {'(all positions)'};    return
    case 'pixel'
        items = {'(whole image)'};
        T = metricsFor(R,'seg');
        if isempty(T) || isempty(getNested(R,'sMap')), return; end
        items = [items, segItems(T)];
        return
    case {'seg','dvs'}
        T = metricsFor(R, d.unit);
        if isempty(T), items = {'(everything)'}; return; end
        if strcmp(tableDomain(T),'dvs')
            items = dvsItems(T);
        else
            items = segItems(T);
        end
end
if isempty(items), items = {'(everything)'}; end
end

function items = segItems(T)
%segItems  The vessel selections a segment table supports, gated on its contents.
items = {};
vn = T.Properties.VariableNames;
cat = []; if ismember('category',vn), cat = double(T.category); end
typ = strings(0,1); if ismember('type',vn), typ = strtrim(string(T.type)); end

isLum = true(height(T),1);
if ~isempty(cat), isLum = (cat==5); end
if any(isLum), items{end+1} = 'All vessels (lumen)'; end
for t = {'Artery','Arteries (lumen)'; 'Vein','Veins (lumen)'; 'Uncertain','Uncertain (lumen)'}'
    if ~isempty(typ) && any(isLum & (typ==t{1})), items{end+1} = t{2}; end %#ok<AGROW>
end
if ~isempty(cat) && any(cat==1), items{end+1} = 'Parenchyma'; end

% The categories the five named selections do not cover.  These are the entries
% that make the list contents-driven rather than constant: a file whose segments
% were never split into walls simply does not offer a wall.
extra = { 2,'Unsegmented'; 3,'Outer wall'; 4,'Inner wall'; 0,'Background' };
for i = 1:size(extra,1)
    if ~isempty(cat) && any(cat==extra{i,1}), items{end+1} = extra{i,2}; end %#ok<AGROW>
end
items = [items, labelItems(T)];
end

function items = dvsItems(T)
items = {'All DVS'};
vn = T.Properties.VariableNames;
if ~ismember('type',vn), items = [items, labelItems(T)]; return; end
typ = strtrim(string(T.type));
for t = {'Artery','Arteries (DVS)'; 'Vein','Veins (DVS)'; 'Uncertain','Uncertain (DVS)'}'
    if any(typ==t{1}), items{end+1} = t{2}; end %#ok<AGROW>
end
items = [items, labelItems(T)];
end

function items = labelItems(T)
items = {};
if ~ismember('label',T.Properties.VariableNames), return; end
lv = strtrim(string(T.label));
lv = unique(lv(strlength(lv)>0));
for i = 1:numel(lv), items{end+1} = ['Label: ' char(lv(i))]; end %#ok<AGROW>
end

% =====================================================================

function dom = tableDomain(T)
if ismember('category',T.Properties.VariableNames), dom='seg'; else, dom='dvs'; end
end

function rows = selectRows(T, sel, dom)
%selectRows  A named selection - or SEVERAL - to a logical row mask into a metrics
%   table.  Several are ORed: a picture takes every entry at once, because arteries
%   and veins combine spatially rather than competing for one axes.
%
%   THE OTHERWISE-ARM FALLS BACK TO "ALL LUMEN", and that is the reason every
%   selection this file offers has to be named here.  A category entry that reached
%   the fallback would turn a request for the outer wall into a request for every
%   vessel: a plot of the wrong tissue, drawn without a word of complaint - and in a
%   UNION it is worse, because the extra vessels would be added to a picture that
%   still looked plausible.  exFetch('named',sel) asks the same switch whether it
%   knows a selection, and testExplorePlan asks it of every entry whereList produces.
rows = false(height(T),1);
sel  = asSelections(sel);
for i = 1:numel(sel)
    rows = rows | oneSelection(T, sel{i}, dom);
end
% drop rows with a non-finite idx (empty table slots)
if ismember('idx',T.Properties.VariableNames)
    rows = rows & isfinite(double(T.idx));
end
end

function rows = oneSelection(T, sel, dom)
%oneSelection  One entry of the Where list, resolved.  Kept apart from the union so
%   that "what does this one name" is a question with an answer - see isNamedSelection.
n = height(T); rows = false(n,1);
if isempty(sel), sel = ''; end
vn = T.Properties.VariableNames;
hasCat = ismember('category',vn); hasType = ismember('type',vn); hasLbl = ismember('label',vn);
if hasCat, isLum = (T.category==5); else, isLum = true(n,1); end   % DVS rows are all lumen
switch true
    case ~isempty(categoryOf(sel))
        % The tissue categories, numbered from zero in the order the launchers
        % declare them - background, parenchyma, unsegmented, outer edge, inner
        % edge, lumen - which is what makes lumen 5 and parenchyma 1 below.
        if hasCat, rows = double(T.category)==categoryOf(sel); end
    case startsWith(sel,'Label:')
        if hasLbl
            want = strtrim(sel(7:end));
            m = strtrim(string(T.label))==want;
            if strcmp(dom,'seg') && hasCat
                r = m & (T.category==5);
                if ~any(r), r = m & (T.category==1); end  % parenchyma-only label
                rows = r;
            else
                rows = m;
            end
        end
    case contains(sel,'Arter'), if hasType, rows = isLum & (string(T.type)=="Artery"); end
    case contains(sel,'Vein'),  if hasType, rows = isLum & (string(T.type)=="Vein"); end
    case contains(sel,'Uncert'),if hasType, rows = isLum & (string(T.type)=="Uncertain"); end
    case contains(sel,'Parenchyma'), if hasCat, rows = (T.category==1); end
    case contains(sel,'All DVS'), rows = true(n,1);
    case contains(sel,'All vessels'), rows = isLum;
    case isEverythingSel(sel) || isBlankSel(sel), rows = isLum;
    otherwise, rows = isLum;
end
end

function tf = isNamedSelection(sel)
%isNamedSelection  DOES THE SWITCH ABOVE KNOW THIS SELECTION, or would it reach the
%   fallback?  Asked here so a test can assert that nothing the Where list offers
%   ends up silently drawn as every vessel.  Every arm of oneSelection is answered,
%   in the same order, and the two are meant to be read side by side.
sel = char(sel);
tf = ~isempty(categoryOf(sel)) || startsWith(sel,'Label:') || ...
     contains(sel,'Arter') || contains(sel,'Vein') || contains(sel,'Uncert') || ...
     contains(sel,'Parenchyma') || contains(sel,'All DVS') || ...
     contains(sel,'All vessels') || isEverythingSel(sel) || isBlankSel(sel);
end

function tf = isEverythingSel(sel)
%isEverythingSel  The entries whereList writes when there is nothing to narrow BY.
%   They mean every unit, which is what the fallback does, and naming them keeps the
%   fallback itself unreachable from anything the tool offers.
tf = any(strcmpi(char(sel), {'(everything)','(whole image)','(all positions)', ...
                             '(the whole vessel)','all'}));
end

function tf = isBlankSel(sel)
%isBlankSel  No selection at all - the default, which also means every unit.
tf = isempty(strtrim(char(sel)));
end

function s = asSelections(sel)
%asSelections  ONE SELECTION OR SEVERAL, always as a cellstr, so nothing downstream
%   has to ask which it was handed.
if isempty(sel), s = {''}; return; end
if ischar(sel) || isstring(sel), s = {char(sel)}; return; end
s = cellfun(@char, reshape(sel,1,[]), 'UniformOutput', false);
if isempty(s), s = {''}; end
end

function c = categoryOf(sel)
%categoryOf  The category number behind a category selection, [] when the selection
%   is not one of them.
c = [];
switch char(sel)
    case 'Background',   c = 0;
    case 'Unsegmented',  c = 2;
    case 'Outer wall',   c = 3;
    case 'Inner wall',   c = 4;
end
end

function w = weightCol(T)
vn = T.Properties.VariableNames;
if ismember('area',vn), w=double(T.area);
elseif ismember('length',vn), w=double(T.length);
else, w=ones(height(T),1);
end
w(~isfinite(w)|w<=0)=eps;
end

function m = wmean(Y,w,dim)
% Weighted mean along dim, ignoring NaNs (equal weights if w empty). The weight
% vector is reshaped onto `dim` so it broadcasts against 2-D or 3-D Y.
if nargin<3, dim=2; end
if isempty(w), m=mean(Y,dim,'omitnan'); return; end
shp=ones(1,max(dim,2)); shp(dim)=numel(w); w=reshape(w(:),shp);
W = w.*isfinite(Y); Yz=Y; Yz(~isfinite(Y))=0;
m = sum(Yz.*W,dim)./max(sum(W,dim),eps);
end

% ================================================== READING THE LEAF

function [V, why] = readLeaf(R, d, sl)
%readLeaf  The numbers, sliced to the indices the chosen plot fixes.
why = '';
subs = leafSubs(d, sl);
V = getNested(R, d.path);
if isempty(V) || (~isnumeric(V) && ~islogical(V))
    V = []; why = sprintf('%s is not in this recording.', d.path); return
end
% THE SLICE IS TAKEN BEFORE THE CAST, not after.  A stride down a 61 208-sample leaf
% throws away 99% of it, and casting first would build the whole 725 MB in doubles to
% do so.  Nothing else changes: the values are the same values.
if numel(subs)>=2 && ~all(cellfun(@(s) ischar(s), subs))
    try
        V = V(subs{:});
    catch
        why = 'The requested slice is outside this leaf.'; V = []; return
    end
end
V = double(V);
end

% =====================================================================
function v = getNested(R, name)
% Fetch results.<a>.<b>... (dotted) from a full struct; [] if missing.
%
% Lifted from guiExplore and generalised in the two ways a descriptor path needs:
% a component of the form name(<digits>) takes that element of a struct array, and
% a final component naming a table VARIABLE returns the column.  A component is an
% index ONLY when it is name(<digits>), so the metrics column 'std(diameter)' is
% never mistaken for one.
v = [];
if isempty(R)||~isstruct(R), return; end
parts = strsplit(name,'.'); v=R;
for i=1:numel(parts)
    [nm, idx] = parseComp(parts{i});
    if istable(v)
        if ismember(nm, v.Properties.VariableNames), v = v.(nm); else, v=[]; return; end
        continue
    end
    if isstruct(v)&&isfield(v,nm), v=v.(nm); else, v=[]; return; end
    if ~isempty(idx)
        if ~isstruct(v) || idx>numel(v), v=[]; return; end
        v = v(idx);
    elseif isstruct(v) && ~isscalar(v)
        v=[]; return
    end
end
end

function [nm, idx] = parseComp(c)
idx = [];
t = regexp(c,'^(.+)\((\d+)\)$','tokens','once');
if isempty(t), nm = c; else, nm = t{1}; idx = str2double(t{2}); end
end

% ================================================== SHAPES

function M = canonical(V, d)
%canonical  Every container brought to [nUnit x n1 x n2], so nothing downstream has
%   to remember that a per-segment time vector is item LAST while a per-pixel map is
%   two-dimensional and item first.  This is the single simplification that keeps
%   the reductions and the slicing readable.
V  = double(V);
nD = numel(d.dims);
sz = padSize(size(V), max(4, nD+2));
switch d.layout
    case 'unitFirst', M = reshape(V, [sz(1), sz(2), sz(3)]);
    case 'unitLast'
        M = permute(V, [nD+1, 1:nD, nD+2]);
        s = padSize(size(M),3);
        M = reshape(M, [s(1), s(2), s(3)]);
    case 'pixel2D',   M = reshape(V, [sz(1)*sz(2), sz(3), sz(4)]);
    otherwise                                   % 'none': one unit, the whole vessel
        M = reshape(V, [1, numel(V), 1]);
        if nD>=2, M = reshape(V, [1, sz(1), sz(2)]); end
end
end

function M = reduceDims(M, d, ax)
%reduceDims  Average away the dimensions the chosen plot does not draw.  A plain
%   mean, not a weighted one: these are coordinates, and a frequency bin does not
%   have an area.
for i = 1:numel(ax.reduce)
    k = dimPos(d, ax.reduce{i});
    if k>0, M = mean(M, k+1, 'omitnan'); end
end
end

function k = dimPos(d, name)
%dimPos  Which of the leaf's dims this is, 0 when it is not one of them.
k = 0;
if isempty(name), return; end
for i = 1:numel(d.dims)
    if strcmp(d.dims{i}, name), k = i; return; end
end
end

function Y = squeezeTo2D(M, k)
%squeezeTo2D  [nUnit x n1 x n2] to [nUnit x nk], the other dim already reduced.
if k<=0, k = 1; end
Y = permute(M, [1, k+1, 4-k]);
Y = reshape(Y, size(Y,1), size(Y,2));
Y = double(Y);
end

% ================================================== SLICING

function sl = sliceIndices(d, ax, opts, A)
%sliceIndices  Which index of each dimension the plot fixes.  Nothing is ever fixed
%   silently at a value the caller did not ask for without being clamped to the axis
%   first: an out-of-range slice would otherwise read past the end of a lazy file.
sl = struct();
for i = 1:numel(ax.slice)
    nm = ax.slice{i};
    n  = axisLength(A, nm, d);
    sl.(nm) = clampIdx(sliceAsked(opts, nm), n);
end
% A PLOT THAT WALKS ITS AXIS AT A STRIDE CUTS THE LEAF ON THE WAY IN.  The subscript
% is the slider's own positions, so gsData's [61208 x 1480] never becomes 725 MB of
% doubles to throw 99% of away - and the payload reads the coordinates at the same
% positions, from the same function, so the two cannot describe different samples.
if ~isempty(ax.scrub.name) && isfinite(scrubCap(ax.plot))
    n = axisLength(A, ax.scrub.name, d);
    % A length nothing could measure is left alone: an empty subscript would read an
    % empty leaf and the payload would come back blank with nothing to say why.
    if isfinite(n) && n>=1
        sl.(ax.scrub.name) = scrubPositions(n, ax.plot);
    end
end
end

function v = sliceAsked(opts, dimName)
v = 1;
key = sliceKey(dimName);
if isstruct(opts.slice) && isfield(opts.slice, key) && ~isempty(opts.slice.(key))
    v = double(opts.slice.(key));
end
end

function k = sliceKey(dimName)
%sliceKey  The four names the cascade speaks, against the axis names the tree uses.
switch char(dimName)
    case 'f',                                  k = 'f';
    case {'pctLevel','pctBin'},                k = 'pct';
    case 'harmonic',                           k = 'harmonic';
    otherwise,                                 k = 'frame';
end
end

function i = clampIdx(i, n)
i = round(double(i(1)));
if ~isfinite(i) || i<1, i = 1; end
if isfinite(n) && n>=1 && i>n, i = n; end
end

function n = axisLength(A, nm, d)
n = Inf;
if isfield(A, nm) && isscalar(A.(nm).n), n = A.(nm).n; end
k = dimPos(d, nm);
if k>0
    sz = padSize(d.size, 6);
    switch d.layout
        case 'unitFirst', n = min(n, sz(k+1));
        case 'unitLast',  n = min(n, sz(k));
        case 'pixel2D',   n = min(n, sz(k+2));
        otherwise,        n = min(n, sz(k));
    end
end
end

function subs = leafSubs(d, sl)
%leafSubs  The slice, expressed in the LEAF's own dimension order - which is where
%   the layout column earns its keep, because the dims sit after the unit for a
%   per-segment scalar and before it for a per-segment time vector.
n = max(2, numel(d.size));
subs = repmat({':'}, 1, n);
if strcmp(d.layout,'xyPairs'), return; end
switch d.layout
    case 'unitFirst', base = 1;
    case 'pixel2D',   base = 2;
    otherwise,        base = 0;              % 'unitLast' and 'none'
end
for i = 1:numel(d.dims)
    nm = d.dims{i};
    if ~isfield(sl, nm), continue; end
    k = base + i;
    if k>=1 && k<=n, subs{k} = sl.(nm); end
end
end

% ================================================== SMALL THINGS

function v = axisValues(a, n)
%axisValues  The recording's own coordinate vector, cut or padded to the length the
%   data actually has.  A bare index is the fallback and it is never preferred.
if isempty(a.values), v = (1:n).'; return; end
v = clipTo(a.values, n);
end

function v = clipTo(values, n)
v = double(values(:));
if numel(v)>=n, v = v(1:n); else, v = [v; ((numel(v)+1):n).']; end
end

function sz = imageSize(d)
%imageSize  [Y X] for a per-pixel leaf, from the leaf itself.
sz = [];
if ~strcmp(d.layout,'pixel2D'), return; end
s = padSize(d.size,2);
sz = s(1:2);
end

function T = metricsFor(R, unit)
T = [];
if strcmp(unit,'seg'), nm = 'sMetrics'; else, nm = 'dvsMetrics'; end
if isempty(R) || ~isstruct(R), return; end
if isfield(R,nm) && istable(R.(nm)), T = R.(nm); end
end

function ids = idColumn(T, rows)
if ismember('idx', T.Properties.VariableNames)
    ids = double(T.idx(rows));
else
    r = find(rows); ids = r(:);
end
end

function names = rowNames(T, rows)
%rowNames  What to call one unit.  Its label when it has one, its type when it has
%   that, and its index otherwise - never the row number, which means nothing to a
%   reader of the plot.
r  = find(rows);
vn = T.Properties.VariableNames;
names = cell(1,numel(r));
lbl = strings(height(T),1); if ismember('label',vn), lbl = strtrim(string(T.label)); end
typ = strings(height(T),1); if ismember('type',vn),  typ = strtrim(string(T.type)); end
ids = idColumn(T, true(height(T),1));
for i = 1:numel(r)
    k = r(i);
    if strlength(lbl(k))>0,      names{i} = char(lbl(k));
    elseif strlength(typ(k))>0,  names{i} = sprintf('%s %d', char(typ(k)), ids(k));
    else,                        names{i} = sprintf('segment %d', ids(k));
    end
end
end

function names = defaultNames(unit, n)
switch char(unit)
    case 'pixel', pre = 'pixel';
    case 'line',  pre = 'position';
    case 'dvs',   pre = 'vessel';
    case 'seg',   pre = 'segment';
    otherwise,    pre = '';
end
if isempty(pre), names = repmat({''},1,n); return; end
names = arrayfun(@(i) sprintf('%s %d',pre,i), 1:n, 'UniformOutput', false);
end

function m = reduceMode(asked, ax)
m = char(asked);
if isempty(m), m = ax.units; end
if ~any(strcmp(m,{'pooled','points','curves'})), m = ax.units; end
if strcmp(ax.kind,'Vector') && strcmp(m,'points'), m = 'curves'; end
end

function tf = intervalMatches(pth, want)
t = regexp(pth,'intervals\((\d+)\)','tokens','once');
if isempty(t), tf = true; return; end
tf = str2double(t{1}) == double(want);
end

function s = scopeOf(pth)
k = strfind(pth,'.');
if isempty(k), s = ''; else, s = pth(1:k(end)-1); end
end

function sz = padSize(sz, n)
sz = reshape(double(sz),1,[]);
if numel(sz)<n, sz(end+1:n) = 1; end
end

function s = trim(s)
s = strtrim(char(s));
end

function opts = fillOpts(opts)
if isempty(opts) || ~isstruct(opts), opts = struct(); end
def = struct('select','', 'units','', 'rows',[], 'interval',[], ...
             'slice',struct(), 'path','', 'scrub','');
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || (isempty(opts.(fn{i})) && ~isstruct(def.(fn{i})))
        opts.(fn{i}) = def.(fn{i});
    end
end
opts.select = asSelections(opts.select);
opts.units  = char(opts.units);
opts.path   = char(opts.path);
opts.scrub  = char(opts.scrub);
end

function P = blankPayload(plotId)
P = struct('ok',false, 'note','', 'plot',char(plotId), 'kind','', ...
    'x',[], 'Y',[], 'w',[], 'names',{{}}, 'img',[], 'cdata',[], 'yvals',[], ...
    'frames',[], 'paint',[], 'fvals',[], 'xlab','', 'ylab','', 'clab','', 'flab','', ...
    'xlog',false, 'xscale','linear', 'nUnits',0, 'companions',{emptyCompanions()});
end
