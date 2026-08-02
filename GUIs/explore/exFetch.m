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
%   is a file too big to load, whose metrics table is opaque to HDF5 - and there the
%   payload says so in .note rather than pretending.
%
%   image.fromSegments IS THE INVERSE OF A ROW SELECTION and belongs here rather
%   than in a renderer.  Each selected unit's own number is painted into results.sMap
%   (dvsMap for tracked vessels) through the segment's idx, and every pixel that
%   belongs to no selected segment - or to one whose value is not finite - is left
%   NaN so it renders transparent.  It is the only way a per-segment scalar reaches
%   an image, and the only reason a segment index never has to become a coordinate.
%
%   A FILE TOO BIG TO LOAD IS READ IN SLABS.  When results is the lazy handle
%   struct('x_lazy_',true,'x_path_',p) the leaf is fetched with h5read, and the
%   indices the chosen plot fixes are pushed down into start/count so a slice of a
%   [65 x 70 x 289] spectrum reads 65 x 70 numbers and not five million.  MATLAB's
%   HDF5 layer preserves MATLAB orientation, so no permute is involved; the returned
%   size is checked against the requested count, and a read that does not match is
%   discarded rather than trusted.  A lazy file has no readable metrics table - a
%   table is an opaque compound to h5info - so unit selection falls back to every
%   unit with equal weights, and .note says which.
%
%   THE SLIDER'S AXIS IS ANSWERED HERE TOO.  exFetch('scrub',...) reports which
%   dimension a scrubbed plot leaves for a slider, its label, its coordinate values
%   and its length.  The window could not work that last one out for itself: the
%   length of a dimension is a fact about the leaf and the axis registry together,
%   and size(...,3) on a lazy file is a question nobody has read the answer to.
%
%   Syntax:
%      payload = exFetch(results, d, plotId)
%      payload = exFetch(results, d, plotId, opts)
%      items   = exFetch('where', results, d)
%      spec    = exFetch('scrub', results, d, plotId, opts)
%      key     = exFetch('slicekey', dimName)
%
%   INPUTS
%     results  the RESULTS member of a triplet, loaded; or the lazy handle
%              struct('x_lazy_',true,'x_path_',<file>).
%     d        one descriptor from exScan.
%     plotId   one of the ids exPlotRules('ids') lists.
%     opts     struct, every field optional:
%                .select    the Where entry, e.g. 'Arteries (lumen)', 'Label: MCA',
%                           '(whole image)', 'Outer wall'.  Default: everything.
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
%                           can be read from the sibling _s.mat.  Taken from the
%                           handle when the file is lazy.
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
%                .fvals   [nFrame x 1] where each frame sits on the scrub axis
%                .xlab .ylab .clab .flab   the axis titles
%                .xlog    true when x is logarithmic
%                .xscale  'linear' | 'log' | 'ordinal' | 'categorical'
%                .nUnits  how many units survived the selection
%     spec     struct: .name .label .values .n .choices - the dimension a slider
%              walks for this plot, empty .name when the plot has none
%     key      the opts.slice field name for one dimension: 'f' | 'pct' |
%              'harmonic' | 'frame'
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
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function out = exFetch(varargin)

a1 = varargin{1};
if ischar(a1) || isstring(a1)
    switch lower(char(a1))
        case 'where',    out = whereList(varargin{2:end});
        case 'scrub',    out = scrubSpecFor(varargin{2:end});
        case 'slicekey', out = sliceKey(varargin{2});
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
opts = fillOpts(opts, R);

% ---- 0. the guards -------------------------------------------------------------
if ~isempty(d.pairedWith)
    P.note = sprintf('%s is an error bar on %s, not a variable of its own.', ...
        d.leaf, d.pairedWith); return
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
if isempty(U.idx)
    P.note = trim([U.note ' Nothing is selected here.']); return
end
M = M(U.idx,:,:);
P.nUnits = numel(U.idx);

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
    otherwise
        P.note = sprintf('No fetch is defined for %s.', plotId); return
end
% A payload that could not build its array said why in .note; it must not come back
% claiming to be drawable just because it reached the end of the switch.
P.ok = hasContent(P);
end

function tf = hasContent(P)
tf = ~isempty(P.Y) || ~isempty(P.img) || ~isempty(P.cdata) || ~isempty(P.frames);
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
opts = fillOpts(opts, R);
A  = exAxes(R, scopeOf(d.path), opts.path);
ax = exPlotRules('axes', d, plotId, A, opts.scrub);
if isempty(ax.scrub.name), return; end
n = axisLength(A, ax.scrub.name, d);
if ~isfinite(n) || n<1, n = 0; end
s.name    = ax.scrub.name;
s.label   = ax.scrub.label;
s.n       = n;
s.values  = axisValues(ax.scrub, n);
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
if strcmp(red,'pooled')
    P.Y = wmean(v, U.w, 1);
    P.w = sum(U.w(:));
    P.names = {''};
else
    P.Y = v(:);
    P.w = U.w(:);
    P.names = U.names;
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
if isempty(ax.x.name) || strcmp(ax.unitAxis,'line')
    Y = double(M(:,1,1));                             % one value per position
    P.Y = Y(:);
    P.x = axisValues(ax.x, numel(Y));
    P.w = U.w(:);
    P.names = {''};
    return
end
Y = squeezeTo2D(M, k);                                % [nUnit x nx]
n = size(Y,2);
if strcmp(red,'pooled')
    P.Y = reshape(wmean(Y, U.w, 1), n, 1);
    P.w = sum(U.w(:));
    P.names = {''};
else
    P.Y = Y.';                                        % [nx x nUnit]
    P.w = U.w(:);
    P.names = U.names;
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
if isempty(ax.family.name) || strcmp(ax.family.name,'line') || dimPos(d,ax.family.name)==0
    Y = squeezeTo2D(M, kx);                           % [nUnit x nx]
    P.Y = Y.';
    P.x = axisValues(ax.x, size(Y,2));
    P.w = U.w(:);
    P.names = U.names;
    return
end
kf = dimPos(d, ax.family.name);
Yu = wmean(M, U.w, 1);                                % [1 x n1 x n2] pooled
Y  = permute(Yu, [kx+1, kf+1, 1]);                    % [nx x nFam]
Y  = reshape(Y, size(Y,1), size(Y,2));
P.Y = double(Y);
P.x = axisValues(ax.x, size(Y,1));
P.w = sum(U.w(:));
fv  = axisValues(ax.family, size(Y,2));
P.names = arrayfun(@(z) sprintf('%g%%', z), fv(:).', 'UniformOutput', false);
end

% =====================================================================
function P = imagePayload(P, M, U, d)
%imagePayload  A picture.  Every pixel outside the selection is NaN, which is what
%   the renderer turns into transparency rather than into a hole of zeros.
sz = imageSize(d);
if isempty(sz), P.note = trim([P.note ' This recording declares no image size.']); return; end
img = nan(sz);
v = double(M(:,1,1));
img(U.idx) = v;
P.img = img;
end

% =====================================================================
function P = paintPayload(P, M, U, R, d)
%paintPayload  image.fromSegments: each selected unit's own number, painted back
%   into the map by its idx.  The lookup is by index rather than by a comparison per
%   segment - 1480 segments against 300000 pixels is not a loop worth writing.
mapName = 'sMap';
if strcmp(d.unit,'dvs'), mapName = 'dvsMap'; end
map = getNested(R, mapName);
if isempty(map) || ~isnumeric(map)
    P.note = trim([P.note ' This recording has no ' mapName ' to paint into.']); return
end
ids = U.ids(:);
v   = double(M(:,1,1));
if isempty(ids) || ~all(isfinite(ids))
    P.note = trim([P.note ' The units of this recording carry no index to paint by.']);
    return
end
top = max(ids);
lut = nan(top+1,1);
lut(round(ids)+1) = v;
m = double(map);
ok = isfinite(m) & m==round(m) & m>=0 & m<=top;
img = nan(size(m));
img(ok) = lut(m(ok)+1);
P.img = img;
end

% =====================================================================
function P = heatPayload(P, M, U, d, ax)
%heatPayload  Rows are the y axis, columns the x - frequency up the side and time or
%   percentile along the bottom, which is the orientation guiExplore already draws.
Yu = wmean(M, U.w, 1);
kx = dimPos(d, ax.x.name); ky = dimPos(d, ax.y.name);
if kx==0 || ky==0, P.note = trim([P.note ' This leaf has no two axes to spread out.']); return; end
C = permute(Yu, [ky+1, kx+1, 1]);
P.cdata = double(reshape(C, size(C,1), size(C,2)));
P.x     = axisValues(ax.x, size(P.cdata,2));
P.yvals = axisValues(ax.y, size(P.cdata,1));
P.w     = sum(U.w(:));
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
P.w     = U.w(:);
end

% =====================================================================
function P = videoPayload(P, M, U, d, ax)
%videoPayload  A cube whose third dimension is the one the slider walks.
sz = imageSize(d);
if isempty(sz), P.note = trim([P.note ' This recording declares no image size.']); return; end
k = dimPos(d, ax.scrub.name);
if k==0, P.note = trim([P.note ' There is nothing to scrub along.']); return; end
S = squeezeTo2D(M, k);                                % [nUnit x nFrame]
nF = size(S,2);
V = nan([sz nF]);
for j = 1:nF
    fr = nan(sz);
    fr(U.idx) = S(:,j);
    V(:,:,j) = fr;
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

% ================================================================== THE UNITS

function U = unitsFor(R, d, opts, nTotal)
%unitsFor  Which items of the leaf the user asked for, what they weigh, and what
%   they are called.  Four kinds of unit and one of them - the pixel - is selected
%   by geometry rather than by a table row.
% Built field by field, not through struct(): a cell passed to struct() makes a
% struct ARRAY of that cell's size, which is how one line here would quietly turn
% one unit list into nUnits of them.
U = struct('idx',[], 'w',[], 'names',{{}}, 'ids',[], 'note','');
U.idx   = (1:nTotal)';
U.w     = ones(nTotal,1);
U.ids   = (1:nTotal)';
U.names = defaultNames(d.unit, nTotal);

switch d.unit
    case {'seg','dvs'}
        T = metricsFor(R, d.unit);
        if isempty(T)
            U.note = 'The metrics table could not be read here, so every unit counts equally.';
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
        U.w     = w(rows);
        U.ids   = idColumn(T, rows);
        U.names = rowNames(T, rows);

    case 'pixel'
        [keep, note] = pixelMask(R, opts, nTotal);
        U.idx   = find(keep);
        U.w     = ones(numel(U.idx),1);
        U.ids   = U.idx;
        U.names = defaultNames('pixel', numel(U.idx));
        U.note  = note;

    case 'line'
        keep = true(nTotal,1);
        if numel(opts.rows)==2
            lo = max(1, floor(opts.rows(1))); hi = min(nTotal, ceil(opts.rows(2)));
            keep = false(nTotal,1);
            if hi>=lo, keep(lo:hi) = true; end
        end
        U.idx = find(keep);
        U.w   = ones(numel(U.idx),1);
        U.ids = U.idx;
        U.names = arrayfun(@(i) sprintf('position %d',i), U.idx(:).', 'UniformOutput', false);
end
if isempty(U.names) || numel(U.names)~=numel(U.idx)
    U.names = defaultNames(d.unit, numel(U.idx));
end
end

% =====================================================================
function [keep, note] = pixelMask(R, opts, nTotal)
%pixelMask  The whole image, or the pixels of the segments a vessel selection picks.
%   For pixels the selection IS the reduction - there is no separate pooling step -
%   so this is where "Arteries (lumen)" turns into a region of the map.
note = '';
keep = true(nTotal,1);
sel  = opts.select;
if isempty(sel) || strcmpi(sel,'(whole image)') || strcmpi(sel,'all'), return; end

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
%selectRows  A named selection to a logical row mask into a metrics table.
%
%   THE OTHERWISE-ARM FALLS BACK TO "ALL LUMEN", and that is the reason every
%   selection this file offers has to be named here.  A category entry that reached
%   the fallback would turn a request for the outer wall into a request for every
%   vessel: a plot of the wrong tissue, drawn without a word of complaint.  The
%   category arm below was a wrapper around this function until the redesign
%   retired guiExplore's copy; it is folded in so there is one list of selections
%   rather than two that have to agree.
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
    otherwise, rows = isLum;
end
% drop rows with a non-finite idx (empty table slots)
if ismember('idx',vn)
    rows = rows & isfinite(double(T.idx));
end
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

% ================================================== READING THE LEAF, EITHER WAY

function [V, why] = readLeaf(R, d, sl)
%readLeaf  The numbers, sliced on the way in when that saves a read.
why = '';
subs = leafSubs(d, sl);
if isfield(R,'x_lazy_')
    [V, why] = readSlab(char(R.x_path_), d, subs);
    return
end
V = getNested(R, d.path);
if isempty(V) || (~isnumeric(V) && ~islogical(V))
    V = []; why = sprintf('%s is not in this recording.', d.path); return
end
V = double(V);
if numel(subs)>=2 && ~all(cellfun(@(s) ischar(s), subs))
    try
        V = V(subs{:});
    catch
        why = 'The requested slice is outside this leaf.'; V = []; return
    end
end
end

function [V, why] = readSlab(pth, d, subs)
%readSlab  One hyperslab of an HDF5 dataset.  MATLAB's HDF5 layer reports and
%   returns MATLAB orientation, so start and count are in the same dimension order
%   as the size exScan measured - but the returned size is CHECKED against the
%   requested count anyway, because a silently transposed read would give plausible
%   numbers from the wrong place.
why = '';
ds = ['/results/' strrep(d.path,'.','/')];
sz = padSize(d.size, numel(subs));
start = ones(1,numel(sz)); count = sz;
for k = 1:numel(subs)
    if ischar(subs{k}), continue; end
    s = subs{k};
    start(k) = min(s); count(k) = max(s)-min(s)+1;
end
try
    if all(start==1) && isequal(count, sz)
        V = h5read(pth, ds);
    else
        V = h5read(pth, ds, start, count);
    end
catch ME
    why = sprintf('This leaf could not be read from the file (%s).', ME.identifier);
    V = []; return
end
if ~isequal(padSize(size(V),numel(count)), padSize(count,numel(count)))
    why = 'The file returned a block of the wrong shape, so it was not used.';
    V = []; return
end
V = double(V);
% re-index inside the block, so a non-contiguous request still lands right
inner = subs;
for k = 1:numel(subs)
    if ischar(subs{k}), continue; end
    inner{k} = subs{k} - start(k) + 1;
end
if ~all(cellfun(@(s) ischar(s), inner))
    % In range by construction: the block spans min(s) to max(s) of every request,
    % so no guard is needed and a silent one would only hide a real mistake.
    V = V(inner{:});
end
end

% =====================================================================
function v = getNested(R, name)
% Fetch results.<a>.<b>... (dotted) from a full struct; [] if missing/lazy.
%
% Lifted from guiExplore and generalised in the two ways a descriptor path needs:
% a component of the form name(<digits>) takes that element of a struct array, and
% a final component naming a table VARIABLE returns the column.  A component is an
% index ONLY when it is name(<digits>), so the metrics column 'std(diameter)' is
% never mistaken for one.
v = [];
if isempty(R)||~isstruct(R)||isfield(R,'x_lazy_'), return; end
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
if isempty(R) || ~isstruct(R) || isfield(R,'x_lazy_'), return; end
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

function opts = fillOpts(opts, R)
if isempty(opts) || ~isstruct(opts), opts = struct(); end
def = struct('select','', 'units','', 'rows',[], 'interval',[], ...
             'slice',struct(), 'path','', 'scrub','');
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || (isempty(opts.(fn{i})) && ~isstruct(def.(fn{i})))
        opts.(fn{i}) = def.(fn{i});
    end
end
opts.select = char(opts.select);
opts.units  = char(opts.units);
opts.path   = char(opts.path);
opts.scrub  = char(opts.scrub);
if isempty(opts.path) && isstruct(R) && isfield(R,'x_path_')
    opts.path = char(R.x_path_);
end
end

function P = blankPayload(plotId)
P = struct('ok',false, 'note','', 'plot',char(plotId), 'kind','', ...
    'x',[], 'Y',[], 'w',[], 'names',{{}}, 'img',[], 'cdata',[], 'yvals',[], ...
    'frames',[], 'fvals',[], 'xlab','', 'ylab','', 'clab','', 'flab','', ...
    'xlog',false, 'xscale','linear', 'nUnits',0);
end
