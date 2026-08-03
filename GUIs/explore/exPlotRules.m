%exPlotRules - Which plots a variable can legally appear on, and on which axes.
%
%   P = exPlotRules(d) is the list of plots that make sense for one descriptor,
%   first entry the default.  K = exPlotRules('kinds',d) is which of Scalar /
%   Vector / Image / Video it appears under, and ax = exPlotRules('axes',d,plotId)
%   is everything a renderer needs to draw it: what goes on x, on y, on the colour
%   bar, which axis indexes a family of curves, which one a scrub slider walks, and
%   which dimensions have to be sliced or averaged away first.
%
%   It reads NOTHING.  A descriptor in, a decision out - no file, no struct, no
%   sizes beyond the ones exScan already measured.  exFetch is what turns the
%   decision into numbers.
%
%   THE LOAD-BEARING RULE: A SEGMENT INDEX IS NOT A COORDINATE.  There is no
%   meaningful "amplitude against segment number" plot, so no plot here ever puts
%   seg or dvs on an axis.  A per-segment scalar is legal in exactly two ways -
%   REDUCED over a selection (box / bar, then compared across strata) or PAINTED
%   back into space through sMap (image.fromSegments) - and in no other.  That is
%   the whole argument of the redesign and it is enforced structurally: a dim whose
%   axis is a label id is stripped from the plot key and averaged away, so it cannot
%   reach an axis even by the inferred fallback path.
%
%   line AND interval ARE COORDINATES, and this is the distinction the old code
%   blurred.  Position along a vessel is a real ordering, so curve.position and the
%   position-time kymograph are legitimate; protocol order is a real ordering, so
%   interval is a legitimate categorical x on a box - it arrives as the COMPARISON
%   axis of a box or bar, which is why box and bar report x with an empty name and
%   a categorical scale rather than naming an axis of their own.
%
%   ONE LEAF, SEVERAL KINDS.  A per-segment scalar is Scalar (a box across strata)
%   AND Image (the same numbers painted back into the map); an [nSeg x nF x nD]
%   spectrum is Image (heat.ft) AND Vector (a mean over one of the two).  The Kind
%   control filters on membership, so a leaf appears under every kind it can be
%   drawn as - which is the question the user is really answering when they pick
%   one.  KINDS ARE DERIVED FROM THE PLOT LIST, never listed separately, so the two
%   answers cannot drift apart.
%
%   TWO PERCENTILE AXES, ONE PLOT ID.  curve.pct serves both axes of the shape
%   model: the percentile LEVELS (11 of them, read from the sibling settings file)
%   carry scalars.<band>.ampPct, and the percentile BIN CENTRES (10, stored in the
%   results) index the family of curves of fVectors.<band>.ampMeanPct.  Same units,
%   different axes, one shorter than the other by construction - so the label and
%   the values are taken from the axis registry and never from a constant here.
%
%   WHERE THIS TABLE EXTENDS SPEC SECTION 4.1.  Five combinations occur in the real
%   files or the fixtures that the frozen table does not name, and each is filled by
%   the nearest analogous row rather than by inventing a plot id - except the last,
%   which IS a new id and is argued for below:
%
%     pixel x f          per-pixel fVectors [Y X nF]        -> as pixel x pct
%     pixel x f,pct      per-pixel ampMeanPct [Y X nF nPct] -> as pixel x f,timeDWT
%     pixel x plane      commonMask [Y X nPlane]            -> ONE FRAME AT A TIME,
%                        and never a video - see the next paragraph
%     line  x <anything> a per-line myograph <VSM>          -> as seg/dvs, plus the
%                        family of curves that only an ORDERED unit makes meaningful
%     seg/dvs x time     a per-segment time series          -> video.fromSegments,
%                        the video an LSCI recording really has
%
%   A combination with no row returns an EMPTY plot list, deliberately: a leaf the
%   tool can see but cannot draw is a hole in the rule table, and testExplorePlan
%   asserts there are none.  A silent default would hide it.
%
%   A MASK IS NOT A VIDEO (D9).  commonMask is [Y x X x nPlane] and plane is a STACK
%   OF LAYERS, not a coordinate anything was measured against, so scrubbing it plays
%   a slideshow of bookkeeping.  On both reference contrast files it is the only
%   [Y X n] leaf there is - pulsatility.ppx holds scalars and no timeVectors, and the
%   tracked branch has no ppx group at all - so Kind = Video collapsed to one family,
%   one signal and one variable, every step of the cascade auto-hid, and the window
%   showed the shared analysed area.  The plane row therefore offers image.frame and
%   nothing else: looking at layer 2 is a legitimate question and it is the whole of
%   what a mask has to answer.  The change is HERE and not in exSchema, because plane
%   really is a dimension of that array and a schema that denied it would be lying
%   about the shape to win an argument about the plots.
%
%   AND THE VIDEO AN LSCI RECORDING DOES HAVE is a per-segment time series painted
%   back through sMap - sData, gsData, vasomotion.<sig>.timeVectors.<band>.amp,
%   pulsatility.<sig>.timeVectors.fData - one frame per sample, every segment showing
%   its own value.  It is image.fromSegments once per sample, which is why it sits on
%   the seg/dvs x time row and why a per-segment leaf is worth keeping resident.
%   THE CUBE IS NEVER MATERIALISED: 2448 x 1496 x 896 doubles is 26 TB, so exFetch
%   holds the [nUnit x nFrame] matrix and paints the frame the slider is on.  AND A
%   LONG RECORDING IS WALKED AT A STRIDE - gsData has 61 208 samples, and a slider
%   with 61 208 positions is a defect even when every frame is correct - so exFetch
%   caps the slider at a few hundred evenly spaced positions, each a real sample
%   labelled in seconds off the leaf's own time base.  Nothing is averaged, no frame
%   is invented, and the cap lives there rather than here because this file reads
%   nothing and a stride is a fact about the leaf's length.
%
%   Syntax:
%      P  = exPlotRules(d)
%      P  = exPlotRules('plots', d)
%      K  = exPlotRules('kinds', d)
%      ax = exPlotRules('axes', d, plotId)
%      ax = exPlotRules('axes', d, plotId, A)
%      k  = exPlotRules('kindof', plotId)
%      ids = exPlotRules('ids')
%
%   INPUTS
%     d       one descriptor from exScan.
%     plotId  one of:
%               box  bar
%               curve.time  curve.f  curve.pct  curve.harmonic  curve.position
%               curves.family
%               image  image.fromSegments  image.timeAverage  image.frame
%               heat.ft  heat.fpct  heat.positionTime
%               video  video.fromSegments
%     A       the axis registry from exAxes, resolved at the leaf's scope.  Optional:
%             with it the axis labels, scales and coordinate VALUES are the
%             recording's own; without it they are generic fallbacks and .values is
%             empty.  exFetch always passes it.
%     want    'axes' only, optional: which dimension the slider should walk, when the
%             caller has offered the choice in .scrubChoices and the user has made
%             it.  Ignored when the leaf does not have that dimension.
%
%   OUTPUTS
%     P   cellstr of plot ids, first = the default.  {} when the descriptor cannot
%         be drawn AS A VARIABLE OF ITS OWN - which today means a COMPANION: a leaf
%         exSchema declares as riding with another one, in a role (an interval on a
%         scalar, the fit of a scatter, the caption of a plot).  exFetch hands those
%         to the renderer with their host rather than as menu entries.
%     K   cellstr of kinds, in the order the plots introduce them.
%     ax  struct describing one plot:
%           .plot .kind
%           .x .y .c .family .scrub   each .name .label .scale .values .log
%           .scrubChoices  the dims a video slider could walk, when there is a choice
%           .slice         dims that must be fixed at one index before drawing
%           .reduce        dims that are averaged away before drawing
%           .unitAxis      the axis the UNIT itself occupies when it is drawn rather
%                          than reduced ('pixel', 'line'); '' when it is reduced
%           .units         the within-file reduction this plot wants by default:
%                          'pooled' | 'points' | 'curves'
%
%   EXAMPLE
%      D = exScan(load('..._c_BFI_r.mat').results);
%      d = D(strcmp({D.path},'pulsatility.sData.scalars.psPI'));
%      exPlotRules(d)                       % {'box','bar','image.fromSegments'}
%      exPlotRules('kinds',d)               % {'Scalar','Image'}
%      exPlotRules('axes',d,'box').units    % 'points'
%      d = D(strcmp({D.path},'sData'));
%      exPlotRules(d)                       % {'curve.time','video.fromSegments'}
%
%   DEPENDS ON
%     nothing.
%
% See also: exFetch, exScan, exAxes, exSchema, exModality, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function varargout = exPlotRules(varargin)

if nargin==0, error('exPlotRules:args','A descriptor or an action is required.'); end

a1 = varargin{1};
if isstruct(a1), varargout{1} = plotsFor(a1); return; end

switch lower(char(a1))
    case 'plots',  varargout{1} = plotsFor(varargin{2});
    case 'kinds',  varargout{1} = kindsFor(varargin{2});
    case 'axes',   varargout{1} = axesFor(varargin{2:end});
    case 'kindof', varargout{1} = kindOf(varargin{2});
    case 'ids',    varargout{1} = allPlotIds();
    otherwise
        error('exPlotRules:action','Unknown action ''%s''.', char(a1));
end
end

% ================================================================= THE RULE TABLE

function P = plotsFor(d)
%plotsFor  Spec section 4.1, one switch per unit class.  The order inside each arm
%   IS the preference order - the first entry is what the cascade selects when the
%   user picks the variable and says nothing else.
P = {};
if ~isempty(d.pairedWith), return; end   % a companion: it rides with its host
key = plotKey(d);

switch unitClass(d.unit)
    case 'item'                                        % seg / dvs
        % THE TIME ROW CARRIES THE RECONSTRUCTION.  A value per segment per sample
        % painted back through sMap, one frame per sample, is the video an LSCI
        % recording actually has - the same operation image.fromSegments does once,
        % done once per sample.  It is not on the f,time row: a spectrum cube would
        % have to be sliced at one frequency first, and that is a second question.
        switch key
            case '',         P = {'box','bar','image.fromSegments'};
            case 'pct',      P = {'curve.pct'};
            case 'harmonic', P = {'curve.harmonic'};
            case 'f',        P = {'curve.f'};
            case 'time',     P = {'curve.time','video.fromSegments'};
            case 'f,time',   P = {'heat.ft','curve.f','curve.time'};
            case 'f,pct',    P = {'curves.family','heat.fpct'};
        end

    case 'pixel'
        % A PER-PIXEL RESULT IS A PICTURE, AND ONLY A PICTURE.  It has no unit a
        % reader can name, so there is nothing for a curve or a number to be OF: an
        % ROI curve over a per-pixel time series used to be offered here and it is
        % not a quantity the file holds, it is an average over whatever happened to
        % be selected.  Averaging pixels is what a SEGMENT is for, and the per-segment
        % branch beside it already carries that answer with the units named.
        % A per-LINE myograph leaf is the opposite case and keeps its curves: a line
        % is a position across the vessel, which is a unit exactly as a segment is.
        % AND plane - the 'frame' key - IS NOT A COORDINATE EITHER (D9).  It indexes
        % the layers stacked into commonMask, so one layer at a time is the whole of
        % what it can be asked; a video of it is a slideshow of bookkeeping, and on
        % the reference files it was the ONLY thing Kind = Video had to offer.
        switch key
            case '',                         P = {'image'};
            case 'frame',                    P = {'image.frame'};
            case {'pct','harmonic','f'}
                                             P = {'video','image.frame'};
            case 'time',                     P = {'video','image.timeAverage'};
            case 'f,time',                   P = {'video','image.timeAverage','heat.ft'};
            case 'f,pct',                    P = {'video','image.frame','heat.fpct'};
        end

    case 'line'
        % Everything seg/dvs can do, plus the two things only an ORDERED unit makes
        % meaningful: a curve along the vessel, and one curve per position.
        switch key
            case '',         P = {'curve.position','box'};
            case 'xy',       P = {'curve.position'};
            case 'position', P = {'curve.position'};
            case 'time',     P = {'curve.time','heat.positionTime','curves.family'};
            case 'f',        P = {'curve.f','curves.family'};
            case 'pct',      P = {'curve.pct','curves.family'};
            case 'harmonic', P = {'curve.harmonic','curves.family'};
            case 'f,time',   P = {'heat.ft','curve.f','curve.time'};
            case 'f,pct',    P = {'curves.family','heat.fpct'};
        end

    case 'whole'
        switch key
            case '',         P = {'box'};
            case 'time',     P = {'curve.time'};
            case 'f',        P = {'curve.f'};
            case 'pct',      P = {'curve.pct'};
            case 'harmonic', P = {'curve.harmonic'};
            case {'xy','position'}
                             P = {'curve.position'};
            case 'f,time',   P = {'heat.ft','curve.f','curve.time'};
            case 'f,pct',    P = {'curves.family','heat.fpct'};
        end
end
end

% =====================================================================
function K = kindsFor(d)
%kindsFor  Derived from the plot list, in the order the plots introduce them - so
%   the kinds column of the rule table cannot drift away from the plots column.
K = {};
P = plotsFor(d);
for i = 1:numel(P)
    k = kindOf(P{i});
    if ~any(strcmp(K,k)), K{end+1} = k; end %#ok<AGROW>
end
end

function k = kindOf(plotId)
%kindOf  What the user sees on the axes, which is not the same as what the array is.
switch char(plotId)
    case {'box','bar'},                                  k = 'Scalar';
    case {'curve.time','curve.f','curve.pct', ...
          'curve.harmonic','curve.position','curves.family'}, k = 'Vector';
    case {'image','image.fromSegments','image.timeAverage', ...
          'image.frame','heat.ft','heat.fpct','heat.positionTime'}, k = 'Image';
    case {'video','video.fromSegments'},                 k = 'Video';
    otherwise
        error('exPlotRules:plotId','Unknown plot ''%s''.', char(plotId));
end
end

function ids = allPlotIds()
ids = {'box','bar','curve.time','curve.f','curve.pct','curve.harmonic', ...
       'curve.position','curves.family','image','image.fromSegments', ...
       'image.timeAverage','image.frame','heat.ft','heat.fpct', ...
       'heat.positionTime','video','video.fromSegments'};
end

% ================================================================= THE PLOT KEY

function key = plotKey(d)
%plotKey  The dims of a leaf, reduced to the handful of CLASSES the rule table is
%   written against.  Two axes of the same class are the same question: timeWT and
%   gsTime are both a time base, pctLevel and pctBin are both a percentile axis.
%
%   A LABEL ID IS DROPPED HERE, and that is what makes the "seg is not a coordinate"
%   rule structural rather than a convention.  Only the inferred fallback path can
%   produce such a dim (a leaf whose second size happens to equal the segment
%   count); dropping it sends the leaf to the no-dims row, where it is reduced.
if strcmp(d.layout,'xyPairs'), key = 'xy'; return; end
c = {};
for i = 1:numel(d.dims)
    a = axisClass(d.dims{i});
    if isempty(a), continue; end
    c{end+1} = a; %#ok<AGROW>
end
key = strjoin(c,',');
end

function c = axisClass(nm)
%axisClass  '' means "not drawable" - a label id, or an axis nothing knows about.
switch char(nm)
    case {'time','gsTime','timeWT','timeDWT'}, c = 'time';
    case 'f',                                  c = 'f';
    case {'pctLevel','pctBin'},                c = 'pct';
    case 'harmonic',                           c = 'harmonic';
    case 'plane',                              c = 'frame';
    case 'line',                               c = 'position';
    otherwise,                                 c = '';
end
end

function u = unitClass(unit)
switch char(unit)
    case {'seg','dvs'}, u = 'item';
    case 'pixel',       u = 'pixel';
    case 'line',        u = 'line';
    otherwise,          u = 'whole';
end
end

% ================================================================= THE AXES

function ax = axesFor(d, plotId, A, want)
%axesFor  What one plot draws, and what has to happen to the array first.
%   WANT names which dimension the slider should walk, when the caller has been
%   offered a choice and made one.  It is honoured only if it is one this leaf
%   actually has, so a stale choice left over from another variable falls back to the
%   default rather than producing an empty picture.
if nargin<3, A = struct(); end
if nargin<4, want = ''; end
plotId = char(plotId);
ax = blankAxes(plotId);

tD = dimOfClass(d,'time');       % the time base this leaf lives on, when it has one
fD = dimOfClass(d,'f');
pD = dimOfClass(d,'pct');
hD = dimOfClass(d,'harmonic');
lD = dimOfClass(d,'position');
% commonMask's plane needs no name of its own: the only plots it reaches are the
% scrubbed ones, which take their dimension from scrubAx.

vLab = varLabel(d);

switch plotId
    % ---- Scalar: one number per unit, compared ACROSS strata -------------------
    case {'box','bar'}
        % x is deliberately unnamed: it is the comparison axis the cascade picks -
        % group, animal, type, recording, or INTERVAL - and interval is a real
        % ordering, which is why the scale is categorical rather than absent.
        ax.x = mkAx('', 'selection', 'categorical', [], A);
        ax.y = mkVar(vLab);
        ax.reduce = drawableDims(d, {});
        if strcmp(plotId,'box'), ax.units = 'points'; else, ax.units = 'pooled'; end

    % ---- Vector: curves --------------------------------------------------------
    case 'curve.time'
        ax.x = mkAx(tD, 'time (s)', 'linear', tD, A);
        ax.y = mkVar(vLab);
        ax.reduce = drawableDims(d, {tD});

    case 'curve.f'
        ax.x = mkAx(fD, 'frequency (Hz)', 'log', fD, A);
        ax.x.log = true;
        ax.y = mkVar(vLab);
        ax.reduce = drawableDims(d, {fD});

    case 'curve.pct'
        ax.x = mkAx(pD, 'percentile (%)', 'linear', pD, A);
        ax.y = mkVar(vLab);
        ax.reduce = drawableDims(d, {pD});

    case 'curve.harmonic'
        ax.x = mkAx(hD, 'harmonic', 'ordinal', hD, A);
        ax.y = mkVar(vLab);
        ax.reduce = drawableDims(d, {hD});

    case 'curve.position'
        % Two shapes reach here: a per-line leaf, whose x IS the line axis, and an
        % xyPairs leaf, which carries its own x in column 1 and is therefore NOT
        % sized against the line axis - see exSchema on why verifying it would mark
        % a correct leaf suspect.
        if strcmp(d.layout,'xyPairs')
            ax.x = mkAx('line','position along the vessel','linear','',A);
            ax.x.values = [];              % exFetch takes it from column 1
        elseif ~isempty(lD)
            ax.x = mkAx(lD,'position along the vessel','linear',lD,A);
            ax.reduce = drawableDims(d, {lD});
        else
            ax.x = mkAx('line','position along the vessel','linear','line',A);
            ax.unitAxis = 'line';
        end
        ax.y = mkVar(vLab);

    case 'curves.family'
        % One axes, several curves.  The family index is the second drawable dim
        % when there is one - the percentile bins of ampMeanPct - and otherwise the
        % UNIT, which is only meaningful because a line is an ordered coordinate.
        if ~isempty(fD) && ~isempty(pD)
            ax.x      = mkAx(fD,'frequency (Hz)','log',fD,A); ax.x.log = true;
            ax.family = mkAx(pD,'percentile (%)','linear',pD,A);
            ax.reduce = drawableDims(d, {fD,pD});
        else
            xn = firstNonEmpty({tD,fD,pD,hD});
            ax.x      = mkAx(xn, defaultLabel(xn), defaultScale(xn), xn, A);
            ax.x.log  = strcmp(defaultScale(xn),'log');
            ax.family = mkAx('line','position along the vessel','linear','line',A);
            ax.unitAxis = 'line';
            ax.units    = 'curves';
            ax.reduce   = drawableDims(d, {xn});
        end
        ax.y = mkVar(vLab);

    % ---- Image: a 2-D field, a picture or a heatmap -----------------------------
    case {'image','image.frame','image.timeAverage','image.fromSegments'}
        ax.x = mkAx('','','linear',[],A);
        ax.y = mkAx('','','linear',[],A);
        ax.c = mkVar(vLab);
        switch plotId
            case 'image'
                ax.unitAxis = 'pixel';
            case 'image.fromSegments'
                % The inverse of a row selection: each unit's own number, painted
                % back into the map.  It is the only way a per-segment scalar
                % reaches an image, and the units are kept, never pooled.
                ax.unitAxis = 'pixel';
                ax.units    = 'points';
            case 'image.frame'
                ax.unitAxis = 'pixel';
                ax.scrub    = scrubAx(d, A, want);
                ax.slice    = drawableDims(d, {});
                ax.scrubChoices = ax.slice;
            case 'image.timeAverage'
                ax.unitAxis = 'pixel';
                ax.reduce   = intersectDims(d, {tD});
                ax.slice    = setdiffDims(drawableDims(d,{}), ax.reduce);
        end

    case 'heat.ft'
        ax.x = mkAx(tD,'time (s)','linear',tD,A);
        ax.y = mkAx(fD,'frequency (Hz)','log',fD,A);
        ax.c = mkVar(vLab);
        ax.reduce = drawableDims(d, {tD,fD});

    case 'heat.fpct'
        ax.x = mkAx(pD,'percentile (%)','linear',pD,A);
        ax.y = mkAx(fD,'frequency (Hz)','log',fD,A);
        ax.c = mkVar(vLab);
        ax.reduce = drawableDims(d, {pD,fD});

    case 'heat.positionTime'
        % The kymograph: the UNIT is the y axis, so nothing is reduced over units.
        ax.x = mkAx(tD,'time (s)','linear',tD,A);
        ax.y = mkAx('line','position along the vessel','linear','line',A);
        ax.c = mkVar(vLab);
        ax.unitAxis = 'line';
        ax.units    = 'curves';
        ax.reduce   = drawableDims(d, {tD});

    % ---- Video: a field with a dimension you scrub ------------------------------
    case 'video'
        ax.x = mkAx('','','linear',[],A);
        ax.y = mkAx('','','linear',[],A);
        ax.c = mkVar(vLab);
        ax.unitAxis = 'pixel';
        ax.scrub    = scrubAx(d, A, want);
        ax.scrubChoices = drawableDims(d, {});
        ax.slice    = setdiffDims(ax.scrubChoices, {ax.scrub.name});

    case 'video.fromSegments'
        % THE INVERSE OF A ROW SELECTION, ONCE PER SAMPLE.  Everything image.
        % fromSegments says holds here - the unit is drawn rather than reduced, so
        % the units are kept and never pooled - and the extra fact is the scrub: the
        % dimension it walks is the leaf's OWN time base, taken from the registry and
        % labelled in seconds, never in frame numbers.
        ax.x = mkAx('','','linear',[],A);
        ax.y = mkAx('','','linear',[],A);
        ax.c = mkVar(vLab);
        ax.unitAxis = 'pixel';
        ax.units    = 'points';
        ax.scrub    = scrubAx(d, A, want);
        ax.scrubChoices = drawableDims(d, {});
        ax.slice    = setdiffDims(ax.scrubChoices, {ax.scrub.name});

    otherwise
        error('exPlotRules:plotId','Unknown plot ''%s''.', plotId);
end

ax.kind = kindOf(plotId);
if isempty(ax.y.label) && ~strcmp(ax.kind,'Image') && ~strcmp(ax.kind,'Video')
    ax.y = mkVar(vLab);
end
end

% =====================================================================
function a = scrubAx(d, A, want)
%scrubAx  Which dimension the slider walks.  With two candidates - the frequency
%   and time of a spectrum cube - TIME is the default, because a movie of a wave
%   passing over the tissue is the thing anyone wants to see first; the other is
%   offered in .scrubChoices, and a caller who has offered that choice passes the
%   answer back as WANT.  A frequency sweep and a time sweep of one cube are
%   different questions, so this is a choice and not a preference.
if nargin<3, want = ''; end
dd = drawableDims(d, {});
if isempty(dd), a = mkAx('','','linear',[],A); return; end
if ~isempty(want) && any(strcmp(dd, want))
    nm = char(want);
else
    t = dimOfClass(d,'time');
    if ~isempty(t), nm = t; else, nm = dd{1}; end
end
a = mkAx(nm, defaultLabel(nm), defaultScale(nm), nm, A);
end

% =====================================================================
function dd = drawableDims(d, exclude)
%drawableDims  The leaf's dims that a plot could put on an axis, minus the ones the
%   caller has already spoken for.  A label id never appears - see plotKey.
dd = {};
for i = 1:numel(d.dims)
    if isempty(axisClass(d.dims{i})), continue; end
    if any(strcmp(exclude, d.dims{i})), continue; end
    dd{end+1} = d.dims{i}; %#ok<AGROW>
end
end

function out = intersectDims(d, want)
out = {};
for i = 1:numel(d.dims)
    if any(strcmp(want, d.dims{i})), out{end+1} = d.dims{i}; end %#ok<AGROW>
end
end

function out = setdiffDims(list, drop)
out = list(~ismember(list, drop));
end

function nm = dimOfClass(d, cls)
%dimOfClass  The leaf's own name for one axis class, '' when it has none.
nm = '';
for i = 1:numel(d.dims)
    if strcmp(axisClass(d.dims{i}), cls), nm = d.dims{i}; return; end
end
end

function s = firstNonEmpty(c)
s = '';
for i = 1:numel(c), if ~isempty(c{i}), s = c{i}; return; end, end
end

% ================================================================= AXIS SPECS

function a = mkAx(name, label, scale, fromA, A)
%mkAx  One axis.  When the recording's own registry is to hand its label, scale and
%   coordinate VALUES win over the generic ones - which is the whole reason exFetch
%   passes it, and the reason the two percentile axes cannot be papered over by a
%   constant here.
a = struct('name',char(name), 'label',char(label), 'scale',char(scale), ...
           'values',[], 'log',strcmp(scale,'log'));
if isempty(fromA) || isempty(A) || ~isstruct(A), return; end
key = char(fromA);
if ~isfield(A, key), return; end
r = A.(key);
if ~isempty(r.label), a.label = r.label; end
if ~isempty(r.scale) && ~strcmp(r.scale,'categorical'), a.scale = r.scale; end
a.values = r.values;
a.log    = strcmp(a.scale,'log');
end

function a = mkVar(lab)
%mkVar  The variable itself - the y of a curve, the colour of an image.
a = struct('name','', 'label',char(lab), 'scale','linear', 'values',[], 'log',false);
end

function s = varLabel(d)
%varLabel  What the variable is called, with its unit.  Written for a biologist:
%   exSchema has already spelled the phrase, this only adds the unit.
s = d.label;
if isempty(s), s = d.leaf; end
if ~isempty(d.yUnit), s = sprintf('%s (%s)', s, d.yUnit); end
end

function s = defaultLabel(nm)
switch char(nm)
    case {'time','gsTime','timeWT','timeDWT'}, s = 'time (s)';
    case 'f',                                  s = 'frequency (Hz)';
    case {'pctLevel','pctBin'},                s = 'percentile (%)';
    case 'harmonic',                           s = 'harmonic';
    case 'plane',                              s = 'mask layer';
    case 'line',                               s = 'position along the vessel';
    otherwise,                                 s = char(nm);
end
end

function s = defaultScale(nm)
switch char(nm)
    case 'f',         s = 'log';
    case 'harmonic',  s = 'ordinal';
    case 'plane',     s = 'ordinal';
    otherwise,        s = 'linear';
end
end

function ax = blankAxes(plotId)
e = struct('name','', 'label','', 'scale','linear', 'values',[], 'log',false);
ax = struct('plot',char(plotId), 'kind','', 'x',e, 'y',e, 'c',e, ...
    'family',e, 'scrub',e, 'scrubChoices',{{}}, 'slice',{{}}, 'reduce',{{}}, ...
    'unitAxis','', 'units','pooled');
end
