%exSchema - What each container in a results tree holds, and what to call it.
%
%   S = exSchema() is the whole table, once.  r = exSchema('match',path) is the row
%   that claims one dotted leaf path, resolved for that path; and
%   exSchema('label',path) is the phrase to put in front of a biologist.  It reads
%   NOTHING: no file, no struct, no sizes.  It is the hand-written half of the shape
%   model, and exScan supplies the other half by measuring.
%
%   THIS FILE MIRRORS Core/Vasomotion/assembleVasomotionTree.m, which is the
%   declared single source of truth for the <VSM> layout.  Every producer -
%   runVasomotion's per-segment and per-pixel paths, getMyographVasomotion - builds
%   its sub-tree through that one function, so the field NAMES and the branch
%   structure below are its field names and its branches.  When a leaf is added
%   there, it is added here.
%
%   WHY THERE IS A 'layout' COLUMN AT ALL.  Blind size-matching against the axis
%   registry gets all three of these wrong, and each of them is a real leaf in the
%   reference data:
%
%     per-segment scalars / fVectors / spectrum are ITEM FIRST     [nSeg x ...]
%     per-segment timeVectors are ITEM LAST                        [nT x nSeg]
%     per-pixel is item first and TWO-DIMENSIONAL                  [Y x X x ...]
%
%   (runVasomotion.m:127 for the accumulator, :714-716 for the per-pixel reshapes.)
%   With nT and nSeg both in the thousands, a scanner guessing from position would
%   label vasomotion.sData.timeVectors.VB.amp - [2024 x 1480] - as 1480 samples of
%   2024 segments.  The item axis cannot be inferred; it has to be declared.
%
%   TWO MORE LAYOUTS EXIST BECAUSE THE DATA IS NOT ALWAYS A GRID.
%     xyPairs   propagation lagByRow is [nRow x 2]: column 1 IS the row coordinate
%               and column 2 the lag.  It carries its own x, which is why it is not
%               verified against the line axis - measured on the reference pressure
%               recording, nY is 480 while lagByRow has 475 rows on the first window
%               and 480 on the other two, because a row with no usable correlation
%               contributes no pair.  Verifying it would mark a correct leaf suspect.
%     pairWith  propagation speedCI is [1x2], an interval ON speed rather than a
%               curve.  It rides with its scalar as an error bar and is never
%               offered as a variable of its own.
%
%   THE UNIT FOLLOWS THE SIGNAL, AND SOMETIMES ONLY THE DATA KNOWS IT.  sData and
%   gsData belong to segments, dvsData and dvsDiameter to tracked vessels, ppx to
%   pixels - so one row covers all four and the signal token in the path decides.
%   A myograph <VSM> is the case the schema deliberately does NOT decide: with
%   s.perLine false the leaf holds ONE item (the reference pressure recording:
%   fVectors.ampMean is [1x70]), with it true, one per row across the vessel.  Such
%   a row is returned with unit 'auto' and exScan reads the count off the leaf.
%
%   results.vasomotion.ppxs GETS NO ROW.  Its producer - the runSegmentAveragingTEMP
%   block in runVasomotion - declares it temporary scaffolding to be deleted with
%   that block.  Leaving it out is a decision, not an oversight: if it is wanted it
%   costs one row here.  Until then it falls to exScan's inferred path and appears
%   under Other, which is visible rather than silent.
%
%   ONE FAMILY IS NOT IN THE SPEC'S LIST.  The metrics tables get family 'Metrics'.
%   sMetrics holds tortuosity beside blood flow index beside pulsatility index; it
%   is not flow, not pulsatility and not a map, and filing it under any of them
%   would put a segment's length in a flow menu.  Flagged for the author.
%
%   Syntax:
%      S   = exSchema()
%      r   = exSchema('match', dottedPath)
%      txt = exSchema('label', dottedPath)
%
%   INPUTS
%     dottedPath  one leaf path as exScan builds it, struct-array elements carrying
%                 their index: 'intervals(2).vasomotion.outer.spectrum.amp'.
%
%   OUTPUTS
%     S   struct array, one element per container row: .id .pattern .family .unitBy
%         .unit .layout .dims .phrase .yUnit .leaves .leafOverride .pairedWith
%     r   the row RESOLVED for that path - the same fields, plus .signal .branch
%         .leaf .label - or [] when nothing claims it.  .unit may be 'auto'.
%     txt the human phrase, '' when nothing claims the path.
%
%   EXAMPLE
%      exSchema('label','vasomotion.sData.timeVectors.VB.amp')
%      % 'Vasomotion band - amplitude envelope'
%      exSchema('match','vasomotion.sData.timeVectors.VB.amp').layout   % 'unitLast'
%
%   DEPENDS ON
%     nothing.
%
% See also: exModality, exAxes, exScan, assembleVasomotionTree, runVasomotion,
%           guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function varargout = exSchema(action, varargin)

if nargin==0, varargout{1} = schemaTable(); return; end
switch lower(char(action))
    case 'match', varargout{1} = matchRow(varargin{:});
    case 'label', varargout{1} = labelFor(varargin{:});
    case 'table', varargout{1} = schemaTable();
    otherwise
        error('exSchema:action','Unknown action ''%s''.', char(action));
end
end

% =====================================================================
function r = matchRow(dotted)
%matchRow  The first row whose pattern claims this path, resolved for it.  ORDER
%   MATTERS: the special propagation leaves are listed before the catch-all one, so
%   lagByRow is xyPairs rather than a scalar with 950 mysterious entries.
r = [];
p = char(dotted);
S = schemaTable();
for i = 1:numel(S)
    t = regexp(p, S(i).pattern, 'names', 'once');
    if isempty(t), continue; end
    r = resolve(S(i), t);
    return
end
end

% =====================================================================
function r = resolve(row, tok)
%resolve  One row instance for one path: the signal decides the unit and, for the
%   per-pixel signal, the layout too; the leaf may override the dims.
r = row;
r.signal = tokOr(tok,'sig');
r.branch = tokOr(tok,'branch');
r.leaf   = tokOr(tok,'leaf');

if strcmp(row.unitBy,'signal')
    switch r.signal
        case {'sData','gsData'},        r.unit = 'seg';
        case {'dvsData','dvsDiameter'}, r.unit = 'dvs';
        case 'ppx',                     r.unit = 'pixel'; r.layout = 'pixel2D';
        case {'sMetrics'},              r.unit = 'seg';
        case {'dvsMetrics'},            r.unit = 'dvs';
        otherwise
            % A myograph measure ('outer') or a wire channel ('channel1'): whether
            % it is one vessel or one row per line is a property of the run, and
            % only the leaf's own item count can say.  exScan resolves it.
            r.unit = 'auto';
    end
end

o = row.leafOverride;
if ~isempty(r.leaf) && isvarname(r.leaf) && isstruct(o) && isfield(o, r.leaf)
    r.dims = o.(r.leaf);
end

[tail, r.yUnit] = leafInfo(r);
r.label = composeLabel(r, tail);
end

% =====================================================================
function txt = labelFor(dotted)
r = matchRow(dotted);
if isempty(r), txt = ''; return; end
txt = r.label;
end

% =====================================================================
function txt = composeLabel(r, tail)
%composeLabel  '<what it belongs to> - <what it is>', written for a biologist: no
%   field names, no wrapper names, nothing in brackets.  A leaf whose phrase is
%   empty is named by its container alone, which is how 'Outer diameter' comes out
%   as itself rather than 'Outer diameter - outer'.
head = r.phrase;
switch r.branch
    case 'VB', head = 'Vasomotion band';
    case 'CB', head = 'Comparison band';
end
head = strrep(head, '<signal>', signalPhrase(r.signal));
if isempty(tail), txt = head;
elseif isempty(head), txt = tail;
else, txt = [head ' - ' tail];
end
end

% =====================================================================
function [phrase, yUnit] = leafInfo(r)
%leafInfo  This row's word for one leaf, and the unit that leaf is measured in.
%   Unknown leaves keep their own name: a tree that grew a field this table has not
%   heard of is still shown, spelled the way the producer spelled it.
phrase = r.leaf; yUnit = r.yUnit;
L = r.leaves;
for i = 1:size(L,1)
    if strcmp(L{i,1}, r.leaf)
        phrase = L{i,2};
        if size(L,2)>=3 && ~isempty(L{i,3}), yUnit = L{i,3}; end
        return
    end
end
end

% =====================================================================
function s = signalPhrase(sig)
switch sig
    case 'sData',       s = 'Segment signal';
    case 'dvsData',     s = 'Vessel signal';
    case 'dvsDiameter', s = 'Vessel diameter';
    case 'gsData',      s = 'Segment signal at full time resolution';
    case 'ppx',         s = 'Per-pixel';
    case 'sMetrics',    s = 'Segment';
    case 'dvsMetrics',  s = 'Vessel';
    case 'outer',       s = 'Outer diameter';
    case 'mid',         s = 'Mid-wall diameter';
    case 'inner',       s = 'Inner diameter';
    otherwise,          s = sig;
end
end

% ================================================================= THE TABLE

function S = schemaTable()
%schemaTable  The ~20 container rows of the shape model.  Built once per call and
%   deliberately not cached: it is a few microseconds of struct assembly, and a
%   persistent copy is one more thing that can be stale while a developer edits it.

% The myograph window prefix.  A pressure recording keeps its windows flat; a wire
% one hangs them off the channel axis, and both spellings reach here.
IV  = '(?:channel\(\d+\)\.)?intervals\(\d+\)\.';
% The negative lookahead is how "ppxs gets no row" is ENFORCED rather than merely
% intended: without it the signal token would swallow ppxs like any other name and
% file the temporary scaffolding as a real per-line vasomotion tree.
VSM = ['^(?:' IV ')?vasomotion\.(?<sig>(?!ppxs\.)[A-Za-z]\w*)\.'];
PUL = '^pulsatility\.(?<sig>sData|dvsData|dvsDiameter|ppx)\.';
MYO = ['^' IV];

S = row('rawSegment', '^(?<sig>sData|dvsData)$', 'Flow', ...
        'signal','','unitLast',{'time'}, '<signal>','a.u.', {});
S(end+1) = row('rawDiameter', '^(?<sig>dvsDiameter)$', 'Diameter', ...
        'signal','','unitLast',{'time'}, '<signal>','px', {});
S(end+1) = row('rawGuided', '^(?<sig>gsData)$', 'Flow', ...
        'signal','','unitLast',{'gsTime'}, '<signal>','a.u.', {});

S(end+1) = row('metrics', '^(?<sig>sMetrics|dvsMetrics)\.(?<leaf>.+)$', 'Metrics', ...
        'signal','','unitFirst',{}, '<signal>','', metricsLeaves());

S(end+1) = row('maps', ['^(?<leaf>imgBFI|imgK|imgI|mask|sMap|pMap|cMask|dvsMap|' ...
        'mapType|regionsMask)$'], 'Maps', 'fixed','pixel','pixel2D',{}, '','', mapLeaves());
S(end+1) = row('extended', '^extendedMetrics\.(?<leaf>\w+)$', 'Maps', ...
        'fixed','pixel','pixel2D',{}, '','', mapLeaves());
S(end+1) = row('commonMask', '^(?<leaf>commonMask)$', 'Maps', ...
        'fixed','pixel','pixel2D',{'plane'}, '','', mapLeaves());

% ---- vasomotion: the four branches assembleVasomotionTree builds ---------------
S(end+1) = row('vsmScalars', [VSM 'scalars\.(?<branch>VB|CB)\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitFirst',{}, '','a.u.', vsmScalarLeaves());
S(end).leafOverride = ov('ampPct',{'pctLevel'});

S(end+1) = row('vsmFVectorsBand', [VSM 'fVectors\.(?<branch>VB|CB)\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitFirst',{'f'}, '','a.u.', vsmFVectorLeaves());
S(end).leafOverride = ov('ampMeanPct',{'f','pctBin'});

S(end+1) = row('vsmFVectors', [VSM 'fVectors\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitFirst',{'f'}, 'Frequency spectrum','a.u.', vsmFVectorLeaves());

S(end+1) = row('vsmTimeVectors', [VSM 'timeVectors\.(?<branch>VB|CB)\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitLast',{'timeWT'}, '','a.u.', vsmTimeLeaves());

S(end+1) = row('vsmSpectrum', [VSM 'spectrum\.(?<leaf>amp|phase)$'], 'Vasomotion', ...
        'signal','','unitFirst',{'f','timeDWT'}, 'Time-frequency map','a.u.', vsmSpectrumLeaves());

% ---- pulsatility ---------------------------------------------------------------
S(end+1) = row('pulsScalars', [PUL 'scalars\.(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitFirst',{}, 'Pulse','', pulsScalarLeaves());
S(end).leafOverride = ov('hAmp',{'harmonic'},'hPhase',{'harmonic'});

S(end+1) = row('pulsTimeVectors', [PUL 'timeVectors\.(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitLast',{'time'}, 'Pulse','a.u.', {'fData','averaged pulse','a.u.'});

% ---- the myograph window: diameter, then propagation ---------------------------
S(end+1) = row('myoDiameter', [MYO 'diameter\.(?<sig>outer|mid|inner)$'], 'Diameter', ...
        'fixed','whole','none',{'time'}, '<signal>','px', {});
S(end+1) = row('myoDiameterStats', [MYO 'diameter\.stats\.(?<sig>\w+)\.(?<leaf>\w+)$'], 'Diameter', ...
        'fixed','whole','none',{}, '<signal>','px', statLeaves());

% The three special propagation leaves come FIRST - the catch-all below would
% otherwise claim them and call a 475-row table of lags a scalar.
S(end+1) = row('myoLagByRow', [MYO 'propagation\.(?<sig>\w+)\.(?<leaf>lagByRow)$'], 'Propagation', ...
        'fixed','line','xyPairs',{}, 'Propagation','s', propLeaves());
S(end+1) = row('myoSpeedCI', [MYO 'propagation\.(?<sig>\w+)\.(?<leaf>speedCI)$'], 'Propagation', ...
        'fixed','whole','none',{}, 'Propagation','px/s', propLeaves());
S(end).pairedWith = 'speed';
S(end+1) = row('myoPropMetrics', [MYO 'propagation\.(?<sig>\w+)\.metrics\.(?<leaf>\w+)$'], 'Propagation', ...
        'fixed','whole','none',{}, 'Propagation','', propMetricLeaves());
S(end+1) = row('myoProp', [MYO 'propagation\.(?<sig>\w+)\.(?<leaf>\w+)$'], 'Propagation', ...
        'fixed','whole','none',{}, 'Propagation','', propLeaves());
end

% =====================================================================
function r = row(id, pattern, family, unitBy, unit, layout, dims, phrase, yUnit, leaves)
r = struct('id',id, 'pattern',pattern, 'family',family, 'unitBy',unitBy, ...
    'unit',unit, 'layout',layout, 'dims',{dims}, 'phrase',phrase, 'yUnit',yUnit, ...
    'leaves',{leaves}, 'leafOverride',struct(), 'pairedWith','', ...
    'signal','', 'branch','', 'leaf','', 'label','');
end

function o = ov(varargin)
%ov  The handful of leaves whose dims differ from their container's.
o = struct();
for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
end

function v = tokOr(tok, f)
if isstruct(tok) && isfield(tok,f), v = char(tok.(f)); else, v = ''; end
end

% ================================================================= LEAF WORDS

function L = vsmScalarLeaves()
L = { 'ampMean',        'mean amplitude',                      'a.u.'
      'ampStd',         'amplitude variability',               'a.u.'
      'ampSkew',        'amplitude skewness',                  ''
      'ampPct',         'amplitude percentiles',               'a.u.'
      'fCentMean',      'mean centre frequency',               'Hz'
      'fCentStd',       'centre frequency variability',        'Hz'
      'fSprdMean',      'mean frequency spread',               'Hz'
      'fSprdStd',       'frequency spread variability',        'Hz'
      'shapePeak',      'peak sharpness',                      ''
      'nPeakMean',      'mean number of peaks',                ''
      'nPeakStd',       'variability in the number of peaks',  ''
      'durFlareMean',   'mean burst duration',                 's'
      'durFlareStd',    'burst duration variability',          's'
      'durSilenceMean', 'mean quiet duration',                 's'
      'durSilenceStd',  'quiet duration variability',          's'
      'ampFlareMean',   'mean amplitude during bursts',        'a.u.'
      'ampFlareStd',    'amplitude variability during bursts', 'a.u.'
      'ampSilenceMean', 'mean amplitude when quiet',           'a.u.'
      'ampSilenceStd',  'amplitude variability when quiet',    'a.u.' };
end

function L = vsmFVectorLeaves()
L = { 'ampMean',    'mean amplitude',              'a.u.'
      'ampStd',     'amplitude variability',       'a.u.'
      'ampSkew',    'amplitude skewness',          ''
      'ampMeanPct', 'mean amplitude by percentile','a.u.'
      'ampFlare',   'amplitude during bursts',     'a.u.'
      'ampSilence', 'amplitude when quiet',        'a.u.' };
end

function L = vsmTimeLeaves()
L = { 'amp',       'amplitude envelope',   'a.u.'
      'fCent',     'centre frequency',     'Hz'
      'fSprd',     'frequency spread',     'Hz'
      'nPeak',     'number of peaks',      ''
      'maskFlare', 'burst periods',        ''
      'rData',     'reconstructed signal', 'a.u.' };
end

function L = vsmSpectrumLeaves()
L = { 'amp',   'amplitude', 'a.u.'
      'phase', 'phase',     'rad' };
end

function L = pulsScalarLeaves()
L = { 'psMin',      'minimum',                      'a.u.'
      'psMax',      'maximum',                      'a.u.'
      'psMean',     'mean',                         'a.u.'
      'psStd',      'variability',                  'a.u.'
      'psTimeMin',  'time of the minimum',          's'
      'psTimeMax',  'time of the maximum',          's'
      'psPI',       'pulsatility index',            ''
      'psRI',       'resistance index',             ''
      'psSymRatio', 'symmetry ratio',               ''
      'psR2',       'fit quality',                  ''
      'pdMin',      'diameter minimum',             'px'
      'pdMax',      'diameter maximum',             'px'
      'pdMean',     'diameter mean',                'px'
      'pdStd',      'diameter variability',         'px'
      'pdTimeMin',  'time of the diameter minimum', 's'
      'pdTimeMax',  'time of the diameter maximum', 's'
      'pdPI',       'diameter pulsatility index',   ''
      'pdRI',       'diameter resistance index',    ''
      'pdSymRatio', 'diameter symmetry ratio',      ''
      'pdR2',       'diameter fit quality',         ''
      'hAmp',       'harmonic amplitudes',          'a.u.'
      'hPhase',     'harmonic phases',              'rad' };
end

function L = statLeaves()
L = { 'mean',             'mean',                        'px'
      'std',              'variability',                 'px'
      'min',              'minimum',                     'px'
      'max',              'maximum',                     'px'
      'amplitude',        'peak-to-peak change',         'px'
      'measuredFraction', 'fraction of frames measured', ''
      'validFraction',    'fraction of frames accepted', '' };
end

function L = propLeaves()
L = { 'domFreq',         'dominant frequency',           'Hz'
      'speed',           'speed',                        'px/s'
      'speedCI',         'speed confidence interval',    'px/s'
      'slope',           'lag slope',                    ''
      'R2',              'fit quality',                  ''
      'pValue',          'significance',                 ''
      'nRows',           'rows used',                    ''
      'confidence',      'confidence',                   ''
      'belowResolution', 'below the timing resolution',  ''
      'lagByRow',        'lag along the vessel',         's' };
end

function L = propMetricLeaves()
L = { 'medianCorr',      'median correlation',           ''
      'R2',              'fit quality',                  ''
      'pValue',          'significance',                 ''
      'rowFraction',     'fraction of rows used',        ''
      'totalLagSamples', 'total lag',                    ''
      'belowResolution', 'below the timing resolution',  ''
      'nRows',           'rows used',                    '' };
end

function L = mapLeaves()
L = { 'imgBFI',      'blood flow map',            'a.u.'
      'imgK',        'contrast map',              'a.u.'
      'imgI',        'intensity image',           'a.u.'
      'imgStdBFI',   'blood flow variability map','a.u.'
      'mask',        'analysed area',             ''
      'commonMask',  'shared analysed area',      ''
      'sMap',        'segment map',               ''
      'pMap',        'parenchyma cell map',       ''
      'cMask',       'tissue category map',       ''
      'dvsMap',      'tracked vessel map',        ''
      'mapType',     'vessel type map',           ''
      'regionsMask', 'selected regions',          '' };
end

function L = metricsLeaves()
%metricsLeaves  The columns of sMetrics / dvsMetrics.  CLR is the segment's path
%   length over its end-to-end distance (runDynamicSegmentation.m:167) - 1 for a
%   straight vessel, 2 for a coil - so it is tortuosity and is named as such.
L = { 'idx',            'index',                    ''
      'category',       'tissue category',          ''
      'length',         'length',                   'px'
      'CLR',            'tortuosity',               ''
      'diameter',       'diameter',                 'px'
      'std(diameter)',  'diameter variability',     'px'
      'area',           'area',                     'px'
      'nearestVesIdx',  'nearest vessel',           ''
      'RegID',          'region',                   ''
      'RegOverlap',     'overlap with the region',  ''
      'R2',             'fit quality',              ''
      'overlapMask',    'overlap with the mask',    ''
      'overlapSelf',    'self overlap',             ''
      'BFI',            'blood flow index',         'a.u.'
      'std(BFI)',       'blood flow variability',   'a.u.'
      'psPI',           'pulsatility index',        ''
      'psRI',           'resistance index',         ''
      'psMean',         'mean over the pulse',      'a.u.'
      'psTimeMin',      'time of the minimum',      's'
      'psTimeMax',      'time of the maximum',      's'
      'psSymRatio',     'symmetry ratio',           ''
      'pdPI',           'diameter pulsatility index',''
      'pdRI',           'diameter resistance index', ''
      'pdMean',         'mean diameter over the pulse','px'
      'pdTimeMin',      'time of the diameter minimum','s'
      'pdTimeMax',      'time of the diameter maximum','s'
      'pdSymRatio',     'diameter symmetry ratio',  ''
      'typeConfidence', 'vessel type confidence',   '' };
end
