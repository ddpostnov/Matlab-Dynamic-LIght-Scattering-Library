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
%               offered as a variable of its own - see the companion table below,
%               which is what says so.
%
%   A LEAF MAY RIDE WITH ANOTHER LEAF, IN A NAMED ROLE.  Four of the propagation
%   leaves are not quantities anyone wants a menu entry for.  They are parts of the
%   answer another leaf already gives:
%
%     speedCI          an INTERVAL on speed       -> an error bar on its box
%     slope            the FIT of lagByRow        -> the line through the scatter
%     confidenceLevel  \  the CAPTION of          -> the sentence beside the plot
%     confidenceText   /  lagByRow
%
%   companionTable() is the ONE place that says so - family, host leaf, rider leaf,
%   role - and both halves of the mechanism are read off it.  A RIDER resolves with
%   .pairedWith naming its host and .pairRole naming what it is for, so a variable
%   list that filters on .pairedWith drops all four; exSchema('companions',...)
%   answers the other direction, so exFetch can hand a host's riders to a renderer
%   and no renderer ever reaches into a results tree for an error bar or a caption.
%
%   TWO MORE ROLES ARE DECLARED ON A ROW INSTEAD, because their host has no leaf
%   name to key a table row on.
%
%     'pooled'  A myograph window stores its diameter twice on purpose -
%               .lines.<measure>, one value per row across the vessel per frame, and
%               .<measure>, those same numbers averaged along the vessel - and the
%               second is the POOLED form of the first.  The leaf table cannot say
%               so, because the measure is the signal token of both paths and neither
%               has a leaf name; so the row that claims the trace names its host with
%               .pairedWith = 'lines.<sig>' and .pairRole = 'pooled'.  ONE quantity,
%               ONE menu entry: the per-line leaf, whose 'Units are: pooled'
%               reproduces the stored trace because a mean over lines with equal
%               weights is exactly what wrote it.
%     'mask'    The same window's .measured and .valid are the QUALITY of that
%               measurement - which points were detected rather than interpolated,
%               which frames had both walls in the field of view.  Nobody wants a
%               menu entry for either, and both belong to the diameter beside them,
%               so they ride on it: .pairedWith = 'lines', .pairRole = 'mask'.  There
%               is no <sig> to interpolate here, and that is not an omission - ONE
%               .measured is shared by all three measures (it is the wall-centre
%               one's), so naming a single measure as its host would be a claim the
%               data does not make.
%
%   ONLY 'pooled' IS CONDITIONAL, and exScan is what resolves it.  A run with
%   s.keepArrays false stores the trace and its statistics alone, and there the trace
%   is the only form of the quantity in the file - so exScan drops the declaration
%   when the host is not there.  The other five roles are not conditional: a producer
%   that writes a fit, a caption or a quality flag writes the leaf it is about in the
%   same breath - myographDiameterBranch adds .lines and .measured together and
%   removes them together.
%
%   THE CAPTION IS NOT A DESCRIPTOR, WHICH IS WHY THE TABLE IS KEYED ON NAMES.
%   confidenceLevel and confidenceText are char, and exScan skips text leaves as
%   labels - so there is no descriptor to hang a declaration on and no size to
%   match.  A leaf NAME reaches both, and the host's own path supplies the
%   container it is a sibling in.  The FAMILY scopes the table so that a leaf called
%   'slope' somewhere else is not silently swallowed by the propagation rule.
%
%   A BAND TOKEN DOES NOT MEAN THE SAME THING IN EVERY CONTAINER, which is why the
%   phrase in front of it belongs to the ROW and not to composeLabel.  Under scalars
%   and timeVectors, VB names the frequency band the number was measured in, and
%   'Vasomotion band - mean amplitude' is exactly right.  Under fVectors it does not:
%
%     getVasomotionMetrics.m:522    spctFlare = mean(wtts(:,mask),2)
%
%   wtts is the WHOLE wavelet spectrum and mask selects the time samples during which
%   that band was bursting - so fVectors.VB.ampFlare is the full spectrum averaged
%   over the moments the vasomotion band was flaring, and 'Vasomotion band - amplitude
%   during bursts' reads as a contradiction of its own [nUnit x nF] shape.  All three
%   of that container's leaves share the property (ampSilence takes ~mask,
%   ampMeanPct :443-444 takes one percentile bin of the same envelope), so the row
%   carries the head 'Whole frequency spectrum' and its own leaf words, in which
%   <band> is what chose the time.
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
%   TWO FAMILIES ARE NOT IN THE SPEC'S LIST.  The metrics tables get family
%   'Metrics'.  sMetrics holds tortuosity beside blood flow index beside pulsatility
%   index; it is not flow, not pulsatility and not a map, and filing it under any of
%   them would put a segment's length in a flow menu.  Flagged for the author.
%
%   And a WIRE myograph's recorded samples get family 'Recording'.  A LabChart
%   channel measures force, so it is neither Flow nor Diameter, and it is the raw
%   trace rather than anything derived from it - which is why the family sorts FIRST
%   in the menu, where a speckle recording puts Flow.  It is the one family whose
%   leaves are the recording itself and not a result of analysing it.
%
%   THE STIMULUS RESPONSE IS A THIRD, AND IT BRINGS AN AXIS WITH IT.  results.nvc
%   gets family 'Response': runNVC cuts every stimulus repetition out of the whole
%   recording and measures each one on its own, so esMetrics.<signal>.<metric> is
%   [nSeg x nEp] - segment against REPETITION.  Every family before it was
%   (segment), (segment, time), (Y, X) or (Y, X, time); this one adds the EPOCH
%   axis, registered by exAxes and drawable exactly as time is, because whether a
%   response fades over a session is a slope along it.  It has TWO containers -
%   esMetrics per segment and ppx per pixel, the same markers measured on two sets
%   of units - and the prefix (ns / nd / ng) is captured as the BRANCH, see the rows
%   for why.  There is no third: the step fits no model, so nothing is measured on
%   the average of the repetitions.
%
%   AND esMetrics IS LEGITIMATELY IN ONE OF TWO STATES, which is why no row here
%   tests for either.  An ordinary product carries results.nvc.epochStart, so the
%   epoch axis exists and every esMetrics leaf is [nSeg x nEp]; a REPRESENTATIVE-
%   repetition product - runNVC's optional collapse - has no epochStart, so there is
%   no epoch axis, and the same leaves resolve as ordinary per-segment scalars
%   through the same rows.  The declaring leaf is what makes that free.
%
%   THE PULSATILITY PREFIX IS A BRANCH IN THE SAME WAY, AND IT NAMES THREE QUANTITIES.
%   ps is a pulsatile FLOW, pd a pulsatile DIAMETER and pv a pulsatile PLASMA VOLUME -
%   the fluorescence branch's, where an intravascular tracer makes the intensity
%   proportional to the labelled plasma in the light path.  No two of them may be
%   pooled, so each has its own row, its own head and its own unit, and one word table
%   serves all three.  The only leaves the producers leave unprefixed are the fitted
%   model's hAmp / hPhase, and they are therefore the one pair this file cannot tell
%   apart between two branches - see the row that claims them.
%
%   THE BOLUS TRANSIT IS A FOURTH, AND IT BRINGS NO AXIS AT ALL.  results.ctth gets
%   family 'Transit'.  It is the stimulus response's shape with one dimension fewer:
%   the same two containers (esMetrics per segment, ppx per pixel), the same captured
%   prefix (bs / bd), the same markers-plus-confidence table - but a recording has ONE
%   bolus, so every leaf is one number per unit and the segment axis the explorer
%   already has carries all of it.  A 'bolus' axis would have been an axis of length
%   one on every product that has it.
%
%   It also has a THIRD kind of leaf the response family has no equivalent of: three
%   numbers that belong to the RECORDING - the smallest spread of transit times it
%   could resolve, whether it could resolve one at all, and how much of every absolute
%   time was the injection.  They are not a summary of anybody's segments; they are
%   what the per-segment markers have to be read against, which is why they are
%   variables rather than a note on a report page.
%
%   THE WALL MOTION IS A FIFTH FAMILY AND IT IS NOT A TREE AT ALL.  results.wallMotion
%   is FLAT - one struct of [nRow x 1] arrays indexed exactly like results.sMetrics -
%   and it sits on the same product as results.pulsatility, which branches by signal
%   and then by container.  They are SIBLINGS with different conventions, not one tree
%   with two shapes, and each is read by its own rows.  What it holds is the EVIDENCE
%   behind the eight wm* columns of the metrics table: the same measurements before the
%   three gates, the floor each was judged against, and which gate refused a row.
%
%   AND THE VASCULAR DENSITY IS A SIXTH, WHOSE UNIT IS AN AREA.  results.topology.
%   metrics is one row per analysed area - the whole field first, then one per drawn
%   region - so its unit is neither a segment nor a pixel, and the row index is a LABEL
%   ID in the same sense a segment index is: there is no "density against region
%   number" plot, so no axis is registered for it.  That is the bolus argument again -
%   a coordinate nothing can be walked along is not a coordinate - and the rows are
%   declared item FIRST so that a box carries one point per area rather than reading
%   the first row and calling it the answer.
%
%   ITS TWO HISTOGRAMS GET NO ROW, and that is a decision.  results.topology.calibre
%   and .evd are [nScope x nBin] counts against EDGES that are one longer by
%   construction - the pctLevel / pctBin trap, in a second place - so drawing them
%   would cost two more coordinates registered for one plot each, which the step's own
%   report page already draws.  They fall to exScan's inferred path and appear under
%   Other, which is visible rather than silent; if they are wanted the axes are the
%   seven-file edit and this is where it starts.
%
%   TWO PRODUCERS NAME THEIR LEAVES FROM A SETTING, and they are the only two: the
%   timing markers are T<p>, one per entry of s.nvcAreaPcts or s.ctthLevelPcts, so a
%   recording processed with [5 25 50 75 95] carries five leaves this file cannot list.
%   A row whose leaf names are derived that way declares a .leafRule - a pattern and the
%   words to build from it - which is consulted when the literal table has no entry, so
%   a level nobody anticipated is still spelled out rather than shown as 'T25'.  THE TWO
%   RULES SAY DIFFERENT THINGS about the same-looking name, which is exactly why each
%   row names its own: an NVC T50 is when half the accumulated change had been
%   delivered, a bolus T50 is when the vessel had filled half way to its plateau.
%
%   Syntax:
%      S   = exSchema()
%      r   = exSchema('match', dottedPath)
%      txt = exSchema('label', dottedPath)
%      K   = exSchema('companions', family, leafName)
%
%   INPUTS
%     dottedPath  one leaf path as exScan builds it, struct-array elements carrying
%                 their index: 'intervals(2).vasomotion.outer.spectrum.amp'.
%     family      one of the family names a row carries, e.g. 'Propagation'.
%     leafName    the HOST leaf's own field name, e.g. 'lagByRow'.
%
%   OUTPUTS
%     S   struct array, one element per container row: .id .pattern .family .unitBy
%         .unit .layout .dims .phrase .tail .yUnit .leaves .leafRule .leafOverride
%         .pairedWith .pairRole
%     r   the row RESOLVED for that path - the same fields, plus .signal .branch
%         .leaf .label - or [] when nothing claims it.  .unit may be 'auto', and
%         .pairedWith / .pairRole are filled here rather than on the raw row,
%         because whether a leaf is a rider is usually a property of the LEAF; a row
%         that declares them says it of every leaf it claims, with <sig> standing for
%         the signal token that matched.
%     txt the human phrase, '' when nothing claims the path.
%     K   struct array .leaf .role - the leaves that ride with that host, in table
%         order.  Empty when it has none.
%
%   EXAMPLE
%      exSchema('label','vasomotion.sData.timeVectors.VB.amp')
%      % 'Vasomotion band - amplitude envelope'
%      exSchema('label','vasomotion.sData.fVectors.VB.ampFlare')
%      % 'Whole frequency spectrum - while the vasomotion band was bursting'
%      exSchema('match','vasomotion.sData.timeVectors.VB.amp').layout   % 'unitLast'
%      {exSchema('companions','Propagation','lagByRow').role}  % {'fit','caption','caption'}
%      exSchema('label','nvc.esMetrics.sData.nsT90')
%      % 'Segment signal - time by which 90 % of the change had accumulated'
%      exSchema('label','pulsatility.sData.scalars.pvPI')
%      % 'Pulsatile plasma volume - pulsatility index'
%      exSchema('label','pulsatility.sData.scalars.psPI')
%      % 'Pulsatile flow - pulsatility index'
%
%   DEPENDS ON
%     nothing.
%
% See also: exModality, exAxes, exScan, assembleVasomotionTree, runVasomotion,
%           runNVC, runIntensityPulsatility, runMotionEnhancement,
%           runTopologyAnalysis, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function varargout = exSchema(action, varargin)

if nargin==0, varargout{1} = schemaTable(); return; end
switch lower(char(action))
    case 'match',      varargout{1} = matchRow(varargin{:});
    case 'label',      varargout{1} = labelFor(varargin{:});
    case 'table',      varargout{1} = schemaTable();
    case 'companions', varargout{1} = companionsOf(varargin{:});
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

% WHETHER THIS LEAF IS A RIDER is usually a property of the leaf, not of the
% container: the speed and its interval are claimed by two different rows of the same
% container, and the caption is claimed by the catch-all.  So it is resolved here,
% off the one companion table, rather than assigned onto a row by hand.
[r.pairedWith, r.pairRole] = companionOf(row.family, r.leaf);

% AND SOMETIMES A WHOLE CONTAINER RIDES.  diameter.<measure> is the pooled-over-lines
% form of diameter.lines.<measure> - one quantity in two containers - and there is no
% leaf name to key that on, because the measure is the SIGNAL token in both paths.  A
% ROW may therefore declare it, naming its host relative to its own container with
% <sig> standing for whichever measure matched.  The leaf table still wins: a row-
% level declaration says what is true of every leaf the row claims, and a named leaf
% saying otherwise is the more specific statement.
if isempty(r.pairedWith) && ~isempty(row.pairedWith)
    r.pairedWith = strrep(row.pairedWith, '<sig>', r.signal);
    r.pairRole   = row.pairRole;
end

[tail, r.yUnit] = leafInfo(r);
r.label = composeLabel(r, tail);
end

% ============================================================== THE COMPANIONS

function C = companionTable()
%companionTable  WHICH LEAF RIDES WITH WHICH, AND IN WHAT ROLE.  Four columns:
%   the family it applies in, the HOST leaf, the RIDER leaf, and the role the rider
%   plays for the host.  This is the whole declaration - the rider's "do not offer
%   me as a variable" flag and the host's "here is what comes with me" list are both
%   read off it, so the two can never disagree.
%
%   THE ROLES ARE WHAT A RENDERER SWITCHES ON, so each is a word about the DRAWING
%   and not about the number:
%     'interval'  a [lo hi] pair around a scalar   -> an error bar
%     'fit'       the slope of a fitted line       -> the line through the scatter
%     'caption'   a finished sentence for a reader -> shown beside the plot, never
%                 paraphrased and never recomputed
%     'mask'      the QUALITY of the number beside it - which points were really
%                 measured, which rule rejected which repetition
%   ('pooled' is the fifth, and the only one no renderer switches on: the pooled form
%   is what the host's own units control already produces, so the rider is simply not
%   offered.  It is declared on a row rather than here - see the file header.)
%
%   A ROLE MAY BE DECLARED EITHER WAY, and which one is not a style choice: the
%   myograph's .measured / .valid ride on diameter.lines, which has no leaf name to
%   key a row on, so they are declared on their ROWS; the propagation speed's
%   confidence interval rides on speedCI, which does have one, so it is declared
%   HERE and the mechanism is the ordinary one.
%
%   THE FAMILY IS PART OF THE KEY.  A leaf called 'slope' in some other family is a
%   quantity of its own until somebody says otherwise here.
%
%   NO RESPONSE LEAF RIDES ON ANOTHER.  The retired step's per-rule rejection
%   breakdown did - it was [rule x repetition] beside a per-repetition flag - and the
%   model-free step has no such pair: how much of a (segment, repetition) can be
%   believed is two ordinary numbers, Conf and ConfMin, each a variable in its own
%   right.
C = { 'Propagation', 'speed',      'speedCI',         'interval'
      'Propagation', 'lagByRow',   'slope',           'fit'
      'Propagation', 'lagByRow',   'confidenceLevel', 'caption'
      'Propagation', 'lagByRow',   'confidenceText',  'caption'  };
end

function [host, role] = companionOf(family, leaf)
%companionOf  Is this leaf a rider, and on what?  '' and '' when it is a quantity in
%   its own right, which is every leaf but the handful named above.
host = ''; role = '';
if isempty(leaf) || isempty(family), return; end
C = companionTable();
for i = 1:size(C,1)
    if strcmp(C{i,1},family) && strcmp(C{i,3},leaf)
        host = C{i,2}; role = C{i,4}; return
    end
end
end

function K = companionsOf(family, hostLeaf)
%companionsOf  The other direction: what rides with this leaf.  exFetch resolves the
%   names against the host's own container, so a file that does not carry one simply
%   contributes nothing.
K = struct('leaf',{},'role',{});
if nargin<2 || isempty(hostLeaf) || isempty(family), return; end
C = companionTable();
for i = 1:size(C,1)
    if strcmp(C{i,1},char(family)) && strcmp(C{i,2},char(hostLeaf))
        K(end+1) = struct('leaf',C{i,3},'role',C{i,4}); %#ok<AGROW>
    end
end
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
%
%   THE BAND FILLS THE HEAD ONLY WHEN THE ROW LEFT IT EMPTY.  A row that spells out
%   its own head has said that the band is not what the leaf belongs to - see the
%   fVectors note in the file header - and names the band in its leaf words instead,
%   through <band>.
head = r.phrase;
if isempty(head)
    head = bandPhrase(r.branch);
end
head = strrep(head, '<signal>', signalPhrase(r.signal));
tail = strrep(tail, '<band>',   bandNoun(r.branch));
if isempty(tail), txt = head;
elseif isempty(head), txt = tail;
else, txt = [head ' - ' tail];
end
end

% =====================================================================
function s = bandPhrase(branch)
%bandPhrase  The band as the head of a label: what the number was measured in.
switch branch
    case 'VB', s = 'Vasomotion band';
    case 'CB', s = 'Comparison band';
    otherwise, s = '';
end
end

function s = bandNoun(branch)
%bandNoun  The band inside a sentence, where it is what chose the time samples
%   rather than what the quantity belongs to.
switch branch
    case 'VB', s = 'the vasomotion band';
    case 'CB', s = 'the comparison band';
    otherwise, s = 'the band';
end
end

% =====================================================================
function [phrase, yUnit] = leafInfo(r)
%leafInfo  This row's word for one leaf, and the unit that leaf is measured in.
%   Unknown leaves keep their own name: a tree that grew a field this table has not
%   heard of is still shown, spelled the way the producer spelled it.  A container
%   whose pattern captures NO leaf at all falls back to the row's own .tail, which
%   is how a whole container that holds one quantity still says what that quantity
%   is instead of repeating its signal.
%
%   THE LITERAL TABLE WINS OVER THE RULE, so a producer that derives most of a family
%   from a setting can still name one member of it explicitly.
phrase = r.leaf; yUnit = r.yUnit;
if isempty(phrase), phrase = r.tail; end
L = r.leaves;
for i = 1:size(L,1)
    if strcmp(L{i,1}, r.leaf)
        phrase = L{i,2};
        if size(L,2)>=3 && ~isempty(L{i,3}), yUnit = L{i,3}; end
        return
    end
end
R = r.leafRule;
if ~isempty(R) && ~isempty(r.leaf) && ~isempty(regexp(r.leaf, R{1}, 'once'))
    phrase = regexprep(r.leaf, R{1}, R{2});
    if numel(R)>=3 && ~isempty(R{3}), yUnit = R{3}; end
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

S(end+1) = row('maps', ['^(?<leaf>imgBFI|imgK|imgI|imgBackground|mask|sMap|pMap|cMask|' ...
        'dvsMap|mapType|mapTree|regionsMask)$'], 'Maps', 'fixed','pixel','pixel2D',{}, ...
        '','', mapLeaves());
S(end+1) = row('extended', '^extendedMetrics\.(?<leaf>\w+)$', 'Maps', ...
        'fixed','pixel','pixel2D',{}, '','', mapLeaves());
S(end+1) = row('commonMask', '^(?<leaf>commonMask)$', 'Maps', ...
        'fixed','pixel','pixel2D',{'plane'}, '','', mapLeaves());

% ---- vasomotion: the four branches assembleVasomotionTree builds ---------------
S(end+1) = row('vsmScalars', [VSM 'scalars\.(?<branch>VB|CB)\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitFirst',{}, '','a.u.', vsmScalarLeaves());
S(end).leafOverride = ov('ampPct',{'pctLevel'});

% THE BAND HERE CHOSE THE TIME, IT DID NOT NARROW THE FREQUENCY (file header).  The
% head says so, and the leaf words name the band as the thing that did the choosing.
S(end+1) = row('vsmFVectorsBand', [VSM 'fVectors\.(?<branch>VB|CB)\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitFirst',{'f'}, 'Whole frequency spectrum','a.u.', vsmFVectorBandLeaves());
S(end).leafOverride = ov('ampMeanPct',{'f','pctBin'});

S(end+1) = row('vsmFVectors', [VSM 'fVectors\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitFirst',{'f'}, 'Frequency spectrum','a.u.', vsmFVectorLeaves());

S(end+1) = row('vsmTimeVectors', [VSM 'timeVectors\.(?<branch>VB|CB)\.(?<leaf>\w+)$'], 'Vasomotion', ...
        'signal','','unitLast',{'timeWT'}, '','a.u.', vsmTimeLeaves());

S(end+1) = row('vsmSpectrum', [VSM 'spectrum\.(?<leaf>amp|phase)$'], 'Vasomotion', ...
        'signal','','unitFirst',{'f','timeDWT'}, 'Time-frequency map','a.u.', vsmSpectrumLeaves());

% ---- pulsatility ---------------------------------------------------------------
% THE PREFIX IS THE BRANCH HERE TOO, and it carries THREE physical quantities rather
% than the response family's two.  ps is a pulsatile FLOW, pd a pulsatile DIAMETER and
% pv a pulsatile PLASMA VOLUME - the fluorescence branch's, where the tracer is
% intravascular and the intensity is proportional to the labelled plasma in the light
% path.  No two of them may be pooled: one is a speed, one a width and one a volume,
% and pd is quantised at a quarter of a pixel while pv and the wall motion beside it
% are not.  So each gets its own row, its own head and its own unit, and one word
% table serves all three - which is the ns/nd argument above, said of a third branch.
%
% THE HEAD IS WHERE THAT REACHES A READER.  'Pulsatile plasma volume - pulsatility
% index' cannot be mistaken for 'Pulsatile flow - pulsatility index', and because the
% branch is captured, guiExplore>varIdOf treats the two as different variables even
% before anybody reads the words.
S(end+1) = row('pulsScalars', [PUL 'scalars\.(?<branch>ps)(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitFirst',{}, 'Pulsatile flow','', pulsScalarLeaves('a.u.'));
S(end+1) = row('pulsScalarsDiameter', [PUL 'scalars\.(?<branch>pd)(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitFirst',{}, 'Pulsatile diameter','', pulsScalarLeaves('px'));
S(end+1) = row('pulsScalarsVolume', [PUL 'scalars\.(?<branch>pv)(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitFirst',{}, 'Pulsatile plasma volume','', pulsScalarLeaves('a.u.'));

% AND THE HARMONIC COEFFICIENTS ARE THE ONE PAIR THE PRODUCERS LEAVE UNPREFIXED, so
% this row cannot say which quantity they are the harmonics OF: both pulsatility
% wrappers write scalars.hAmp / scalars.hPhase whatever the markers beside them are
% called.  It is the last row of the four on purpose - the three above claim their own
% leaves first - and it keeps the neutral head.  A speckle hAmp and a fluorescence one
% therefore resolve as ONE variable, which is a hole in the prefix rule and is in the
% PRODUCER rather than here: naming them pshAmp / pvhAmp would close it, at the cost of
% a re-processing of every product that carries a fitted model.  The export closes its
% half by tagging the column with the prefix its siblings carry.
S(end+1) = row('pulsHarmonics', [PUL 'scalars\.(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitFirst',{}, 'Pulse','', pulsHarmonicLeaves());
S(end).leafOverride = ov('hAmp',{'harmonic'},'hPhase',{'harmonic'});

S(end+1) = row('pulsTimeVectors', [PUL 'timeVectors\.(?<leaf>\w+)$'], 'Pulsatility', ...
        'signal','','unitLast',{'time'}, 'Pulse','a.u.', {'fData','averaged pulse','a.u.'});

% ---- the stimulus response: one number per segment PER REPETITION ---------------
% THE PREFIX IS THE BRANCH, and capturing it as one is what keeps a diameter response
% out of a flow menu.  runNVC writes the prefixed name into the tree - nsPeak on a
% flow trace, ndPeak on the tracked diameter, ngPeak on the guided contrast - and the
% data model's rule is that an nd* number must never be pooled with an ns* one,
% because one is a length in pixels and the other a flow index.  guiExplore>varIdOf
% already treats .branch as part of what makes two leaves the same variable, so
% capturing the prefix there enforces the rule structurally AND leaves <leaf> as the
% bare quantity - which is why one word table serves all three branches instead of
% three copies of it differing only in a two-letter prefix.  (bandPhrase / bandNoun
% never see these: both are consulted only when a row leaves its head empty or writes
% <band> into a leaf phrase, and these rows do neither.)
%
% THE HEAD IS THE SIGNAL AND NOTHING MORE, and that is a change of 06-Aug-2026 with a
% reason.  It used to end in 'repetition by repetition', to tell esMetrics apart from
% the fit of the repetition-average that sat beside it; there is no fit and no second
% container now, so the phrase distinguishes nothing - and it would be a false
% statement on a representative-repetition product, where the same leaves are one
% number per segment.  The family says these are stimulus responses and the leaf
% words say what each one is.
%
% THE TIMING LEAVES ARE NAMED FROM s.nvcAreaPcts, so they are matched by a rule
% rather than listed: see .leafRule below and the file header.
S(end+1) = row('nvcTrials', ...
        '^nvc\.esMetrics\.(?<sig>sData|gsData|dvsData)\.(?<branch>ns|ng)(?<leaf>\w+)$', ...
        'Response', 'signal','','unitFirst',{'epoch'}, ...
        '<signal>','', nvcLeaves('a.u.'));
S(end).leafRule = nvcTimeRule();
S(end+1) = row('nvcTrialsDiameter', ...
        '^nvc\.esMetrics\.(?<sig>dvsDiameter)\.(?<branch>nd)(?<leaf>\w+)$', ...
        'Response', 'signal','','unitFirst',{'epoch'}, ...
        '<signal>','', nvcLeaves('px'));
S(end).leafRule = nvcTimeRule();

% Per-pixel is the SAME MARKERS on a second set of units, so it takes the same word
% table.  It carries the two confidence numbers beside whichever markers were asked
% for and no Trust cube: the numbers are what a reader looks at, and a logical map
% adds nothing a threshold on the confidence does not already say.
S(end+1) = row('nvcPpx', '^nvc\.(?<sig>ppx)\.(?<branch>ns)(?<leaf>\w+)$', 'Response', ...
        'signal','','pixel2D',{'epoch'}, ...
        '<signal>','', nvcLeaves('a.u.'));
S(end).leafRule = nvcTimeRule();

% ---- what happened to each repetition, at the level of the whole recording -------
% These are the numbers the recording-wide decision was made on, one per repetition:
% how much of the segmented field responded coherently, and how many seconds of
% frames the camera did not deliver.  Neither is per segment - the per-segment answer
% is Conf / ConfMin / Trust inside esMetrics - which is why they are a row of their
% own with the whole recording as the unit.
S(end+1) = row('nvcQuality', '^nvc\.(?<leaf>areaFrac|timeLoss)$', 'Response', ...
        'fixed','whole','none',{'epoch'}, 'Repetition quality','', nvcQualityLeaves());

S(end+1) = row('nvcEpochTrust', '^nvc\.(?<leaf>epochTrust)$', 'Response', ...
        'fixed','whole','none',{'epoch'}, '','', nvcEpochLeaves());

S(end+1) = row('nvcStimulus', '^nvc\.(?<leaf>stimulus)$', 'Response', ...
        'fixed','whole','none',{'time'}, '','', nvcEpochLeaves());

% ---- the bolus transit: results.ctth ------------------------------------------
% A BOLUS RECORDING HAS ONE BOLUS, so unlike the stimulus response there is no epoch
% axis and every leaf below is ONE NUMBER PER UNIT.  That is why this family brought no
% new coordinate with it: the segment axis the explorer already has carries all of it.
%
% THE PREFIX IS CAPTURED FOR THE REASON ns/nd IS.  bs is a transit measured on an
% intensity trace and bd one measured on the tracked diameter, and the two must never
% pool - one is a time read off a dye curve and the other a time read off a width.
% Capturing it in <branch> makes guiExplore>varIdOf treat them as different variables
% structurally, and leaves <leaf> as the bare quantity so one word table serves both.
%
% THE TIMING LEAVES ARE NAMED FROM s.ctthLevelPcts, so they are matched by a rule -
% see ctthTimeRule, and note that it says something DIFFERENT from the NVC one: here
% the level is a share of the plateau the tracer reached, not of an accumulated change.
S(end+1) = row('ctthSegments', ...
        '^ctth\.esMetrics\.(?<sig>sData|dvsData)\.(?<branch>bs)(?<leaf>\w+)$', ...
        'Transit', 'signal','','unitFirst',{}, '<signal>','', ctthLeaves('a.u.'));
S(end).leafRule = ctthTimeRule();
S(end+1) = row('ctthSegmentsDiameter', ...
        '^ctth\.esMetrics\.(?<sig>dvsDiameter)\.(?<branch>bd)(?<leaf>\w+)$', ...
        'Transit', 'signal','','unitFirst',{}, '<signal>','', ctthLeaves('px'));
S(end).leafRule = ctthTimeRule();

% Per-pixel is the SAME MARKERS on a second set of units, so it takes the same word
% table - exactly as the stimulus response's ppx row does.
S(end+1) = row('ctthPpx', '^ctth\.(?<sig>ppx)\.(?<branch>bs)(?<leaf>\w+)$', 'Transit', ...
        'signal','','pixel2D',{}, '<signal>','', ctthLeaves('a.u.'));
S(end).leafRule = ctthTimeRule();

% ---- what belongs to the RECORDING rather than to a vessel --------------------
% THE FLOOR IS THE ONE THAT MATTERS AND IT IS WHY THESE ARE HERE.  A spread of transit
% times means nothing without the smallest spread this recording could resolve, and the
% infusion term says how much of every absolute time was the injection rather than the
% circulation.  Both are one number per recording, so they are the whole recording's
% variables and not a summary of anybody's segments.
S(end+1) = row('ctthRecording', '^ctth\.(?<leaf>cthFloor|cthUsable|inputM1)$', ...
        'Transit', 'fixed','whole','none',{}, 'The recording','', ctthRecordingLeaves());

% The measured input itself, on the recording clock, and which vessels it came from.
S(end+1) = row('ctthInput', '^ctth\.(?<leaf>aif)$', 'Transit', ...
        'fixed','whole','none',{'time'}, '','', ctthInputLeaves());
S(end+1) = row('ctthInputUnits', '^ctth\.(?<leaf>aifUnits)$', 'Transit', ...
        'fixed','seg','unitFirst',{}, '','', ctthInputLeaves());

% ---- the wall motion: results.wallMotion ---------------------------------------
% A FIFTH FAMILY, AND IT IS NOT A TREE.  results.wallMotion is FLAT - one struct of
% [nRow x 1] arrays indexed exactly like results.sMetrics, wallMotion.idx ==
% sMetrics.idx - so a join anywhere on the product is by row and every leaf below is an
% ordinary per-segment scalar.  It sits beside results.pulsatility on the same product
% and the two are SIBLINGS with different conventions: the pulsatility tree branches by
% signal and then by container, this one does not branch at all.
%
% WHAT IT IS FOR IS THE EVIDENCE.  The eight wm* columns in sMetrics are the answer
% after three gates, NaN in all eight wherever any of them refused; these are the same
% measurements BEFORE the gates, with the floor they were judged against and the flag
% saying which gate refused - so a segment that reports nothing still says why.
%
% EVERY DISPLACEMENT IS PEAK TO PEAK, per SEGMENT, and Core/MotionEnhancement's README
% quotes zero-to-peak per CUT throughout: these numbers are about twice that README's
% and are averaged over a whole vessel rather than found at its hottest cut.  The two
% sets are not comparable, which is why no word here invites the comparison.
S(end+1) = row('wallMotionRows', ['^wallMotion\.(?<leaf>idx|wallPxRaw|centrePxRaw|nullPx|' ...
        'nullSpreadPx|centreNullPx|conf|centreConf|cohere|centreCohere|nCut|nCutOffered|' ...
        'diameterPx|diameterMaskUm|taperPx|edgeUm|pathPx|clr|spurFrac|' ...
        'gateSize|gateCuts|gateConf|cleared)$'], ...
        'Wall motion', 'fixed','seg','unitFirst',{}, 'Wall motion','', wallMotionLeaves());

% AND THREE NUMBERS BELONG TO THE MEASUREMENT RATHER THAN TO A VESSEL: how many
% matched controls every confidence was divided by, and - on a continuous recording -
% which band was magnified and which one was the control.  They are what the per-vessel
% confidences have to be read against, exactly as the transit floor is.
S(end+1) = row('wallMotionRun', '^wallMotion\.(?<leaf>nControls|bandHz|offBandHz)$', ...
        'Wall motion', 'fixed','whole','none',{}, 'The measurement','', wallMotionRunLeaves());

% ---- vascular density: results.topology ----------------------------------------
% THE UNIT IS A SCOPE, AND A SCOPE IS NOT A COORDINATE.  results.topology.metrics is
% one row per analysed AREA - row 1 the whole analysed field, then one per drawn region
% - so a leaf here is [nScope x 1] and the row index is a label id in exactly the sense
% a segment index is: there is no "density against region number" plot and the tool
% must not offer one.  So no axis was registered for it, for the reason S5 refused a
% bolus axis: a coordinate nothing can be walked along is not a coordinate.
%
% ITEM FIRST IS WHAT KEEPS EVERY SCOPE A POINT.  With layout 'none' a box would read
% M(:,1,1) - the whole field's row - and silently drop every drawn region while looking
% completely right.  Declared item-first, the box carries one point per analysed area
% and the whole field is the first of them.
%
% A DENSITY IS A PROPERTY OF AN AREA, so none of this is in sMetrics and none of it can
% be averaged over segments: runTopologyAnalysis deliberately writes no per-segment
% column, and putting one number on a thousand rows is exactly what would invite it.
S(end+1) = row('topologyMetrics', '^topology\.metrics\.(?<leaf>.+)$', ...
        'Vascular density', 'fixed','whole','unitFirst',{}, 'The analysed area','', ...
        topologyLeaves());

% The distance transform itself, kept only when the step was asked for it.
S(end+1) = row('topologyMap', '^topology\.(?<leaf>evdMap)$', ...
        'Vascular density', 'fixed','pixel','pixel2D',{}, '','', topologyMapLeaves());

% ---- the myograph window: diameter, then propagation ---------------------------
% THE MEASUREMENT IS PER LINE, and it is the exact analogue of a speckle recording's
% sData: [samples x items], item last, one item per row across the vessel, with the
% item count declared beside it as .nY.  Same layout, same rule, one row apart.
S(end+1) = row('myoDiameterLines', [MYO 'diameter\.lines\.(?<sig>outer|mid|inner)$'], ...
        'Diameter', 'fixed','line','unitLast',{'time'}, '<signal>','px', {});
S(end).tail = 'along the vessel';

% THE WINDOW'S OWN TRACE IS THE SAME QUANTITY POOLED OVER THOSE LINES
% (runMyographDiameter's diameterBranch - the mean over the rows that were really
% measured), so it is declared as riding with them: one quantity, one menu entry,
% and 'Units are: pooled' on the per-line leaf reproduces this array.  A file that
% carries no .lines is the case exScan resolves - see its note on the pooled role -
% and there the trace is a variable of its own, spelled as it always was.
S(end+1) = row('myoDiameter', [MYO 'diameter\.(?<sig>outer|mid|inner)$'], 'Diameter', ...
        'fixed','whole','none',{'time'}, '<signal>','px', {});
S(end).tail       = 'averaged along the vessel';
S(end).pairedWith = 'lines.<sig>';
S(end).pairRole   = 'pooled';
S(end+1) = row('myoDiameterStats', [MYO 'diameter\.stats\.(?<sig>\w+)\.(?<leaf>\w+)$'], 'Diameter', ...
        'fixed','whole','none',{}, '<signal>','px', statLeaves());

% THE TWO WALL POSITIONS ARE VARIABLES, and worth being ones: the kymograph of a wall
% position is the picture that shows whether the detection wandered, which the
% diameter alone cannot - both walls drifting together leave the difference flat.
% They are the per-line diameter's shape exactly (the diameter IS their difference),
% so they take the same row, and no plot id had to be invented for them.
S(end+1) = row('myoDiameterWalls', ...
        [MYO 'diameter\.(?<leaf>wallL|wallR)\.(?<sig>outer|mid|inner)$'], ...
        'Diameter', 'fixed','line','unitLast',{'time'}, '<signal>','px', wallLeaves());

% AND THE TWO QUALITY FLAGS RIDE WITH THEM, in a role of their own.  Neither is a
% quantity anyone wants a menu entry for: .measured says which points were detected
% rather than filled, .valid which frames had both walls in the field of view, and
% both are the QUALITY of the diameter beside them.  They are declared on the ROW,
% like the pooled trace and for the same reason - their host is diameter.lines, which
% has no leaf name because the measure is its signal token.  Unlike the pooled role
% this one is NOT conditional: myographDiameterBranch writes the arrays and the flags
% in one breath and removes them in one breath, so a file carrying a flag and no
% array does not exist.
S(end+1) = row('myoMeasured', [MYO 'diameter\.(?<leaf>measured)$'], 'Diameter', ...
        'fixed','line','unitLast',{'time'}, '','', qualityLeaves());
S(end).pairedWith = 'lines';  S(end).pairRole = 'mask';
S(end+1) = row('myoValid', [MYO 'diameter\.(?<leaf>valid)$'], 'Diameter', ...
        'fixed','whole','none',{'time'}, '','', qualityLeaves());
S(end).pairedWith = 'lines';  S(end).pairRole = 'mask';

% ---- the wire myograph: the recorded signal itself -----------------------------
% A WIRE RECORDING'S SAMPLES ARE ITS MEASUREMENT, and they are stored in whichever of
% two places the recording has reached: whole on the channel until the intervals step
% cuts it, then inside each window.  Both are the same quantity and take the same
% declaration; two rows, because the two paths cannot be one pattern without making
% the channel index optional in a place it is not.
%
% THE CHANNEL INDEX IS THE SIGNAL TOKEN, which is the whole of why it is captured.
% Four chambers on one LabChart file produce four leaves spelled identically, and the
% Signal control keys on .signal - so without it all four would collapse into one
% entry and three chambers would be unreachable.  The index is never SHOWN: exScan
% picks the channel's real name up on the way past and the menu reads that.
%
% SO THE CHANNEL PREFIX IS REQUIRED HERE, unlike in IV above, and for two reasons.
% Only a WIRE window ever carries samples - a pressure window leaves .samples empty,
% which exScan skips as an absent branch - so nothing is excluded by demanding it.
% And it cannot be made optional even if one wanted to: MATLAB's regexp drops a NAMED
% token that sits inside a non-capturing optional group, silently and without leaving
% the field in the struct, so '(?:channel\((?<sig>\d+)\)\.)?' matches the path and
% then reports no signal at all - and four chambers collapse into one menu entry
% again, which is the bug this token exists to prevent.
S(end+1) = row('wireWindowSamples', ...
        '^channel\((?<sig>\d+)\)\.intervals\(\d+\)\.samples\.(?<leaf>data)$', ...
        'Recording', 'fixed','whole','none',{'time'}, '','', sampleLeaves());
S(end+1) = row('wireChannelSamples', '^channel\((?<sig>\d+)\)\.(?<leaf>data)$', ...
        'Recording', 'fixed','whole','none',{'time'}, '','', sampleLeaves());

% The three special propagation leaves come FIRST - the catch-all below would
% otherwise claim them and call a 475-row table of lags a scalar.
S(end+1) = row('myoLagByRow', [MYO 'propagation\.(?<sig>\w+)\.(?<leaf>lagByRow)$'], 'Propagation', ...
        'fixed','line','xyPairs',{}, 'Propagation','s', propLeaves());
% speedCI keeps a row of its own for its SHAPE - [1x2] is not a curve - while what it
% rides on is declared once, in companionTable, beside the three that ride on
% lagByRow.  Setting it here as well would be a second place to keep in step.
S(end+1) = row('myoSpeedCI', [MYO 'propagation\.(?<sig>\w+)\.(?<leaf>speedCI)$'], 'Propagation', ...
        'fixed','whole','none',{}, 'Propagation','px/s', propLeaves());
S(end+1) = row('myoPropMetrics', [MYO 'propagation\.(?<sig>\w+)\.metrics\.(?<leaf>\w+)$'], 'Propagation', ...
        'fixed','whole','none',{}, 'Propagation','', propMetricLeaves());
S(end+1) = row('myoProp', [MYO 'propagation\.(?<sig>\w+)\.(?<leaf>\w+)$'], 'Propagation', ...
        'fixed','whole','none',{}, 'Propagation','', propLeaves());
end

% =====================================================================
function r = row(id, pattern, family, unitBy, unit, layout, dims, phrase, yUnit, leaves)
%row  One container row.  .leafRule is empty on all but the handful of rows whose
%   producer DERIVES its leaf names from a setting - see nvcTimeRule - and is a
%   {pattern, words, unit} triple consulted only when the literal table has no entry.
r = struct('id',id, 'pattern',pattern, 'family',family, 'unitBy',unitBy, ...
    'unit',unit, 'layout',layout, 'dims',{dims}, 'phrase',phrase, 'tail','', ...
    'yUnit',yUnit, 'leaves',{leaves}, 'leafRule',{{}}, 'leafOverride',struct(), ...
    'pairedWith','', 'pairRole','', 'signal','', 'branch','', 'leaf','', 'label','');
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
%vsmFVectorLeaves  The BAND-LESS container: per-frequency moments over the whole
%   recording, which belong to no band and are named as themselves.
L = { 'ampMean',    'mean amplitude',              'a.u.'
      'ampStd',     'amplitude variability',       'a.u.'
      'ampSkew',    'amplitude skewness',          ''
      'ampMeanPct', 'mean amplitude by percentile','a.u.'
      'ampFlare',   'amplitude during bursts',     'a.u.'
      'ampSilence', 'amplitude when quiet',        'a.u.' };
end

function L = vsmFVectorBandLeaves()
%vsmFVectorBandLeaves  The same six leaves UNDER A BAND, where the band picked the
%   time samples and not the frequencies.  Each phrase finishes the sentence the
%   row's head starts, so the whole label reads 'Whole frequency spectrum - while
%   the vasomotion band was bursting' rather than claiming to be a band's spectrum.
%   The three moments keep their own words: they are never written under a band, and
%   a tree that started to would still be named for what it is.
L = { 'ampMean',    'mean amplitude',                             'a.u.'
      'ampStd',     'amplitude variability',                      'a.u.'
      'ampSkew',    'amplitude skewness',                         ''
      'ampMeanPct', 'grouped by how strong <band> was',           'a.u.'
      'ampFlare',   'while <band> was bursting',                  'a.u.'
      'ampSilence', 'while <band> was quiet',                     'a.u.' };
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

function L = pulsScalarLeaves(u)
%pulsScalarLeaves  ONE WORD TABLE FOR ALL THREE PULSATILITY BRANCHES.  u is the
%   trace's OWN unit - a.u. for a flow index and for a plasma-volume waveform, px for
%   the tracked diameter - and it is an argument rather than a column here because the
%   same nine markers are taken on all three and only the unit and the head differ.
%   That is the arrangement nvcLeaves and ctthLeaves already use, and the reason is the
%   same: three copies of one table differing by a two-letter prefix is three places to
%   keep in step.
%
%   THE TIMES ARE ON THE BEAT'S OWN CLOCK, which is one cardiac cycle and not the
%   recording, so they are seconds from the foot of the beat and never absolute times.
L = { 'Min',      'minimum',              u
      'Max',      'maximum',              u
      'Mean',     'mean over the beat',   u
      'Std',      'variability',          u
      'TimeMin',  'time of the minimum',  's'
      'TimeMax',  'time of the maximum',  's'
      'PI',       'pulsatility index',    ''
      'RI',       'resistance index',     ''
      'SymRatio', 'symmetry ratio',       ''
      'R2',       'fit quality',          '' };
end

function L = pulsHarmonicLeaves()
%pulsHarmonicLeaves  The fitted model's coefficients.  They are written only when the
%   model was asked for, and they are the two leaves whose name does not say which
%   quantity they belong to - so the unit here is the flow branch's and a diameter
%   model's coefficients are in pixels despite it.  See the row above.
L = { 'hAmp',   'harmonic amplitudes', 'a.u.'
      'hPhase', 'harmonic phases',     'rad' };
end

function L = nvcLeaves(u)
%nvcLeaves  ONE WORD TABLE FOR THE STIMULUS RESPONSE AND ALL THREE BRANCHES.  u is
%   the trace's OWN unit - a.u. for a flow index, px for the tracked diameter - and it
%   is an argument rather than the row's default because the same quantity is measured
%   on both and only its unit differs.  Everything dimensionless says so explicitly,
%   so an empty third column here means "this quantity really has no unit" and never
%   "nobody filled it in".
%
%   THE THREE RELATIVES ARE THE MARKERS TO READ, and the reason is measured: dividing
%   by the baseline is what makes a number comparable between two segments and between
%   two recordings, and on the reference recording it is the relatives that survive a
%   second estimator of the same vessel where the absolute amplitudes do not.
%
%   AucRel IS IN SECONDS AND THAT IS NOT A SLIP.  It is the accumulated change divided
%   by the baseline level, so the amplitude divides out and the time does not: read it
%   as "the response delivered as much extra flow as this many seconds of resting
%   flow".  A reader who sees no unit on it would take it for a fraction, which is why
%   the second is written down here rather than left blank.
%
%   THE TIMES ARE ALL MEASURED FROM THE STIMULUS MARK, so they may be negative and a
%   negative one means the trace was already rising before the stimulus.  The windows
%   in the SETTINGS are on the repetition clock, which is the clock the protocol is
%   written in; these are not.  They are NOT LISTED HERE because their names come from
%   a setting - see nvcTimeRule.
%
%   AND EVERY ROW CARRIES THE TWO CONFIDENCE NUMBERS whatever else was asked for.  A
%   reader who selects only the amplitudes still gets to see how far to believe them,
%   which is the whole reason the producer never gates them.
L = { % --- the levels the response is measured against ----------------------------
      'Bl',       'baseline level',                                  u
      'BlStd',    'baseline noise',                                  u
      'Fn',       'level at the end of the repetition',              u
      'FnStd',    'noise at the end of the repetition',              u
      'Ep',       'average level over the whole repetition',         u
      % --- the amplitudes, absolute then relative to baseline ---------------------
      'St',       'average change while the stimulus was on',        u
      'StRel',    'change while the stimulus was on, against baseline', ''
      'Pk',       'peak change from baseline',                       u
      'PkRel',    'peak change, against baseline',                   ''
      'Auc',      'change accumulated over the repetition',          [u '*s']
      'AucRel',   'accumulated change, in seconds of baseline',      's'
      % --- how far the repetition is to be believed, for THIS unit ----------------
      'Conf',     'confidence in this repetition',                   ''
      'ConfMin',  'the weakest single check on this repetition',     ''
      'Trust',    'the repetition was kept for this unit',           '' };
end

function R = nvcTimeRule()
%nvcTimeRule  THE ONE LEAF FAMILY IN THIS FILE THAT IS NAMED FROM A SETTING.  The
%   timing markers are one per entry of s.nvcAreaPcts - T10, T50 and T90 by default,
%   but a recording processed with [5 25 50 75 95] carries five - so a literal table
%   would leave whichever levels it had not anticipated shown as their own field
%   names.  The pattern gives every level its words, and the level itself is what the
%   words are about.
%
%   THE LEVEL IS A SHARE OF THE ACCUMULATED CHANGE, not of the peak: the marker is the
%   moment that much of the repetition's total change had been delivered.  That is why
%   T90 is late in an epoch whose response does not come back - the tail is still
%   accumulating - and it is the property that makes these times comparable between
%   segments of very different amplitude.
R = { '^T(\d+)$', 'time by which $1 % of the change had accumulated', 's' };
end

function L = ctthLeaves(u)
%ctthLeaves  ONE WORD TABLE FOR THE BOLUS TRANSIT AND BOTH ITS BRANCHES.  u is the
%   trace's OWN unit - a.u. for a dye curve, px for the tracked diameter - and it is an
%   argument rather than the row's default because the same construction is applied to
%   both and only the unit differs.  Everything dimensionless says so explicitly, so an
%   empty third column means "this quantity really has no unit" and never "nobody filled
%   it in".
%
%   THREE MARKERS COMPARE BETWEEN RECORDINGS AND THE REST DO NOT, and the words say so
%   rather than leaving a reader to find out.  The delay, the mean transit time and the
%   spread are DIFFERENCES against the recording's own arterial curve, so the injection -
%   which lasts seconds and adds itself to every absolute time - is in both terms and
%   cancels.  T0, the level times and the time of the peak are absolute, and two
%   recordings injected differently cannot be compared on them.
%
%   THE SPREAD IS BLANK BELOW THE RECORDING'S FLOOR, not small.  ctth.cthFloor is the
%   smallest spread that recording could resolve, and a value under it is left empty
%   because the arithmetic that produces it is a difference of two numbers that agree to
%   one per cent - a plausible answer there is not a small answer, it is no answer.
%
%   THE BLOOD VOLUME IS RELATIVE AND IT IS EXPOSED TO HAZE.  It is the plateau this unit
%   reached against the plateau the arterial curve reached, which is the volume-like
%   quantity for a tracer that does not clear - but scattered light from the vessels
%   floods the tissue around them, so until the glow has been removed it orders units
%   within one recording rather than measuring a fraction.
L = { % --- the levels the transit is measured against ------------------------------
      'Bl',       'level before the injection',                      u
      'BlStd',    'noise before the injection',                      u
      'Fn',       'level at the end of the recording',               u
      'FnStd',    'noise at the end of the recording',               u
      'Pk',       'peak above the pre-injection level',              u
      % --- the amplitudes, absolute then relative -----------------------------------
      'Step',     'plateau above the pre-injection level',           u
      'PkRel',    'how far the peak overshot the plateau',           ''
      'BvRel',    'relative blood volume',                           ''
      % --- absolute times, NOT comparable between recordings ------------------------
      'T0',       'onset of the rise',                               's'
      'TPeak',    'time of the peak',                                's'
      % --- the transit markers, referenced to the arterial curve --------------------
      'Delay',    'delay after the artery',                          's'
      'Mtt',      'mean transit time',                               's'
      'Cth',      'spread of transit times',                         's'
      % --- the shape of the leading edge -------------------------------------------
      'EdgeWid',  'width of the leading edge',                       's'
      'RiseUp',   'steepest rise, in plateaus per second',           '1/s'
      'FallDn',   'steepest fall, in plateaus per second',           '1/s'
      % --- how far this unit is to be believed --------------------------------------
      'Conf',     'confidence in this vessel',                       ''
      'ConfMin',  'the weakest single check on this vessel',         ''
      'Trust',    'this vessel passed both thresholds',              '' };
end

function R = ctthTimeRule()
%ctthTimeRule  The bolus timing markers, one per entry of s.ctthLevelPcts.  A literal
%   table would leave whichever levels it had not anticipated shown as their own field
%   names, so the pattern gives every level its words.
%
%   THE LEVEL IS A SHARE OF THE PLATEAU, and that is where this differs from the
%   stimulus response's rule of the same shape.  A bolus of a non-clearing tracer never
%   comes back, so the CURVE ITSELF is the cumulative and there is nothing to accumulate
%   first - the marker is the moment the unit had filled that far towards the level it
%   ends at.  Reading it as "that much of the change had happened" would be right; as
%   "that much of an area" would not.
R = { '^T(\d+)$', 'time by which it had filled to $1 % of its plateau', 's' };
end

function L = ctthRecordingLeaves()
%ctthRecordingLeaves  The three numbers that belong to the recording and not to a
%   vessel.  All three are what the per-vessel markers have to be read against.
L = { 'cthFloor',  'smallest spread this recording can resolve',   's'
      'cthUsable', 'the recording can resolve a spread at all',    ''
      'inputM1',   'how much of every time was the injection',     's' };
end

function L = ctthInputLeaves()
%ctthInputLeaves  The arterial curve the whole measurement is quoted against, and the
%   vessels it was derived from.  They carry the WHOLE label rather than the tail of
%   one, because their rows have no head to hang it on - there is no signal these belong
%   to, they belong to the recording.
L = { 'aif',       'The arterial curve',                            ''
      'aifUnits',  'Vessels the arterial curve came from',          '' };
end

function L = wallMotionLeaves()
%wallMotionLeaves  THE EVIDENCE BEHIND THE EIGHT wm* COLUMNS, one row per segment.
%
%   THE MEASUREMENT AND ITS FLOOR ARE THE PAIR TO READ TOGETHER.  A wall that moves
%   0.01 px is not a small movement, it is a movement the instrument cannot separate
%   from its own noise - so the same measurement made on the same vessel with the
%   heartbeat's timing scrambled away is stored beside it, and the confidence is the
%   ratio of the two.  Below its own floor a segment reports nothing rather than a
%   small number, and these arrays are what say so.
%
%   EVERY DISPLACEMENT IS PEAK TO PEAK, over the averaged beat, and averaged over the
%   places along the vessel it was measured at.
%
%   THE THREE GATE FLAGS ARE WHY, NOT WHETHER.  A row that cleared has all three true;
%   a row that did not names the one that refused it - too narrow to measure, too few
%   places along it, or no more movement than its own scrambled copy.
L = { 'idx',            'segment it belongs to',                              ''
      'wallPxRaw',      'the wall, before the gates',                         'px'
      'centrePxRaw',    'the whole vessel sliding, before the gates',         'px'
      'nullPx',         'the wall, with the beat scrambled away',             'px'
      'nullSpreadPx',   'how much the scrambled copies differed from each other',   'px'
      'centreNullPx',   'the whole vessel, with the beat scrambled away',     'px'
      'conf',           'the wall against its scrambled copy',                ''
      'centreConf',     'the sliding against its scrambled copy',             ''
      'cohere',         'how well the places along the vessel agreed',         ''
      'centreCohere',   'how well they agreed about the sliding',              ''
      'nCut',           'places along the vessel it was measured at',          ''
      'nCutOffered',    'places along the vessel that were offered',           ''
      'diameterPx',     'width the wall fit settled on',                       'px'
      'diameterMaskUm', 'width from the segmentation, which is what the gate read', 'um'
      'taperPx',        'how much the width varied along the vessel',          'px'
      'edgeUm',         'how far the wall sat from the centre line',           'um'
      'pathPx',         'length of the centre line',                           'px'
      'clr',            'tortuosity',                                          ''
      'spurFrac',       'share of the centre line that was a spur',            ''
      'gateSize',       'it was wide enough to measure',                       ''
      'gateCuts',       'it offered enough places to measure',                 ''
      'gateConf',       'it moved more than its scrambled copy',               ''
      'cleared',        'it passed every check',                               '' };
end

function L = wallMotionRunLeaves()
%wallMotionRunLeaves  What the whole measurement was made with.  The count of matched
%   controls is what every confidence above was divided by; the two frequencies exist
%   only on a continuous recording, where a band is magnified and a neighbouring band
%   is the control.
L = { 'nControls',  'matched controls each vessel was compared with', ''
      'bandHz',     'frequencies that were magnified',                'Hz'
      'offBandHz',  'frequencies used as the control',                'Hz' };
end

function L = topologyLeaves()
%topologyLeaves  The columns of results.topology.metrics, one row per analysed area.
%
%   THE DENOMINATOR IS THE ANALYSED AREA AND NOT THE FRAME, which is why the first
%   column is here at all: every fraction below is of that area, so a recording whose
%   trust mask threw half the field away is still described honestly.
%
%   TWO OF THEM ARE COMPARATORS AND NOT ANATOMY, and the words say so rather than
%   leaving a reader to find out.  This is a flat picture of a vasculature that is not
%   flat, so a crossing at two depths cannot be told from a branching and the junction
%   count is not a count of branchings; and the capillary bed is not resolved at this
%   pixel size, so the distance to the nearest vessel is not a diffusion distance.
%   Both compare fields within one study and neither is an absolute.
%
%   THE UNITS FOLLOW THE CALIBRATION.  With micrometres per pixel given, lengths are
%   micrometres, areas mm2 and densities mm/mm2 or mm-2; without it everything is in
%   pixels.  results.topology.units records which was in force, and no word here
%   claims one - which is why the unit column below is empty for every density.
L = { 'scope',               'which area this row is about',                    ''
      'areaAnalysed',        'area that was analysed',                          ''
      'areaFraction',        'share of that area that is vessel',               ''
      'unsegmentedFraction', 'share that looked like vessel but could not be resolved', ''
      'lengthDensity',       'vessel length per area',                          ''
      'junctionDensity',     'vessel meetings per area',                        ''
      'endpointDensity',     'vessel ends per area',                            ''
      'calibreMedian',       'median vessel width',                             ''
      'calibreIQR',          'spread of vessel widths',                         ''
      'tortuosityMedian',    'median tortuosity',                               ''
      'evdMean',             'mean distance from tissue to the nearest vessel', ''
      'evdMedian',           'median distance from tissue to the nearest vessel',''
      'evdP95',              'distance the furthest twentieth of the tissue is beyond', ''
      'evdFractionBeyond',   'share of tissue further from a vessel than the chosen distance', ''
      'evdCoverage',         'share of tissue whose nearest vessel was inside the field',''
      'nSegments',           'vessels counted',                                 '' };
end

function L = topologyMapLeaves()
L = { 'evdMap', 'distance to the nearest vessel', '' };
end

function L = nvcQualityLeaves()
%nvcQualityLeaves  What was true of a REPETITION, at the level of the whole
%   recording.  Neither is a summary of the per-segment confidence and neither should
%   be read as one: the responding area is the share of the segmented field that
%   responded coherently, which is what decides whether a repetition is kept, and the
%   lost time is what the camera did not deliver during it.
%
%   A REPETITION IN WHICH A TENTH OF THE FIELD RESPONDED IS A GOOD REPETITION.  Only
%   part of the field is ever a responder, so the share is not expected to approach
%   one and a low value is not by itself a bad recording.
L = { 'areaFrac', 'how much of the field responded',        ''
      'timeLoss', 'frames the camera did not deliver',      's' };
end

function L = nvcEpochLeaves()
%nvcEpochLeaves  The two leaves that describe the repetitions themselves rather than
%   a quantity measured in one.  They carry the WHOLE label rather than the tail of
%   one, because their rows have no head to hang it on - there is no signal these
%   belong to.
L = { 'epochTrust',  'Repetitions that were kept',            ''
      'stimulus',    'The stimulus',                          '' };
end

function L = wallLeaves()
%wallLeaves  The two edges of the vessel, in the words a biologist reads them in.
%   'left' and 'right' are positions IN THE IMAGE, which is what they are: the
%   detector works along an image row and reports where the dark band starts and
%   where it ends.
L = { 'wallL', 'left wall position',  'px'
      'wallR', 'right wall position', 'px' };
end

function L = qualityLeaves()
%qualityLeaves  The two flags that say how much of a window was really seen.  They
%   are never offered - they ride with the per-line diameter - but a rider keeps its
%   words, because exFetch labels a companion from this table and a label that fell
%   back to the field name would put 'measuredFraction' in front of a reader.
%   They carry the WHOLE label rather than the tail of one, because their row has no
%   head to hang it on: there is no single measure these belong to.
L = { 'measured', 'Where the walls were detected rather than filled', ''
      'valid',    'Frames with both walls in view',                   '' };
end

function L = sampleLeaves()
%sampleLeaves  A wire recording's own signal.  No unit is claimed: what a LabChart
%   channel is measured in is a property of that channel (mN here, but not always),
%   the channel carries it, and this table reads nothing.
L = { 'data', 'Recorded signal', '' };
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
%propLeaves  THE FOUR RIDERS ARE NAMED HERE TOO.  They are never offered as
%   variables, but exFetch labels a companion from this table, and a label that fell
%   back to the field name would put 'confidenceText' in front of a biologist.
L = { 'domFreq',         'dominant frequency',           'Hz'
      'speed',           'speed',                        'px/s'
      'speedCI',         'speed confidence interval',    'px/s'
      'slope',           'lag slope',                    ''
      'R2',              'fit quality',                  ''
      'pValue',          'significance',                 ''
      'nRows',           'rows used',                    ''
      'confidence',      'confidence',                   ''
      'confidenceLevel', 'confidence',                   ''
      'confidenceText',  'what the fit does and does not support', ''
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
L = { 'imgBFI',        'blood flow map',            'a.u.'
      'imgK',          'contrast map',              'a.u.'
      'imgI',          'intensity image',           'a.u.'
      'imgBackground', 'background that was removed from the intensity image', 'a.u.'
      'imgStdBFI',     'blood flow variability map','a.u.'
      'mask',          'analysed area',             ''
      'commonMask',    'shared analysed area',      ''
      'sMap',          'segment map',               ''
      'pMap',          'parenchyma cell map',       ''
      'cMask',         'tissue category map',       ''
      'dvsMap',        'tracked vessel map',        ''
      'mapType',       'vessel type map',           ''
      'mapTree',       'vascular tree map',         ''
      'regionsMask',   'selected regions',          '' };
end

function L = metricsLeaves()
%metricsLeaves  The columns of sMetrics / dvsMetrics.  CLR is the segment's path
%   length over its end-to-end distance (runDynamicSegmentation.m:167) - 1 for a
%   straight vessel, 2 for a coil - so it is tortuosity and is named as such.
%
%   THIS IS THE ONE TABLE WHERE ALL THREE PULSATILITY FAMILIES CAN SIT SIDE BY SIDE,
%   and each of the three says what it is a pulsatility OF rather than leaving the
%   unmarked one to be guessed at.  A ps* column is a pulsatile FLOW; a pd* column is a
%   pulsatile DIAMETER, read off a per-frame width that cannot change by less than a
%   quarter of a pixel - the right instrument for vasomotion and vasoreactivity, where
%   a vessel moves by whole pixels, and the wrong one for a heartbeat, where the change
%   is a tenth of that step; a pv* column is a pulsatile PLASMA VOLUME, how much the
%   labelled plasma in a vessel rises and falls over one beat as a fraction of its own
%   mean.  NO TWO OF THEM MAY BE POOLED.
%
%   AND THE wm* COLUMNS ARE A FOURTH QUANTITY ON THE SAME ROWS: how far the wall itself
%   moved, measured under a pixel against the same vessel with the heartbeat's timing
%   scrambled away.  They answer the question pd* is too coarse for, they are PEAK TO
%   PEAK and per segment, and a row that is NaN in all eight was refused by a gate
%   while its pv* neighbours on the same row remain perfectly good - the two families
%   share no gates.
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
      'psPI',           'flow pulsatility index',   ''
      'psRI',           'flow resistance index',    ''
      'psMean',         'mean flow over the beat',  'a.u.'
      'psTimeMin',      'time of the lowest flow',  's'
      'psTimeMax',      'time of the highest flow', 's'
      'psSymRatio',     'flow symmetry ratio',      ''
      'pdPI',           'diameter pulsatility index',''
      'pdRI',           'diameter resistance index', ''
      'pdMean',         'mean diameter over the beat','px'
      'pdTimeMin',      'time of the narrowest diameter','s'
      'pdTimeMax',      'time of the widest diameter','s'
      'pdSymRatio',     'diameter symmetry ratio',  ''
      % the six the fluorescence beat duplicates out of results.pulsatility.  A
      % plasma volume is neither of the two above and pools with neither.
      'pvPI',           'plasma volume pulsatility index', ''
      'pvRI',           'plasma volume resistance index',  ''
      'pvMean',         'mean plasma volume over the beat','a.u.'
      'pvTimeMin',      'time of the smallest plasma volume','s'
      'pvTimeMax',      'time of the largest plasma volume', 's'
      'pvSymRatio',     'plasma volume symmetry ratio',    ''
      % the eight the wall-motion step writes.  Peak to peak, per segment, and NaN in
      % all eight wherever a gate refused the row.
      'wmWallPx',       'wall movement, peak to peak',     'px'
      'wmWallUm',       'wall movement, peak to peak',     'um'
      'wmDilRel',       'width change over the beat',      '%'
      'wmCentrePx',     'movement of the whole vessel, peak to peak', 'px'
      'wmCentreUm',     'movement of the whole vessel, peak to peak', 'um'
      'wmDilShare',     'how much of the movement was width change rather than sliding', ''
      'wmConf',         'wall movement against its scrambled copy',   ''
      'wmCuts',         'places along the vessel it was measured at', ''
      'typeConfidence', 'vessel type confidence',   ''
      % what setVascularTree writes, identical in name and shape on both branches
      'treeID',         'tree it belongs to',       ''
      'generation',     'branchings from the root', ''
      'strahlerOrder',  'branching order',          ''
      'hierarchyConfidence', 'confidence in its place in the tree', ''
      'flowPotential',  'how far along the flow it sits', ''
      'isRoot',         'it is where a tree begins',      ''
      'isOutlet',       'it is where a tree ends',        ''
      % the five the transit step duplicates out of results.ctth, so a reader who works
      % from the table alone still gets the markers that compare between recordings and
      % the number that says how far to believe them.  bd* is the same construction on
      % the tracked diameter and must never be pooled with bs*.
      'bsDelay',        'delay after the artery',   's'
      'bsMtt',          'mean transit time',        's'
      'bsCth',          'spread of transit times',  's'
      'bsBvRel',        'relative blood volume',    ''
      'bsConf',         'transit confidence',       ''
      'bdDelay',        'diameter delay after the artery',  's'
      'bdMtt',          'diameter mean transit time',       's'
      'bdCth',          'diameter spread of transit times', 's'
      'bdBvRel',        'relative diameter change',         ''
      'bdConf',         'diameter transit confidence',      '' };
end
