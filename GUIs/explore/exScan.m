%exScan - Every plottable number in a results tree, and what it is a function of.
%
%   D = exScan(results) walks a whole results tree and returns one DESCRIPTOR per
%   plottable leaf: where it lives, whose it is, which axes it varies along, and
%   what to call it.  It is the reader that replaces the hand-written plot
%   catalogue - adding a leaf to the vasomotion tree makes it plottable, rather than
%   requiring someone to edit a switch statement.
%
%   D = exScan(results, path) does the same and additionally knows which file it
%   came from, which lets exAxes reach the sibling settings file for the percentile
%   levels.  When `results` is the lazy handle struct('x_lazy_',true,'x_path_',p)
%   the walk is done through h5info instead, and nothing is loaded.
%
%   THE THREE SKIP RULES ARE NOT TIDINESS - each is forced by the real data, and
%   getting one wrong offers the user a plot that cannot be drawn.
%
%     AN EMPTY NUMERIC IS AN ABSENT BRANCH.  Every wire-myograph file in the
%     reference set carries channel(i).intervals(k).diameter and .propagation as
%     fields holding double [0x0] - a wire myograph measures neither.  isfield is
%     true and the field is meaningless.  Without this rule every wire recording in
%     the library gets a diameter plot it has no data for.
%
%     TEXT AND CELL LEAVES ARE LABELS, NEVER VARIABLES.  method, direction,
%     speedUnit, confidenceLevel, confidenceText, qualityFlags, channelName,
%     gsType, name, and the string columns of the metrics tables.  One of them
%     matters enough to keep: channelName is the REAL name of a wire channel
%     ('channel 1') whose struct field has been sanitised ('channel1'), and the
%     user must see the former - so it rides on the descriptor as .signalLabel.
%
%     BOOKKEEPING IS SKIPPED.  meta, comments, timeCrop, frames, channels,
%     timeStamp, nY - and, added here for the same reason the spec skips `frames`,
%     tStart and tEnd, which are that same window boundary written in seconds
%     instead of frame numbers.  (comments is a 51-to-74 element event log on the
%     wire files, with time / text / channel / record.  It is not a variable, but it
%     IS the natural annotation overlay for any time-axis plot.  Deliberately out of
%     scope; noted so it stays a decision.)
%
%   AND A FOURTH, WHICH THE SHAPE MODEL IMPLIES RATHER THAN STATES: A COORDINATE IS
%   NOT A VARIABLE.  time, gsTime, timeWT, timeDWT, f, pctCenters and harmonics are
%   the axis registry's raw material.  Left in, every real file in the library would
%   report half a dozen leaves the schema does not cover, and the `inferred` signal -
%   which is meant to mean "this tree grew something new, go look" - would be noise
%   on every scan.
%
%   NOTHING IS EVER SILENTLY DROPPED FOR BEING WRONG.  After the schema names a
%   leaf's dims, the trailing sizes are checked against the resolved axis lengths.
%   A leaf that does not verify keeps its descriptor, gets .suspect true and a .note
%   saying which axis disagreed and by how much, and is offered last.  A leaf no row
%   claims is classified by the fallback - unit from the signal token in its path,
%   dims by matching sizes against the axis registry - and marked .inferred under
%   family 'Other'.  That is what lets the tool survive a new branch without an
%   edit, and what makes a missing schema row visible instead of invisible.
%
%   WHOLE OR LINE IS READ OFF THE DATA.  A myograph <VSM> is one unit when the run
%   had s.perLine false and one per row across the vessel when it did not; the
%   schema returns unit 'auto' and the item count decides.  On the reference
%   pressure recording every <VSM> leaf comes back with one unit, hence 'whole'.
%
%   THE LAZY PATH IS COMPLETE FOR THE BRANCHES THAT NEED IT.  Measured on the real
%   3.63 GB LSCI_...._t_BFI_r.mat: h5info(f,'/results') returns in 45 ms and the walk
%   reaches every leaf of its vasomotion tree - all four signals, scalars /
%   fVectors / timeVectors / spectrum, down to spectrum/amp [1480 x 70 x 289].  It
%   STOPS at struct arrays, which appear as opaque object-reference datasets, and
%   tables are opaque the same way (sMetrics reports as a [6x1] object - the table's
%   six internal fields, not its 1480 rows).  Neither limitation overlaps with the
%   big files: struct arrays are a myograph shape and myograph files are megabytes,
%   while the multi-gigabyte files are per-pixel speckle results whose heavy
%   branches are plain scalar-struct trees.  Everything the h5 walk produces carries
%   .lazy true, so exFetch knows to go through h5read.
%
%   Syntax:
%      D = exScan(results)
%      D = exScan(results, path)
%
%   INPUTS
%     results  the RESULTS member of a triplet, loaded; or the lazy handle
%              struct('x_lazy_',true,'x_path_',<file>).
%     path     the _r.mat it came from.  Optional for a loaded tree, and only used
%              to find the sibling _s.mat; taken from the handle when lazy.
%
%   OUTPUTS
%     D  1xN struct array, sorted by path so a loaded scan and a lazy one of the
%        same tree are comparable:
%          .path        dotted and resolvable: 'intervals(2).vasomotion.outer.f'.
%                       A component is a struct-array INDEX only when it is
%                       name(<digits>), which is why the metrics column
%                       'std(diameter)' can never be mistaken for one
%          .family      Flow | Pulsatility | Vasomotion | Diameter | Propagation |
%                       Maps | Metrics | Other
%          .signal      sData | dvsData | dvsDiameter | gsData | ppx | outer | mid |
%                       inner | <sanitised channel name> | ''
%          .signalLabel the channel's real name, for the legend; '' when there is one
%          .leaf        the field name itself
%          .branch      'VB' | 'CB' | '' when the container has no band branch
%          .unit        seg | dvs | pixel | line | whole
%          .layout      unitFirst | unitLast | pixel2D | xyPairs | none
%          .dims        cellstr of axis names, in order
%          .size        size() of the leaf
%          .nUnits      how many items it holds, read through .layout
%          .kinds       {} here - exPlotRules fills it in session 2
%          .label       the human phrase, from exSchema
%          .yUnit       'a.u.' | 'Hz' | 's' | 'px' | 'px/s' | 'rad' | ''
%          .suspect     true when a size did not agree with its axis
%          .note        why, when it did not
%          .inferred    true when no schema row claimed it
%          .lazy        true when it must be fetched through h5read
%          .pairedWith  'speed' for a confidence interval; '' otherwise.  A
%                       descriptor with this set is an error bar on another
%                       variable, so a variable LIST filters on isempty(.pairedWith)
%
%   EXAMPLE
%      S = load('LSCI_20240809_1ADCF08BP_c_BFI_r.mat');
%      D = exScan(S.results);
%      numel(D)                                  % 105
%      D(find(strcmp({D.leaf},'hAmp'),1)).dims   % {'harmonic'}
%
%   DEPENDS ON
%     exAxes, exSchema (GUIs/explore).
%
% See also: exModality, exAxes, exSchema, guiExplore, assembleVasomotionTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function D = exScan(results, resultsPath)

if nargin<2, resultsPath = ''; end
D = emptyDescriptors();
if isempty(results) || ~isstruct(results) || ~isscalar(results), return; end

lazy = isfield(results,'x_lazy_');
if lazy && isempty(resultsPath) && isfield(results,'x_path_')
    resultsPath = char(results.x_path_);
end

ctx = struct('R',results, 'path',char(resultsPath), 'lazy',lazy, ...
             'ax',containers.Map('KeyType','char','ValueType','any'));

if lazy
    D = walkH5(D, ctx);
else
    D = walkStruct(D, results, '', ctx, '');
end
D = sortByPath(D);
end

% ================================================================== THE WALK

function D = walkStruct(D, node, pth, ctx, sigLabel)
%walkStruct  One scalar struct: its fields, in the order the producer wrote them.

% A channel, or a <VSM> written by the wire path, carries the REAL channel name
% beside the sanitised field name.  Pick it up on the way past.
sigLabel = labelOf(node, pth, sigLabel);

fn = fieldnames(node);
for i = 1:numel(fn)
    f = fn{i};
    if isSkippedName(f), continue; end
    v = node.(f);
    p = joinPath(pth, f);

    if isstruct(v)
        if isWindowContainer(f) || ~isscalar(v)
            for k = 1:numel(v)
                D = walkStruct(D, v(k), sprintf('%s(%d)',p,k), ctx, sigLabel);
            end
        else
            D = walkStruct(D, v, p, ctx, sigLabel);
        end
        continue
    end
    if istable(v),  D = walkTable(D, v, p, ctx, sigLabel); continue; end
    if ~isnumeric(v) && ~islogical(v), continue; end   % text and cells are labels
    if isempty(v),  continue; end                      % an absent branch
    D(end+1) = describe(p, size(v), ctx, sigLabel); %#ok<AGROW>
end
end

% =====================================================================
function D = walkTable(D, T, pth, ctx, sigLabel)
%walkTable  A metrics table, one descriptor per NUMERIC column.  The string columns
%   (type, label) are the vessel's classification and belong in a legend, not on a
%   y-axis.  The column name goes into the path as its last component, which stays
%   unambiguous because an index component is always name(<digits>).
vn = T.Properties.VariableNames;
for j = 1:numel(vn)
    v = T.(vn{j});
    if ~isnumeric(v) && ~islogical(v), continue; end
    if isempty(v), continue; end
    D(end+1) = describe(joinPath(pth,vn{j}), size(v), ctx, sigLabel); %#ok<AGROW>
end
end

% =====================================================================
function D = walkH5(D, ctx)
%walkH5  The same walk over an HDF5 header, for a file too big to load.
try
    info = h5info(ctx.path, '/results');
catch
    return
end
D = walkH5Group(D, info, '', ctx);
end

function D = walkH5Group(D, g, pth, ctx)
for i = 1:numel(g.Datasets)
    d = g.Datasets(i);
    if isSkippedName(d.Name), continue; end
    if ~isPlainNumeric(d), continue; end
    sz = reshape(d.Dataspace.Size,1,[]);
    if isempty(sz) || any(sz==0), continue; end
    D(end+1) = describe(joinPath(pth,d.Name), padSize(sz,2), ctx, ''); %#ok<AGROW>
end
for i = 1:numel(g.Groups)
    nm = g.Groups(i).Name;
    k = strfind(nm,'/'); nm = nm(k(end)+1:end);
    if isSkippedName(nm), continue; end
    D = walkH5Group(D, g.Groups(i), joinPath(pth,nm), ctx);
end
end

function tf = isPlainNumeric(d)
%isPlainNumeric  Is this dataset a number array, or one of the opaque things?  A
%   struct array is stored as object references; a table, a cell and a string are
%   stored as composites that h5info describes but cannot flatten.
tf = false;
if isfield(d,'Datatype') && isfield(d.Datatype,'Class') && ...
        any(strcmp(d.Datatype.Class,{'H5T_REFERENCE','H5T_COMPOUND','H5T_VLEN'}))
    return
end
cls = 'double';
for a = 1:numel(d.Attributes)
    if strcmp(d.Attributes(a).Name,'MATLAB_class'), cls = char(d.Attributes(a).Value); end
end
tf = ~any(strcmp(cls,{'struct','table','cell','char','string','categorical','datetime','function_handle'}));
end

% ================================================================== ONE LEAF

function d = describe(pth, sz, ctx, sigLabel)
%describe  One leaf: matched, resolved, measured and checked.
d = blankDescriptor();
d.path        = pth;
d.size        = sz;
d.lazy        = ctx.lazy;
d.signalLabel = sigLabel;

A = axesFor(ctx, scopeOf(pth));
r = exSchema('match', pth);

if isempty(r)
    d = inferLeaf(d, A);
    return
end

d.family     = r.family;
d.signal     = r.signal;
d.branch     = r.branch;
d.leaf       = r.leaf;
d.layout     = r.layout;
d.dims       = r.dims;
d.label      = r.label;
d.yUnit      = r.yUnit;
d.pairedWith = r.pairedWith;

d.nUnits = unitCount(sz, d.layout, numel(d.dims));
if strcmp(r.unit,'auto')
    % Whether a myograph <VSM> holds the whole vessel or one item per row across
    % it is a property of the RUN (s.perLine), so only the data can say.
    if d.nUnits<=1, d.unit = 'whole'; else, d.unit = 'line'; end
else
    d.unit = r.unit;
end

d = verify(d, A);
end

% =====================================================================
function n = unitCount(sz, layout, nDims)
%unitCount  How many items a leaf holds, read through its declared layout.  This is
%   the one number that cannot be guessed from position - see exSchema's header.
sz = padSize(sz, 4);
switch layout
    case 'unitFirst', n = sz(1);
    case 'unitLast',  n = sz(nDims+1);
    case 'pixel2D',   n = sz(1)*sz(2);
    case 'xyPairs',   n = sz(1);
    otherwise,        n = 1;                 % 'none': the whole vessel, one item
end
end

% =====================================================================
function d = verify(d, A)
%verify  Do the sizes agree with the axes?  A leaf that does not verify is KEPT and
%   named, because a silent drop is how a real result goes missing.
%
%   The UNIT count is deliberately not checked: seg and dvs have no coordinate
%   vector to disagree with, only a count taken from the same file, so the check
%   would compare a number against itself.  What is checked is the trailing
%   coordinates, and that nothing is left over afterwards.
if ~isempty(d.pairedWith), return; end       % an interval on a scalar, not a curve

sz  = padSize(d.size, 6);
nd  = max(2, numel(d.size));
nDim = numel(d.dims);
switch d.layout
    case 'unitFirst', first = 2;             extra = nDim + 2;
    case 'unitLast',  first = 1;             extra = nDim + 2;
    case 'pixel2D'
        first = 3;                           extra = nDim + 3;
        if isfield(A,'pixel')
            if ~isequal(sz(1:2), reshape(A.pixel.n,1,[]))
                d = flag(d, sprintf('image is %s, the recording''s images are %s', ...
                    mat2str(sz(1:2)), mat2str(reshape(A.pixel.n,1,[]))));
            end
        else
            d = flag(d, 'the recording declares no image size');
        end
    case 'xyPairs'
        % It carries its own x in column 1, so there is nothing to check it
        % against - and checking it against the line axis would be WRONG: on the
        % reference pressure recording nY is 480 while the first window's lagByRow
        % has 475 rows, one per row that produced a usable correlation.
        if sz(2)~=2, d = flag(d, sprintf('expected two columns, found %d', sz(2))); end
        return
    otherwise, first = 1;                    extra = nDim + 1;   % 'none'
end

for i = 1:nDim
    nm = d.dims{i};
    if ~isfield(A, nm)
        d = flag(d, sprintf('this recording declares no %s axis', nm));
        continue
    end
    want = A.(nm).n;
    got  = sz(first+i-1);
    if got ~= want
        d = flag(d, sprintf('%s is %d long here but the %s axis has %d', ...
            d.leaf, got, nm, want));
    end
end

for k = extra:nd
    if sz(k)~=1
        d = flag(d, sprintf('has a dimension of length %d that no axis explains', sz(k)));
    end
end
end

function d = flag(d, why)
d.suspect = true;
if isempty(d.note), d.note = why; else, d.note = [d.note '; ' why]; end
end

% =====================================================================
function d = inferLeaf(d, A)
%inferLeaf  A leaf no row claims.  The unit comes from the nearest ancestor that
%   has one - which in practice is the SIGNAL token in the path, since the signal is
%   the ancestor that owns the units - and the dims by matching the remaining sizes
%   against the axis registry.  It goes under 'Other' and is marked, so a branch
%   that grew without a schema row is visible rather than absent.
d.inferred = true;
d.family   = 'Other';
c = pathComponents(d.path);
d.leaf     = c{end};
d.label    = c{end};

sz = padSize(d.size, 4);
[d.unit, d.signal] = inferUnit(c, sz, A);

switch d.unit
    case 'pixel', d.layout = 'pixel2D'; first = 3;
    case 'whole', d.layout = 'none';    first = 1;
    otherwise,    d.layout = 'unitFirst'; first = 2;
end
d.nUnits = unitCount(sz, d.layout, 0);

nd = max(2, numel(d.size));
for k = first:nd
    nm = axisOfLength(A, sz(k));
    if ~isempty(nm), d.dims{end+1} = nm; end
end
end

function [unit, sig] = inferUnit(comps, sz, A)
%inferUnit  The signal token names the unit whenever the path has one.
unit = 'whole'; sig = '';
for i = 1:numel(comps)
    switch comps{i}
        case {'sData','gsData'},        unit = 'seg';   sig = comps{i}; return
        case {'dvsData','dvsDiameter'}, unit = 'dvs';   sig = comps{i}; return
        case 'ppx',                     unit = 'pixel'; sig = comps{i}; return
    end
end
% No signal in the path: fall back to the shape itself.
if isfield(A,'pixel') && numel(sz)>=2 && isequal(sz(1:2),reshape(A.pixel.n,1,[]))
    unit = 'pixel'; return
end
for nm = {'seg','dvs','line'}
    if isfield(A,nm{1}) && isscalar(A.(nm{1}).n) && sz(1)==A.(nm{1}).n
        unit = nm{1}; return
    end
end
end

function nm = axisOfLength(A, n)
%axisOfLength  Which coordinate is this dimension?  Ordered axes are preferred:
%   a length that matches both a time base and a segment count is far more likely
%   to be the time base, and a label id is not something to plot against anyway.
nm = '';
if n<=1, return; end
order = {'time','timeWT','timeDWT','gsTime','f','pctLevel','pctBin','harmonic', ...
         'plane','line','seg','dvs'};
for i = 1:numel(order)
    if isfield(A,order{i}) && isscalar(A.(order{i}).n) && A.(order{i}).n==n
        nm = order{i}; return
    end
end
end

% ================================================================== AXES

function A = axesFor(ctx, scope)
%axesFor  The axis registry at one scope, remembered.  exScan asks for the same
%   dozen scopes over and over - once per leaf - and on a lazy file each miss is an
%   HDF5 header read.
M = ctx.ax;                                 % containers.Map is a handle
if M.isKey(scope), A = M(scope); return; end
A = exAxes(ctx.R, scope, ctx.path);
M(scope) = A;                               %#ok<NASGU> - it mutates the map
end

% ================================================================== NAMES

function tf = isSkippedName(f)
%isSkippedName  Bookkeeping, and the coordinates themselves.  See the header for
%   why each group is here; tStart/tEnd are this module's addition, for the same
%   reason the shape model skips `frames`.
persistent S
if isempty(S)
    S = { 'meta','comments','timeCrop','frames','channels','timeStamp','nY', ...
          'tStart','tEnd', ...
          'time','gsTime','timeWT','timeDWT','f','pctCenters','harmonics' };
end
tf = any(strcmp(f, S));
end

function tf = isWindowContainer(f)
%isWindowContainer  The two struct arrays that hold analysed windows.  They are
%   ALWAYS indexed in the path, even when a recording happens to have exactly one,
%   so a one-window file and a three-window file produce the same shape of path -
%   which is what lets a consumer key on it.
tf = any(strcmp(f, {'intervals','channel'}));
end

function lab = labelOf(node, pth, lab)
%labelOf  The channel's real name.  A wire myograph sanitises it into a field name
%   ('channel 1' -> channel1, 'Channel 4' -> Channel4) and keeps the original both
%   on the channel element and inside the <VSM> it wrote, so either can supply it.
if isfield(node,'channelName') && (ischar(node.channelName) || isstring(node.channelName))
    s = char(string(node.channelName));
    if ~isempty(s), lab = s; end
elseif ~isempty(regexp(pth,'channel\(\d+\)$','once')) && isfield(node,'name')
    s = char(string(node.name));
    if ~isempty(s), lab = s; end
end
end

% ================================================================== PATHS

function p = joinPath(a,b)
if isempty(a), p = b; else, p = [a '.' b]; end
end

function s = scopeOf(pth)
k = strfind(pth,'.');
if isempty(k), s = ''; else, s = pth(1:k(end)-1); end
end

function c = pathComponents(pth)
c = strsplit(char(pth),'.');
for i = 1:numel(c)
    t = regexp(c{i},'^(.+)\(\d+\)$','tokens','once');
    if ~isempty(t), c{i} = t{1}; end
end
end

function sz = padSize(sz, n)
sz = reshape(double(sz),1,[]);
if numel(sz)<n, sz(end+1:n) = 1; end
end

% =====================================================================
function D = sortByPath(D)
%sortByPath  A stable order, so the loaded walk and the h5 walk of one tree line up
%   and a pinned inventory means something.  Indices are zero-padded in the sort KEY
%   only, so intervals(10) does not land between intervals(1) and intervals(2).
if numel(D)<2, return; end
keys = cell(1,numel(D));
for i = 1:numel(D), keys{i} = padIndices(D(i).path); end
[~,ord] = sort(keys);
D = D(ord);
end

function s = padIndices(p)
s = char(p);
[tok,st,en] = regexp(s,'\((\d+)\)','tokens','start','end');
for k = numel(tok):-1:1
    s = [s(1:st(k)-1) sprintf('(%06d)',str2double(tok{k}{1})) s(en(k)+1:end)];
end
end

% =====================================================================
function d = blankDescriptor()
d = struct('path','', 'family','', 'signal','', 'signalLabel','', 'leaf','', ...
    'branch','', 'unit','', 'layout','', 'dims',{{}}, 'size',[], 'nUnits',0, ...
    'kinds',{{}}, 'label','', 'yUnit','', 'suspect',false, 'note','', ...
    'inferred',false, 'lazy',false, 'pairedWith','');
end

function D = emptyDescriptors()
D = blankDescriptor(); D(:) = [];
end
