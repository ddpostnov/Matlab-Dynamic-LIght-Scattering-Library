%exAxes - The coordinate vectors one results tree measures things against.
%
%   A = exAxes(results) returns the axis REGISTRY of a whole results struct, and
%   a = exAxes(results, scopePath) returns it as seen from inside a sub-tree.  An
%   axis is a named, shared coordinate vector: the thing that goes on the x of a
%   plot, or the thing a leaf's trailing size has to agree with before the leaf can
%   be believed.  exScan resolves every leaf through this registry, so everything
%   the tool knows about "how long is this dimension and what is it" is here.
%
%   NEAREST ENCLOSING SCOPE WINS.  This is the rule the whole module exists for,
%   and it is load-bearing INSIDE A SPECKLE FILE, not only in a myograph one.  In
%   the reference contrast recording:
%
%       results.vasomotion.timeWT           [2024x1]    timeDWT [289x1]
%       results.vasomotion.gsData.timeWT   [50576x1]    timeDWT [256x1]
%
%   both in the same file, because the guided-contrast signal is sampled on its own
%   time base.  A scanner that resolved axes at the root alone would mis-size every
%   gsData leaf - and gsData's timeVectors really are [50576 x 1480].  A myograph
%   window is the same rule again: intervals(k).vasomotion.<meas> carries its own
%   f / timeWT / timeDWT / pctCenters, and each window has a different one.
%
%   A WINDOW'S OWN TIME BASE AND ROW COUNT ARE DECLARED ON ITS DIAMETER BLOCK, and
%   this function absorbs them at the WINDOW scope rather than at the block's.
%   intervals(k).diameter carries .time (the window's frame times) and .nY (how
%   many rows across the vessel were measured), and both describe the window, not
%   the diameters: propagation's lagByRow is indexed by the same rows and lives one
%   branch over.  A wire recording has diameter as an empty placeholder, so nothing
%   is absorbed and neither axis exists - which is correct, since a wire myograph
%   measures no rows.
%
%   THERE ARE TWO PERCENTILE AXES AND ONLY ONE OF THEM IS IN THE RESULTS FILE.
%   This is the trap that made the first draft of the shape model wrong:
%
%       scalars.<band>.ampPct     [nUnit x 11]        percentile LEVELS
%       fVectors.VB.ampMeanPct    [nUnit x nF x 10]   percentile BIN CENTRES
%       vasomotion.pctCenters     [1 x 10]            only the BINS are stored
%
%   pctCenters is the midpoints of the levels (runVasomotion.m:696), so it is one
%   shorter BY CONSTRUCTION.  Treating them as one axis marks scalars.<band>.ampPct
%   suspect on every file in the library.  The LENGTH of pctLevel therefore needs no
%   file access at all - it is numel(pctCenters)+1 - and only the VALUES have to be
%   fetched, from the sibling settings file the _d/_r/_s triplet guarantees is
%   there:
%
%       h5read(<stem>_s.mat, '/settings/runVasomotion/pcts')  ->  [0 10 ... 100]
%
%   measured 0.011 s against the 53.8 MB _s.mat of the reference set, with no full
%   load.  That read is DEFERRED until a scope actually has a percentile axis, so a
%   file with no vasomotion tree never opens its settings at all, and it is memoised
%   per file so a scan of a thousand leaves pays for it once.  A myograph writes the
%   same setting under runMyographVasomotion, so both names are tried.  When the
%   sibling is missing, or its pcts no longer has the length the data was written
%   with, the axis falls back to linspace(0,100,n) and is marked .inferred - a
%   non-default s.pcts must never be silently misread.
%
%   THIS IS THE ONE FILE THE EXPLORER OPENS BESIDE THE _r, AND IT IS MEANT TO STAY.
%   D7 makes _r the whole of what the tool reads, and every other sibling read went
%   with it - the _d.mat slices and the VideoReader frames alike.  This one is not
%   of that kind: pcts is a SETTING, not data.  It is the x-axis of
%   scalars.<band>.ampPct and it is not written into the results at all, so without
%   it a percentile curve is drawn on an assumed 0:10:100 and a non-default run is
%   silently misread.  The _d/_r/_s triplet guarantees the file is there, one
%   dataset is read out of it rather than the file, and the fallback above marks the
%   axis inferred when it is not.  Do not delete it under the _r-only rule.
%
%   SEG AND DVS ARE NOT COORDINATES.  A segment index is a label id: there is no
%   meaningful "amplitude vs segment number" plot, and the tool must not offer one.
%   They are registered with .ordered false and no .values, which is what tells
%   exPlotRules that a per-segment scalar has to be reduced over a selection or
%   painted back into space through sMap.  line and interval are the opposite case -
%   position along a vessel and protocol order are both real orderings - and keeping
%   the two classes apart is the point of registering them separately.
%
%   Syntax:
%      A = exAxes(results)
%      A = exAxes(results, scopePath)
%      A = exAxes(results, scopePath, resultsPath)
%      exAxes('flush')
%
%   INPUTS
%     results     the RESULTS member of a triplet, loaded.
%     scopePath   dotted path of the sub-tree to resolve from, '' for the root.
%                 Struct-array elements carry their index: 'intervals(2).vasomotion'.
%     resultsPath the _r.mat this tree came from, so the sibling _s.mat can be
%                 found for the percentile LEVELS.  Optional; without it pctLevel is
%                 inferred.
%
%   OUTPUTS
%     A  struct whose FIELDS ARE THE AXIS NAMES present at that scope.  An axis
%        that this tree does not define is simply absent, so isfield(A,'f') is the
%        test for "does this recording have a frequency axis".  Each field holds:
%          .values    the coordinate vector, [] when only the length is known
%          .n         its length; [Y X] for the pixel axis, which is two-dimensional
%          .ordered   true when the index means something (time, position),
%                     false for a label id (segment, vessel)
%          .scale     'linear' | 'log' | 'categorical' | 'ordinal'
%          .label     the axis title, written for a biologist
%          .inferred  true when the values were guessed rather than read
%
%        The axis names: time gsTime timeWT timeDWT f pctLevel pctBin harmonic
%        plane seg dvs pixel line interval.
%
%   EXAMPLE
%      S = load('LSCI_20240809_1ADCF08BP_t_BFI_r.mat');
%      exAxes(S.results,'vasomotion.gsData').timeWT.n     % 50576
%      exAxes(S.results,'vasomotion.sData').timeWT.n      %  2024
%
%   DEPENDS ON
%     myographIntervals (Core/Myograph) for the interval axis.
%
% See also: exModality, exSchema, exScan, guiExplore, assembleVasomotionTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function A = exAxes(results, scopePath, resultsPath)

if nargin==1 && (ischar(results)||isstring(results)) && strcmpi(char(results),'flush')
    fileCache('flush'); A = struct(); return
end
if nargin<2 || isempty(scopePath),   scopePath   = ''; end
if nargin<3 || isempty(resultsPath), resultsPath = ''; end

A = struct();
if isempty(results) || ~isstruct(results) || ~isscalar(results), return; end

% ---- 1. the axes that are defined once, for the whole recording ---------------
A = rootOnlyAxes(A, results);

% ---- 2. the shared coordinate vectors, root first then every enclosing scope --
scopes = [{''}, prefixesOf(scopePath)];
for i = 1:numel(scopes)
    A = absorbAt(A, results, scopes{i});
end

% ---- 3. the percentile LEVELS, only now that we know whether they are wanted --
A = fillPctLevel(A, resultsPath);
end

% =====================================================================
function A = rootOnlyAxes(A, R)
%rootOnlyAxes  The axes a recording has ONE of: its units, its image, its windows.
%   None of them is redefined deeper in the tree, so resolving them at the root is
%   not a shortcut - a sub-tree that held its own segment count would be a
%   different recording.

% seg / dvs: the metrics table is the authority, because it is the thing that
% NAMES the units; the data matrix is the fallback, for a recording that carries
% the signal without the table.
n = tableHeight(R,'sMetrics');
if isempty(n), n = unitDimOf(R,'sData'); end
if ~isempty(n), A.seg = mkAxis([],n,false,'categorical','segment'); end

n = tableHeight(R,'dvsMetrics');
if isempty(n), n = unitDimOf(R,'dvsData'); end
if isempty(n), n = unitDimOf(R,'dvsDiameter'); end
if ~isempty(n), A.dvs = mkAxis([],n,false,'categorical','vessel'); end

% pixel: the image geometry, [Y X].  sMap first - it is the map every per-segment
% quantity is painted back through - then whatever else is a picture.
for c = {'sMap','mask','imgK','imgBFI','imgI','cMask','pMap','dvsMap'}
    sz = sizeOf(R,c{1});
    if numel(sz)>=2 && all(sz(1:2)>1), A.pixel = mkAxis([],sz(1:2),false,'categorical','pixel'); break; end
end

% plane: how many masks were intersected into commonMask
sz = sizeOf(R,'commonMask');
if numel(sz)>=3 && sz(3)>0, A.plane = mkAxis((1:sz(3))',sz(3),true,'ordinal','mask layer'); end

% gsTime: the guided-contrast signal's own time base, a root-level vector
v = valueOf(R,'gsTime');
if ~isempty(v), A.gsTime = mkAxis(v(:),numel(v),true,'linear','time (s)'); end

% interval: the analysed windows, in protocol order.  A speckle recording has none.
k = numel(myographIntervals(R));
if k>0, A.interval = mkAxis((1:k)',k,true,'categorical','interval'); end
end

% =====================================================================
function A = absorbAt(A, R, scope)
%absorbAt  The shared coordinate vectors declared AT this scope, overriding
%   whatever an enclosing scope said.  Nothing here is conditional on where in the
%   tree we are: a node that carries an 'f' declares the frequency axis, whether it
%   is results.vasomotion or intervals(3).vasomotion.inner.
v = valueOf(R,join2(scope,'time'));
if ~isempty(v), A.time = mkAxis(v(:),numel(v),true,'linear','time (s)'); end

v = valueOf(R,join2(scope,'timeWT'));
if ~isempty(v), A.timeWT = mkAxis(v(:),numel(v),true,'linear','time (s)'); end

v = valueOf(R,join2(scope,'timeDWT'));
if ~isempty(v), A.timeDWT = mkAxis(v(:),numel(v),true,'linear','time (s)'); end

v = valueOf(R,join2(scope,'f'));
if ~isempty(v), A.f = mkAxis(v(:),numel(v),true,'log','frequency (Hz)'); end

v = valueOf(R,join2(scope,'harmonics'));
if ~isempty(v), A.harmonic = mkAxis(v(:),numel(v),true,'ordinal','harmonic'); end

v = valueOf(R,join2(scope,'pctCenters'));
if ~isempty(v)
    A.pctBin = mkAxis(v(:),numel(v),true,'linear','percentile (%)');
    % the LEVELS are one longer by construction; their values come later, or not
    A.pctLevel = mkAxis([],numel(v)+1,true,'linear','percentile (%)');
end

% A WINDOW'S OWN COORDINATES, declared one branch down on its diameter block.
v = valueOf(R,join2(scope,'diameter.time'));
if ~isempty(v), A.time = mkAxis(v(:),numel(v),true,'linear','time (s)'); end
v = valueOf(R,join2(scope,'diameter.nY'));
if isscalar(v) && v>0
    A.line = mkAxis((1:double(v))',double(v),true,'linear','position along the vessel');
end
end

% =====================================================================
function A = fillPctLevel(A, resultsPath)
%fillPctLevel  The percentile levels, read from the sibling settings file - the one
%   place this module deliberately reaches outside the results tree it was handed.
%   Deferred to here so a recording with no percentile axis never opens a second
%   file, and memoised per file so a whole scan pays once.
if ~isfield(A,'pctLevel'), return; end
n = A.pctLevel.n;
p = fileCache('pcts', resultsPath);
if numel(p)==n
    A.pctLevel.values = p(:);
else
    % No sibling, or it no longer describes the run that wrote this tree.  Assume
    % the library default spacing and SAY SO, so a non-default s.pcts cannot be
    % read as if it were 0:10:100.
    A.pctLevel.values   = linspace(0,100,n)';
    A.pctLevel.inferred = true;
end
end

% =====================================================================
function a = mkAxis(values, n, ordered, scale, label)
a = struct('values',values,'n',n,'ordered',logical(ordered), ...
           'scale',char(scale),'label',char(label),'inferred',false);
end

% ================================================== READING A NODE

function v = valueOf(R, dotted)
%valueOf  The VALUE of one leaf, or [] when the path does not lead to one.
v = resolveEager(R, dotted);
if ~isnumeric(v) && ~islogical(v), v = []; end
end

function sz = sizeOf(R, dotted)
%sizeOf  The SIZE of one leaf.
sz = [];
v = resolveEager(R, dotted);
if isnumeric(v) || islogical(v), sz = size(v); end
end

function n = unitDimOf(R, name)
%unitDimOf  How many UNITS a raw signal matrix holds.  sData / dvsData / gsData are
%   written [nT x nUnit] - item LAST - so this is the second dimension and not the
%   first, which is the single most common way to get this tree wrong.
n = []; sz = sizeOf(R,name);
if numel(sz)>=2 && sz(2)>0 && sz(1)>0, n = sz(2); end
end

function n = tableHeight(R, name)
%tableHeight  Rows of a metrics table, or [] when the recording has none.
n = [];
if ~isfield(R,name), return; end
v = R.(name);
if istable(v), n = height(v); end
end

function v = resolveEager(R, dotted)
%resolveEager  Walk a dotted path through a loaded struct.  'intervals(2).x' takes
%   element 2; a component is an INDEX only when it is name(<digits>), which is why
%   a table column called 'std(diameter)' can never be mistaken for one.
v = R;
c = strsplit(char(dotted), '.');
for i = 1:numel(c)
    [nm, idx] = parseComp(c{i});
    if ~isstruct(v) || ~isfield(v,nm), v = []; return; end
    v = v.(nm);
    if ~isempty(idx)
        if ~isstruct(v) || idx>numel(v), v = []; return; end
        v = v(idx);
    elseif isstruct(v) && ~isscalar(v)
        v = []; return
    end
end
end

function [nm, idx] = parseComp(c)
idx = [];
t = regexp(c,'^(.+)\((\d+)\)$','tokens','once');
if isempty(t), nm = c; else, nm = t{1}; idx = str2double(t{2}); end
end

% ================================================== PATHS

function p = prefixesOf(dotted)
%prefixesOf  Every enclosing scope of a path, outermost first, so absorbing them in
%   order leaves the NEAREST one in place.
p = {};
if isempty(dotted), return; end
c = strsplit(char(dotted),'.');
for i = 1:numel(c), p{end+1} = strjoin(c(1:i),'.'); end %#ok<AGROW>
end

function s = join2(a,b)
if isempty(a), s = b; else, s = [a '.' b]; end
end

% ================================================== THE PER-FILE MEMO

function out = fileCache(what, pth)
%fileCache  The percentile levels, remembered per file.
%
%   They are a pure function of a file that nothing in this session writes, and they
%   are asked for once per SCOPE by exScan - a dozen times or more on one tree.  The
%   key carries the file's size and timestamp, so a fixture rewritten between two
%   runs of a test is re-read rather than remembered wrong; exAxes('flush') empties
%   it outright.
persistent pcts
if isempty(pcts),  pcts  = containers.Map('KeyType','char','ValueType','any'); end

if strcmp(what,'flush'), pcts = containers.Map('KeyType','char','ValueType','any');
                         out = []; return
end
out = [];
if isempty(pth), return; end
key = stampOf(pth);
if isempty(key), return; end

if strcmp(what,'pcts')
    if ~pcts.isKey(key), pcts(key) = readSiblingPcts(pth); end
    out = pcts(key);
end
end

function p = readSiblingPcts(resultsPath)
%readSiblingPcts  The percentile LEVELS, from the settings half of the triplet.
%   Two producers write the same setting under two names - runVasomotion for a
%   speckle recording and runMyographVasomotion for a myograph one - so both are
%   tried before giving up.  h5read reads the one dataset, never the file.
p = [];
[d,n,e] = fileparts(char(resultsPath));
if ~strcmpi(e,'.mat') || ~endsWith(n,'_r'), return; end
sPath = fullfile(d,[n(1:end-2) '_s' e]);
if ~isfile(sPath), return; end
for g = {'runVasomotion','runMyographVasomotion'}
    try
        p = h5read(sPath, ['/settings/' g{1} '/pcts']);
        p = double(p(:));
        if ~isempty(p), return; end
    catch
        p = [];
    end
end
end

function key = stampOf(pth)
%stampOf  A cache key that changes when the file does.
key = '';
d = dir(char(pth));
if isempty(d) || d(1).isdir, return; end
key = sprintf('%s|%d|%.6f', lower(char(pth)), d(1).bytes, d(1).datenum);
end
