%exMyograph - What a myograph recording is called, and what its SOURCE can show.
%
%   The variable index reads a results TREE.  Two things a myograph plot needs are
%   not in one, and this module owns both.
%
%   THE NAMES.  A window's own name does not identify it.  myographIntervals
%   flattens a wire recording's channel(i).intervals(k) into one list, and on the
%   four-channel reference file every element of that list is called 'interval1' -
%   four chambers, one name.  Anything keying on the bare name averages them
%   together, so a name here is always qualified by its channel ('interval1 -
%   channel 1'), and the flat index behind a descriptor path is recomposed from
%   both the channel and the window rather than read off the window alone.
%
%   THE THREE VIEWS THAT ARE NOT IN THE RESULTS AT ALL.  A myograph window stores
%   the line-AVERAGED diameter; the per-line detail and the detected walls live once,
%   in the SOURCE beside the results.  They are the only evidence that edge detection
%   found the right edges and that the lines moved together - the assumption a
%   propagation speed is fitted under - and no myograph step writes a report page, so
%   guiExplore is the only place they can be looked at.  They are read here and
%   handed back in the shape the payload wants:
%
%       'diametermap'   the window's per-line diameter as [time x line], decimated in
%                       time so a long recording draws at the size of the axes rather
%                       than at the size of the file
%       'walls'         one frame of the recording with that measure's detected walls
%                       over it, flagged when the frame was invalid
%
%   THE SOURCE IS READ THE SAME WAY THE RESULTS ARE - loaded when it is small, opened
%   as an HDF5 handle when it is not, and sliced rather than loaded either way.  A
%   results file that has been moved away from its triplet simply reports nothing,
%   which is why guiExplore offers these three views only when the _d.mat is there.
%
%   Syntax:
%      n     = exMyograph('count', results)
%      names = exMyograph('names', results)
%      nm    = exMyograph('name', results, k)
%      k     = exMyograph('index', results, dottedPath)
%      S     = exMyograph('signals', results)
%      nm    = exMyograph('measure', fieldName)
%      [M,t,y,unit] = exMyograph('diametermap', results, resultsPath, k, signalName)
%      W     = exMyograph('walls', results, resultsPath, k, signalName)
%
%   INPUTS
%     results      the RESULTS member of a myograph triplet, loaded.
%     k            an index into the FLAT window list, as 'names' returns it.
%     dottedPath   a descriptor path, e.g. 'channel(2).intervals(1).vasomotion.Near'.
%     resultsPath  the *_MYO_r.mat the tree came from; the source is its _d sibling.
%     signalName   one of the names 'signals' reports, e.g. 'outer diameter'.
%     fieldName    a diameter measure's struct field: 'outer' | 'mid' | 'inner'.
%
%   OUTPUTS
%     n      how many analysed windows the recording holds, 0 when it is not a
%            myograph result at all.
%     names  {1 x n} the windows in recording order, each qualified by its channel.
%     nm     one window's name, or one measure's name for a reader.
%     k      the flat window index behind a path, 0 when the path names no window.
%     S      struct array .name .field .kind - what was analysed, as the user should
%            see it ('outer diameter', or a wire channel under the real name
%            LabChart gave it) and as the tree spells it.  .kind is 'measure' or
%            'channel'.
%     M,t,y,unit  the per-line diameter [nTime x nLine], its time base, the image
%            rows it was measured on, and the unit those numbers are in.  M empty
%            when the recording is not beside the result.
%     W      struct .frame .left .right .rows .valid .time - the walls of the middle
%            frame of the window.  .frame empty when the video is not there.
%
%   EXAMPLE
%      S = load('WT794 uden endothel_MYO_r.mat');
%      exMyograph('names', S.results)          % {'interval1','interval2','interval3'}
%      W = exMyograph('walls', S.results, p, 1, 'outer diameter');
%
%   DEPENDS ON
%     myographIntervals (Core/Myograph).
%
% See also: exModality, exScan, exFetch, guiExplore, myographIntervals
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function varargout = exMyograph(what, varargin)

switch lower(char(what))
    case 'count',       varargout{1} = numel(flatWindows(varargin{1}));
    case 'names',       varargout{1} = windowNames(varargin{:});
    case 'name',        varargout{1} = windowName(varargin{:});
    case 'index',       varargout{1} = flatIndexOf(varargin{:});
    case 'signals',     varargout{1} = signalList(varargin{:});
    case 'measure',     varargout{1} = measureLabel(varargin{:});
    case 'diametermap', [varargout{1:4}] = diameterMap(varargin{:});
    case 'walls',       varargout{1} = wallFrame(varargin{:});
    otherwise
        error('exMyograph:what','''%s'' is not something this module answers.', char(what));
end
end

% ===================================================================== THE WINDOWS

function iv = flatWindows(R)
%flatWindows  THE WINDOWS OF A MYOGRAPH RECORDING, in the one shape this module
%   reads.  A pressure myograph keeps them flat in results.intervals; a wire myograph
%   splits them by CHANNEL, because one LabChart file is several chambers.
%   myographIntervals knows about both and hands back one list with .channelName
%   written onto every element, so every k here is an index into THAT list and
%   nothing below has to ask which myograph it is looking at.
iv = [];
if isempty(R) || ~isstruct(R) || isfield(R,'x_lazy_'), return; end
if ~isfield(R,'intervals') && ~isfield(R,'channel'), return; end
iv = myographIntervals(R);
end

function iv = windowAt(R,k)
%windowAt  One window of the flat list, or [] when the index is not one.
iv = [];
IV = flatWindows(R);
if k>=1 && k<=numel(IV), iv = IV(k); end
end

function nm = windowNames(R)
%windowNames  What the analysed windows are called, in recording order.  A wire
%   myograph names them PER CHANNEL: two chambers whose baseline window is called
%   'baseline' are two different windows, and a filter that could not tell them apart
%   would silently average the chambers together.
nm = {};
iv = flatWindows(R);
for k = 1:numel(iv), nm{end+1} = windowName(R,k); end %#ok<AGROW>
end

function nm = windowName(R,k)
%windowName  One window's name, or what it is when nobody named it, qualified by the
%   channel it belongs to when there is one.
iv = flatWindows(R);
nm = '';
if k<1 || k>numel(iv), return; end
if isfield(iv(k),'name'), nm = char(string(iv(k).name)); end
if isempty(nm), nm = sprintf('interval %d',k); end
ch = '';
if isfield(iv(k),'channelName'), ch = char(string(iv(k).channelName)); end
if ~isempty(ch), nm = [nm ' - ' ch]; end
end

function k = flatIndexOf(R, pth)
%flatIndexOf  Which of the FLAT windows a descriptor path belongs to.  A pressure
%   myograph keeps them in results.intervals(k) and the answer is k; a wire myograph
%   hangs them off results.channel(i).intervals(k), and myographIntervals walks the
%   channels in tree order - so the flat index is every earlier channel's windows plus
%   k.  Two chambers both calling their first window 'interval1' are different
%   windows, and anything keying on the bare name would average them together.
k = 0;
tI = regexp(char(pth),'intervals\((\d+)\)','tokens','once');
if isempty(tI), return; end
iv = str2double(tI{1});
tC = regexp(char(pth),'channel\((\d+)\)','tokens','once');
if isempty(tC), k = iv; return; end
ch = str2double(tC{1});
off = 0;
if isstruct(R) && isfield(R,'channel') && isstruct(R.channel)
    for i = 1:min(ch-1, numel(R.channel))
        if isfield(R.channel(i),'intervals'), off = off + numel(R.channel(i).intervals); end
    end
end
k = off + iv;
end

% ===================================================================== THE SIGNALS

function S = signalList(R)
%signalList  WHAT WAS ANALYSED in this recording, as the user should see it and as
%   the tree stores it.  A pressure myograph analyses diameter MEASURES; a wire
%   myograph analyses CHANNELS, whose real names ('Force 1 (mN)') cannot be struct
%   field names and are therefore kept beside their trees.  Unioned over the windows
%   in the order they first appear.
S = struct('name',{},'field',{},'kind',{});
IV = flatWindows(R);
for k = 1:numel(IV)
    iv = IV(k);
    for b = {'diameter','propagation','vasomotion'}
        if ~isfield(iv,b{1}), continue; end
        B = iv.(b{1});
        if ~isstruct(B) || ~isscalar(B), continue; end
        if strcmp(b{1},'diameter')
            if ~isfield(B,'stats') || ~isstruct(B.stats), continue; end
            flds = fieldnames(B.stats)';
        else
            flds = fieldnames(B)';
        end
        for j = 1:numel(flds)
            nm = measureLabel(flds{j}); kind = 'measure';
            if strcmp(b{1},'vasomotion') && isstruct(B.(flds{j})) && ...
                    isfield(B.(flds{j}),'channelName') && ~isempty(B.(flds{j}).channelName)
                nm = char(string(B.(flds{j}).channelName)); kind = 'channel';
            end
            if any(strcmp(nm,{S.name})), continue; end
            S(end+1) = struct('name',nm,'field',flds{j},'kind',kind); %#ok<AGROW>
        end
    end
end
end

function nm = measureLabel(field)
%measureLabel  The diameter measures, named for a reader rather than for the code.
switch char(field)
    case 'outer', nm = 'outer diameter';
    case 'mid',   nm = 'wall-centre diameter';
    case 'inner', nm = 'luminal diameter';
    otherwise,    nm = char(field);
end
end

function f = signalField(R,name)
%signalField  The tree field behind a selection, '' when this recording has none.
%   A CALLER MAY HAND EITHER SPELLING.  The cascade keys its Signal step on the tree
%   field ('mid'), so that the diameter the tree stores and the per-line views read
%   from the recording land on ONE entry rather than two; a caller working from a
%   list the user saw has the display name ('wall-centre diameter').  Both resolve.
f = '';
S = signalList(R);
k = find(strcmp({S.name}, name),1);
if ~isempty(k), f = S(k).field; return; end
k = find(strcmp({S.field}, name),1);
if ~isempty(k), f = S(k).field; end
end

function meas = measureNames(R)
%measureNames  The diameter measures this recording carries, in the order they were
%   measured - which is the order of the third dimension of source.data.
meas = {};
iv = flatWindows(R);
for k = 1:numel(iv)
    d = iv(k).diameter;
    if isstruct(d) && isscalar(d) && isfield(d,'stats') && isstruct(d.stats)
        meas = fieldnames(d.stats)'; return
    end
end
end

function m = measureIndex(R,name)
%measureIndex  Which plane of the per-line arrays a selection is, read off the tree's
%   own order rather than assumed.  Falls back to the wall centre.
m = min(2,max(1,numel(measureNames(R))));
k = find(strcmp(measureNames(R), signalField(R,name)),1);
if ~isempty(k), m = k; end
end

% ====================================================================== THE SOURCE

function S = sourceHandle(rPath)
%sourceHandle  The SOURCE beside a myograph RESULTS file: the whole struct when it is
%   small, an HDF5 handle when it is not - the same rule the results themselves
%   already follow.
S = [];
d = strrep(char(rPath),'_r.mat','_d.mat');
if ~isfile(d), return; end
info = dir(d);
if ~isempty(info) && info(1).bytes > 1.5e9
    S = struct('x_lazy_',true,'x_path_',d);
else
    L = load(d,'source');
    if isfield(L,'source'), S = L.source; end
end
end

function v = sourceField(S,name)
%sourceField  One small SOURCE field (a time base, a row range, a flag).
v = [];
if isempty(S), return; end
if isfield(S,'x_lazy_')
    try, v = h5read(S.x_path_, ['/source/' name]); catch, end
else
    if isfield(S,name), v = S.(name); end
end
end

function sz = sourceSize(S,name)
%sourceSize  The shape of one big SOURCE array, without reading it.
sz = [];
if isempty(S), return; end
if isfield(S,'x_lazy_')
    try
        I = h5info(S.x_path_, ['/source/' name]); sz = I.Dataspace.Size;
    catch
    end
else
    if isfield(S,name), sz = size(S.(name)); end
end
end

function A = sourceSlab(S,name,start,count,stride)
%sourceSlab  One block of a big SOURCE array: `count` elements per dimension, spaced
%   `stride` apart, starting at `start`.  Read straight out of the file when the
%   SOURCE is a handle, indexed out of the struct when it is not.
A = [];
if isempty(S), return; end
if isfield(S,'x_lazy_')
    try, A = h5read(S.x_path_, ['/source/' name], start, count, stride); catch, end
    return
end
if ~isfield(S,name), return; end
subs = cell(1,numel(start));
for i = 1:numel(start)
    subs{i} = start(i) : stride(i) : start(i)+(count(i)-1)*stride(i);
end
A = S.(name)(subs{:});
end

function p = recordingPath(S,rPath)
%recordingPath  THE recording this product was made from: the path the entry step
%   wrote down, else the video sitting beside the product - the same order the
%   diameter step uses, so this finds one exactly when the pipeline would.
p = '';
v = sourceField(S,'fName');
if ~isempty(v)
    if isnumeric(v), v = char(v(:)'); end
    p = char(string(v));
end
if ~isempty(p) && isfile(p), return; end
[fPath,stem] = fileparts(regexprep(char(rPath),'_MYO_r\.mat$',''));
for ext = {'.avi','.mp4','.mov','.mkv'}
    cand = fullfile(fPath,[stem ext{1}]);
    if isfile(cand), p = cand; return; end
end
p = '';
end

function rows = measuredRows(S,nY)
%measuredRows  The image rows the vessel was actually measured on.  Rows outside them
%   carry an interpolated fill and are not a measurement, so nothing draws them.
rr = double(sourceField(S,'rowRange'));
rows = 1:nY;
if numel(rr)==2
    r0 = max(1,round(rr(1))); r1 = min(nY,round(rr(2)));
    if r1>=r0, rows = r0:r1; end
end
end

function [M,t,y,unit] = diameterMap(R,rPath,k,name)
%diameterMap  The per-line diameter of one window as [time x line], decimated in time
%   so a long recording draws at the size of the axes rather than of the file.
M = []; t = []; y = []; unit = 'px';
S = sourceHandle(rPath);
if isempty(S), return; end
iv = windowAt(R,k);
if isempty(iv) || ~isfield(iv,'frames') || numel(iv.frames)~=2, return; end
i0 = double(iv.frames(1)); i1 = double(iv.frames(2));
sz = sourceSize(S,'data');
if numel(sz)<3 || i1>sz(1) || i1<i0, return; end
rows = measuredRows(S,sz(2));
m = min(measureIndex(R,name), sz(3));
nT = i1-i0+1;
dec = max(1, ceil(nT/2000));
n = floor((nT-1)/dec)+1;
M = double(sourceSlab(S,'data',[i0 rows(1) m],[n numel(rows) 1],[dec 1 1]));
if isempty(M), return; end
tAll = double(sourceField(S,'time'));
if numel(tAll)>=i1, t = tAll(i0:dec:i0+(n-1)*dec); else, t = (0:n-1)'; end
t = t(:)'; y = rows;
if ~isempty(sourceField(S,'pixelSize')), unit = 'px'; end
end

function W = wallFrame(R,rPath,k,name)
%wallFrame  One frame of the recording, with the detected walls of the selected
%   measure over it.  The MIDDLE frame of the window, as the one most likely to be
%   representative of it, and the walls come back flagged when that frame was invalid
%   - a wall had left the field of view, so its diameter is a lower bound.
W = struct('frame',[],'left',[],'right',[],'rows',[],'valid',true,'time',NaN);
S = sourceHandle(rPath);
if isempty(S), return; end
iv = windowAt(R,k);
if isempty(iv) || ~isfield(iv,'frames') || numel(iv.frames)~=2, return; end
sz = sourceSize(S,'wallL');
if numel(sz)<3, return; end
fi = min(max(1,round(mean(double(iv.frames)))), sz(1));
rows = measuredRows(S,sz(2));
m = min(measureIndex(R,name), sz(3));
W.left  = double(sourceSlab(S,'wallL',[fi rows(1) m],[1 numel(rows) 1],[1 1 1])); W.left = W.left(:);
W.right = double(sourceSlab(S,'wallR',[fi rows(1) m],[1 numel(rows) 1],[1 1 1])); W.right = W.right(:);
W.rows  = rows(:);
if numel(W.left)~=numel(W.rows) || isempty(W.left), W = resetWalls(W); return; end
tAll = double(sourceField(S,'time'));
if numel(tAll)>=fi, W.time = tAll(fi); end
vv = sourceField(S,'valid');
if ~isempty(vv) && numel(vv)>=fi, W.valid = logical(vv(fi)); end
video = recordingPath(S,rPath);
if isempty(video), return; end
try
    vr = VideoReader(video);
    vr.CurrentTime = min(max(W.time,0), max(0, vr.Duration - 1/max(vr.FrameRate,eps)));
    W.frame = readFrame(vr);
catch
    W.frame = [];
end
end

function W = resetWalls(W)
W.left = []; W.right = []; W.rows = []; W.frame = [];
end
