%exMyograph - What a myograph recording's analysed windows are called, and which one
%a descriptor path belongs to.
%
%   The variable index reads a results TREE, and a myograph keeps its results inside
%   the WINDOWS it was analysed in.  A window is not a leaf, so nothing in the index
%   can say what one is called, which chamber it was recorded on, or which of them a
%   path belongs to - and those three, plus the reader's name for a diameter measure,
%   are the whole of what this module answers.
%
%   IT READS THE RESULTS AND NOTHING ELSE.  This module used to reach into the _d.mat
%   beside the result for the per-line diameter, and into the recording itself for a
%   frame with the detected walls drawn on it.  Both are gone.  The per-line diameter
%   is now a RESULT - runMyographDiameter writes results.intervals(k).diameter.lines,
%   one value per row across the vessel per frame - so the index finds it like any
%   other leaf and nothing has to be read from a sibling.  The wall frame cannot be
%   rescued the same way and is not: it is a frame of the RECORDING, and a video
%   frame is source by definition.  What is lost HERE is the "did it find the right
%   edges" picture; what survives here is the position-time kymograph, which still
%   shows a wall that wandered as a stripe, and .stats.<measure>.measuredFraction /
%   .validFraction, which carry the same fact as numbers.
%
%   THE PICTURE ITSELF STILL EXISTS - getMyographWallFrame (Core/Myograph) re-opens
%   the recording and hands back one frame with that frame's own walls on it - and
%   this module deliberately does not offer it.  The argument above was never that
%   the picture is worthless; it is that THIS tool reads a results file and a video
%   frame is not in one.  A step that is already holding the recording open is in a
%   different position, so the interval editor draws it while the operator decides
%   where the windows go, which is where that check belongs anyway.  Adding a
%   variable here that quietly opened an '.avi' would undo the rule this whole
%   module is built on.
%
%   THE NAMES.  ADDRESSING A WINDOW AND NAMING IT ARE TWO DIFFERENT JOBS, and this
%   module keeps them apart.  myographIntervals flattens a wire recording's
%   channel(i).intervals(k) into one list, and on the four-channel reference file
%   every element of that list is called 'interval1' - four chambers, one name.
%   Anything keying on the bare name to ADDRESS a window would average four chambers
%   together, so the flat index behind a descriptor path is recomposed from both the
%   channel and the window rather than read off the window alone.
%
%   BUT A NAME IS NOT AN ADDRESS.  What a window's name is FOR is the protocol step -
%   baseline, drug, washout - and that is an axis a reader wants on an x, shared
%   across chambers, recordings and animals.  The chamber is already separated by the
%   SIGNAL axis, because a wire recording's signal IS its channel under the real name
%   LabChart gave it; qualifying the name with it as well separated the chamber twice
%   and left the protocol axis with nothing to be.  So 'name' is the window's own
%   BARE name, 'channel' is a question of its own, and a recording that has no window
%   called 'washout' contributes nothing to that x-category rather than being padded.
%
%   WHY THIS IS STILL A FILE OF ITS OWN, now that the source reads are gone and five
%   questions are left.  Two of them could not live anywhere else: 'index' parses a
%   DESCRIPTOR PATH, which is the variable index's idea and not the recording's, and
%   'measure' spells a struct field for a biologist, which is a GUI voice and not
%   something Core should own.  The other three are one-line reads off
%   myographIntervals, and they are here so that the split this header argues for -
%   a bare name to label with, a channel to separate by, a flat index to address with
%   - is visible in ONE place rather than re-derived at each call site.  Two questions
%   that used to be here are gone instead of moved: 'signals' and 'count' were read
%   only by the three invented variables, and 'signals' was a SECOND answer to the
%   question exScan's labelOf already answers off channelName - the drift this whole
%   redesign exists to prevent.
%
%   Syntax:
%      names = exMyograph('names', results)
%      nm    = exMyograph('name', results, k)
%      ch    = exMyograph('channel', results, k)
%      k     = exMyograph('index', results, dottedPath)
%      nm    = exMyograph('measure', fieldName)
%
%   INPUTS
%     results      the RESULTS member of a myograph PAIR, loaded.  There is no source
%                  member: a myograph product is _MYO_r + _MYO_s and nothing else.
%     k            an index into the FLAT window list, as 'names' returns it.
%     dottedPath   a descriptor path, e.g. 'channel(2).intervals(1).vasomotion.Near'.
%     fieldName    a diameter measure's struct field: 'outer' | 'mid' | 'inner'.
%
%   OUTPUTS
%     names  {1 x n} the windows in recording order, one entry per window and each
%            its BARE name - so on a wire recording the same protocol step appears
%            once per channel, spelled identically.  A caller building an axis out
%            of them unions; a caller addressing one uses the index.
%     nm     one window's bare name, or one measure's name for a reader.
%     ch     the channel a window was analysed on, '' on a pressure myograph.
%     k      the flat window index behind a path, 0 when the path names no window.
%
%   EXAMPLE
%      S = load('WT794 uden endothel_MYO_r.mat');
%      exMyograph('names', S.results)          % {'interval1','interval2','interval3'}
%      exMyograph('channel', S.results, 1)     % '' - a pressure myograph has none
%
%   DEPENDS ON
%     myographIntervals (Core/Myograph).
%
% See also: exModality, exScan, exFetch, guiExplore, myographIntervals
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function varargout = exMyograph(what, varargin)

switch lower(char(what))
    case 'names',       varargout{1} = windowNames(varargin{:});
    case 'name',        varargout{1} = windowName(varargin{:});
    case 'channel',     varargout{1} = windowChannel(varargin{:});
    case 'index',       varargout{1} = flatIndexOf(varargin{:});
    case 'measure',     varargout{1} = measureLabel(varargin{:});
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
if isempty(R) || ~isstruct(R), return; end
if ~isfield(R,'intervals') && ~isfield(R,'channel'), return; end
iv = myographIntervals(R);
end

function nm = windowNames(R)
%windowNames  What the analysed windows are called, in recording order - ONE ENTRY
%   PER WINDOW, so the k'th entry is the k'th window and an addressing caller can
%   still index it.  On a wire recording the same protocol step therefore appears
%   once per channel under one identical spelling, which is the point: a caller
%   building the protocol axis unions them and gets the step back, once.
nm = {};
iv = flatWindows(R);
for k = 1:numel(iv), nm{end+1} = windowName(R,k); end %#ok<AGROW>
end

function nm = windowName(R,k)
%windowName  ONE WINDOW'S BARE NAME, or what it is when nobody named it.  It names
%   the PROTOCOL STEP and nothing else - see the file header on why the channel was
%   taken out of it, and windowChannel for where the channel went.
iv = flatWindows(R);
nm = '';
if k<1 || k>numel(iv), return; end
if isfield(iv(k),'name'), nm = char(string(iv(k).name)); end
if isempty(nm), nm = sprintf('interval %d',k); end
end

function ch = windowChannel(R,k)
%windowChannel  WHICH CHAMBER a window was analysed on, under the real name LabChart
%   gave it; '' on a pressure myograph, which records one vessel and has no channel
%   axis at all.  myographIntervals writes it onto every flat window, so this is one
%   lookup and not a walk of the channel tree.
ch = '';
iv = flatWindows(R);
if k<1 || k>numel(iv), return; end
if isfield(iv(k),'channelName'), ch = char(string(iv(k).channelName)); end
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

% ====================================================================== THE MEASURES

function nm = measureLabel(field)
%measureLabel  The diameter measures, named for a reader rather than for the code.
switch char(field)
    case 'outer', nm = 'outer diameter';
    case 'mid',   nm = 'wall-centre diameter';
    case 'inner', nm = 'luminal diameter';
    otherwise,    nm = char(field);
end
end

