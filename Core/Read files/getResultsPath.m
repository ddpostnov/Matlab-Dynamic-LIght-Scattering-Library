%getResultsPath - Where a file goes when the results live apart from the recordings.
%
%   A project may keep its RESULTS somewhere other than its RECORDINGS.  The user
%   names a results folder beside the root folder, and everything this library
%   writes - the _d/_r/_s triplets, the report pages, the assembled documents -
%   goes there instead of beside the raw data.  THE SUBPATH IS MIRRORED WHOLE, so
%   <root>\A\B\Mouse1.rls becomes <results>\A\B\Mouse1.rls and only the root part
%   ever changes.  This is the one place that rule is written.
%
%   IT DEFAULTS TO NOTHING HAPPENING.  A project that never names a results folder
%   - which is every project on disk today - has an empty resultsFolder, or one
%   equal to rootFolder, and gets its own path straight back.  The mapping is
%   opt-in, and no call site branches on whether it was taken.
%
%   THE RESULTS FOLDER IS TESTED BEFORE THE ROOT FOLDER, ALWAYS.  The ordinary
%   choice is a results folder INSIDE the root folder (<root>\Results), and then
%   every path under the results folder is also under the root folder.  Asking
%   "is it in the results tree already?" first is what stops <root>\Results\A
%   becoming <root>\Results\Results\A, and it is what makes the map idempotent:
%   a downstream step, whose fName is already a product path, maps to itself, so
%   nothing has to know which kind of step it is.
%
%   A PATH UNDER NEITHER FOLDER IS LEFT WHERE IT IS.  A hand-listed absolute
%   fNames{1} writes beside itself, exactly as it did before this existed.  The
%   alternative - inventing a subpath for it - would collide two same-stemmed
%   recordings from two different disks into one results folder.
%
%   AND IT CREATES THE FOLDER IT NAMES, when the mapping actually moved the path
%   and only then.  A get* with a side effect is a smell; the alternative is the
%   same two lines of mkdir at every call site, and the only reason to ask this
%   question is to write there.  fName is a FILE path, so the folder made is the
%   one that would contain it.  THE INVERSE CREATES NOTHING: 'root' is asked in
%   order to READ a recording that already exists, and a raw-data folder conjured
%   by a lookup is litter.
%
%   TWO CALLERS ASK THIS, AND ONLY ONE OF THEM IS ABOUT TO WRITE.  A wrapper asks
%   where its product GOES, and wants the room made and a throw if it cannot be -
%   that is 'make', the default, and it is every writing call site.  The workbench
%   asks where a recording's products WOULD be so it can go and LOOK for them, and
%   it asks once per file on every repaint of the Files tab: making a folder there
%   would litter the results tree with one empty folder per scanned recording, cost
%   a disk stat per file per repaint, and let an unreachable results folder throw
%   out of an edit box.  That is 'name'.  It is the same rule and the same
%   arithmetic - the only difference is whether the answer is prepared for.
%
% Syntax:
%    p = getResultsPath(fName, s)
%    p = getResultsPath(fName, s, direction)
%    p = getResultsPath(fName, s, direction, mode)
%    p = getResultsPath(fName, rootFolder, resultsFolder)
%    p = getResultsPath(fName, rootFolder, resultsFolder, direction)
%    p = getResultsPath(fName, rootFolder, resultsFolder, direction, mode)
%
% Inputs:
%    fName         - char/string FILE path.  A cell array is mapped elementwise
%                    and keeps its size and shape; empty entries stay empty.
%    s             - a struct carrying the two folders as .rootFolder and
%                    .resultsFolder: the wrapper's parameter struct, or the
%                    reporting context from reportOpen, which resolves the same
%                    two fields.  Either may be absent or empty.
%    rootFolder    - where the RECORDINGS are.
%    resultsFolder - where the RESULTS go.  Empty, or equal to rootFolder, means
%                    no mapping at all.
%    direction     - 'results' (default) maps recordings -> results;
%                    'root' is the inverse, results -> recordings, for the two
%                    places that go looking for a product's raw recording.  In
%                    the two-folder form the direction is the 4th argument, so
%                    that a results folder literally named 'root' is unambiguous.
%    mode          - 'make' (default) creates the folder the mapped file would
%                    need, and throws when it cannot; 'name' only names it, for a
%                    caller deciding where to LOOK rather than where to write.
%                    Ignored going 'root', which never creates anything.
%
% Outputs:
%    p - the mapped path, or the input unchanged, or a cell array of them.
%
% Examples:
%    outBase = getResultsPath(s.fName, s);        % <root>\A\M.rls -> <results>\A\M.rls
%    save(strrep(outBase,'.rls','_t_K_d.mat'),'source','-v7.3');
%    rawName = getResultsPath(productName, s, 'root');
%    where   = getResultsPath(p, root, results, 'results', 'name');  % just tell me
%
% See also: getProductPath, getFileNamesList, reportSave, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function p = getResultsPath(fName, arg2, arg3, arg4, arg5)

if nargin < 2
    error('getResultsPath:noFolders', ...
        ['Both folders are required: getResultsPath(fName,s) or ' ...
         'getResultsPath(fName,rootFolder,resultsFolder).']);
end

% ---- the two call forms ----------------------------------------------------
direction = 'results';
mode      = 'make';
if isstruct(arg2)
    if nargin > 4
        error('getResultsPath:tooManyInputs', ...
            'The struct form takes a direction and a mode and nothing else after them.');
    end
    rootFolder    = folderField(arg2,'rootFolder');
    resultsFolder = folderField(arg2,'resultsFolder');
    if nargin > 2, direction = arg3; end
    if nargin > 3, mode      = arg4; end
else
    if nargin < 3
        error('getResultsPath:noResultsFolder', ...
            'The two-folder form needs both rootFolder and resultsFolder.');
    end
    rootFolder    = asChar(arg2);
    resultsFolder = asChar(arg3);
    if nargin > 3, direction = arg4; end
    if nargin > 4, mode      = arg5; end
end

direction = lower(asChar(direction));
if ~any(strcmp(direction,{'results','root'}))
    error('getResultsPath:badDirection', ...
        'Direction must be ''results'' or ''root'', not ''%s''.', direction);
end
mode = lower(asChar(mode));
if ~any(strcmp(mode,{'make','name'}))
    error('getResultsPath:badMode', ...
        'Mode must be ''make'' or ''name'', not ''%s''.', mode);
end

% ---- rule 1: no results folder, no mapping ---------------------------------
% The default, and the whole reason an existing project sees no change: the input
% comes back verbatim, whatever its type, and no folder is created.
rootFolder    = tidyRoot(rootFolder);
resultsFolder = tidyRoot(resultsFolder);
if isempty(rootFolder) || isempty(resultsFolder) || samePath(rootFolder,resultsFolder)
    p = fName;
    return
end

if iscell(fName)
    p = cell(size(fName));
    for idx = 1:1:numel(fName)
        if isempty(fName{idx}), p{idx} = fName{idx}; continue; end
        p{idx} = oneName(fName{idx}, rootFolder, resultsFolder, direction, mode);
    end
    return
end
p = oneName(fName, rootFolder, resultsFolder, direction, mode);
end

% =====================================================================
function p = oneName(fName, rootFolder, resultsFolder, direction, mode)
%oneName  Rules 2-4, on ONE path.  THE RESULTS FOLDER IS TESTED FIRST in both
%   directions, which is what makes a results folder nested inside the root folder
%   work: under that arrangement every results path is also a root path, and the
%   deeper folder has to win or the map would either nest twice going out or never
%   come back going in.
p = char(fName);
q = tidyPath(p);

if isUnder(q, resultsFolder)
    if strcmp(direction,'root')
        p = swapRoot(q, resultsFolder, rootFolder);
    end
    return                          % rule 2: already there, unchanged
end

if strcmp(direction,'results') && isUnder(q, rootFolder)
    p = ensureFolder(swapRoot(q, rootFolder, resultsFolder), mode);   % rule 3
end
                                    % rule 4: under neither, unchanged
end

% =====================================================================
function p = ensureFolder(p, mode)
%ensureFolder  The folder the mapped FILE goes in (see the header).  A failure here
%   is raised rather than swallowed: it means the results folder cannot be written
%   to at all, which is worth hearing at the naming site instead of as whatever the
%   save says three lines later.  reportSave, which must never kill a run, makes its
%   call inside the try it already had, and a caller that is only deciding where to
%   LOOK asks for 'name' and never reaches the disk here at all.
if strcmp(mode,'name'), return, end
d = fileparts(p);
if isempty(d) || isfolder(d), return, end
[ok,msg] = mkdir(d);
if ~ok
    error('getResultsPath:cannotCreate', ...
        'Could not create the results folder "%s": %s', d, msg);
end
end

% =====================================================================
function tf = isUnder(p, root)
%isUnder  p is root itself, or something inside it.  THE BOUNDARY IS A SEPARATOR,
%   not a character count: '\\server\shareX\y' is not under '\\server\share', and a
%   bare prefix test would say it is.
n  = numel(root);
tf = false;
if isempty(root) || numel(p) < n, return, end
if ~samePath(p(1:n), root),       return, end
tf = numel(p)==n || p(n+1)==filesep;
end

% =====================================================================
function p = swapRoot(p, fromRoot, toRoot)
%swapRoot  Only the root part changes; the subpath is carried over whole.
rest = p(numel(fromRoot)+1:end);
while ~isempty(rest) && rest(1)==filesep
    rest(1) = [];
end
if isempty(rest)
    p = toRoot;
else
    p = fullfile(toRoot, rest);
end
end

% =====================================================================
function q = tidyPath(p)
%tidyPath  One spelling of the separators, so a hand-typed 'C:/Data' compares
%   against a composed 'C:\Data'.  ONLY on Windows: a backslash is a legal
%   character in a POSIX file name, so translating it there would rewrite names.
q = char(p);
if ispc, q = strrep(q,'/',filesep); end
end

% =====================================================================
function q = tidyRoot(p)
%tidyRoot  A folder without its trailing separator(s) - '<root>\' and '<root>' are
%   the same folder, and a user types either.
q = tidyPath(asChar(p));
while ~isempty(q) && q(end)==filesep
    q(end) = [];
end
end

% =====================================================================
function tf = samePath(a,b)
%samePath  Windows file names differ in case without differing.
if ispc
    tf = strcmpi(a,b);
else
    tf = strcmp(a,b);
end
end

% =====================================================================
function v = folderField(s, name)
v = '';
if isfield(s,name) && ~isempty(s.(name))
    v = asChar(s.(name));
end
end

% =====================================================================
function c = asChar(v)
if isempty(v)
    c = '';
    return
end
c = char(v);
end
