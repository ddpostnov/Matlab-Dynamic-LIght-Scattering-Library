%wbRename - Rename a WHOLE recording on disk: plan the moves, check them, do them.
%
%   A recording is not one file.  It is the _d/_r/_s triplet of every branch it has
%   produced ('_t_K', '_c_K', '_t_BFI', ...), the report images written beside them,
%   the raw recording they came from, and any workbook exported from them - all of
%   them named after ONE recording identity.  Renaming a single member is therefore
%   not on offer: every wrapper finds its siblings with
%   strrep(fName,'_d.mat','_s.mat'), so a lone _d rename breaks the set the moment
%   anything runs.  This module renames the identity, and everything named after it
%   moves together (author decision D4, 2026-07-31).
%
%   THE NAME GRAMMAR IS wbFileModel'S, and it is never re-parsed here.  The identity
%   is read off the model (.folder, .roiPrefix, .stem), the recording's folder is
%   listed ONCE, and each new name is the old one with that identity prefix replaced
%   - no name is invented any other way.  What belongs to the recording is decided
%   twice over:
%     * a '.mat' belongs when wbFileModel parses it to the SAME identity, so a
%       rename of 'Foo' can never carry away 'Foo_b2_t_K_d.mat' (identity 'Foo_b2')
%       or 'Foo2_t_K_d.mat' (a different recording entirely);
%     * a RAW RECORDING (.rls / .mraw / .cxd / a video) is a recording in its own
%       right, so it belongs only when its name IS the identity - renaming 'Foo'
%       must never move 'Foo_a1.rls';
%     * anything else - a report image, a workbook - belongs when its name is the
%       identity followed by '_'.  Those have no grammar to check, so the boundary
%       character is the whole fence, and the caller is expected to SHOW THE LIST
%       before applying it.
%
%   ALL-OR-NOTHING WHERE IT CAN BE.  'check' refuses the whole plan if any target
%   already exists, if two entries would land on the same name, if a source has gone
%   missing, or if a composed name is not a usable file name.  'apply' re-checks
%   before it moves anything, and if a move still fails part-way it puts the ones
%   that already moved back and says so.
%
% Syntax:
%    list     = wbRename('plan',  model, newStem)
%    [tf,why] = wbRename('check', list)
%    out      = wbRename('apply', list)
%
% Inputs:
%    model   - a wbFileModel struct of ANY file of the recording; it supplies the
%              folder and the identity ([roiPrefix stem]) the plan is built from.
%    newStem - char, the new stem.  The 'RoiN_' crop prefix is part of the identity
%              and is KEPT - what is being replaced is the stem alone.
%    list    - the struct array returned by 'plan'.
%
% Outputs:
%    list - 1xK struct array, one entry per file that would move, folder order:
%             .old .new         full paths
%             .oldName .newName bare names (what a confirmation dialog lists)
%           1x0 when the recording has no files on disk, or when the stem is
%           unchanged (there is then nothing to do).
%    tf   - 'check': true when the plan is safe to apply.
%    why  - 'check': one line saying what is wrong ('' when tf is true).
%    out  - 'apply': struct with
%             .ok         logical, every file moved;
%             .moved      1xK cellstr of the new paths (empty unless .ok);
%             .failed     1xK cellstr of the bare names that would not move;
%             .rolledBack logical, the already-moved files were put back;
%             .why        one line, '' when .ok.
%
% Example:
%    m    = wbFileModel('D:\data\PSY01_a1_t_K_d.mat');
%    list = wbRename('plan', m, 'PSY01_a2');
%    if wbRename('check', list), out = wbRename('apply', list); end
%
% See also: wbFileModel, wbArtifacts, guiWorkbench, wbSession
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function varargout = wbRename(action, varargin)

switch action
    case 'plan',  varargout{1} = planOf(varargin{:});
    case 'check', [varargout{1:nargout}] = checkOf(varargin{:});
    case 'apply', varargout{1} = applyOf(varargin{:});
    otherwise
        error('wbRename:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function list = emptyPlan()
%emptyPlan  The 1x0 shape of a plan, declared once so a caller can index it blind.
list = struct('old',{},'new',{},'oldName',{},'newName',{});
end

% =====================================================================
function list = planOf(model, newStem)
%planOf  What would move, and where to.  ONE directory listing of everything named
%   after the recording, filtered by belongsTo, with the identity prefix swapped.
list = emptyPlan();
if isempty(model) || ~isstruct(model) || ~isfield(model,'folder'), return; end
folder = char(model.folder);
if isempty(folder) || ~isfolder(folder), return; end

base    = [char(model.roiPrefix) char(model.stem)];
newBase = [char(model.roiPrefix) char(newStem)];
if isempty(base) || strcmp(base, newBase), return; end   % nothing to rename

identity = char(model.identity);
d = dir(fullfile(folder,[base '*']));
for i = 1:numel(d)
    if d(i).isdir, continue; end
    [tf, rest, ext] = belongsTo(folder, d(i).name, base, identity);
    if ~tf, continue; end
    newName = [newBase rest ext];
    list(end+1) = struct('old', fullfile(folder,d(i).name), ...
                         'new', fullfile(folder,newName), ...
                         'oldName', d(i).name, 'newName', newName); %#ok<AGROW>
end
end

% =====================================================================
function [tf, rest, ext] = belongsTo(folder, name, base, identity)
%belongsTo  Is this file part of the recording?  Three answers, in order of how much
%   is known about the name:
%     a '.mat'              the name GRAMMAR decides - wbFileModel must parse it to
%                           the same identity, so a rename of 'Foo' can carry away
%                           neither 'Foo2_t_K_d.mat' nor 'Foo_b2_t_K_d.mat';
%     a RAW RECORDING       a recording of its own, so only an exact name match -
%                           renaming 'Foo' must never move 'Foo_a1.rls';
%     anything else         a report image or a workbook written beside the data, so
%                           the identity prefix plus its boundary character is the
%                           whole fence, and the caller shows the list first.
tf = false; rest = '';
[~, bare, ext] = fileparts(name);
if numel(bare) < numel(base) || ~strncmp(bare, base, numel(base)), return; end
rest = bare(numel(base)+1:end);
if ~isempty(rest) && rest(1) ~= '_', return; end
if strcmpi(ext,'.mat')
    m  = wbFileModel(fullfile(folder,name));
    tf = strcmp(m.identity, identity);
elseif isRecordingContainer(ext)
    tf = isempty(rest);
else
    tf = true;
end
end

% =====================================================================
function tf = isRecordingContainer(ext)
%isRecordingContainer  Does this extension hold a RAW RECORDING?  Asked of
%   wbFileModel rather than answered from a list written down here: it narrows the
%   modality vocabulary for every container it recognises (.rls, .mraw, .cxd, a
%   video) and leaves it whole for anything it does not (a .jpg report page, an
%   .xlsx).  So the container list has one definition, and it is the parser's.
tf = numel(wbFileModel('modalities',ext)) < numel(wbFileModel('modalities'));
end

% =====================================================================
function [tf, why] = checkOf(list)
%checkOf  Is this plan safe to apply?  Every reason to refuse is checked BEFORE any
%   file moves, so a refused rename leaves the folder exactly as it was.
tf = false; why = '';
if isempty(list)
    why = 'there is nothing to rename - no file of that recording is on disk.'; return
end
for i = 1:numel(list)
    if ~isUsableName(list(i).newName)
        why = sprintf(['"%s" is not a usable file name - a name cannot contain ' ...
            '\\ / : * ? " < > | or end in a space or a dot.'], list(i).newName); return
    end
end
for i = 1:numel(list)
    if ~isfile(list(i).old)
        why = sprintf('"%s" is no longer on disk - refresh the list and try again.', ...
            list(i).oldName); return
    end
end
if numel(unique({list.new})) < numel(list)
    why = 'two of those files would end up with the same name.'; return
end
sources = {list.old};
for i = 1:numel(list)
    if any(strcmpi(list(i).new, sources)), continue; end   % a case-only rename
    if isfile(list(i).new) || isfolder(list(i).new)
        why = sprintf('"%s" already exists in that folder.', list(i).newName); return
    end
end
tf = true;
end

% =====================================================================
function tf = isUsableName(name)
%isUsableName  A plain file name: no path separator, no character Windows rejects,
%   no trailing space or dot (which Windows silently drops), and not empty.
tf = false;
name = char(name);
if isempty(strtrim(name)), return; end
if any(ismember(name, '\/:*?"<>|')) || any(double(name) < 32), return; end
if name(end)==' ' || name(end)=='.', return; end
tf = true;
end

% =====================================================================
function out = applyOf(list)
%applyOf  Do the moves.  Re-checks first (a plan can be handed here from anywhere),
%   then moves in order; a failure part-way puts the ones that already moved back,
%   so the folder is left either wholly renamed or wholly untouched.
out = struct('ok',false,'moved',{cell(1,0)},'failed',{cell(1,0)}, ...
             'rolledBack',false,'why','');
[tf, why] = checkOf(list);
if ~tf, out.why = why; return; end

done = zeros(1,0);
for i = 1:numel(list)
    [ok, msg] = tryMove(list(i).old, list(i).new);
    if ~ok
        out.failed{end+1} = list(i).oldName;
        out.why = sprintf('could not rename "%s"%s', list(i).oldName, reason(msg));
        out.rolledBack = rollBack(list, done);
        return
    end
    done(end+1) = i; %#ok<AGROW>
end
out.moved = {list.new};
out.ok    = true;
end

% =====================================================================
function [ok, msg] = tryMove(src, dst)
%tryMove  movefile, with its error turned into the same status/message pair it
%   returns for an ordinary failure (a locked file can throw rather than report).
try
    [ok, msg] = movefile(src, dst);
    ok = logical(ok);
catch ME
    ok = false; msg = ME.message;
end
end

% =====================================================================
function tf = rollBack(list, done)
%rollBack  Put the files that already moved back where they came from, newest first.
tf = true;
for k = numel(done):-1:1
    if ~tryMove(list(done(k)).new, list(done(k)).old), tf = false; end
end
end

% =====================================================================
function s = reason(msg)
%reason  A movefile message folded onto one line, in brackets ('' when there is none).
s = strtrim(char(msg));
if isempty(s), return; end
s = [' (' regexprep(s,'\s+',' ') ')'];
end
