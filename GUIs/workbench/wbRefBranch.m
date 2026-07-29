%wbRefBranch - Which FILE of an animal's reference recording a per-animal step uses.
%
%   The workbench pins ONE reference RECORDING per animal (an identity, valid for
%   every file of that animal whatever its type or experimental group - spec D5).
%   A recording owns several co-registered branch products ('_t_K', '_c_K',
%   '_t_BFI', '_c_BFI', ...), so each per-animal step still has to say WHICH of
%   them it wants: registration templates on the CONTRAST side, vessel typing
%   paints on the CARDIAC side.  That declaration is the registry's refBranch
%   field; this function resolves it against what is actually on disk.
%
%   The rule mirrors wbExecutor>desiredStage for the step INPUTS - prefer the
%   declared branch, and FALL BACK to any other branch product of that same
%   recording rather than failing, since a half-processed project legitimately has
%   only one branch yet.  The caller decides what to do with a fall-back; the
%   Constructor tab warns, it never blocks.
%
%   Candidates are located by base name + the step's own inGlob and filtered to the
%   reference recording's IDENTITY, so a sibling recording sharing a name prefix
%   can never be picked up.  Within the wanted branch the contrast STAGE the
%   project actually uses (t or s) is preferred, which is why cstage is an input.
%
%   Pure lookup: no state, no settings, no side effects - only dir().
%
% Syntax:
%    out = wbRefBranch(step, refModel)
%    out = wbRefBranch(step, refModel, cstage)
%
% Inputs:
%    step     - one step struct from wbStepRegistry (reads .refBranch and .inGlob).
%    refModel - the reference recording: a wbFileModel struct, or its IDENTITY as
%               a char (what guiWorkbench stores per animal), or '' for none.
%    cstage   - 't' | 's', the contrast flag this project's files carry (default
%               't').  Only orders the candidates INSIDE the contrast branch.
%
% Outputs:
%    out - struct with fields:
%       path       full path of the file the step will use ('' when none exists).
%       name       its bare file name ('' when none).
%       wanted     the branch asked for ('contrast'|'cardiac'|'').
%       branch     the branch actually resolved ('' when none).
%       status     'ok'       - a file of the wanted branch was found;
%                  'fallback' - only another branch exists, and it was taken;
%                  'noref'    - no reference recording was given;
%                  'missing'  - the reference owns no file this step can read.
%       note       one short human-readable line for the log / warning banner.
%       candidates cellstr of every file considered, wanted branch first.
%
% See also: wbStepRegistry, wbExecutor, guiWorkbench, wbFileModel
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbRefBranch(step, refModel, cstage)

if nargin<3 || isempty(cstage), cstage = 't'; end
cstage = char(cstage);

wanted = '';
if isstruct(step) && isfield(step,'refBranch'), wanted = char(step.refBranch); end
out = struct('path','','name','','wanted',wanted,'branch','', ...
             'status','noref','note','no reference recording pinned for this animal', ...
             'candidates',{{}});

model = modelOf(refModel);
if isempty(model), return; end

cand = candidatesFor(model, step);
if isempty(cand)
    out.status = 'missing';
    out.note   = sprintf('the reference %s has no %s file on disk yet', ...
        shortName(model.identity), globTail(step));
    return
end

% ---- order: the wanted branch first (contrast stage preferred inside it) ------
[hit, rest] = splitByBranch(cand, wanted, cstage);
out.candidates = [hit, rest];
if ~isempty(hit)
    out.path = hit{1};  out.status = 'ok';  out.note = '';
else
    out.path = rest{1}; out.status = 'fallback';
    out.note = sprintf('the reference %s has no %s branch - falling back to %s', ...
        shortName(model.identity), branchWord(wanted), shortName(rest{1}));
end
out.name   = shortName(out.path);
bm         = wbFileModel(out.path);
out.branch = bm.branch;
if isempty(wanted), out.status = 'ok'; out.note = ''; end   % no preference = never a fallback
end

% =====================================================================
function m = modelOf(refModel)
%modelOf  Accept a wbFileModel, a recording identity, or '' (no reference).
m = [];
if isempty(refModel), return; end
if isstruct(refModel), m = refModel; return; end
m = wbFileModel(char(refModel));    % an identity parses to folder + stem + identity
end

% =====================================================================
function c = candidatesFor(model, step)
%candidatesFor  Every file of THIS recording matching the step's input glob.
c = {};
if isempty(model.folder) || ~isfolder(model.folder), return; end
tail = globTail(step);
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.folder,[base '*' tail]));
for i = 1:numel(d)
    fp = fullfile(d(i).folder, d(i).name);
    cm = wbFileModel(fp);
    if ~strcmp(cm.identity, model.identity), continue; end   % a name-prefix sibling
    c{end+1} = fp; %#ok<AGROW>
end
end

% =====================================================================
function t = globTail(step)
%globTail  '_K_d.mat' / '_BFI_d.mat' - the step's input glob without its wildcard.
t = '';
if isstruct(step) && isfield(step,'inGlob'), t = regexprep(char(step.inGlob),'^\*',''); end
end

% =====================================================================
function [hit, rest] = splitByBranch(cand, wanted, cstage)
%splitByBranch  Candidates of the wanted branch (preferred stage first), then the rest.
rest = {};
if isempty(wanted)
    hit = cand; return                        % no preference: every candidate qualifies
end
pref = {}; other = {};
for i = 1:numel(cand)
    m = wbFileModel(cand{i});
    if ~strcmp(m.branch, wanted)
        rest{end+1} = cand{i}; %#ok<AGROW>
    elseif strcmp(m.stage, cstage)
        pref{end+1} = cand{i}; %#ok<AGROW>   % e.g. _t when the project is temporal
    else
        other{end+1} = cand{i}; %#ok<AGROW>  % same branch, the other flag (_s)
    end
end
hit = [pref, other];
end

% =====================================================================
function w = branchWord(b)
%branchWord  How a branch is named in a message.
switch b
    case 'contrast', w = 'contrast (_t/_s)';
    case 'cardiac',  w = 'cardiac (_c)';
    case '',         w = 'any';
    otherwise,       w = b;
end
end

% =====================================================================
function n = shortName(p)
[~,nm,ex] = fileparts(char(p));
n = [nm ex];
end
