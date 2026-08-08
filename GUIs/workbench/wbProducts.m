%wbProducts - Which files of a recording belong to THIS session, and to what.
%
%   A recording's folder accumulates.  A project processed last month leaves a
%   '_c_K' triplet and its report pages beside the '_t_K' one this session
%   configured, a branch that was tried once and abandoned leaves a whole
%   pipeline behind, and nothing on disk says which of them the session in front
%   of you is about.  The workbench TABLES already know - they are built from
%   wbTypeSelection('rows',...), i.e. from the configuration - but two disk globs
%   were not: the input glob that fills a wrapper's fNames (wbBatchPlan>
%   resolveStepInputs) and the artifact glob that fills the result list and the
%   per-column PDFs (wbArtifacts).  This module is the one place that answers the
%   question both of them have to ask, so "listed in the table", "allowed into the
%   run" and "listed in the reports" are ONE fact rather than three.
%
%   It is pure: it reads the registry, the type selection and the directory, and
%   it changes nothing.  It owns no state and touches no figure.  DELETING is the
%   caller's act - 'below' only NAMES the files a re-run supersedes.
%
%   THE DELETION SET COMES FROM THE REGISTRY, NOT FROM THE FILE NAME.  This is the
%   non-obvious decision here and the next reader will otherwise "simplify" it
%   back to a name-prefix test, which silently deletes nothing.  A derived product
%   SUBSTITUTES the product token rather than extending the name: runBFI.m:103 is
%   strrep(fName,'_K_r.mat','_BFI_r.mat'), so 'Mouse1_t_K_r.mat' becomes
%   'Mouse1_t_BFI_r.mat' and there is no '_K_BFI' for a prefix test to find.  The
%   flag chain cannot separate ancestor from descendant either - '_t_K' and '_t_BFI'
%   both parse to flags {t}.  So WHICH STEPS are downstream is asked of the
%   registry's requires graph (wbInvalidate), WHAT each of them writes is its own
%   outSuffix parsed into a (stage, product) SIGNATURE, and only WHICH FILES ARE
%   THOSE is answered by the names on disk.  The registry orders; the file name
%   identifies.
%
%   A SIGNATURE, NEVER A COMPOSED NAME.  outTransform is a strrep with a hardcoded
%   from-token, so it is not a safe way to compose an output name: a rule written
%   '_t_K_r.mat' -> '_x_K_r.mat' is a NO-OP on an '_s' or '_c' input and hands back
%   THE INPUT FILE ITSELF, which a deletion set would then propose for deletion.  The
%   step that made that concrete - the external cycle, whose declared '_e_K_r.mat'
%   and whose actually-written 'Mouse1_t_e_K_r.mat' never agreed - has been retired
%   with the whole epoch-averaged product (2026-08-05), and every outTransform left
%   in the registry rewrites an extension or the product token alone.  The rule stays
%   because it is what makes this module survive the NEXT step whose composed name
%   and written name disagree: matching a PARSED signature asks what a file IS rather
%   than what a rule would have called it.
%
%   WHEN ANCESTRY IS NOT RECOVERABLE FROM THE NAME, UNDER-DELETE (spec D9a).
%   Over-deleting destroys the author's data; under-deleting leaves one narrow
%   case of the bug.  A file matching a downstream signature whose flag chain
%   cannot be tied to the stage being rewritten - a 'Mouse1_c_K_r.mat' when '_t_K'
%   is the one being recomputed - is returned SEPARATELY as unattributable, for
%   the caller to name in the log, and is never deleted on a guess.
%
%   A REPORT PAGE IS ATTRIBUTED THROUGH THE FILE IT DESCRIBES.  wbFileModel parses
%   '.mat' products only: handed a '.jpg' it returns the whole name as the identity
%   and leaves stage / branch / product / role empty.  A page is named
%   '<stem>_rep_<tag>.jpg' after the product it belongs to (reportSave.m:86), so
%   'Mouse1_t_K_r_rep_segments.jpg' is stripped at its last '_rep_' and read as
%   'Mouse1_t_K_r.mat' -> stage 't', and 'Mouse1_rep_contrast.jpg' as 'Mouse1' ->
%   no flag, i.e. the recording itself.  A PRE-RENAME page ('Mouse1_t_K_cm.jpg')
%   carries no '_rep_' marker and cannot be attributed at all: 'stageOf' returns
%   the UNKNOWN token for it and 'admits' FAILS OPEN.  That is a policy choice -
%   hiding a real report the user has is worse than listing a stale one.  To
%   reverse it, make admits return false for the unknown token (one line, in
%   admitsPath below).
%
% Syntax:
%    flags        = wbProducts('flags', reg, typeSel, type, contrastStage)
%    stage        = wbProducts('stageOf', path)
%    tf           = wbProducts('admits', flags, path)
%    stage        = wbProducts('writes', step, contrastStage)
%    [files,unk]  = wbProducts('below', reg, model, stepId, stage)
%    tok          = wbProducts('unknown')
%
% Inputs:
%    reg           - the struct array from wbStepRegistry.
%    typeSel       - the wbTypeSelection map (the per-type configuration).
%    type          - one recording type token ('' = nothing configured).
%    contrastStage - 't' | 's', the flag the settings-driven producer writes for
%                    this type (the host owns that answer - guiWorkbench>
%                    contrastStageForType / contrastStageForModel).
%    path          - one file path or bare name (a '.mat' product or a report page).
%    flags         - the admissible stage-flag list from 'flags'.  EMPTY = no
%                    answer, which is NO FENCE, never "nothing is allowed".
%    step          - one wbStepRegistry element.
%    model         - the wbFileModel of the recording (supplies folder + identity).
%    stepId        - the step being re-run.
%    stage         - the stage flag that re-run is about to rewrite ('t', 'c', ...).
%
% Outputs:
%    flags - 1xK cellstr of admissible stage flags, always including '' (the raw
%            recording carries no flag).
%    stage - a single stage flag, '' for a flagless name, or the UNKNOWN token.
%    tf    - whether this file belongs to the session those flags describe.
%    files - 1xK cellstr of full paths a re-run supersedes: the downstream
%            products AND their own report pages.  Empty when there is nothing to
%            remove, which is the normal case.
%    unk   - 1xM cellstr of paths that match a downstream signature but cannot be
%            tied to the stage being rewritten.  NOT deleted (spec D9a).
%    tok   - the token 'stageOf' returns when a name cannot be attributed.
%
% See also: wbBatchPlan, wbArtifacts, wbExecutor, wbInvalidate, wbFileModel,
%           wbTypeSelection, wbStepRegistry, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function varargout = wbProducts(action, varargin)

switch lower(char(action))
    case 'flags',   varargout{1} = admissibleFlags(varargin{:});
    case 'stageof', varargout{1} = stageOfPath(varargin{:});
    case 'admits',  varargout{1} = admitsPath(varargin{:});
    case 'writes',  varargout{1} = stageWritten(varargin{:});
    case 'below',   [varargout{1:max(1,nargout)}] = filesBelow(varargin{:});
    case 'unknown', varargout{1} = unknownStage();
    otherwise
        error('wbProducts:action','Unknown action ''%s''.', char(action));
end
end

% =====================================================================
function tok = unknownStage()
%unknownStage  The answer for a name whose stage cannot be read.  It is not a
%   stage flag (wbFileModel's set is t|s|c|b), so it can never collide with a
%   real one, and '' is already taken by the flagless recording itself.
tok = '?';
end

% =====================================================================
function flags = admissibleFlags(reg, sel, type, cstage)
%admissibleFlags  The stage flags THIS RECORDING'S CONFIGURATION produces (spec
%   D7) - the same question the derived monitor table asks, so a file that is not
%   in the table cannot be in the run either.
%
%   Every branch row the type has contributes the flag its own steps stamp: the
%   row's raw producer ('_t'/'_s' for contrast, '_c' for the internal cycle) and
%   any derived step ticked on it that writes a new stage of its own.  A step that
%   appends in place, or that keeps its input's stage (BFI: '_t_K' -> '_t_BFI'),
%   adds nothing - its products already carry a flag this list holds.  NOTHING HERE
%   NAMES A FLAG: every one of them is read off the registry's outSuffix, so a new
%   entry step for a new modality is a registry edit and this function does not
%   change.  (No DERIVED step stamps a stage of its own today - the external cycle
%   was the only one and it is retired - so in practice the list is the raw
%   producers'.  The loop is still over every ticked step, because which kind a step
%   is, is the registry's answer and not this file's.)
if nargin<4, cstage = ''; end
flags = {''};                                   % the raw recording carries no flag
if nargin<3 || isempty(type) || isempty(sel) || isempty(reg), return; end

brs = wbTypeSelection('rows', sel, reg, char(type));
for b = 1:numel(brs)
    ids = [wbTypeSelection('steps',     sel, reg, char(type), brs{b}), ...
           wbTypeSelection('inherited', sel, reg, char(type), brs{b})];
    for i = 1:numel(ids)
        st = stageWritten(stepById(reg, ids{i}), cstage);
        if ~isempty(st) && ~any(strcmp(st, flags)), flags{end+1} = st; end %#ok<AGROW>
    end
end
end

% =====================================================================
function st = stageWritten(step, cstage)
%stageWritten  The stage flag a step STAMPS on the triplet it writes, '' when it
%   writes none of its own (an in-place step, or one like BFI that keeps whatever
%   stage its input carried).  Read from the step's own outSuffix through
%   wbFileModel - the idiom guiWorkbench>rowFlagFor:997 already uses.
%
%   ONE PRODUCER'S STAGE IS A SETTING, not a constant: the contrast step writes
%   '_t' or '_s' per s.contrastType, and the host resolves that per type.  Which
%   producer that is, is DERIVED rather than named: the given contrast flag is
%   re-parsed against this step's own product/role, and it wins only when it lands
%   in the same analysis BRANCH the step's declared suffix does.  A cardiac producer
%   is therefore untouched by it.
st = '';
if nargin<2, cstage = ''; end
if isempty(step) || ~isstruct(step) || ~isfield(step,'outSuffix'), return; end
if isempty(step.outSuffix), return; end

m  = wbFileModel(['x' char(step.outSuffix{1}) '.mat']);      % '_t_K_d' -> stage t
st = m.stage;
if isempty(st) || isempty(cstage) || isempty(m.product), return; end
alt = wbFileModel(['x_' char(cstage) '_' m.product '_' m.role '.mat']);
if ~isempty(alt.branch) && strcmp(alt.branch, m.branch), st = char(cstage); end
end

% =====================================================================
function st = stageOfPath(pth)
%stageOfPath  WHICH STAGE DOES THIS FILE BELONG TO - for a '.mat' product and for
%   a report page alike, because wbFileModel answers only the first (see header).
p = char(pth);
[~, nm, ext] = fileparts(p);

k = strfind(nm, '_rep_');
if ~isempty(k)
    st = stageOfName([nm(1:k(end)-1) '.mat']);   % the product the page describes
    return
end
if strcmpi(ext, '.mat')
    st = stageOfName(p);
    return
end
st = unknownStage();          % a raw recording or a pre-rename page: fails open
end

function st = stageOfName(nm)
m = wbFileModel(nm); st = m.stage;
end

% =====================================================================
function tf = admitsPath(flags, pth)
%admitsPath  Is this file part of the session those flags describe?  THE rule that
%   makes "listed in the table" and "allowed into the run" one fact.  An empty
%   flag list is NO ANSWER and therefore no fence - a caller that cannot say what
%   its session contains may never be the reason a file disappears from it.
if isempty(flags), tf = true; return; end
if ischar(flags), flags = {flags}; end
st = stageOfPath(pth);
tf = strcmp(st, unknownStage()) || any(strcmp(st, flags));
end

% =====================================================================
function [files, unattributable] = filesBelow(reg, model, stepId, stage)
%filesBelow  WHAT A RE-RUN OF THIS STEP MUST REMOVE (spec D9): every product of
%   this recording that was DERIVED from the one about to be rewritten, plus the
%   report pages those products own.
%
%   Three moves, and each is here because the other two cannot do its job:
%     1. which STEPS are downstream        - wbInvalidate's requires graph;
%     2. what each of them WRITES          - its outSuffix, parsed into a
%                                            (stage, product) signature;
%     3. which FILES on disk are those     - this recording's own names, matched
%                                            against those signatures and scoped
%                                            to the stage being rewritten.
%   The step's OWN triplet is deliberately not in the list: save() replaces a file
%   whole (no -append anywhere in the library, and every producer clears its
%   variables per file), so the thing being recomputed overwrites itself.
%
%   Nothing outside the recording's own results folder is ever named.
files = cell(1,0); unattributable = cell(1,0);
if isempty(model) || ~isstruct(model), return; end
if isempty(model.resultsFolder) || ~isfolder(model.resultsFolder), return; end
stage = char(stage);
if isempty(stage), return; end     % no stage to scope by: nothing is attributable

% ---- 1+2. the signatures of every step downstream of this one -----------------
sigs = downstreamSignatures(reg, stepId);
if isempty(sigs), return; end

% ---- 3. this recording's own files, matched against them ----------------------
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.resultsFolder,[base '*']));
victims = cell(1,0);
for i = 1:numel(d)
    if d(i).isdir, continue; end
    fp = fullfile(d(i).folder, d(i).name);
    m  = wbFileModel(fp);
    % the BASE, not the identity - the folder half of an identity comparison is a
    % tautology within one listing, and a false negative once that listing is the
    % results folder and the model came from the raw tree
    if ~strcmp([m.roiPrefix m.stem], base), continue; end
    if isempty(m.product), continue; end        % the raw recording, or a page
    [hit, sure] = matchesSignature(sigs, m, stage);
    if hit && sure
        victims{end+1} = fp; %#ok<AGROW>
    elseif hit
        unattributable{end+1} = fp; %#ok<AGROW>
    end
end
if isempty(victims), return; end

% ---- 4. a doomed product takes its own report pages with it -------------------
files = [victims, pagesOf(model.resultsFolder, victims)];
end

% =====================================================================
function sigs = downstreamSignatures(reg, stepId)
%downstreamSignatures  The (stage, product) pairs written by every step BELOW this
%   one.  wbInvalidate's step-id form IS that graph (self + requires-descendants);
%   the self is dropped, and a step that appends in place declares no outSuffix and
%   so contributes nothing - which is why re-running contrast removes a BFI and an
%   external cycle but never a segmentation.
sigs = struct('stage',{},'product',{});
ids = wbInvalidate(reg, stepId);
ids = ids(~strcmp(ids, char(stepId)));
for i = 1:numel(ids)
    st = stepById(reg, ids{i});
    if isempty(st) || isempty(st.outSuffix), continue; end
    for k = 1:numel(st.outSuffix)
        m = wbFileModel(['x' char(st.outSuffix{k}) '.mat']);
        if isempty(m.product), continue; end
        if any(strcmp(m.stage,{sigs.stage}) & strcmp(m.product,{sigs.product})), continue; end
        sigs(end+1) = struct('stage',m.stage,'product',m.product); %#ok<AGROW>
    end
end
end

% =====================================================================
function [hit, sure] = matchesSignature(sigs, m, stage)
%matchesSignature  Is this file one of the downstream products, and can it be tied
%   to the stage being rewritten?  Two shapes of signature, and they read the file
%   differently:
%
%     the step stamps its OWN stage - the file's stage then says which STEP wrote it,
%       so the base it came from has to be looked for in the rest of the flag chain,
%       and a name carrying only the step's own flag is UNATTRIBUTABLE.  No step in
%       the registry is of this kind today: the external cycle was, and its stacked
%       'Mouse1_t_e_K_d' is what the two-token flag chain existed to read.
%     the step KEEPS its input's stage (BFI '_BFI') - the file's own stage names
%       the pipeline it belongs to, so '_t_BFI' is ours and '_c_BFI' is simply the
%       cardiac side's, not an ambiguity worth reporting.
hit = false; sure = false;
for i = 1:numel(sigs)
    if ~strcmp(sigs(i).product, m.product), continue; end
    if isempty(sigs(i).stage)
        if strcmp(m.stage, stage), hit = true; sure = true; return; end
    else
        if ~strcmp(sigs(i).stage, m.stage), continue; end
        hit = true;
        sure = any(strcmp(stage, m.flags));
        if sure, return; end
    end
end
end

% =====================================================================
function pgs = pagesOf(folder, victims)
%pagesOf  The report pages a doomed product owns.  They are '<thatStem>_rep_*' by
%   construction (reportSave.m:86), so no tail list and no registry lookup is
%   needed here - which is also why a page written by a step this module has never
%   heard of still goes with its data.
pgs = cell(1,0);
for i = 1:numel(victims)
    [~, stem] = fileparts(victims{i});
    d = dir(fullfile(folder,[stem '_rep_*']));
    for k = 1:numel(d)
        if d(k).isdir, continue; end
        pgs{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
    end
end
pgs = unique(pgs,'stable');
end

% =====================================================================
function s = stepById(reg, id)
s = [];
if isempty(reg), return; end
s = reg(strcmp({reg.id}, char(id)));
if ~isempty(s), s = s(1); end
end
