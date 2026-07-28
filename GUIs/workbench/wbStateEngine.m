%wbStateEngine - Per-(file,step) state for the Processing Workbench matrix.
%
%   For one recording (a wbFileModel handle) and the step registry, compute the
%   state of every step from three inputs - mirroring removeProcessedFiles but
%   richer:
%     1. Disk markers.  A step is DONE when its settings gating field is present
%        in the recording's SETTINGS files.  The field is read from the UNION of
%        every '<identity>*_s.mat' the recording owns, so detection survives the
%        BFI rename (_K_->_BFI_) and any deleted originals - each step appends
%        its gating field to the settings, and downstream _s files carry it
%        forward.  The REAL gating field per step is used (01 A1).
%     2. Prerequisites.  Every id in step.requires must be done, else unavailable.
%     3. Settings fingerprint.  The current step settings are compared to the
%        stored settings.(gatingField); a mismatch flips a done step to stale.
%
%   States: 'unavailable' | 'ready' | 'done' | 'stale' | 'error'.  ('error' is
%   set by the executor at run time, never here.)
%
% Syntax:
%    st = wbStateEngine(model, reg)
%    st = wbStateEngine(model, reg, curSettings)
%    st = wbStateEngine(model, reg, curSettings, opts)
%
% Inputs:
%    model       - a struct from wbFileModel (identifies the recording).
%    reg         - the struct array from wbStepRegistry.
%    curSettings - (optional) struct whose fields are step ids, each a resolved
%                  settings struct s (from wbSettingsModel).  Used only for the
%                  staleness fingerprint; omit/[] to skip (done stays done).
%    opts        - (optional) struct with:
%                    .sData  containers.Map name->settings struct, injected in
%                            place of a disk scan (for pure unit tests).
%
% Outputs:
%    st - 1xN struct array aligned to reg, with fields:
%           id, applicable (logical), state (char), reason (char).
%
% See also: wbFileModel, wbStepRegistry, wbInvalidate, removeProcessedFiles
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function st = wbStateEngine(model, reg, curSettings, opts)

if nargin<3, curSettings = []; end
if nargin<4, opts = struct(); end

% ---- gather the recording's accumulated settings fields ---------------------
[doneFields, storedMap] = gatherSettings(model, opts);

% id -> gating field, for prerequisite testing
ids   = {reg.id};
gate  = {reg.gatingField};
gateOf = containers.Map(ids, gate);

st = repmat(struct('id','','applicable',false,'state','','reason',''),1,numel(reg));

for k = 1:numel(reg)
    step = reg(k);
    st(k).id = step.id;
    st(k).applicable = any(strcmp(model.modality, step.modalities));

    if ~st(k).applicable
        st(k).state  = 'unavailable';
        st(k).reason = 'modality';
        continue
    end

    doneOnDisk = ~isempty(step.gatingField) && ismember(step.gatingField, doneFields);
    requiresDone = prereqsDone(step.requires, gateOf, doneFields);
    hasInput = inputAvailable(step, model, requiresDone, doneOnDisk);

    if doneOnDisk
        stored = [];
        if isKey(storedMap, step.gatingField), stored = storedMap(step.gatingField); end
        if fingerprintMatches(step, curSettings, stored)
            st(k).state = 'done';
        else
            st(k).state  = 'stale';
            st(k).reason = 'settings changed';
        end
    elseif ~requiresDone
        st(k).state  = 'unavailable';
        st(k).reason = 'prerequisites not met';
    elseif ~hasInput
        st(k).state  = 'unavailable';
        st(k).reason = 'no input file';
    else
        st(k).state = 'ready';
    end
end
end

% =====================================================================
function [doneFields, storedMap] = gatherSettings(model, opts)
%gatherSettings  Union of settings fields across the recording's _s.mat files.
doneFields = {};
storedMap  = containers.Map('KeyType','char','ValueType','any');

if isfield(opts,'sData') && ~isempty(opts.sData)
    % injected: containers.Map name -> settings struct
    names = keys(opts.sData);
    for i = 1:numel(names)
        settings = opts.sData(names{i});
        [doneFields, storedMap] = absorb(settings, doneFields, storedMap);
    end
    return
end

if isempty(model.folder) || ~isfolder(model.folder), return; end
d = dir(fullfile(model.folder,'*_s.mat'));
for i = 1:numel(d)
    cand = wbFileModel(fullfile(d(i).folder, d(i).name));
    if ~strcmp(cand.role,'s'), continue; end
    if ~strcmp(cand.identity, model.identity), continue; end   % same recording only
    try
        S = load(fullfile(d(i).folder, d(i).name),'settings');
    catch
        continue
    end
    if ~isfield(S,'settings'), continue; end
    [doneFields, storedMap] = absorb(S.settings, doneFields, storedMap);
end
end

% =====================================================================
function [doneFields, storedMap] = absorb(settings, doneFields, storedMap)
%absorb  Fold one settings struct's gating fields into the running union.
if ~isstruct(settings), return; end
fn = fieldnames(settings);
for j = 1:numel(fn)
    if ~ismember(fn{j}, doneFields), doneFields{end+1} = fn{j}; end %#ok<AGROW>
    if ~isKey(storedMap, fn{j}), storedMap(fn{j}) = settings.(fn{j}); end
end
end

% =====================================================================
function tf = prereqsDone(requires, gateOf, doneFields)
tf = true;
for i = 1:numel(requires)
    r = requires{i};
    if ~isKey(gateOf, r), tf = false; return; end
    gf = gateOf(r);
    if isempty(gf) || ~ismember(gf, doneFields), tf = false; return; end
end
end

% =====================================================================
function tf = inputAvailable(step, model, requiresDone, doneOnDisk)
%inputAvailable  Whether the step could run now (its input exists).
if doneOnDisk, tf = true; return; end            % already ran => input existed
if isempty(step.requires)
    % entry step (contrast/internalCycle): needs the raw recording on disk
    ext = rawExtFromGlob(step.inGlob);
    if isempty(ext)
        tf = requiresDone;
    else
        tf = isfile(fullfile(model.folder,[model.stem ext]));
    end
else
    tf = requiresDone;                            % the prereq's output IS the input
end
% guided-style steps additionally need the raw recording alongside the product
if tf && step.needsRaw
    tf = isfile(fullfile(model.folder,[model.stem '.rls'])) || ...
         isfile(fullfile(model.folder,[model.stem '.cxd']));
end
end

% =====================================================================
function ext = rawExtFromGlob(glob)
%rawExtFromGlob  '.rls'/'.cxd' from a raw inGlob like '*.rls'; '' for .mat globs.
ext = '';
if contains(glob,'.mat'), return; end
[~,~,e] = fileparts(glob);
ext = e;
end

% =====================================================================
function tf = fingerprintMatches(step, curSettings, stored)
%fingerprintMatches  Current step settings vs the stored run-time settings.
%   Compares only the step's own tunable fields present in BOTH structs; missing
%   comparisons (no current snapshot / older run) count as a match.
tf = true;
if isempty(curSettings) || isempty(stored) || ~isstruct(stored), return; end
if ~isfield(curSettings, step.id), return; end
cur = curSettings.(step.id);
if ~isstruct(cur), return; end

fields = stepFields(step);
ignore = {'libraryFolder','fNamesCopyTo','optimizer','metric','categories', ...
          'refFName','fName','progressFcn','stageFcn','cancelFcn'};
for i = 1:numel(fields)
    f = fields{i};
    if ismember(f, ignore), continue; end
    if isfield(cur,f) && isfield(stored,f)
        if ~isequaln(cur.(f), stored.(f)), tf = false; return; end
    end
end
end

% =====================================================================
function fields = stepFields(step)
%stepFields  The tunable field names a step reads (from its settingGroups + shared keys).
fields = {};
if ~isempty(step.settingGroups)
    fields = [step.settingGroups{:,2}];
end
fields = unique([fields, step.sharedKeys]);
end
