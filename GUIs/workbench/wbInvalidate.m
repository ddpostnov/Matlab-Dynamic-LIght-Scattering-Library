%wbInvalidate - Forward-cascade set of steps to reset after a settings edit.
%
%   When the user edits a setting belonging to step X (or a shared key that some
%   steps read), step X and every step that transitively DEPENDS on X (its
%   requires-descendants) must be reset forward - upstream steps stay done.
%   The same cascade fires when a step is actually re-run (its outputs changed).
%   This is a pure function of the registry's requires graph.
%
% Syntax:
%    stepIds = wbInvalidate(reg, edited)
%    cells   = wbInvalidate(reg, edited, files)
%
% Inputs:
%    reg    - the struct array from wbStepRegistry.
%    edited - either a step id ('segmentation') or a shared-key name
%             ('trustLimitsK').  A step id seeds the cascade at that step; a
%             shared key seeds it at EVERY step that reads that key.
%    files  - (optional) cell array of file keys (or a wbFileModel array); when
%             given, the result is the (file,step) grid to reset.
%
% Outputs:
%    stepIds - forward step ids in registry order (self + descendants).
%    cells   - Mx2 cell {fileKey, stepId} when files is supplied.
%
% See also: wbStepRegistry, wbStateEngine
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbInvalidate(reg, edited, files)

ids = {reg.id};

% ---- seed steps: the edited step, or every step reading the shared key ------
edited = char(edited);
if any(strcmp(edited, ids))
    seed = {edited};
else
    seed = ids(cellfun(@(st) ismember(edited, stepFields(st)), num2cell(reg)));
end

% ---- transitive closure over "who requires whom" ----------------------------
hit = false(1,numel(reg));
for i = 1:numel(seed)
    hit = hit | reachable(reg, ids, seed{i});
end
stepIds = ids(hit);                              % already in registry order

if nargin < 3 || isempty(files)
    out = stepIds;
    return
end

% ---- expand to (file, step) cells -------------------------------------------
fileKeys = fileKeysOf(files);
out = cell(numel(fileKeys)*numel(stepIds), 2);
row = 0;
for fi = 1:numel(fileKeys)
    for si = 1:numel(stepIds)
        row = row + 1;
        out{row,1} = fileKeys{fi};
        out{row,2} = stepIds{si};
    end
end
end

% =====================================================================
function mask = reachable(reg, ids, startId)
%reachable  Logical mask over reg of startId plus all its requires-descendants.
mask = false(1,numel(reg));
mask(strcmp(ids,startId)) = true;
changed = true;
while changed
    changed = false;
    for k = 1:numel(reg)
        if mask(k), continue; end
        % staleness is deliberately BROAD: any declared producer counts, including
        % an alternative the step did not happen to consume (wbPrereqs 'all')
        if any(ismember(wbPrereqs('all', reg(k)), ids(mask)))
            mask(k) = true;
            changed = true;                      % a new hit may pull in more
        end
    end
end
end

% =====================================================================
function fields = stepFields(step)
fields = {};
if ~isempty(step.settingGroups)
    fields = [step.settingGroups{:,2}];
end
fields = unique([fields, step.sharedKeys]);
end

% =====================================================================
function keys = fileKeysOf(files)
if isstruct(files)                                % wbFileModel array
    keys = arrayfun(@(m) m.identity, files, 'UniformOutput', false);
elseif ischar(files)
    keys = {files};
else
    keys = files;
end
keys = keys(:).';
end
