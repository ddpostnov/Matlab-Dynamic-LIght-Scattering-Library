%wbArtifacts - Resolve a step's on-disk report images for one recording.
%
%   The pipeline wrappers emit preview JPG/PNGs next to the data triplet
%   (_c.jpg from contrast, _cm.jpg/_vs.jpg from segmentation, _registration.png
%   from registration, ...).  A step's expected report globs are declared in its
%   wbStepRegistry spec (step.artifacts, a set of file-name TAILS).  This helper
%   resolves those tails against a recording's base name and returns the ones
%   that actually exist on disk, so the workbench can surface a finished cell's
%   reports as clickable links.  It never creates a new artifact format - it only
%   locates the images the wrappers already write (01-pipeline-map.md §7).
%
%   Matching is by base name + TAIL (endsWith), so the stage/product flags that
%   sit between them ('_t_K' in 'Foo_t_K_cm.jpg') do not have to be spelled out:
%   a '_cm.jpg' tail finds 'Foo_t_K_cm.jpg' and a 'Roi2_' crop's variant alike.
%
% Syntax:
%    files = wbArtifacts(model, step)      % existing artifacts for one step
%    files = wbArtifacts(model, reg)       % union over a whole registry array
%
% Inputs:
%    model - a wbFileModel struct (supplies folder + roiPrefix + stem = the
%            recording identity the images are named after).
%    step  - one wbStepRegistry element (uses step.artifacts), OR the whole
%            registry array (every element's artifacts are unioned).
%
% Outputs:
%    files - 1xK cellstr of full paths to the artifact images that exist on disk
%            (empty 1x0 when the step declares none or none are present yet).
%
% See also: wbStepRegistry, wbFileModel, wbExecutor, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function files = wbArtifacts(model, step)

files = cell(1,0);
if isempty(model) || ~isstruct(model), return; end
if isempty(model.folder) || ~isfolder(model.folder), return; end

% collect the requested tails (one step, or every step in a registry array)
tails = cell(1,0);
for k = 1:numel(step)
    a = step(k).artifacts;
    if ~isempty(a), tails = [tails, a(:)']; end %#ok<AGROW>
end
if isempty(tails), return; end
tails = unique(tails,'stable');

% one directory listing of everything named after this recording's base
base = [model.roiPrefix model.stem];
d = dir(fullfile(model.folder,[base '*']));
if isempty(d), return; end
names = {d.name};

seen = containers.Map('KeyType','char','ValueType','logical');
for t = 1:numel(tails)
    tail = tails{t};
    for i = 1:numel(names)
        nm = names{i};
        if endsWith(nm, tail) && ~isKey(seen, nm)
            seen(nm) = true;
            files{end+1} = fullfile(model.folder, nm); %#ok<AGROW>
        end
    end
end
end
