%wbArtifacts - Resolve a step's on-disk report images for one recording.
%
%   The pipeline wrappers emit a report page next to the data triplet, named
%   '<stem>_rep_<stage>.jpg' by Core/Reporting (_rep_contrast from contrast,
%   _rep_categories/_rep_segments from segmentation, ...).  A step's expected
%   report tails are declared in its wbStepRegistry spec (step.artifacts).  This
%   helper resolves those tails against a recording's base name and returns the
%   ones that actually exist on disk, so the workbench can surface a finished
%   cell's reports as clickable links.  It never creates a new artifact format -
%   it only locates the images the wrappers already write (01-pipeline-map.md §7).
%
%   BOTH GENERATIONS ARE LOOKED FOR.  Before unified reporting the same pages had
%   cryptic tails ('_c.jpg', '_cm.jpg', '_vs.jpg', '_registration.png'); those live
%   on in step.legacyArtifacts and are asked for in the same pass, so a folder
%   processed last month still lists its reports.  Both lists are matched over ONE
%   directory listing, so the second generation costs an endsWith pass and no I/O.
%
%   Matching is by base name + TAIL (endsWith), so the stage/product flags that
%   sit between them ('_t_K' in 'Foo_t_K_rep_segments.jpg') do not have to be
%   spelled out: a '_rep_segments.jpg' tail finds 'Foo_t_K_rep_segments.jpg' and a
%   'Roi2_' crop's variant alike.
%
%   THE WORKING SET IS A FENCE (round-2 item 8).  That base-name listing is the
%   whole folder, so a page left by an EARLIER session - the '_c_K' pipeline of a
%   project that now configures only contrast - matched a tail like any other and
%   joined the result list and the per-column PDF.  The optional third argument is
%   the session's admissible stage flags (wbProducts 'flags'), and a page whose
%   stage is not one of them is dropped.  ABSENT ARGUMENT = NO FENCE, so the
%   headless callers and the tests keep the behaviour they were written against,
%   and a page whose stage cannot be read at all is kept either way (wbProducts
%   fails open - hiding a real report is worse than listing a stale one).
%
% Syntax:
%    files = wbArtifacts(model, step)         % existing artifacts for one step
%    files = wbArtifacts(model, reg)          % union over a whole registry array
%    files = wbArtifacts(model, step, flags)  % ... limited to the session's branches
%
% Inputs:
%    model - a wbFileModel struct (supplies folder + roiPrefix + stem = the
%            recording identity the images are named after).
%    step  - one wbStepRegistry element (uses step.artifacts and
%            step.legacyArtifacts), OR the whole registry array (every element's
%            tails are unioned).
%    flags - (optional) the admissible stage flags of this recording's session,
%            from wbProducts('flags',...).  {} (or omitted) = no fence.
%
% Outputs:
%    files - 1xK cellstr of full paths to the artifact images that exist on disk
%            (empty 1x0 when the step declares none or none are present yet).
%
% See also: wbStepRegistry, wbFileModel, wbProducts, wbExecutor, makeReportPdf,
%           guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

%------------- BEGIN CODE --------------
function files = wbArtifacts(model, step, flags)

if nargin < 3, flags = {}; end

files = cell(1,0);
if isempty(model) || ~isstruct(model), return; end
if isempty(model.folder) || ~isfolder(model.folder), return; end

% collect the requested tails (one step, or every step in a registry array).  The
% tails a run WRITES and the ones it used to write are asked for together, so a
% folder processed before the reporting rename still lists its reports.
tails = cell(1,0);
for k = 1:numel(step)
    tails = [tails, tailsOf(step(k),'artifacts'), tailsOf(step(k),'legacyArtifacts')]; %#ok<AGROW>
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
        if ~endsWith(nm, tail) || isKey(seen, nm), continue; end
        seen(nm) = true;                              % judged once, whatever the verdict
        if ~isempty(flags) && ~wbProducts('admits', flags, nm), continue; end
        files{end+1} = fullfile(model.folder, nm); %#ok<AGROW>
    end
end
end

% =====================================================================
function t = tailsOf(step, field)
%tailsOf  One tail list of a step spec as a row cellstr, tolerating a spec that
%   predates the field (an old registry struct read back from somewhere else).
t = cell(1,0);
if ~isfield(step,field), return; end
a = step.(field);
if isempty(a), return; end
if ischar(a), a = {a}; end
t = a(:)';
end
