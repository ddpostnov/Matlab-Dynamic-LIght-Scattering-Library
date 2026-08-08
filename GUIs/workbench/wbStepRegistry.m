%wbStepRegistry - Declarative spec of the Processing-Workbench steps.
%
%   Returns the ordered array of pipeline-step specifications that drives every
%   part of the workbench: the matrix columns, the on-disk gating, the settings
%   panels, execution, and the report links.  Adding a modality or a step later
%   is a data edit here, not a code rewrite.  The specs are transcribed FROM THE
%   CODE (the Wrappers/ functions and the Launchers/ %Example-of-s blocks) as of
%   branch main - see claude-docs/processing-workbench/01-pipeline-map.md.
%
%   STEPS ARE GROUPED BY MODALITY FAMILY, and listed in dependency order WITHIN each
%   family - which is why an entry step (myoVideo, intensity) can sit below a consumer
%   of another family:
%     contrast · internalCycle · setRegions · segmentation · dynamicSegmentation ·
%     guided · registration · BFI · vasomotion · pulsatility · nvc ·
%     fitVasoreactivity · vesselTypes · vascularTree ·
%     myoVideo · labChart · myoPresetIntervals · myoCrop · myoDiameter ·
%     myoIntervals · myoPropagation · myoVasomotion ·
%     intensity · intensityCycle · setRegionsI · segmentationI ·
%     dynamicSegmentationI · backgroundRemoval · intensityBolus · ctth · topology ·
%     motionEnhancement · vasomotionI · pulsatilityI · vesselTypesI · vascularTreeI
%
% Syntax:
%    reg = wbStepRegistry()
%    reg = wbStepRegistry(modality)      % filter to steps exposing that modality
%    reg = wbStepRegistry(modalities)    % cellstr: the UNION over several, in
%                                        %   registry order (a mixed working set)
%    reg = wbStepRegistry('filter', reg, modality)   % the same filter + prune,
%                                        %   applied to a step array you hold
%    wbStepRegistry('validate', reg)     % throw on a malformed step array
%    twin = wbStepRegistry('intensityTwin', step)    % the EPFL twin of one step
%                                        %   (see INTENSITY TWINS below)
%
% Outputs:
%    reg - 1xN struct array; each element has fields:
%       id           stable key                      ('contrast')
%       label        column header                   ('Contrast')
%       wrapper      seam function handle            (@runContrastFromRLS)
%       arity        'perFile' | 'perAnimal'          (perAnimal = 2-D fNames row)
%       inGlob       dir glob it consumes            ('*.rls')
%       outSuffix    new triplet suffixes            ({'_t_K_d','_t_K_r','_t_K_s'})
%       outKind      'new'|'inplace'|'prefix'|'none' (how outputs relate to input)
%       outTransform strrep rule on the _r name      (struct from/to, or [])
%       gatingField  settings.<field> written        ('runContrastFromRLS')
%       requires     upstream step ids, ALL needed  ({} | {'setRegions'} ...)
%       requiresAny  upstream step ids, ANY one is enough - the branch-agnostic
%                    middle of the pipeline (see below); wbPrereqs owns the rule
%       conflictsWith  step ids that CANNOT be selected together with this one
%                    ({} for every v1 step); declared symmetrically on BOTH sides,
%                    enforced in wbTypeSelection>setOne (see below)
%       produces     logical product tokens          ({'contrast_t'})
%       interactive  false | @(s)tf                  (blocks for user input?)
%       needsRaw     true for runGuided*             (also passes the raw file)
%       artifacts    report-image tails a run WRITES ({'_rep_contrast.jpg'})
%       legacyArtifacts  the tails the same step wrote BEFORE the unified reporting
%                    rename, kept so images from earlier runs are still listed
%                    ({'_c.jpg'}); see the note below
%       settingGroups Nx2 cell {label,{fields}}      (param-editor groups)
%       basicFields  the fields a protocol actually tunes  ({'contrastType',...});
%                    everything else in settingGroups is ADVANCED (see below)
%       sharedKeys   fields propagated between steps ({'trustLimitsK',...})
%       presets      struct of named default bundles (.default = launcher values)
%       tips         struct field->tooltip           (from the %Example comments)
%       enums        struct field->allowed values    (dropdown items)
%       labels       struct field->row label         (OPT-IN; absent = the field
%                    name, which is what a protocol is written in anyway)
%       note         one sentence a person has to know BEFORE ticking this step,
%                    appended to its tooltip in the Constructor ('' for most steps)
%       modalities   which modalities expose it       ({'LSCI'})
%       branch       'contrast' | 'cardiac' | 'angiogram' | '' (which file branch
%                    this step's column belongs to; wbFileModel also derives 'bolus'
%                    and 'myograph')
%       branchScope  'one' | 'all' | 'copy'  - HOW MANY of a recording's branch
%                    products the step consumes (see below)
%       fanOut       'flat' | 'animal' - the SHAPE of the fNames a call gets when
%                    several recordings are batched into it (see below)
%       fileOrder    'independent' | 'ordered' - whether the wrapper READS ACROSS
%                    its file list, i.e. whether the list may be cut up (see below)
%       refBranch    'contrast' | 'cardiac' | '' - which BRANCH of the animal's
%                    REFERENCE RECORDING this step prefers (see below)
%
% Notes (for the author):
%   * FOUR MODALITY FAMILIES ARE HERE: speckle (LSCI), the PRESSURE MYOGRAPH
%     (PMYO, a video), the WIRE MYOGRAPH (WMYO, a LabChart '.adicht' recording) and
%     FLUORESCENCE (EPFL, a '.cxd' stack).  The last of those is arriving a step at a
%     time: 'intensity' (the angiogram), 'intensityCycle' (the cardiac cycle) and
%     'intensityBolus' (the injection) are the three entry steps, the four shape-finding
%     and trace twins ('setRegionsI', 'segmentationI', 'dynamicSegmentationI',
%     'vasomotionI') sit after them, 'backgroundRemoval', 'ctth', 'topology',
%     'motionEnhancement' and 'pulsatilityI' after those, and the vascular pair
%     ('vesselTypesI', 'vascularTreeI') closes it.  THE TWELVE PIPELINE STEPS ARE ALL
%     HERE; what is still open is 'registrationI' and 'guidedI' (master-plan Q12).
%   * ONE STEP KEEPS A TWIN'S ID AND RUNS A DIFFERENT WRAPPER: 'pulsatilityI' points at
%     runIntensityPulsatility, not at runPulsatility [D13, author 07-Aug-2026].  The id
%     stays because a step declares ONE inGlob and that is the whole reason a twin needs
%     its own id; the WRAPPER differs because runPulsatility names the physical quantity
%     in its COLUMN PREFIX ('ps' = pulsatile flow) and would write psPI for a pulsatile
%     intensity, while runVasomotion names it in a TREE NODE and therefore twins safely.
%     Where a quantity is named is what decides whether a step can be twinned at all.
%   * FIVE EPFL STEPS ARE NOT TWINS, and each entry says why on its own line.
%     'backgroundRemoval' takes the dye's own scattered light off the picture, and there
%     is no dye and therefore no such haze on the speckle side.  'ctth' measures a
%     tracer, and there is no tracer.  'topology' WOULD run on a contrast product and
%     must not: its "vessel mask" would come from a flow image, so an occluded vessel
%     would be counted as tissue.  'motionEnhancement' has nothing on the speckle side
%     to twin FROM - what is there is a per-frame diameter quantised at a quarter of a
%     pixel, which is a different instrument rather than an earlier version of this one.
%     'pulsatilityI' has a speckle counterpart and cannot use it, for the reason in the
%     note above - so it is the one entry here whose ID is a twin's and whose WRAPPER is
%     its own.  A step being modality-blind is not the same as its ANSWER being
%     modality-blind, which is why none of the five is a one-line twin.
%     TWO STEPS CAN DECLARE ONE BRANCH - 'internalCycle' and 'intensityCycle' both
%     say 'cardiac' - because a recording is a '.rls' OR a '.cxd' and never both, so no
%     recording ever carries two products on that branch and the modality filter leaves
%     exactly one of the pair standing.
%   * INTENSITY TWINS.  A registry step declares ONE inGlob, and five workbench files
%     consume it, so a step that must serve both '_K'/'_BFI' and '_I' needs SEPARATE
%     IDS rather than a per-modality override - which is also the only shape that is
%     safe under the UNION registry a mixed working set uses.  intensityTwinOf below is
%     the one place that derivation lives; it is reachable as
%     wbStepRegistry('intensityTwin',step) and its own header says exactly what it
%     carries over and what it refuses to guess.  Four twins are at the end of this file
%     and are BUILT THROUGH IT rather than written out, so the shipped registry
%     exercises the derivation instead of only a test doing so.  Each then sets what the
%     helper deliberately clears, and three of them additionally set fields the helper
%     CANNOT derive: 'vasomotionI' the glob and the branch (a '*_BFI_r.mat' glob yields
%     the branch-agnostic '*_I_r.mat', which is right for a shape-finding step and wrong
%     for one that needs a clock, and branch 'contrast' is cleared to '' when the answer
%     is 'angiogram'); 'vesselTypesI' its refBranch; 'vascularTreeI' its glob, both
%     branch fields and one preset.  THE PATTERN IN ALL THREE IS THE SAME: the helper
%     rewrites the PRODUCT TOKEN and can say nothing about which BRANCH of that product
%     a step wants, because that is a fact about where the measurement was made.
%   * WHICH WAY ROUND THE VESSELS ARE IS NOT A SETTING ANYWHERE (author, 2026-08-08).
%     A contrast image's vessels are dark as a matter of physics, and a fluorescence
%     recording processed by this library has bright ones - that is the branch's
%     stated assumption, not a question put to the operator.  getVesselPolarity
%     resolves both arms from the product token, so no step here carries a field for
%     it and no step infers it from a file name of its own.
%   * NO MYOGRAPH STEP WRITES A REPORT IMAGE (author, 2026-08-01, re-affirmed
%     2026-08-02 when the wire myograph landed): the myograph narrates the three
%     ordinary lines per recording and nothing else, so every myograph step declares
%     no artifacts and the per-column PDF has nothing to assemble for them.  This
%     holds for the wire myograph too - what was recorded is looked at in the
%     interval editor, which opens on the channels with the comments marked.
%   * THE MYOGRAPH BLOCK IS SHAPED DIFFERENTLY, and deliberately.  ONE ENTRY STEP
%     PER MODALITY (myoVideo for the pressure myograph, labChart for the wire one)
%     creates the recording's '_MYO' triplet and every later step appends to that
%     same triplet in place - so the product carries a PRODUCT token and no stage
%     flag, one recording has exactly one myograph result set, and the 'below'
%     deletion set of a re-run is empty.  Both modalities write the SAME token on
%     purpose: the interval editor and the vasomotion step are literally the same
%     registry entries for each, which needs one inGlob.  They all declare branch
%     'myograph' rather than '': a branch-agnostic derived step is treated as
%     RECORDING-level by wbTypeSelection (ticked once on the anchor row and
%     inherited), which these are not.
%   * REAL gating fields (differ from the function name, 01 A1): runBFI->calculateBFI.
%     ONE, and it used to be two.  THE LIST DID NOT GROW when runInternalCycle was
%     renamed to runContrastInternalCycle (2026-08-07), and it SHRANK when the transit
%     step was rebuilt (2026-08-08): runCTTH's field was 'ctthCalculation' and is now
%     'runCTTH'.  Both times the sidecars on disk were migrated instead, because a
%     settings field is named after its writer and a deliberate mismatch costs more
%     than the migration does.
%   * BRANCH: 'contrast' is the temporal OR spatial contrast side (t|s) - the user
%     picks per s.contrastType and the analysis (segmentation/BFI/vasomotion/...)
%     runs on either; 'cardiac' is the internal-cycle side (c); '' is branch-
%     agnostic.  vasomotion and nvc are 'contrast' (t|s only, NOT c).  State stays
%     at the recording level here; per-file column filtering by branch/stage is a
%     Phase-3 concern (this field is the hint it uses).
%   * STAGE flag is a SINGLE token: the internal cycle is '_c_K', the angiogram is
%     '_a_I'.  A cycle's contrast base (t|s) is recorded in the SETTINGS, not the
%     name, since a project is not expected to mix bases - so the suffix stays simple,
%     showing only what the next step needs (t_K / s_K / c_K / a_I).
%     THERE IS NO EPOCH-AVERAGED '_e' STAGE ANY MORE (2026-08-05).  The NVC step
%     measures every stimulus repetition in place on the contrast product, so the
%     whole triplet it used to write - 4.7 GB per recording - stopped existing along
%     with the five pipeline steps that fed it.  wbFileModel no longer parses an '_e'
%     name either: the flag left the stage set with the stage.
%   * 'guided' maps to runGuidedContrast (intensity variant deferred).
%     setRegions is modelled perFile (only registration/vesselTypes are perAnimal).
%   * EXCEL EXPORT IS NOT A STEP (author, 2026-07-28).  It is neither configured in
%     the Constructor nor run by the Processor: it is a STANDALONE tool that reads a
%     finished session (spec §7/D9).  Everything here writes data products; the
%     export writes a report, and mixing the two put a workbook column in a
%     processing matrix that nobody wanted to tick.
%   * BRANCHSCOPE - one raw recording yields SEVERAL co-registered products ('_t_K',
%     '_c_K', ...), and the launchers are explicit about how many of them a
%     step touches: '*_t_K_r.mat' means one branch, '*_K_r.mat' means all of them.
%     A workbench ROW is the recording, not the file, so each step declares the
%     fan-out its launcher cell implies (wbExecutor>buildFNames does the resolving):
%       'one'  - a single file, disambiguated by branch/desiredStage.  The default,
%                and correct for a step that is meaningful on ONE branch only
%                (vasomotion on _t, pulsatility/vascularTree on _c, the entry steps
%                that read the raw recording).
%       'all'  - every branch product of the recording, as an Nx1 fNames column, the
%                way the launcher passes '*_K_r.mat' (BFI, dynamicSegmentation).
%                The wrapper's own per-file loop then covers the branches, one
%                workbench cell for the lot.
%       'copy' - the step RUNS on the contrast-side file and the result is inherited
%                by the other branches through the wrapper's s.fNamesCopyTo
%                (setRegions, segmentation).  This is what makes the interactive /
%                expensive work happen once per recording instead of once per branch,
%                and it mirrors launcher STEPs 3+5.  wbExecutor derives the target
%                list from the recording's own siblings, so the field never has to be
%                typed by hand.
%     branchScope applies to a perAnimal step too, where it fans out WITHIN each
%     member of the animal rather than across rows: 'one' gives the animal one file
%     per recording, 'all' gives it every branch product of every recording, still
%     as a single Nx1 column with the reference first.  It never multiplies the
%     animal's columns - a perAnimal step is one work item however many files it
%     hands the wrapper.
%     RESOLVED (Phase 6, 2026-07-29): registration and vesselTypes are 'all', and the
%     fluorescence twin inherits it - a recording's angiogram, cardiac cycle and bolus
%     are painted or inherited in one call, with refBranch saying which of the three
%     opens the editor ('bolus', where the tracer arrival that drives the guess is).
%     Their launcher cells are branch-wide ('Roi*_K_r.mat' and 'Roi*_BFI_r.mat',
%     both of which match _t AND _c), and leaving them at 'one' did not just under-
%     cover - it silently picked the WRONG file, because resolveStepInputs falls
%     back to the first dir() match when the step declares no branch, and
%     '_c_K_r.mat' sorts before '_t_K_r.mat'.  A two-pipeline recording was
%     therefore registered on its cardiac product only, and never got vessel types
%     onto '_t_BFI' - which is exactly where vasomotion's per-segment table lives.
%     This is the same trap the 'copy' note above records for setRegions.
%   * FANOUT - branchScope says how many files ONE recording contributes; fanOut
%     says how the files of SEVERAL recordings are laid out once the workbench
%     batches them into a single call (wbBatchPlan does the shaping):
%       'flat'   - the default: a COLUMN of every file of every recording in the
%                  call, and s.fNamesCopyTo in the runSegmentation convention (one
%                  ROW of targets per source file).  It is what every wrapper that
%                  simply loops its file list wants.
%       'animal' - a 2-D list whose ROWS ARE ANIMALS, and s.fNamesCopyTo in the
%                  setRegions ELEMENTWISE convention (the same size as fNames).
%     Only setRegions declares 'animal', and the reason is its carry-forward: the
%     ROIs drawn on a file are re-offered as editable objects on the NEXT file of
%     the same row and reset at the next row, so the row has to BE the animal.
%     That is exactly the launcher's list -
%     getFileNamesList(rootFolder,'*_t_K_r.mat','[A-Z]+\d+') - whose grouping
%     pattern is the animal token.  Left flat, every recording would arrive as its
%     own row and the drawn region could never persist from one recording of an
%     animal to the next, which is the whole point of the step (author, 2026-07-29).
%   * FILEORDER - branchScope says how many files one recording brings and fanOut
%     how several recordings are laid out; fileOrder says whether the resulting list
%     may be CUT UP.  It is what lets one recording fail without taking the rest of
%     the call with it: for an 'independent' step the executor invokes the wrapper
%     once per recording, so a throw reddens that one cell and the remaining
%     recordings still run; for an 'ordered' step it makes the single call it always
%     made, because the wrapper's cross-file state is exactly what would be lost.
%       'independent' - the wrapper LOOPS its file list and each iteration stands on
%                       its own.  The default, and true of everything that is not
%                       named below.
%       'ordered'     - the wrapper READS ACROSS the list: a template, a carry-
%                       forward, or a per-file index into something else.  Each of
%                       the four says WHY on its own line.
%     IT CANNOT BE DERIVED from the other fields, which is why it is declared:
%     nvc is perFile / flat / branchScope 'one' / no refBranch - the exact profile of
%     an independent step - yet it indexes a per-file stimulus list, so every rule
%     that guesses gets it wrong (author, 2026-07-31, of the step nvc replaced).
%     It also decides where a run STOPS after a failure: errors never abort mid-step,
%     but the run ends when the expansion reaches the first 'ordered' step while a
%     recording is still red, since a cross-file step given a half-processed set is
%     the one thing worse than stopping (spec D7/D8).
%   * RAW PRODUCER vs DERIVED CONSUMER (author, 2026-07-28).  A step whose inGlob
%     is NOT a '*.mat' glob reads the raw recording and writes a NEW, independent
%     triplet: contrast writes '_t_K' (or '_s_K' when the type's contrastType is
%     spatial), the internal cycle '_c_K'.  Every other step APPENDS to one of
%     those products.  So one recording can carry TWO independent result sets, and
%     the Constructor asks its question per (type, BRANCH): the raw producers are
%     ticked once per type and each of them brings its branch ROW into existence;
%     the derived consumers are ticked per row.  The split is DERIVED from inGlob
%     and the row membership from 'branch' (empty = branch-agnostic, i.e. offered
%     on every row) - wbTypeSelection owns both rules, and nothing anywhere lists
%     the members by name.
%   * REFBRANCH - the perAnimal steps work against the animal's REFERENCE
%     RECORDING, which is pinned once per animal (a recording IDENTITY, valid for
%     every file of that animal whatever its type or experimental group).  That
%     recording owns several branch products, so each step declares WHICH branch of
%     it it needs: registration templates on the contrast side ('_t'/'_s'), vessel
%     typing paints on the cardiac side ('_c'), and the vascular tree - which is
%     perFile but is derived on the cardiac product and inherited by its partners -
%     prefers '_c' too.  THE FLUORESCENCE PAIR SAYS 'bolus' FOR BOTH, and for the same
%     reason rather than a different one: a step that reads the timing runs where the
%     timing was measured, and on that branch the transit-time step measures it on the
%     injection.  '' = no preference.  It is the reference twin of
%     branchScope, resolved by wbRefBranch with the same fall-back-to-any-branch
%     rule wbExecutor>desiredStage uses for the inputs; the Constructor tab shows
%     the resolved file per animal and the executor consumes it.  NOTHING here is a
%     UI choice - the user pins a recording, not a branch.
%
%   * THE MODALITY FILTER PRUNES WHAT IT REMOVES (2026-08-01).  Filtering to one
%     modality - or to the union of the modalities a working set actually holds -
%     drops steps, and a surviving step may still NAME one of them: a step serving
%     two modalities lists both of its possible producers in requiresAny, and only
%     one of them exists once the filter has run.  So a second pass strikes every
%     id the filter removed out of the survivors' requires / requiresAny /
%     conflictsWith.  A cascade, a tooltip or a greyed-cell reason must never name
%     a step the user cannot see, and a requiresAny whose default producer belongs
%     to the OTHER modality would otherwise tick a step this one does not run.
%     A cellstr argument is the union (a mixed folder), a char one modality; the
%     order is always the registry's own.
%   * CONFLICTSWITH.  Two steps that are alternative ways to do one thing cannot
%     both be selected - the myograph's pre-set intervals versus its interval
%     editor, for instance.  Each names the other, and ticking one unticks its
%     conflicts on that row (wbTypeSelection>setOne), reporting them so the
%     Constructor logs the move; unticking triggers nothing.  Declaring a conflict
%     that is also a prerequisite is a REGISTRY BUG and is thrown at construction
%     rather than resolved silently - the two rules would pull the same box in
%     opposite directions and whichever won would be an accident of ordering.
%   * REQUIRES vs REQUIRESANY.  '*_K_r.mat' steps (setRegions, segmentation, BFI,
%     registration) read ANY branch product of a recording, so their real
%     prerequisite is "some entry step has run", not "contrast has run".  A purely
%     pulsatile protocol starts at internalCycle and never computes a contrast
%     product - listing 'contrast' as a hard requirement made it tick a step it
%     does not run.  They therefore declare requiresAny={'contrast','internalCycle'};
%     the FIRST entry is the default producer a one-click chain pulls in.  Hard
%     chains (dynamicSegmentation after segmentation, pulsatility after
%     BFI+internalCycle) stay in 'requires'.  wbPrereqs is the single definition
%     of "satisfied".
%
%   * ARTIFACTS, AND WHY THE TAIL LIST SURVIVED THE RENAME (2026-07-29).  Unified
%     reporting gave every report image a legible name, '<stem>_rep_<stage>.jpg',
%     so 'artifacts' now reads like the pages it finds ('_rep_segments.jpg' rather
%     than '_vs.jpg').  legacyArtifacts carries the pre-rename tails alongside, so a
%     folder processed before the rename still lists its reports; wbArtifacts unions
%     the two over ONE directory listing, so recognising the old names costs an
%     endsWith pass and nothing else.
%     The tempting next move - drop the per-step lists for a single '_rep_*' glob,
%     after which a wrapper that starts writing a new page would need no registry
%     edit at all - DOES NOT WORK, and the reason is worth writing down.
%     wbArtifacts resolves against the RECORDING's base name ('Mouse1*'), not
%     against a product, so every report of a recording matches that glob whatever
%     step wrote it: 'Mouse1_rep_contrast.jpg' and 'Mouse1_t_K_rep_segments.jpg'
%     alike.  The stage token in the tail is the ONLY thing that attributes an image
%     to a column.  With a bare '_rep_*' every column would claim every report, the
%     result list would name the wrong step, and each per-column PDF would be the
%     same document.  The list stays; the rename bought legibility, not brevity.
%   * BASIC vs ADVANCED (author, 2026-07-28).  Most users never touch most of a
%     step's settings.  basicFields names the handful a protocol is actually
%     written around - the ones the launchers' %Example blocks comment on - and the
%     settings panel renders those first, under 'Basic', with everything else
%     collapsed below under 'Advanced'.  It is a DISPLAY split only: both halves
%     resolve, save and invalidate identically, and a step with no basicFields
%     shows all of its fields as Basic.
%   * PARALLEL SWITCHES, AND WHY THEY SAY 'parfor' (author, 2026-07-30).  Three steps
%     carry a 'Parallel' group of logical fields that bound their own parfor: one
%     field per distinct parallel LOOP, not per step, which is why vasomotion has
%     three of them.  They are ADVANCED (no step lists them in basicFields) because a
%     protocol is never written around them - they are a machine decision.
%     Their labels ('Use parfor: segments', 'Use parfor: per pixel', ...) come from
%     the author verbatim and are a DELIBERATE exception to the house rule that a
%     user-visible string carries no code vocabulary: 'parfor' is jargon, and it is
%     the right jargon for an audience that runs MATLAB and knows what a pool costs.
%     Do not "helpfully" rewrite them to 'Parallel: per pixel'.
%     LABELS is the mechanism that made this possible - the settings panel used the
%     raw field name as the row label, and 'parforVasomotionSegments' is not a label.
%     It is opt-in per field, so nothing else in the panel changed.
%     parforMyographLines joined them with the myograph block (the per-line
%     vasomotion loop).  Two more axes exist in the library and are NOT here,
%     because their steps are not in the registry: parforCTTHPixels (.cxd/CTTH) and
%     fitDLSI's 'parforFit' name/value (DLSI).  They are reachable from their
%     launchers and inherit the same true-when-absent default.
%
% See also: wbFileModel, wbStateEngine, wbSettingsModel, wbInvalidate,
%           wbPrereqs, wbRefBranch, wbTypeSelection, wbTypePresets
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function reg = wbStepRegistry(modality, arg2, arg3)

% the two non-registry entry points: hand either one a step array and it gets the
% same rules, which is how a hand-built fragment (a test, a trimmed registry that
% a caller assembled itself) is checked and filtered by the definitions here
if nargin>=2 && ischar(modality) && strcmp(modality,'validate')
    validateRegistry(arg2);  reg = arg2;  return
end
if nargin>=3 && ischar(modality) && strcmp(modality,'filter')
    validateRegistry(arg2);  reg = filterTo(arg2, arg3);  return
end
% the twin derivation is reachable from outside for the same reason validate and
% filter are: it is a RULE, and a rule nothing can call is a rule nothing can check
if nargin>=2 && ischar(modality) && strcmp(modality,'intensityTwin')
    if nargin>=3, reg = intensityTwinOf(arg2, arg3); else, reg = intensityTwinOf(arg2); end
    return
end

reg = base();  reg(1) = [];   % empty 1x0 struct array with the right fields

% ---- 1. contrast (temporal) --------------------------------------------------
s = base();
s.id='contrast'; s.label='Contrast'; s.wrapper=@runContrastFromRLS;
% s.contrastType chooses the flag: temporal -> _t_K, spatial -> _s_K (same branch)
s.inGlob='*.rls'; s.outSuffix={'_t_K_d','_t_K_r','_t_K_s'}; s.outKind='new';
s.outTransform=struct('from','.rls','to','_t_K_r.mat');   % spatial: '_s_K_r'
s.gatingField='runContrastFromRLS'; s.requires={}; s.produces={'contrast'};
s.interactive=@(ss) isfield(ss,'manualMask') && isequal(ss.manualMask,1);
s.artifacts={'_rep_contrast.jpg'}; s.legacyArtifacts={'_c.jpg'}; s.branch='contrast';
s.settingGroups={ 'Contrast calculation',{'contrastType','contrastKernel','decimFactor','decimMethod'};
                  'Performance',{'procType'};
                  'Initial masking',{'trustLimitsK','trustLimitsI','minTrust','manualMask'} };
s.basicFields={'contrastType','contrastKernel','decimFactor','manualMask'};
s.sharedKeys={'trustLimitsK','trustLimitsI','libraryFolder'};
s.enums=struct('contrastType',{{'temporal','spatial'}},'decimMethod',{{'sharp','leaking'}}, ...
               'procType',{{'gpu','cpu'}});
s.presets=struct('default',struct('contrastType','temporal','contrastKernel',25, ...
    'decimFactor',25,'decimMethod','sharp','procType','gpu','trustLimitsK',[0.001 0.99], ...
    'trustLimitsI',[1 254],'minTrust',[0.99 0.99],'manualMask',0));
s.tips=struct('contrastType','''temporal'' or ''spatial''', ...
    'contrastKernel','must be odd: 25 for temporal, 5 or 7 for spatial', ...
    'decimFactor','output framerate = original / decimFactor', ...
    'decimMethod','''sharp'' (temporal only) or ''leaking''', ...
    'procType','''gpu'' for spatial if a high-end GPU is available, else ''cpu''', ...
    'trustLimitsK','[min max] expected contrast (fastest,slowest flow)', ...
    'trustLimitsI','[min max] expected intensity', ...
    'minTrust','per-pixel trust in relation to dark/saturated frame fraction', ...
    'manualMask','1 = interactive ROI selection of the analysed area');
reg(end+1)=s;

% ---- 2. internalCycle (cardiac) ---------------------------------------------
s = base();
s.id='internalCycle'; s.label='Internal cycle'; s.wrapper=@runContrastInternalCycle;
% THE WRAPPER WAS RENAMED WHEN THE SECOND INTERNAL-CYCLE STEP ARRIVED (2026-08-07), and
% THE GATING FIELD WENT WITH IT.  'intensityCycle' below collapses a .cxd onto the same
% 'cardiac' branch, so an unqualified 'runInternalCycle' no longer said which of the two
% had written a product.  Leaving the field spelt the old way was the alternative and was
% rejected: this file already documents a gating field that differs from its function
% name (calculateBFI) as a surprise to look out for, and a SECOND one invented on
% purpose would be worse than the migration it saves - which was 20 sidecars
% in claude-tests and nothing at all on the author's data disk.  The STEP ID is unchanged:
% it names the step, not the writer, and every session file keyed on it still resolves.
s.inGlob='*.rls'; s.outSuffix={'_c_K_d','_c_K_r','_c_K_s'}; s.outKind='new';
s.outTransform=struct('from','.rls','to','_c_K_r.mat');
s.gatingField='runContrastInternalCycle'; s.requires={}; s.produces={'cardiac'};
s.artifacts={'_rep_cycle-detect.jpg','_rep_cycle-average.jpg'};
s.legacyArtifacts={'_ic1.jpg','_ic2.jpg'}; s.branch='cardiac';
s.settingGroups={ 'Contrast calculation',{'contrastKernelS'};
                  'Pulse detection area',{'maskLimitsK','maskLimitsI'};
                  'Output masking',{'trustLimitsK','trustLimitsI','minTrust'};
                  'Frequency band',{'maxFrqIni','minFrqIni','rangeFrq'};
                  'Exclusion criteria',{'excludeFirstNCycles','coeffsSTD','coeffsRel','coeffsAbs'};
                  'Cycle calculation',{'method','decimationSpace','framesToAverage', ...
                     'contrastKernelT','contrastKernelPreproc','interpFactor','smoothCoef1','minPromCoef'} };
s.basicFields={'method','maxFrqIni','minFrqIni','excludeFirstNCycles'};
% trustLimitsK/I stay shared with the contrast step and now mean the SAME thing there
% and here - the propagated mask - which is what the link always claimed.  The pulse
% detection area is this step's alone and is deliberately NOT shared; minTrust is not
% shared on the contrast step either.
s.sharedKeys={'trustLimitsK','trustLimitsI','libraryFolder'};
s.enums=struct('method',{{'sLSCIMM','tLSCIMM','ltLSCIMM','sLSCIMMM'}});
s.presets=struct('default',struct('maskLimitsK',[0.01 0.3],'maskLimitsI',[5 250], ...
    'trustLimitsK',[0.001 0.99],'trustLimitsI',[1 254],'minTrust',[0.99 0.99], ...
    'contrastKernelS',5,'maxFrqIni',20,'minFrqIni',1,'rangeFrq',1,'excludeFirstNCycles',0, ...
    'coeffsSTD',[3 2 2 2 2 3 3 2 2],'coeffsRel',[0.5 0.1],'coeffsAbs',2,'method','sLSCIMM', ...
    'decimationSpace',4,'framesToAverage',1,'contrastKernelT',25,'contrastKernelPreproc',5, ...
    'interpFactor',10,'smoothCoef1',1/3,'minPromCoef',1/4));
s.tips=struct('method','sLSCIMM recommended; ltLSCIMM for high-quality data', ...
    'contrastKernelS','spatial contrast kernel, must be odd: 5 or 7', ...
    'contrastKernelPreproc','spatial kernel used to detect the pulse, must be odd: 5 or 7', ...
    'maxFrqIni','initial max frequency of interest, Hz', ...
    'minFrqIni','initial min frequency of interest, Hz', ...
    'coeffsSTD','pulse-rejection coefficients relative to feature std', ...
    'maskLimitsK','contrast range of the pixels averaged to detect the pulse', ...
    'maskLimitsI','intensity range of the pixels averaged to detect the pulse', ...
    'trustLimitsK','[min max] expected contrast (fastest,slowest flow)', ...
    'trustLimitsI','[min max] expected intensity', ...
    'minTrust','per-pixel trust in relation to dark/saturated frame fraction');
reg(end+1)=s;

% ---- 3. setRegions ----------------------------------------------------------
s = base();
s.id='setRegions'; s.label='Regions'; s.wrapper=@setRegions;
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setRegions'; s.requires={}; s.requiresAny={'contrast','internalCycle'};
s.produces={'regionsMask'};
s.artifacts={'_rep_regions.jpg'};                     % NEW in Session 3: the drawn
                                                      % regions, drawn or copied - no
                                                      % legacy tail, it wrote none
s.interactive=true;                                   % always opens the ROI editor
s.branch=''; s.branchScope='copy';                    % draw on _t, inherit onto _c/_e
s.fanOut='animal';                                    % rows = animals: ROIs carry along a row
s.fileOrder='ordered';                                % the ROIs drawn on one file are re-offered
                                                      % on the NEXT file of the row: cut the row up
                                                      % and the carry-forward is gone
s.settingGroups={ 'Regions',{'nRegions'} };           % the drawing is interactive; how MANY may
                                                      % be drawn is the one thing a protocol fixes
s.basicFields={'nRegions'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('nRegions',1));     % the wrapper's own default is unlimited (Inf)
s.tips=struct('nRegions','how many regions you may draw on each recording');
reg(end+1)=s;

% ---- 4. segmentation --------------------------------------------------------
s = base();
s.id='segmentation'; s.label='Segmentation'; s.wrapper=@runSegmentation;
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runSegmentation'; s.requires={}; s.requiresAny={'contrast','internalCycle'};
s.produces={'segmentation'};
s.artifacts={'_rep_categories.jpg','_rep_segments.jpg'};
s.legacyArtifacts={'_cm.jpg','_vs.jpg'};              % '_vs' came from showSegmentsPreview
s.branch=''; s.branchScope='copy';                % computed on contrast side, copied to cardiac
s.settingGroups={ 'Contrast',{'trustLimitsK'};
                  'Categorization',{'lSizeN','sSizeN','sens','sSizeScale','deSens','lThinN','imOpen','iEdge','eEdge','diffusionSchedule'};
                  'Labelling & traces',{'sStat','sMinL','prchNSize','correctNodes','simR','difR'};
                  'Parallel',{'parforSegmentationLabels'} };
                  % no 'Copy to siblings' panel: branchScope 'copy' means wbExecutor
                  % derives s.fNamesCopyTo from the recording's own branch products
s.basicFields={'lSizeN','sSizeN','sens','sMinL'};
% parforSegmentationLabels is SHARED with dynamicSegmentation: both steps drive the
% same core (getSegmentationLabels), so the machine-level choice is one tick, not two.
s.sharedKeys={'trustLimitsK','prchNSize','sMinL','correctNodes','simR','difR','sStat', ...
    'parforSegmentationLabels','libraryFolder'};
s.enums=struct('sStat',{{'median','mean'}},'diffusionSchedule',{{'multiscale','single'}});
% diffusionSchedule was NOT a setting until 2026-08-08: getPixelCategories chose the
% schedule off the product token, so a tuning decision was welded to a file name.  It
% is 'multiscale' here and 'single' on the intensity twin, which is exactly what each
% branch did before - nobody has measured which suits a fluorescence image, and now
% they can without a code change.
s.presets=struct('default',struct('trustLimitsK',[0.001 0.99],'lSizeN',141,'sSizeN',9, ...
    'sens',0.1,'sSizeScale',1,'deSens',1,'lThinN',2,'imOpen',0,'iEdge',2,'eEdge',2, ...
    'diffusionSchedule','multiscale', ...
    'categories',{{'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}}, ...
    'sStat','median','sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3,'difR',0.4, ...
    'parforSegmentationLabels',true));
s.labels=struct('parforSegmentationLabels','Use parfor: segment labels', ...
                'diffusionSchedule','Vessel smoothing');
s.tips=struct('lSizeN','odd, ~2x the largest vessel', ...
    'sSizeN','odd, ~2x the small-vessel diameter', ...
    'sens','segmentation sensitivity (raise to catch faint vessels)', ...
    'sMinL','minimum segment length', ...
    'prchNSize','parenchymal pixel neighbourhood', ...
    'diffusionSchedule',['how hard the picture is smoothed before vessels are looked ' ...
       'for in it.  Several passes from coarse to fine cope with a grainy image; one ' ...
       'pass keeps more of the small vessels'], ...
    'parforSegmentationLabels',parforTip('the per-segment label growing'));
reg(end+1)=s;

% ---- 5. dynamicSegmentation -------------------------------------------------
s = base();
s.id='dynamicSegmentation'; s.label='Dynamic segmentation'; s.wrapper=@runDynamicSegmentation;
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runDynamicSegmentation'; s.requires={'segmentation'}; s.produces={'dynamicSeg'};
s.artifacts={'_rep_segments.jpg'};                    % re-emitted, overwriting the static one
s.legacyArtifacts={'_vs.jpg'}; s.branch=''; s.branchScope='all';
s.settingGroups={ 'Labelling (match segmentation)',{'sMinL','prchNSize','correctNodes','simR','difR'};
                  'Dynamic segmentation',{'sMinP2R2','sMaxLBI','sMaxCLR','sMaxKK','iniNSize','sMaxP2D'};
                  'Quality & interpolation',{'gSizeN','minOverlapMask','minOverlapSelf','pInterpF'};
                  'Parallel',{'parforSegmentationLabels'} };
s.basicFields={'sMinL','minOverlapMask'};
s.sharedKeys={'sMinL','prchNSize','correctNodes','simR','difR','parforSegmentationLabels', ...
    'libraryFolder'};
s.presets=struct('default',struct('sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3, ...
    'difR',0.4,'sMinP2R2',0.95,'sMaxLBI',(1/7)/15,'sMaxCLR',1.3,'sMaxKK',0.3,'iniNSize',7, ...
    'sMaxP2D',3,'gSizeN',3,'minOverlapMask',0.6,'minOverlapSelf',0.2,'pInterpF',4, ...
    'parforSegmentationLabels',true));
s.labels=struct('parforSegmentationLabels','Use parfor: segment labels');
s.tips=struct('sMinP2R2','min accepted R^2 of the 3-degree polynomial fit', ...
    'sMaxCLR','max chord-length ratio (1 straight, 2 coiled)', ...
    'minOverlapMask','min overlap of centre line with the per-frame mask', ...
    'parforSegmentationLabels',parforTip('the per-segment label growing'));
reg(end+1)=s;

% ---- 6. guided --------------------------------------------------------------
s = base();
s.id='guided'; s.label='Guided'; s.wrapper=@runGuidedContrast;
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runGuidedContrast'; s.requires={'segmentation'}; s.produces={'guidedTraces'};
s.needsRaw=true; s.branch='contrast';
s.settingGroups={};                                   % uses sMap + raw; no tunable params
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct());
reg(end+1)=s;

% ---- 7. registration (per animal) --------------------------------------------
s = base();
s.id='registration'; s.label='Registration'; s.wrapper=@runRegistration;
s.arity='perAnimal';                                   % 2-D fNames row; col 1 = template
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runRegistration';                      % code writes runRegistration (header drift A4)
s.requires={}; s.requiresAny={'contrast','internalCycle'}; s.produces={'registered'};
s.interactive=@(ss) ~(isfield(ss,'silent') && isscalar(ss.silent) && isequal(ss.silent,true));
s.artifacts={'_rep_registration.jpg'};                % was a .png; now a report page
s.legacyArtifacts={'_registration.png'};
s.branch=''; s.refBranch='contrast';                  % template = the reference's _t|_s
s.branchScope='all';                                   % EVERY product of every member (launcher: 'Roi*_K_r.mat')
s.fileOrder='ordered';                                 % column 1 is the TEMPLATE every other
                                                       % file is registered onto
s.settingGroups={ 'Registration',{'tFormType','matchSegmentation','prchNSize','silent','forceMethod','rotationLimit'} };
s.basicFields={'tFormType','silent','forceMethod'};
s.sharedKeys={'prchNSize','libraryFolder'};
s.enums=struct('tFormType',{{'affine','rigid','similarity','translation'}}, ...
               'forceMethod',{{'auto','correlation','intensity'}});
s.presets=struct('default',struct('tFormType','affine','matchSegmentation',true,'prchNSize',50, ...
    'silent',true,'forceMethod','auto','rotationLimit',45));
s.tips=struct('silent','true = pick the best transform automatically (no manual landmarks)', ...
    'forceMethod','''auto'' tries both and keeps the best; force one only if auto misbehaves', ...
    'rotationLimit','degrees; reject registrations rotating beyond this ([] = none)');
reg(end+1)=s;

% ---- 8. BFI -----------------------------------------------------------------
s = base();
s.id='BFI'; s.label='BFI'; s.wrapper=@runBFI;
s.inGlob='*_K_r.mat'; s.outSuffix={'_BFI_d','_BFI_r','_BFI_s'}; s.outKind='new';
s.outTransform=struct('from','_K_r.mat','to','_BFI_r.mat');   % branch-preserving strrep
s.gatingField='calculateBFI';                         % REAL field (differs from fn name)
s.requires={}; s.requiresAny={'contrast','internalCycle'}; s.produces={'bfi'};
s.branch=''; s.branchScope='all';                     % _t AND _c both need a BFI product
s.settingGroups={ 'Conversion',{'deleteOriginal','method'} };
s.basicFields={'deleteOriginal'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('method',{{'basic'}});
s.presets=struct('default',struct('deleteOriginal',false,'method','basic'));
s.tips=struct('deleteOriginal','delete the original _K_ triplet after conversion', ...
    'method','only "basic" (=1/K^2) is available');
reg(end+1)=s;

% ---- 9. vasomotion ---------------------------------------------------------
s = base();
s.id='vasomotion'; s.label='Vasomotion'; s.wrapper=@runVasomotion;
s.inGlob='*_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';   % contrast-side BFI (_t_BFI or _s_BFI)
s.gatingField='runVasomotion'; s.requires={'BFI'}; s.produces={'vasomotion'};
s.branch='contrast';
s.settingGroups={ 'Bands',{'vFR','cFR','wFR','wVPO'};
                  'Normalisation',{'normalisation','normsize','tgtFS'};
                  'Peaks & percentiles',{'pcts','otsuMaxN','otsuElbow','nPeakProm'};
                  'Signals & levels',{'vsmSignals','segVsmReturn','ppxVsmReturn'};
                  'Parallel',{'parforVasomotionSegments','parforVasomotionPixels', ...
                              'parforVasomotionAveraging'} };
s.basicFields={'vFR','cFR','tgtFS','vsmSignals'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('normalisation',{{'mean','median','mmean','mmedian'}});
s.presets=struct('default',struct('vFR',[0.05 0.25],'cFR',[0.4 0.6],'wFR',[0.01 1],'wVPO',10, ...
    'normalisation','median','normsize',101,'tgtFS',1,'pcts',0:10:100,'otsuMaxN',5, ...
    'otsuElbow',0.05,'nPeakProm',0.10, ...
    'vsmSignals',{{'sData','dvsData','dvsDiameter','gsData'}}, ...
    'segVsmReturn',{{'bands','moments','series','clustering','spectrum'}},'ppxVsmReturn',[], ...
    'parforVasomotionSegments',true,'parforVasomotionPixels',true, ...
    'parforVasomotionAveraging',true));
% This step has THREE different parallel loops, so it gets three switches; the wording
% is the author's own and deliberately says 'parfor' (see the note under LABELS below).
s.labels=struct('parforVasomotionSegments','Use parfor: segments', ...
    'parforVasomotionPixels','Use parfor: per pixel', ...
    'parforVasomotionAveraging','Use parfor: segment averaging');
s.tips=struct('vFR','vasomotion frequency band [lo hi], Hz', ...
    'cFR','control (cardiac) frequency band [lo hi], Hz', ...
    'nPeakProm','VB peak-count prominence as a fraction of the band range', ...
    'segVsmReturn','which per-segment levels to store (set of tokens)', ...
    'ppxVsmReturn','per-pixel analysis ([] = off; e.g. {''bands''} to enable)', ...
    'parforVasomotionSegments',parforTip('the per-segment analysis'), ...
    'parforVasomotionPixels',parforTip('the per-pixel analysis'), ...
    'parforVasomotionAveraging',parforTip('the per-segment averaging'));
reg(end+1)=s;

% ---- 10. pulsatility --------------------------------------------------------
s = base();
s.id='pulsatility'; s.label='Pulsatility'; s.wrapper=@runPulsatility;
s.inGlob='*_c_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runPulsatility'; s.requires={'BFI','internalCycle'}; s.produces={'pulsatility'};
s.branch='cardiac';
s.settingGroups={ 'Harmonic model',{'nHarm'};
                  'Analysis levels',{'segPulsReturn','ppxPulsReturn'};
                  'Parallel',{'parforPulsatilityPixels'} };
s.basicFields={'nHarm'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('nHarm',5, ...
    'segPulsReturn',{{'markers','model','reconstruction'}},'ppxPulsReturn',{{'markers'}}, ...
    'parforPulsatilityPixels',true));
s.labels=struct('parforPulsatilityPixels','Use parfor: per pixel');
s.tips=struct('nHarm','number of harmonics in the sinusoidal cardiac model', ...
    'segPulsReturn','per-segment levels: markers / model / reconstruction', ...
    'ppxPulsReturn','per-pixel maps ([] = off; non-empty enables marker maps)', ...
    'parforPulsatilityPixels',parforTip('the per-pixel fit'));
reg(end+1)=s;

% ---- 11. nvc (stimulus-locked response, every repetition on its own) ---------
s = base();
s.id='nvc'; s.label='NVC response'; s.wrapper=@runNVC;
% THE WHOLE RECORDING, not an average of it.  This step cuts every stimulus
% repetition out of the contrast-side BFI product and measures each one separately,
% which is why there is no epoch-averaged product anywhere in the pipeline any more -
% the step it replaces built a second triplet whose only content was that average.
% The input is exactly what vasomotion and the vasoreactivity fit take.
s.inGlob='*_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runNVC';
% REQUIRES NAMES STEPS, NOT CAPABILITIES.  wbPrereqs matches these against step IDS;
% 'produces' is a separate vocabulary nothing in the dependency graph reads.
s.requires={'BFI'}; s.requiresAny={'contrast'}; s.produces={'nvcResponse'};
% A ROW, NOT A FILE BRANCH.  Rows exist only where a RAW PRODUCER declares one
% (wbTypeSelection>branchesOf): contrast, cardiac, myograph.  The step this replaces
% was first declared on a branch of its own, the one wbFileModel>branchOf used to
% derive from an '_e' file name, and was untickable on every row that exists while
% every test passed.  Neither the branch nor the file exists now, so the trap is gone;
% the value is what it is for the reason that outlived it.
s.branch='contrast';
% THE STIMULUS TIMING IS A PER-FILE LIST indexed by position in fNames, so the file
% order is part of the meaning: cut the list and the files up differently and every
% epoch is cut at the wrong time.
s.fileOrder='ordered';
% THE ONE INTERACTIVE SETTING IS THE EPOCH EDITOR.  It replaced the old picker and its
% misspelt switch; a step is interactive exactly when a person has to look at a window
% before the run can go on.
s.interactive=@(ss) isfield(ss,'nvcEditEpochs') && isequal(ss.nvcEditEpochs,true);
s.artifacts={'_rep_nvccuts.jpg','_rep_nvcconfidence.jpg','_rep_nvcresponse.jpg', ...
    '_rep_nvctrials.jpg'};
% THERE IS NO FITTING GROUP ANY MORE (author, 06-Aug-2026: "at this stage I do not see
% value in fitting the data").  What replaced it is two groups of trust settings - one
% that judges a repetition OF ONE SEGMENT and one that judges the repetition for the
% whole recording - because those are two different decisions made at two levels.
s.settingGroups={ 'Stimulation',{'stimStartType','stimOffset','epochsN','epochDurationSec', ...
                     'epochBaselineSec','epochStimStartSec','stimDurationSec','epochFinaleSec'};
                  'How a response is measured',{'nvcPeakGraceSec','nvcAreaPcts'};
                  'Signals & levels',{'nvcSignals','segNvcReturn','ppxNvcReturn'};
                  'How much to trust one repetition of one segment', ...
                     {'nvcConfThreshold','nvcConfMinThreshold','nvcReturnScale', ...
                      'nvcDevRules','nvcDevScale','rejectFirstEpoch'};
                  'Which repetitions the recording trusts', ...
                     {'nvcEpochAreaFrac','nvcEpochConnected','nvcEditEpochs'};
                  'Representative repetition',{'nvcRepresentative'};
                  'Parallel',{'parforNvcSegments','parforNvcPixels'} };
s.basicFields={'stimStartType','stimOffset','epochsN','epochDurationSec', ...
    'epochBaselineSec','epochStimStartSec','stimDurationSec','epochFinaleSec', ...
    'nvcEditEpochs'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('stimStartType',{{'offset','manual'}});
% nvcEditEpochs IS TRUE HERE AND FALSE IN THE WRAPPER, and that split is deliberate:
% the wrapper must be scriptable, the workbench is where a person is choosing.
s.presets=struct('default',struct('stimStartType','offset','stimOffset',0,'epochsN',20, ...
    'epochDurationSec',30,'epochBaselineSec',[0 10],'epochStimStartSec',10, ...
    'stimDurationSec',5,'epochFinaleSec',[-5 0], ...
    'nvcPeakGraceSec',2,'nvcAreaPcts',[10 50 90], ...
    'nvcSignals',{{'sData','gsData','dvsData','dvsDiameter'}}, ...
    'segNvcReturn',{{'levels','amplitudes','times'}},'ppxNvcReturn',[], ...
    'nvcConfThreshold',0.6,'nvcConfMinThreshold',0.2,'nvcReturnScale',5, ...
    'nvcDevRules',false,'nvcDevScale',3,'rejectFirstEpoch',1, ...
    'nvcEpochAreaFrac',0.10,'nvcEpochConnected',true,'nvcEditEpochs',true, ...
    'nvcRepresentative',false, ...
    'parforNvcSegments',true,'parforNvcPixels',true));
s.labels=struct('nvcEditEpochs','Check the repetitions myself', ...
                'nvcRepresentative','Replace the recording with the average', ...
                'parforNvcSegments','Use parfor: segments', ...
                'parforNvcPixels','Use parfor: per pixel');
s.tips=struct( ...
    'stimStartType','a fixed delay after the recording starts, or a list of clock times', ...
    'stimOffset','seconds from the start of the recording to the start of the FIRST repetition', ...
    'epochsN','how many times the stimulus was repeated', ...
    'epochDurationSec','how long one repetition is, seconds - stimulus and recovery together', ...
    'epochBaselineSec','the quiet stretch at the start of each repetition the response is measured against, seconds', ...
    'epochStimStartSec','when the stimulus starts within the repetition, seconds', ...
    'stimDurationSec','how long the stimulus lasts, seconds', ...
    'epochFinaleSec','the stretch at the END of the repetition where flow should be back to baseline, counted back from its end', ...
    'nvcPeakGraceSec',['how long after the stimulus stops the peak may still arrive, ' ...
       'seconds.  On the reference recording the peak is inside two seconds for 95 % ' ...
       'of segments'], ...
    'nvcAreaPcts',['the response is timed by when these percentages of the flow ' ...
       'increase had been delivered.  Whole percentages, and each one becomes a ' ...
       'measurement of its own'], ...
    'nvcSignals','which traces to analyse', ...
    'segNvcReturn',['which measurements to store: levels, amplitudes and times.  All ' ...
       'of them are always measured - this only decides what is kept'], ...
    'ppxNvcReturn',['per-pixel maps, one per repetition.  Leave it empty to skip them: ' ...
       'this is the slowest and largest thing the step can be asked for'], ...
    'nvcConfThreshold',['the overall confidence a repetition of a segment has to reach ' ...
       'to be used, from 0 to 1.  It is the geometric mean of every check, so one ' ...
       'failed check pulls it to zero'], ...
    'nvcConfMinThreshold',['the lowest any single check may score.  This is what ' ...
       'cancels a repetition for one clear reason rather than for general mediocrity'], ...
    'nvcReturnScale',['how far flow may still be from its own baseline at the end of a ' ...
       'repetition, in noise levels, before the score starts to fall.  A short rest ' ...
       'between stimuli lowers this score for every segment alike'], ...
    'nvcDevRules',['also compare every repetition of a segment against the other ' ...
       'repetitions of the same segment.  Off by default: the checks that are always ' ...
       'on judge a repetition without reference to any other, so they behave the same ' ...
       'on four repetitions and on forty'], ...
    'nvcDevScale',['how far a repetition may sit from the others of its segment before ' ...
       'it scores zero, in robust noise levels.  With only a handful of repetitions ' ...
       'that spread is itself noisy: at five repetitions about one ordinary repetition ' ...
       'in five already scores below a half'], ...
    'rejectFirstEpoch','never trust the first repetition', ...
    'nvcEpochAreaFrac',['how much of the imaged area has to respond before the ' ...
       'recording trusts a repetition.  Only a part of the field is ever a responder, ' ...
       'so a tenth of it responding coherently can still be a good repetition'], ...
    'nvcEpochConnected',['count only the largest connected responding region, so ' ...
       'scattered segments that passed the checks by chance do not add up to a ' ...
       'response'], ...
    'nvcEditEpochs','look at every repetition and change the decisions yourself', ...
    'nvcRepresentative',['replace the recording with the average of the trusted ' ...
       'repetitions.  THIS CANNOT BE UNDONE: the individual repetitions and the ' ...
       'recording clock are not kept, only the report page records which ' ...
       'repetitions went in, and this step cannot be run a second time on what ' ...
       'it produced'], ...
    'parforNvcSegments',parforTip('the per-segment measurement'), ...
    'parforNvcPixels',parforTip('the per-pixel analysis'));
reg(end+1)=s;

% ---- 12. fitVasoreactivity (drug / gas challenge response fit) ---------------
s = base();
s.id='fitVasoreactivity'; s.label='Vasoreactivity fit'; s.wrapper=@runFitVasoreactivity;
% THE WHOLE RECORDING, not an average of it: a drug challenge is one continuous trace
% with one injection, so this reads the contrast-side BFI product (_t_BFI or _s_BFI)
% exactly as the vasomotion step does.  An epoch or cardiac average carries no
% recording clock to place the injection on, and the wrapper refuses one by name.
s.inGlob='*_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runFitVasoreactivity';
% REQUIRES NAMES STEPS, NOT CAPABILITIES.  wbPrereqs matches these against step IDS;
% 'produces' is a separate vocabulary nothing in the dependency graph reads.
s.requires={'BFI'}; s.requiresAny={'contrast'}; s.produces={'vasoreactivityFit'};
% A ROW, NOT A FILE BRANCH.  Rows exist only where a RAW PRODUCER declares one
% (wbTypeSelection>branchesOf): contrast, cardiac, myograph.  '' would be wrong too -
% it offers the step on the cardiac and myograph rows, where a contrast BFI product
% cannot exist, and it would also make the step recording-level.
s.branch='contrast';
% s.injectionSec IS A PER-FILE LIST, so the file order is part of the meaning: cut the
% list and the files up differently and every fit anchors to the wrong time.
s.fileOrder='ordered';
s.artifacts={'_rep_vasoreactivityfit.jpg'};
s.settingGroups={ 'Challenge',{'injectionSec','baselineSec'};
                  'Model',{'vrcModel','tgtFS'};
                  'Response limits',{'vrcAucSec','vrcStealK','vrcStealFrac','vrcStealSec','vrcWashFrac'};
                  'Analysis levels',{'segVrcReturn','ppxVrcReturn'};
                  'Fitting',{'vrcStarts'};
                  'Parallel',{'parforVrcSegments','parforVrcPixels'} };
s.basicFields={'injectionSec','baselineSec','tgtFS'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('vrcModel',{{'bateman'}});
s.presets=struct('default',struct('vrcModel','bateman','tgtFS',1, ...
    'injectionSec',[],'baselineSec',[],'vrcAucSec',2700,'vrcStealK',2, ...
    'vrcStealFrac',0.15,'vrcStealSec',60,'vrcWashFrac',0.1, ...
    'segVrcReturn',{{'markers','model','reconstruction'}},'ppxVrcReturn',[], ...
    'vrcStarts',16,'parforVrcSegments',true,'parforVrcPixels',true));
s.labels=struct('parforVrcSegments','Use parfor: segments', ...
                'parforVrcPixels','Use parfor: per pixel');
s.tips=struct( ...
    'injectionSec',['when the drug was given, seconds from the start of the ' ...
       'recording.  One number per recording, in the order the files are listed'], ...
    'baselineSec',['the quiet stretch before the injection the response is measured ' ...
       'against, [start end] in seconds.  It has to end before the injection'], ...
    'vrcModel',['a rising and a falling phase on a drifting baseline.  The reported ' ...
       'amplitude is the peak change in the trace''s own units'], ...
    'tgtFS',['the traces are averaged down to this rate before fitting.  The response ' ...
       'takes minutes, so 1 Hz is plenty'], ...
    'vrcAucSec','how long after the injection the area under the response is taken over, seconds', ...
    'vrcStealK','how far below baseline, in baseline noise levels, counts as a dip', ...
    'vrcStealFrac','or this fraction of the response, whichever is the larger drop', ...
    'vrcStealSec','and how long it has to stay there before the recording is flagged, seconds', ...
    'vrcWashFrac','washout time = when the response is back to this fraction of its peak', ...
    'segVrcReturn','per-segment levels: markers / model / reconstruction', ...
    'ppxVrcReturn',['per-pixel maps ([] = off).  This fits every pixel over the whole ' ...
       'recording and takes many minutes'], ...
    'vrcStarts',['how many start points the fit is tried from.  Fewer is faster and ' ...
       'more likely to settle on a wrong answer'], ...
    'parforVrcSegments',parforTip('the per-segment fit'), ...
    'parforVrcPixels',parforTip('the per-pixel analysis'));
reg(end+1)=s;

% ---- 13. vesselTypes (per animal) --------------------------------------------
s = base();
s.id='vesselTypes'; s.label='Vessel types'; s.wrapper=@setVesselTypes;
s.arity='perAnimal';                                   % 2-D fNames row; col 1 = reference
s.inGlob='*_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setVesselTypes'; s.requires={'BFI','segmentation'}; s.produces={'vesselTypes'};
s.artifacts={'_rep_vesseltypes.jpg'};                 % NEW in Session 3: the artery/vein
                                                      % map, painted or inherited - no
                                                      % legacy tail, it wrote none
s.interactive=true;                                   % paint GUI (per-file skipped under useReference)
s.branch=''; s.refBranch='cardiac';                   % paint target = the reference's _c
s.branchScope='all';                                  % EVERY product of every member (launcher: 'Roi*_BFI_r.mat')
s.fileOrder='ordered';                                % column 1 is the PAINT TARGET the rest
                                                      % inherit under useReference
s.settingGroups={ 'Reference',{'useReference'} };
s.basicFields={'useReference'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('useReference',false));
s.tips=struct('useReference','true = paint the first (reference) file only, inherit to the rest');
reg(end+1)=s;

% ---- 14. vascularTree -------------------------------------------------------
s = base();
s.id='vascularTree'; s.label='Vascular tree'; s.wrapper=@setVascularTree;
s.inGlob='*_c_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setVascularTree'; s.requires={'vesselTypes','pulsatility'}; s.produces={'hierarchy'};
s.interactive=@(ss) ~(isfield(ss,'autoOnly') && isscalar(ss.autoOnly) && isequal(ss.autoOnly,true));
s.branch='cardiac'; s.refBranch='cardiac';   % the hierarchy is derived on the _c side
% TWO SETTINGS LEFT THIS PANEL AT S10, and they had not been connected to anything for
% some time: 'phiWeights' ([foot peak -PI]) and 'useHarmonicPhase' were the flow
% potential's controls BEFORE s.flowParams replaced them with a per-parameter weight and
% enable, and setVascularTree reads neither.  A panel offering two settings nothing reads
% is the stub this branch has been removing, and the fluorescence twin below would have
% inherited both.  Launcher_pulsatility_vasomotion's two lines went with them.
% AND FOUR ARRIVED, WHICH THE CORE HAS HONOURED ALL ALONG AND THIS PANEL NEVER OFFERED.
% getVascularTree reads connectivity/minBorder to decide which segments touch, and the
% two bridge radii to rejoin a vessel that a crossing one splits in the 2D projection;
% all four land in H.params, so all four change the derived hierarchy.  Until now only
% a launcher line could reach them, which meant a protocol could move one and a
% workbench run of the same step had no way to follow.  The defaults below are the
% wrapper's own - the values Launcher_intensity's STEP 14 already sets - so no
% hierarchy derived so far changes.  ADVANCED rather than basic: 8/1/50/0 have suited
% every recording processed here, and bridgeTipRadius is the one most likely to move,
% because it is a distance in PIXELS and so depends on the magnification.
% The tree editor's own remaining controls - the parenchyma degree caps and the tip
% band - stay editor-only on purpose: they are read off a hierarchy already on screen
% and re-derived in place, which is not a decision a protocol makes in advance.
s.settingGroups={ 'Hierarchy derivation',{'autoOnly','propagatePartners'}; ...
                  'Segment adjacency',{'connectivity','minBorder'}; ...
                  'FOV bridging',{'bridgeTipRadius','bridgeWallRadius'} };
s.basicFields={'autoOnly'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('autoOnly',false,'propagatePartners',{{'t','s'}}, ...
    'connectivity',8,'minBorder',1,'bridgeTipRadius',50,'bridgeWallRadius',0));
% THE TOOLTIPS ARE THE LAUNCHER'S OWN COMMENTS, WORD FOR WORD, so the panel row and the
% protocol line a user copies from say the same thing rather than two paraphrases.
s.tips=struct('autoOnly','true = derive & save without opening the tree editor', ...
    'propagatePartners','after a _c file is derived, copy the hierarchy to these partners', ...
    'connectivity','4 or 8 pixel adjacency between segments', ...
    'minBorder','minimum shared-border pixels to treat two segments as touching', ...
    'bridgeTipRadius','px gap searched from the vessel ends (0 = off)', ...
    'bridgeWallRadius','px gap searched from the vessel sides (usually 0)');
reg(end+1)=s;

% ---- 15. myoVideo (pressure myograph, entry step) ----------------------------
s = base();
s.id='myoVideo'; s.label='Video'; s.wrapper=@runMyographVideo;
% THE ENTRY STEP OF THE PRESSURE MYOGRAPH, and the only function that creates the
% recording's '_MYO' PAIR; every myograph step after it appends to that one pair in
% place, which is why the product carries no stage flag.  There is no SOURCE member:
% the recording is described in results.recording and re-opened when a frame is
% wanted, never copied.
s.inGlob='*.avi'; s.outSuffix={'_MYO_r','_MYO_s'}; s.outKind='new';
s.outTransform=struct('from','.avi','to','_MYO_r.mat');
s.gatingField='runMyographVideo'; s.requires={}; s.produces={'myograph'};
s.modalities={'PMYO'}; s.branch='myograph';
% 'myograph' rather than '': a branch-agnostic derived step is treated as
% RECORDING-level by wbTypeSelection (ticked once on the anchor row and inherited),
% and these are per-row steps of one pipeline.
s.settingGroups={ 'Calibration',{'pixelSize','rowRange'} };
s.basicFields={'pixelSize','rowRange'};
s.sharedKeys={'rowRange','libraryFolder'};   % the row band is one decision, set once
s.presets=struct('default',struct('pixelSize',[],'rowRange',[1 Inf]));
s.tips=struct('pixelSize','µm per px; leave empty (or 0) to report results in px', ...
    'rowRange','[first last] image row the vessel occupies');
reg(end+1)=s;

% ---- 16. labChart (wire myograph, entry step) --------------------------------
s = base();
s.id='labChart'; s.label='LabChart'; s.wrapper=@runLabChart;
% THE ENTRY STEP OF THE WIRE MYOGRAPH, and the twin of myoVideo: it creates the
% same '_MYO' pair, so the interval editor and the vasomotion step downstream are
% literally the same registry entries for both myographs.  Reading a LabChart file
% IS the step - there is nothing to detect and nothing to measure - so unlike the
% video entry step it fills the product completely (the channels themselves land in
% results.channel) and the '.adicht' is never opened again.
s.inGlob='*.adicht'; s.outSuffix={'_MYO_r','_MYO_s'}; s.outKind='new';
s.outTransform=struct('from','.adicht','to','_MYO_r.mat');
s.gatingField='runLabChart'; s.requires={}; s.produces={'myograph'};
s.modalities={'WMYO'}; s.branch='myograph';
s.settingGroups={ 'What is read',{'records','channels'};
                  'Sampling rate',{'decimate','decimateFS','decimateFilter'} };
s.basicFields={'records','channels','decimate','decimateFS','decimateFilter'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('decimateFilter',{{'mean','median'}});
s.presets=struct('default',struct('records',[],'channels',{{}}, ...
    'decimate',false,'decimateFS',100,'decimateFilter','mean'));
s.tips=struct('records','which LabChart records to read; leave empty for all of them', ...
    'channels','names of the channels to keep; leave empty for all of them', ...
    'decimate',['on: reduce the sampling rate of every channel while the recording ' ...
    'is read.  Everything after this step then works on the reduced recording'], ...
    'decimateFS',['the rate to reduce to.  A whole number of samples goes into each ' ...
    'kept point, so the rate you get is close to this one rather than exactly it, ' ...
    'and a channel recorded at or below it is kept as it was'], ...
    'decimateFilter',['each kept point is the mean or the median of the samples it ' ...
    'replaces.  Median ignores spikes, mean does not']);
s.labels=struct('records','Records to read','channels','Channels to keep', ...
    'decimate','Decimate','decimateFS','Requested frequency, Hz', ...
    'decimateFilter','Averaging');
s.note=['Reading LabChart files needs the ADInstruments SDK in the library''s ' ...
        '"3rd party" folder, and works on Windows only.'];
reg(end+1)=s;

% ---- 17. myoPresetIntervals (before anything is measured) --------------------
s = base();
s.id='myoPresetIntervals'; s.label='Pre-set intervals'; s.wrapper=@setMyographPresetIntervals;
% CHOSEN ON THE VIDEO, BEFORE ANY DIAMETER EXISTS: the diameter step then measures
% only inside these windows and the product is only that long, which is what makes
% measuring 20 minutes of a two-hour recording cost 20 minutes of memory.
s.inGlob='*_MYO_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setMyographPresetIntervals'; s.requires={'myoVideo'};
s.produces={'intervals'};
s.interactive=true; s.needsRaw=true;                  % the recording IS the window
s.conflictsWith={'myoCrop','myoIntervals'};           % one set of windows per recording
s.modalities={'PMYO'}; s.branch='myograph';
s.settingGroups={ 'Profile',{'profileSamples'} };
s.basicFields={'profileSamples'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('profileSamples',1200));
s.tips=struct('profileSamples', ...
    ['how many points the brightness curve has that the windows are picked on - ' ...
     'more is a finer curve and a longer wait before the window opens']);
s.labels=struct('profileSamples','Points in the preview curve');
reg(end+1)=s;

% ---- 18. myoCrop -------------------------------------------------------------
s = base();
s.id='myoCrop'; s.label='Time crop'; s.wrapper=@setMyographCrop;
% ONE window instead of several - the alternative to pre-set intervals.  It records
% the decision; the diameter step is what reads only those frames.  CHANGING IT
% AFTER A DIAMETER EXISTS MEANS MEASURING AGAIN, which the tooltip says out loud.
s.inGlob='*_MYO_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setMyographCrop'; s.requires={'myoVideo'};
s.produces={'timeCrop'};
s.interactive=true; s.needsRaw=true;
s.conflictsWith={'myoPresetIntervals'};
s.note=['Change the crop after the diameter has been measured and the Diameter step ' ...
        'has to run again - the measured arrays cover the old window.'];
s.modalities={'PMYO'}; s.branch='myograph';
s.settingGroups={ 'Profile',{'profileSamples'} };
s.basicFields={'profileSamples'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('profileSamples',1200));
s.tips=struct('profileSamples', ...
    ['how many points the brightness curve has that the window is picked on - ' ...
     'more is a finer curve and a longer wait before the window opens']);
s.labels=struct('profileSamples','Points in the preview curve');
reg(end+1)=s;

% ---- 19. myoDiameter ---------------------------------------------------------
s = base();
s.id='myoDiameter'; s.label='Diameter'; s.wrapper=@runMyographDiameter;
s.inGlob='*_MYO_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runMyographDiameter'; s.requires={'myoVideo'}; s.produces={'diameter'};
s.needsRaw=true;                                      % it reads the video itself
s.modalities={'PMYO'}; s.branch='myograph';
s.settingGroups={ 'Vessel location',{'rowRange','wallContrast','minWallGap','wallProm'};
                  'Image cleaning',{'smoothSigma','dustRadius','dustContrast','brightPrctile'};
                  'Wall measurement',{'edgeMode','subpixel','searchWin'};
                  'Smoothing & outliers',{'tSmoothHz','smoothSpan','tSpan','ySpan','outlierK'};
                  'What is kept',{'keepArrays'} };
s.basicFields={'rowRange','wallContrast','smoothSigma','dustRadius','tSmoothHz','minWallGap'};
s.sharedKeys={'rowRange','libraryFolder'};
s.enums=struct('edgeMode',{{'mid','outer','inner','min'}});
s.presets=struct('default',struct('rowRange',[1 Inf],'wallContrast',0.05,'minWallGap',3, ...
    'wallProm',0.25,'smoothSigma',1.2,'dustRadius',8,'dustContrast',0.06, ...
    'brightPrctile',90,'edgeMode','mid','subpixel',true,'searchWin',[], ...
    'tSmoothHz',1,'smoothSpan',15,'tSpan',25,'ySpan',31,'outlierK',3,'keepArrays',true));
s.tips=struct('rowRange','[first last] image row the vessel occupies', ...
    'wallContrast','how dark a wall must be, relative to the lumen, to count as one', ...
    'minWallGap','the smallest diameter that can be real, px', ...
    'wallProm','below this, one wall has dilated out of view and the frame is flagged', ...
    'smoothSigma','blur applied before detection, px (raise it for noisy recordings)', ...
    'dustRadius','dark specks up to this size are removed before detection (0 = off)', ...
    'tSmoothHz','the fastest the diameter is allowed to change, Hz', ...
    'edgeMode','which of the three diameters is plotted and analysed by default', ...
    'tSpan','window used to spot a diameter that jumps in time, frames', ...
    'ySpan','window used to spot a diameter that jumps along the vessel, rows', ...
    'keepArrays',['on: keep every line''s own diameter and the wall positions ' ...
     'behind it, so the diameter can be looked at position by position and along ' ...
     'the vessel over time.  Off keeps only the averaged trace and makes the ' ...
     'results file far smaller - but the travelling-wave step and the per-line ' ...
     'vasomotion then have nothing to read, and measuring again is the only way ' ...
     'to get them back']);
s.labels=struct('keepArrays','Keep every line''s diameter');
reg(end+1)=s;

% ---- 20. myoIntervals --------------------------------------------------------
s = base();
s.id='myoIntervals'; s.label='Intervals'; s.wrapper=@setMyographIntervals;
% THE WINDOWS THE ANALYSES RUN IN, defined on the diameter that has been measured -
% which is what makes 'baseline', 'drug' and 'washout' three answers from one
% recording instead of one average over all three.
s.inGlob='*_MYO_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setMyographIntervals'; s.requires={};
% ONE STEP, TWO MYOGRAPHS, and requiresAny is what lets it be one: the producer
% differs per modality, and because the Constructor works on a modality-FILTERED
% registry only one of the two survives in any given type - so this collapses to
% whichever entry step that type actually runs.
s.requiresAny={'myoDiameter','labChart'}; s.produces={'intervals'};
s.interactive=true;
s.needsRaw=true;                                      % the pressure myograph's preview
                                                      % plays the recording; the wire
                                                      % myograph gets the path and has
                                                      % no use for it - the channels are
                                                      % already in the product
s.conflictsWith={'myoPresetIntervals'};               % windows are chosen once, not twice
s.modalities={'PMYO','WMYO'}; s.branch='myograph';
s.settingGroups={ 'Signal',{'edgeMode'}; ...
                  'Drawing',{'drawPoints'} };
s.basicFields={'edgeMode'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('edgeMode',{{'mid','outer','inner'}});
s.presets=struct('default',struct('edgeMode','mid','drawPoints',2000));
s.tips=struct('edgeMode',['which of the three measured diameters the trace shows.  ' ...
    'A wire myograph has no diameter: its channel is chosen in the window''s table'], ...
    'drawPoints',['how much of the recording is drawn at once.  A wire myograph ' ...
    'file holds millions of samples per channel and drawing them all makes the ' ...
    'window slow to respond; lower this on a slow machine.  It changes only what ' ...
    'is on screen - the windows are cut out of the full recording either way']);
s.labels=struct('edgeMode','Diameter shown','drawPoints','Points drawn per channel');
s.note=['Moving an interval re-cuts its diameter and clears its propagation and ' ...
        'vasomotion - run those again afterwards.'];
reg(end+1)=s;

% ---- 21. myoPropagation ------------------------------------------------------
s = base();
s.id='myoPropagation'; s.label='Propagation'; s.wrapper=@runMyographPropagation;
s.inGlob='*_MYO_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runMyographPropagation'; s.requires={'myoDiameter'}; s.produces={'propagation'};
s.modalities={'PMYO'}; s.branch='myograph';
s.settingGroups={ 'Signal',{'diameterMeasures','vFR','detrendSec'};
                  'Location quality',{'propMinMeasured','propMinCoh','propArtifactK','propMinRows'};
                  'Delay search',{'propMaxLagFrac','propResolutionSamples'};
                  'Significance',{'propNShuffle'} };
s.basicFields={'diameterMeasures','vFR'};
% diameterMeasures and the band are ONE protocol decision, not one per step: the
% vasomotion step analyses the same diameter over the same oscillation.
s.sharedKeys={'diameterMeasures','vFR','libraryFolder'};
s.presets=struct('default',struct('diameterMeasures',{{'mid'}},'vFR',[0.05 0.25], ...
    'detrendSec',30,'propMinMeasured',0.5,'propMinCoh',0.3,'propArtifactK',4, ...
    'propMinRows',20,'propMaxLagFrac',0.5,'propResolutionSamples',1.0,'propNShuffle',200));
s.tips=struct('diameterMeasures','which diameters to analyse: outer wall, wall centre (mid) or lumen', ...
    'vFR','vasomotion frequency band [lo hi], Hz - it sets how far back a delay is looked for', ...
    'detrendSec','slow drift over this many seconds is removed before comparing locations', ...
    'propMinMeasured','a location is used only if this fraction of it was really measured', ...
    'propMinCoh','how well a location must follow the vessel''s oscillation to be used', ...
    'propMinRows','fewest locations an estimate may rest on', ...
    'propNShuffle','how many times the locations are shuffled to test the result against chance');
reg(end+1)=s;

% ---- 22. myoVasomotion -------------------------------------------------------
s = base();
s.id='myoVasomotion'; s.label='Vasomotion'; s.wrapper=@runMyographVasomotion;
s.inGlob='*_MYO_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runMyographVasomotion'; s.requires={};
% ONE STEP, TWO MYOGRAPHS - see the note on myoIntervals above.  The wrapper reads
% the recording's own results.recording.modality and iterates diameter measures or
% selected CHANNELS accordingly; the analysis behind both is the same core.
s.requiresAny={'myoDiameter','labChart'}; s.produces={'vasomotion'};
s.modalities={'PMYO','WMYO'}; s.branch='myograph';
s.settingGroups={ 'Signal',{'diameterMeasures','perLine'};
                  'Bands',{'vFR','cFR','wFR','wVPO'};
                  'Normalisation',{'normalisation','normsize','tgtFS'};
                  'Peaks & percentiles',{'pcts','otsuMaxN','otsuElbow','nPeakProm'};
                  'Analysis levels',{'segVsmReturn'};
                  'Parallel',{'parforMyographLines'} };
s.basicFields={'diameterMeasures','perLine','vFR','cFR','tgtFS'};
s.sharedKeys={'diameterMeasures','vFR','libraryFolder'};
s.enums=struct('normalisation',{{'mean','median','mmean','mmedian'}});
% the analysis defaults ARE the speckle vasomotion step's, verbatim: the two steps
% call the same core, so a myograph result and an LSCI result are comparable only
% while these stay identical.
s.presets=struct('default',struct('diameterMeasures',{{'mid'}},'perLine',false, ...
    'vFR',[0.05 0.25],'cFR',[0.4 0.6],'wFR',[0.01 1],'wVPO',10, ...
    'normalisation','median','normsize',101,'tgtFS',1,'pcts',0:10:100,'otsuMaxN',5, ...
    'otsuElbow',0.05,'nPeakProm',0.10, ...
    'segVsmReturn',{{'bands','moments','series','clustering','spectrum'}}, ...
    'parforMyographLines',true));
s.labels=struct('parforMyographLines','Use parfor: image lines');
s.tips=struct('diameterMeasures',['which diameters to analyse: outer wall, wall centre (mid) ' ...
    'or lumen.  A wire myograph analyses the channels chosen for each interval instead'], ...
    'perLine',['on: analyse every measured image row separately.  off: the vessel''s ' ...
    'averaged diameter.  A wire myograph ignores it - a channel is one trace'], ...
    'vFR','vasomotion frequency band [lo hi], Hz', ...
    'cFR','control (cardiac) frequency band [lo hi], Hz', ...
    'nPeakProm','VB peak-count prominence as a fraction of the band range', ...
    'segVsmReturn','which analysis levels to store (set of tokens)', ...
    'parforMyographLines',parforTip('the per-line analysis'));
reg(end+1)=s;

% ---- 23. intensity (fluorescence angiogram, entry step) ----------------------
s = base();
s.id='intensity'; s.label='Angiogram'; s.wrapper=@runIntensity;
% THE ENTRY STEP OF THE FLUORESCENCE BRANCH, and a RAW PRODUCER: it reads the .cxd and
% writes a new, independent '_a_I' triplet, the way contrast writes '_t_K'.  It is one
% of three entry steps a .cxd can have - the angiogram here, the cardiac cycle and the
% bolus still to come - and which of them a recording is for is the protocol's answer.
% The step reads the recording ONE PLANE AT A TIME and never holds more than the output
% cube, which is what makes a 60 GB recording processable at all.
s.inGlob='*.cxd'; s.outSuffix={'_a_I_d','_a_I_r','_a_I_s'}; s.outKind='new';
s.outTransform=struct('from','.cxd','to','_a_I_r.mat');
s.gatingField='runIntensity'; s.requires={}; s.produces={'angiogram'};
s.artifacts={'_rep_intensity.jpg'};        % no legacy tail: nothing on this branch
                                           % predates the unified reporting rename
s.modalities={'EPFL'}; s.branch='angiogram';
s.settingGroups={ 'Frame averaging',{'decimFactor'};
                  'What is kept',{'dataTypeOut','saveSource'} };
s.basicFields={'decimFactor','dataTypeOut','saveSource'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('decimFactor',1,'dataTypeOut','single','saveSource',true));
s.labels=struct('decimFactor','Frames averaged into one', ...
                'dataTypeOut','Frames stored as', ...
                'saveSource','Keep the frames');
s.tips=struct('decimFactor','output framerate = original / decimFactor', ...
    'dataTypeOut',['''single'' keeps an averaged frame exactly; leave it empty to ' ...
       'store the frames the way the camera recorded them, which is half the size'], ...
    'saveSource',['on: keep every frame beside the angiogram, which is what the flow, ' ...
       'motion and per-pixel analyses read.  Off keeps only the picture of the ' ...
       'vasculature and the intensity over time - far smaller, and it is the only way ' ...
       'to get an angiogram out of a recording too long to hold at once, but nothing ' ...
       'that reads the frames can then be run on it and reading the recording again ' ...
       'is the only way back']);
reg(end+1)=s;

% ---- 24. intensityCycle (fluorescence cardiac cycle, entry step) -------------
s = base();
s.id='intensityCycle'; s.label='Cardiac cycle'; s.wrapper=@runIntensityInternalCycle;
% THE SECOND ENTRY STEP OF THE FLUORESCENCE BRANCH, and the twin of 'internalCycle' on
% the other modality: it reads the .cxd and writes a new, independent '_c_I' triplet the
% way that step writes '_c_K'.  BOTH DECLARE branch 'cardiac', deliberately - a recording
% is a .cxd OR an .rls and never both, so no recording ever carries two products on that
% branch, and the modality filter is what keeps the two steps apart.  It streams: pass 1
% finds the beats from a spatially decimated reference cube, pass 2 re-reads the accepted
% ones through the reader that is already open, and the recording is never resident.
s.inGlob='*.cxd'; s.outSuffix={'_c_I_d','_c_I_r','_c_I_s'}; s.outKind='new';
s.outTransform=struct('from','.cxd','to','_c_I_r.mat');
s.gatingField='runIntensityInternalCycle'; s.requires={}; s.produces={'cardiac'};
s.artifacts={'_rep_cardiac-detect.jpg','_rep_cardiac-average.jpg'};
% no legacy tail: nothing on this branch predates the unified reporting rename.  The
% tails are the step's OWN and not the contrast cycle's ('_rep_cycle-*'), because
% wbArtifacts resolves against the RECORDING's base name and the stage token in the tail
% is the only thing that attributes an image to a column.
s.modalities={'EPFL'}; s.branch='cardiac';
s.settingGroups={ 'Pulse detection area',{'maskPctI'};
                  'Frequency band',{'maxFrqIni','minFrqIni','rangeFrq'};
                  'Exclusion criteria',{'excludeFirstNCycles','coeffsSTD','coeffsRel', ...
                     'coeffsAbs','gateDetrend','minPromCoef'};
                  'Cycle calculation',{'nPhaseBins','nControls','controlSeed', ...
                     'decimationSpace'} };
s.basicFields={'maskPctI','maxFrqIni','minFrqIni','nPhaseBins','nControls', ...
               'excludeFirstNCycles'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('maskPctI',85,'maxFrqIni',20,'minFrqIni',1, ...
    'rangeFrq',0.3,'excludeFirstNCycles',0,'coeffsSTD',[3 2 2 2 2 3 3 2 2], ...
    'coeffsRel',[0.5 0.1],'coeffsAbs',2,'gateDetrend',3,'minPromCoef',1/4, ...
    'nPhaseBins',25,'nControls',3,'controlSeed',20260807,'decimationSpace',4));
s.labels=struct('maskPctI','Brightest pixels used, %', ...
                'nPhaseBins','Points in the averaged cycle', ...
                'nControls','Comparison beats stored', ...
                'controlSeed','Comparison random seed', ...
                'gateDetrend','Bleaching removed as');
s.tips=struct('maskPctI', ...
       ['the heartbeat is looked for in the brightest pixels, because the dye is in ' ...
        'the plasma.  85 keeps the brightest 15 % of the picture, which measured best ' ...
        'on the reference recordings - taking fewer than that loses more to noise ' ...
        'than it gains'], ...
    'maxFrqIni','initial max frequency of interest, Hz', ...
    'minFrqIni','initial min frequency of interest, Hz', ...
    'rangeFrq',['how far a single beat''s own rate may sit from the recording''s, as a ' ...
        'fraction of it, before that beat is dropped'], ...
    'coeffsSTD','beat-rejection coefficients relative to feature std', ...
    'coeffsRel','how uneven a beat''s two ends, and its level, may be', ...
    'coeffsAbs','how far either side of a dropped beat its neighbours are dropped too', ...
    'excludeFirstNCycles','drop this many beats at the start of the recording', ...
    'gateDetrend',['the dye fades over a recording, and that fading is removed as a ' ...
        'polynomial of this order before the heart rate is looked for'], ...
    'minPromCoef',['how far a beat has to dip below its neighbours to count as one, ' ...
        'as a fraction of the heartbeat''s own size'], ...
    'nPhaseBins',['how many points the averaged beat is described by.  About ten of ' ...
        'them are measured at these frame rates and the rest place those measurements ' ...
        'on one common axis so beats of different length can be averaged together'], ...
    'nControls',['averaged beats built from the same heartbeats with their timing ' ...
        'scrambled, and stored beside the real one.  They are what says how much of an ' ...
        'averaged beat is the heart and how much is the recording, and Wall motion ' ...
        'cannot run without them.  Three lets it take a middle value rather than a ' ...
        'single draw; each one adds as much again to the file'], ...
    'controlSeed',['fixes which scrambling is used, so re-running a recording gives the ' ...
        'same comparison beats.  Change it only to see whether an answer depended on ' ...
        'the draw'], ...
    'decimationSpace',['pixels grouped together while the heartbeat is being found.  ' ...
        'It costs memory only, and the averaged beat is always at full resolution']);
s.note=['The averaged beat is stored in the recording''s own intensity, so how much a ' ...
        'vessel pulses can be compared with a speckle recording directly.  It is ' ...
        'stored with three scrambled-timing copies of itself, which is what Wall ' ...
        'motion measures against.'];
reg(end+1)=s;

% ---- 25-27. the shape-finding twins (regions, segmentation, dynamic segmentation) --
% THE FIRST THREE TWINS THERE HAVE EVER BEEN, and they are built through
% intensityTwinOf rather than written out, so the derivation is exercised by the
% shipped registry and not only by a test.  Each one then sets the three fields the
% helper deliberately CLEARS - requires / requiresAny / conflictsWith name STEP IDS,
% and a twin's producers are different steps.
%
% requiresAny is the THREE fluorescence entry steps with 'intensity' FIRST, because the
% first entry is the default producer a one-click chain pulls in and the angiogram is
% what a shape-finding step normally runs on.  'intensityBolus' is appended LAST for
% exactly that reason - a bolus product needs regions and segmentation like any other,
% but nobody's one-click chain should start by asking for an injection.
%
% THE CHAIN IS THE SPECKLE SIDE'S, EXACTLY.  Regions was briefly a HARD prerequisite
% of Segmentation on this branch, because Regions was where the operator answered
% which way round the vessels were and getVesselPolarity refused to guess.  That
% question is gone (author, 2026-08-08: fluorescence recordings processed here have
% bright vessels), so the reason for the asymmetry is gone with it - drawing nothing
% is a complete answer on '_I' as it always was on '_K', and Segmentation asks for
% nothing Regions writes.
t = intensityTwinOf(stepById(reg,'setRegions'));
t.requiresAny = {'intensity','intensityCycle','intensityBolus'};
% branchScope 'copy' is carried, and it means MORE here than it does on the speckle
% side: one '.cxd' can carry three co-registered products ('_a_I', '_c_I', '_b_I')
% against the speckle side's two, and all three come from one recording through one
% objective - so the regions are drawn once and inherited by every one of them.
% wbExecutor derives the target list from the recording's own siblings, so the extra
% branch costs no edit anywhere.
reg(end+1) = t;

t = intensityTwinOf(stepById(reg,'segmentation'));
t.requiresAny = {'intensity','intensityCycle','intensityBolus'};
% 'single', where the speckle step says 'multiscale' - which is what each branch was
% already doing when the schedule was chosen by the file's product token instead of by
% a setting.  Nobody has measured which suits a fluorescence mean image; this is the
% preset that makes that measurable rather than the one that claims to know.
t.presets.default.diffusionSchedule = 'single';
reg(end+1) = t;

t = intensityTwinOf(stepById(reg,'dynamicSegmentation'));
t.requires = {'segmentationI'};
% NO requiresAny, exactly as on the speckle side: the chain reaches an entry step
% through segmentationI and the tick cascade is transitive, so naming the producers
% again here would only be a second place to edit when the bolus entry step lands.
reg(end+1) = t;

% ---- 28. backgroundRemoval (fluorescence) ------------------------------------
s = base();
s.id='backgroundRemoval'; s.label='Background removal'; s.wrapper=@runBackgroundRemoval;
% WRITTEN OUT BY HAND rather than derived through intensityTwinOf, and that is not an
% oversight: it has no speckle counterpart to be a twin OF.  A contrast image has no
% halo to take off it - the haze this removes is the dye's own light, spread sideways
% by the tissue, and there is no dye on the speckle side.
%
% IT CLEANS results.imgI AND, ONLY IF ASKED, THE FRAMES [D15].  That split is the
% whole design: segmentation and topology read the picture, CTTH and the per-pixel
% passes read the cube, and they are different members of the same product - so
% cleaning the first cannot change the second, and the safety property is a property
% of the data flow rather than a rule in a document.
%
% applyToSource IS NOT INTERLOCKED AGAINST CTTH [D16], by author's ruling of
% 2026-08-07.  It costs a bolus transit time about two seconds of arrival and
% scatters the first-pass area from one vessel to the next, and the result looks
% completely ordinary - so what ships is the record (the settings sidecar says which
% was chosen) and one warning from runCTTH, NOT a conflictsWith and not a greyed
% cell.  Do not add a registry relation for it.
s.inGlob='*_I_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runBackgroundRemoval';
s.requires={}; s.requiresAny={'intensity','intensityCycle','intensityBolus'};
s.produces={'background'};
s.artifacts={'_rep_background.jpg'};
s.modalities={'EPFL'};
s.branch=''; s.branchScope='all';       % offered on every intensity row, and every
                                        % branch product of a recording has its own
                                        % picture to clean - nothing is inherited
s.settingGroups={ 'Calibration',{'pixelSize'};
                  'Glow',{'radiusUm','medianUm'};
                  'Where to apply it',{'applyToSource'} };
s.basicFields={'pixelSize','radiusUm'};
s.sharedKeys={'pixelSize','libraryFolder'};
% pixelSize HAS NO DEFAULT VALUE, deliberately.  2.5 um/px is measured for one rig and
% every length in this step scales with it, so a number here would be a silent wrong
% answer on any other magnification; empty is refused by name instead
% (runBackgroundRemoval:pixelSizeNotSet).  It is in a settings group because a
% required setting no panel can reach is a step nobody can run, and it is shared
% because it is a fact about the microscope rather than about this step.
s.presets=struct('default',struct('pixelSize',[],'radiusUm',200,'medianUm',37.5, ...
    'applyToSource',false));
s.labels=struct('pixelSize','Micrometres per pixel', ...
                'applyToSource','Also clean the recording itself');
s.tips=struct( ...
    'pixelSize','how many micrometres across one pixel is on this microscope', ...
    'radiusUm','the width of the glow around the vessels; set it above the widest vessel you want to keep', ...
    'medianUm','smooths speckle before the glow is measured; keep it below the narrowest vessel you care about', ...
    'applyToSource','clean every frame too, not just the picture. Needed for per-frame diameter and wall motion, and it costs bolus transit times seconds of accuracy - so leave it off if you want timings from this recording');
reg(end+1)=s;

% ---- 29. intensityBolus (fluorescence bolus, entry step) ---------------------
s = base();
s.id='intensityBolus'; s.label='Bolus'; s.wrapper=@runIntensityBolus;
% THE THIRD AND LAST ENTRY STEP OF THE FLUORESCENCE BRANCH.  It splits one '.cxd' into
% a full-rate BOLUS cube and an averaged angiogram and writes a new '_b_I' triplet, the
% way 'intensity' writes '_a_I' and 'intensityCycle' writes '_c_I'.  Which of the three
% a recording is for is the protocol's answer: a 'BB' in the recording code names a
% bolus recording (author, 07-Aug-2026), which is what the launcher globs on.
%
% IT IS INTERACTIVE WHEN THE SPANS ARE EMPTY, and they ship empty.  Marking the bolus on
% the streamed trace is the normal way to use it - the injection is not at the same
% frame in two recordings - so this declares itself interactive rather than pretending a
% batch of them runs unattended.  Filling both spans in the panel makes it silent.
s.inGlob='*.cxd'; s.outSuffix={'_b_I_d','_b_I_r','_b_I_s'}; s.outKind='new';
s.outTransform=struct('from','.cxd','to','_b_I_r.mat');
s.gatingField='runIntensityBolus'; s.requires={}; s.produces={'bolus'};
s.branch='bolus'; s.modalities={'EPFL'};
s.artifacts={'_rep_bolus.jpg'};    % no legacy tail: nothing on this branch predates
                                   % the unified reporting rename
s.interactive=true;
s.settingGroups={ 'Which frames',{'fBolus','fAngio'};
                  'Calibration',{'pixelSize'} };
s.basicFields={'fBolus','fAngio','pixelSize'};
% pixelSize IS SHARED, and this is the step that answers it.  Background removal works
% entirely in micrometres and refuses without it by name; the entry step is the only one
% that sees the recording, so the number typed HERE is the number that reaches every
% micrometre step of the session.  It has no default for the same reason it has none
% there: 2.5 um/px is measured for one rig and a guess is a silent wrong answer.
s.sharedKeys={'pixelSize','libraryFolder'};
s.presets=struct('default',struct('fBolus',[],'fAngio',[],'pixelSize',[]));
s.labels=struct('fBolus','Frames kept at full speed', ...
                'fAngio','Frames averaged into the picture', ...
                'pixelSize','Micrometres per pixel');
s.tips=struct( ...
    'fBolus',['the stretch around the injection, kept frame by frame.  Leave it empty ' ...
       'to mark it on the recording.  Keep about 25 to 30 seconds of it and start at ' ...
       'least a second and a half before the dye arrives - a shorter stretch still ' ...
       'gives the transit time but not the spread of transit times, and one that ' ...
       'starts too late is refused'], ...
    'fAngio',['the stretch averaged into the picture the vessels are found on.  Leave ' ...
       'it empty to use everything after the injection, which is when the vessels are ' ...
       'brightest'], ...
    'pixelSize','how many micrometres across one pixel is on this microscope');
s.note=['The picture is averaged over a different stretch of the recording from the ' ...
        'one kept frame by frame, so the vessels can be found on a filled field and ' ...
        'still be measured from before the dye arrived.'];
reg(end+1)=s;

% ---- 30. ctth (transit time on a bolus product) ------------------------------
s = base();
s.id='ctth'; s.label='Transit time'; s.wrapper=@runCTTH;
% THE LABEL IS FOR A BIOLOGIST.  'CTTH' is the literature's name for one of the numbers
% this produces and not for what the step does, and a registry label is what a user
% picks from - so the column says what is measured and the header of the wrapper says
% what the acronym means.
%
% requires WITHOUT requiresAny, and the precedent is two entries above.
% 'dynamicSegmentationI' declares exactly this shape with the comment "the chain reaches
% an entry step": naming 'segmentationI' is enough, because segmentationI names the
% three entry steps itself and the tick cascade is transitive.
%
% NOTHING IS INTERLOCKED AGAINST BACKGROUND REMOVAL [D16].  Cleaning the frames costs a
% transit time about two seconds of arrival, and the author ruled on 2026-08-07 that the
% combination is nevertheless the user's to make - so what ships is the record (the
% settings sidecar says which was chosen) and one warning from the wrapper, NOT a
% conflictsWith on either step.  Do not add one.
s.inGlob='*_b_I_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runCTTH'; s.requires={'segmentationI'};
s.branch='bolus'; s.modalities={'EPFL'};
s.artifacts={'_rep_ctth-input.jpg','_rep_ctth-markers.jpg', ...
             '_rep_ctth-confidence.jpg','_rep_ctth-maps.jpg'};
s.settingGroups={ 'How a curve is read',{'ctthLevelPcts','ctthPlateauSec', ...
                     'ctthGuardSec','ctthSlopeSec'};
                  'The arterial curve',{'ctthInputAmpPct','ctthInputPct'};
                  'How much to trust a vessel',{'ctthConfThreshold', ...
                     'ctthConfMinThreshold','ctthSettleFrac','ctthInputScale', ...
                     'ctthCthTol'};
                  'What to keep',{'segCtthReturn','ppxCtthReturn','parforCTTHPixels'} };
s.basicFields={'ctthLevelPcts','ctthConfThreshold','ppxCtthReturn'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('ctthLevelPcts',[10 25 50 75 90], ...
    'ctthPlateauSec',1,'ctthGuardSec',0.5,'ctthSlopeSec',0.1, ...
    'ctthInputAmpPct',75,'ctthInputPct',5, ...
    'ctthConfThreshold',0.6,'ctthConfMinThreshold',0.2,'ctthSettleFrac',1, ...
    'ctthInputScale',0.5,'ctthCthTol',0.25, ...
    'segCtthReturn',{{'levels','amplitudes','times','transit','shape'}}, ...
    'ppxCtthReturn',{{}},'parforCTTHPixels',true));
s.labels=struct('ctthLevelPcts','Filled to, %', ...
                'ctthPlateauSec','Final level averaged over, s', ...
                'ctthGuardSec','Quiet time before the arrival, s', ...
                'ctthSlopeSec','Steepest slope measured over, s', ...
                'ctthInputAmpPct','A vessel counts if it is this bright, %', ...
                'ctthInputPct','Earliest-filling vessels used, %', ...
                'ctthConfThreshold','All the checks together', ...
                'ctthConfMinThreshold','The single worst check', ...
                'ctthSettleFrac','Still rising at the end, at most', ...
                'ctthInputScale','May fill ahead of the artery by, s', ...
                'ctthCthTol','Spread check has to come within', ...
                'segCtthReturn','Which numbers to keep', ...
                'ppxCtthReturn','Which to picture per pixel', ...
                'parforCTTHPixels','Use parfor: the per-pixel pictures');
s.tips=struct( ...
    'ctthLevelPcts','how far the tracer has filled, in per cent, at each time that gets reported', ...
    'ctthPlateauSec','how much of the end of the recording is averaged to read the final level', ...
    'ctthGuardSec','how much quiet time to leave between the pre-injection period and the arrival', ...
    'ctthSlopeSec','how long a stretch the steepest rise and fall are measured over', ...
    'ctthInputAmpPct',['how bright a vessel has to be to count towards the arterial ' ...
       'curve.  It is not optional: without it the earliest-filling vessels are the ' ...
       'dimmest ones, which cross their own threshold on noise'], ...
    'ctthInputPct','how many of the earliest-filling vessels make up the arterial curve', ...
    'ctthConfThreshold','how good all the checks together have to be', ...
    'ctthConfMinThreshold','how good the single worst check has to be', ...
    'ctthSettleFrac','how much the curve may still be rising at the end and still be used', ...
    'ctthInputScale','how far ahead of the artery a vessel may fill before it is not trusted', ...
    'ctthCthTol',['how close the spread check has to come before the recording is ' ...
       'judged long enough for it.  Below that the spread is left blank rather than ' ...
       'reported small'], ...
    'segCtthReturn','which groups of numbers to keep', ...
    'ppxCtthReturn',['which of them to also make a picture of, one value per pixel.  ' ...
       'Leave it empty for the vessels alone, which is what a table is read from'], ...
    'parforCTTHPixels',parforTip('the per-pixel pictures'));
s.note=['Every time is measured against the recording''s own arterial curve, which the ' ...
        'step works out from the vessels that fill first - so the delay, the mean ' ...
        'transit time and the spread compare between recordings and the plain times do ' ...
        'not.  The spread is left blank wherever the recording was too short to ' ...
        'resolve it.'];
reg(end+1)=s;

% ---- 31. topology (vascular density on any intensity product) ----------------
s = base();
s.id='topology'; s.label='Vascular density'; s.wrapper=@runTopologyAnalysis;
% THE LABEL IS 'Vascular density', NOT 'Topology', for two reasons.  It is what a
% biologist calls the thing they are asking for, and 'topology' in this library already
% means the flow DAG setVascularTree builds out of pulse-arrival timing - a second
% meaning for the word in the step column would be the kind of collision nobody notices
% until two steps are on screen together.
%
% HAND-WRITTEN AND NOT AN intensityTwin, like backgroundRemoval two entries above: there
% is no speckle counterpart to be a twin OF, and there should not be one.  The core is
% modality-blind and WOULD run on a '_K' product - which is exactly the danger.  There
% the "vessel mask" comes out of a FLOW image, so a vessel with no flow in it is
% invisible in contrast and would be counted as tissue; a vascular density that silently
% drops occluded vessels is worse than no vascular density.  If it is ever wanted it
% needs its own justification, not a one-line twin.
%
% requires WITHOUT requiresAny, the shape 'dynamicSegmentationI' and 'ctth' both use: the
% chain reaches an entry step through segmentationI, which names the three of them
% itself, and the tick cascade is transitive.  Naming a step id that does not exist would
% fail SILENTLY here - validateRegistry checks only conflictsWith, pruneReferences
% strikes the unknown id, and an empty requires reads as an ENTRY step - so the id below
% is one S3 actually landed.
s.inGlob='*_I_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runTopologyAnalysis';
s.requires={'segmentationI'}; s.produces={'topology'};
s.artifacts={'_rep_topology.jpg'};
s.modalities={'EPFL'};
% branch-agnostic and offered on EVERY intensity row - '_a_I', '_c_I' and '_b_I' alike.
% The step never reads results.imgI's provenance, only the mask runSegmentation left, so
% it is one step further removed from the entry product than the segmentation is.  The
% bolus angiogram is not second-best: it is an average over thousands of post-injection
% frames and it segments at least as well.
s.branch=''; s.branchScope='all';
s.settingGroups={ 'Labelling (match segmentation)',{'sMinL','prchNSize','correctNodes','simR','difR'};
                  'Calibration',{'pixelSize'};
                  'What is reported',{'evdBeyond','calibreEdges','evdEdges','keepMaps'};
                  'Parallel',{'parforSegmentationLabels'} };
s.basicFields={'pixelSize','evdBeyond'};
% pixelSize is SHARED so one calibration serves the whole working set - it is a fact
% about the microscope, and the entry step is where it is answered.
% parforSegmentationLabels is shared with both segmentation steps for the reason given
% there: all three drive getSegmentationLabels, so the machine-level choice is one tick.
s.sharedKeys={'sMinL','prchNSize','correctNodes','simR','difR','pixelSize', ...
    'parforSegmentationLabels','libraryFolder'};
% THE HISTOGRAM RANGES AND evdBeyond ARE MEASURED, NOT GUESSED, and they are narrower
% than the research prototype's on purpose: at 0:10:200 um and 0:20:600 um both figures
% came out as a single bar with most of the axis empty, because on real fields all the
% calibre mass is below 60 um and all the tissue-distance mass is below 60 um.  Do not
% widen them back.  evdBeyond is PROVISIONAL at 25 um - at 50 um the fraction of tissue
% beyond it is 5e-5 and 1e-5 on the two example fields, i.e. identically zero.
s.presets=struct('default',struct('sMinL',15,'prchNSize',50,'correctNodes',true, ...
    'simR',0.3,'difR',0.4,'pixelSize',[],'evdBeyond',25, ...
    'calibreEdges',0:5:100,'evdEdges',0:2:60,'keepMaps',false, ...
    'parforSegmentationLabels',true));
s.labels=struct('pixelSize','Micrometres per pixel', ...
                'evdBeyond','Tissue further than', ...
                'calibreEdges','Vessel width bins', ...
                'evdEdges','Distance bins', ...
                'keepMaps','Keep the distance picture', ...
                'parforSegmentationLabels','Use parfor: segment labels');
s.tips=struct( ...
    'sMinL','minimum segment length', ...
    'prchNSize','parenchymal pixel neighbourhood', ...
    'pixelSize','leave empty to report distances and densities in pixels', ...
    'evdBeyond','the share of tissue further than this from the nearest vessel is reported', ...
    'calibreEdges','edges of the vessel-width histogram', ...
    'evdEdges','edges of the tissue-distance histogram', ...
    'keepMaps','on: also store the picture of how far each point is from the nearest vessel', ...
    'parforSegmentationLabels',parforTip('the per-segment label growing'));
% THE THREE SENTENCES THAT SAY WHAT THE NUMBERS ARE NOT.  They are on the report page as
% well, because a settings tooltip is not where the reader of a figure looks.
s.note=['Junction density counts where centre-lines meet. Two vessels crossing at ' ...
    'different depths look the same as one vessel branching, so compare fields, not ' ...
    'animals.  Extravascular distance is how far tissue is from the nearest vessel it ' ...
    'can see. Vessels below the surface are counted as if they were in the plane, and ' ...
    'the smallest ones are not resolved at all, so this is not an oxygen supply ' ...
    'distance.  Length density is total vessel length per unit area. Vessels that ' ...
    'overlap in the picture merge into one line, so this falls short where the ' ...
    'vasculature is densest.'];
reg(end+1)=s;

% ---- 32. motionEnhancement (wall motion on any intensity product) ------------
s = base();
s.id='motionEnhancement'; s.label='Wall motion'; s.wrapper=@runMotionEnhancement;
% HAND-WRITTEN AND NOT AN intensityTwin, the third entry in this family that is: there
% is no 'motionEnhancement' on the speckle side to be a twin OF.  What exists there is
% runDynamicSegmentation's per-frame diameter, and that is a DIFFERENT instrument, not
% an earlier version of this one - it walks integer indices on an interpolated profile,
% so its trace is quantised at a quarter of a pixel while the cardiac wall movement on
% the reference recording is 1.6 to 18 times SMALLER than one such step.  Both should
% exist; neither should be compared with the other.
%
% requires WITHOUT requiresAny, the shape 'dynamicSegmentationI', 'ctth' and 'topology'
% all use: the chain reaches an entry step through segmentationI, which names the three
% of them itself, and the tick cascade is transitive.  THE ID BELOW IS THE EPFL ONE.
% Naming 'segmentation' instead would fail SILENTLY - validateRegistry checks only
% conflictsWith, the modality filter's pruneReferences strikes the unknown id, and an
% empty requires reads as an ENTRY step whose availability test is vacuously true - so
% the step would be offered on every intensity product including never-segmented ones,
% and the wrapper would then throw on a missing results.cMask.
s.inGlob='*_I_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runMotionEnhancement';
s.requires={'segmentationI'}; s.produces={'wallMotion'};
s.artifacts={'_rep_wall-motion.jpg'};
s.modalities={'EPFL'};
% branch-agnostic and offered on EVERY intensity row, because it does something
% DIFFERENT on each rather than the same thing three times: on '_c_I' it measures the
% fundamental of one averaged beat against the scrambled-timing copies stored beside
% it, and on '_a_I' / '_b_I' it measures the same cuts inside a band against an empty
% band of the same width.  That is not a 'copy' and it is not a 'one'.
s.branch=''; s.branchScope='all';
s.settingGroups={ 'Calibration',{'pixelSize'};
                  'Labelling (match segmentation)',{'sMinL','prchNSize','correctNodes','simR','difR'};
                  'What is allowed to carry a number',{'minDiameterUm','confMin','minCohere','minCuts'};
                  'Where the cuts go',{'cutSpacing','endMarginRadii','spanRadii', ...
                     'spanPadUm','smoothRadii','smoothMin','interpF','maxCLR'};
                  'How a cut is measured',{'cutStepPx','widthTol','edgeTol'};
                  'Continuous recordings',{'bandHz'};
                  'Movie',{'writeVideo','videoFolder','alpha','levels'};
                  'Parallel',{'parforSegmentationLabels'} };
s.basicFields={'pixelSize','minDiameterUm','confMin','writeVideo'};
% pixelSize is SHARED so one calibration serves the whole working set, and the labelling
% keys are shared with both segmentation steps and with Vascular density for the reason
% given there: all of them drive getSegmentationLabels and must agree.
s.sharedKeys={'pixelSize','sMinL','prchNSize','correctNodes','simR','difR', ...
    'parforSegmentationLabels','libraryFolder'};
% pixelSize and bandHz ship EMPTY on purpose.  The step refuses both by name rather
% than guessing: 2.5 um/px is measured for one rig, and a heart rate is a fact about a
% preparation.  alpha is the one number here that is quoted from a single recording and
% survives as a default anyway, because it decides how a MOVIE looks and not what any
% column says - and the calibre it is valid for is burnt on every frame.
s.presets=struct('default',struct('pixelSize',[],'minDiameterUm',37,'confMin',2, ...
    'minCohere',0,'minCuts',3,'cutSpacing',3,'endMarginRadii',1.5,'spanRadii',3, ...
    'spanPadUm',30,'smoothRadii',1,'smoothMin',3,'interpF',4,'maxCLR',Inf, ...
    'cutStepPx',0.25,'widthTol',[0.5 2],'edgeTol',1.0,'bandHz',[], ...
    'writeVideo',false,'videoFolder','','alpha',[53.1 106.2 212.5],'levels',3:5, ...
    'sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3,'difR',0.4, ...
    'parforSegmentationLabels',true));
s.labels=struct('pixelSize','Micrometres per pixel', ...
                'minDiameterUm','Smallest vessel to attempt', ...
                'confMin','Times its own comparison', ...
                'minCohere','Least agreement between cuts', ...
                'minCuts','Fewest cuts per vessel', ...
                'cutSpacing','Spacing between cuts', ...
                'spanPadUm','Tissue included beyond the wall', ...
                'bandHz','Heartbeat band, beats per second', ...
                'writeVideo','Also write a movie', ...
                'videoFolder','Movie folder', ...
                'alpha','Movie amplification', ...
                'levels','Movie detail scales', ...
                'parforSegmentationLabels','Use parfor: segment labels');
s.tips=struct( ...
    'pixelSize',['there is no default: every limit in this step is a real distance, ' ...
        'so it cannot run until it knows how big a pixel is'], ...
    'minDiameterUm',['how small a vessel to attempt, in micrometres.  Narrower ones ' ...
        'report nothing at all.  37 is where the smallest vessel with a movement ' ...
        'above its own comparison was found on the reference recording; it is a limit ' ...
        'of that microscope and that preparation, measured once, and a finer pixel ' ...
        'does not lower it'], ...
    'confMin',['how many times larger than its own comparison a vessel''s movement ' ...
        'has to be before it is reported.  Below this the vessel reports nothing ' ...
        'rather than a small number'], ...
    'minCohere',['how much the cuts along one vessel have to agree, from 0 to 1.  ' ...
        'Zero lets the other two rules decide, which is what they already do'], ...
    'minCuts',['a vessel measured across fewer places than this is one place wearing ' ...
        'a vessel''s name'], ...
    'cutSpacing',['how far apart the places a vessel is measured at are, in pixels.  ' ...
        'Neighbouring places see almost the same thing, so more of them buy very little'], ...
    'endMarginRadii',['how far from each end of a vessel to stop measuring, in ' ...
        'multiples of its own width.  The ends sit against branches, and a measurement ' ...
        'there crosses the vessel that branches off'], ...
    'spanRadii','how far across the vessel each measurement reaches, in vessel widths', ...
    'spanPadUm',['how much tissue beyond the wall each measurement includes.  It sets ' ...
        'the level the wall is found at, so shortening it moves every answer'], ...
    'smoothRadii','how much the vessel''s path is smoothed, in multiples of its width', ...
    'smoothMin','the smallest that smoothing may become, in pixels', ...
    'interpF','points per pixel along the vessel''s path', ...
    'maxCLR','refuse a vessel more curved than this.  Blank lets every shape through', ...
    'cutStepPx','how finely the brightness is sampled across the vessel, in pixels', ...
    'widthTol',['how far the measured width may sit from the width the outlines ' ...
        'promised, as a range of ratios, before that place is dropped'], ...
    'edgeTol',['how blurred the wall may be, as a share of the vessel''s own width, ' ...
        'before the place is treated as not crossing a vessel at all'], ...
    'bandHz',['on a continuous recording there is no averaged beat, so this step has ' ...
        'to be told which rhythm to look for.  The Cardiac cycle step reports the ' ...
        'heart rate it found on the same animal'], ...
    'writeVideo',['also write a magnified movie and, beside it, the same movie of the ' ...
        'scrambled-timing comparison.  The two are written together or not at all, ' ...
        'because a magnified movie looks convincing whether or not anything moved.  ' ...
        'They are a few hundred megabytes each'], ...
    'videoFolder',['where the movies go.  Blank puts them in a folder beside the ' ...
        'results, never inside it.  Do not choose a folder that is being synced to ' ...
        'the cloud'], ...
    'alpha',['how much the movie exaggerates the movement, one value per detail ' ...
        'scale.  One setting suits one vessel width; the width it is right for is ' ...
        'written on every frame'], ...
    'levels','which detail scales the movie exaggerates', ...
    'parforSegmentationLabels',parforTip('the per-segment label growing'));
% THE THREE SENTENCES THAT SAY WHAT THE NUMBERS ARE NOT.  They are on the report page
% as well, because a settings tooltip is not where the reader of a figure looks.
s.note=['Every vessel is compared only with itself: the same measurement on the same ' ...
    'vessel with the heartbeat''s timing scrambled away.  A vessel that does not ' ...
    'clear that reports nothing, not a small number, and on the reference recording ' ...
    'that is most of them.  How far the whole vessel slid is reported beside how far ' ...
    'its walls moved, because the two are the same size here and a width change ' ...
    'quoted without it can be the preparation moving.  The movie is an illustration ' ...
    'and is never the measurement.'];
reg(end+1)=s;

% ---- 33. vasomotionI (the fluorescence twin of the vasomotion step) ----------
% A GENUINE TWIN, and the reason it twins while its neighbour below does not is a fact
% about the two wrappers rather than a preference.  runVasomotion names the physical
% quantity in a TREE NODE - results.vasomotion.sData / .dvsData / .dvsDiameter /
% .gsData, chosen by s.vsmSignals - so a trace is a trace and one wrapper serves both
% modalities.  runPulsatility names it in the COLUMN PREFIX, which is why 'pulsatilityI'
% below points at a different wrapper.
%
% NOTHING IS DROPPED, and that is the finding rather than an omission: every setting
% this step has describes the WAVELET TRANSFORM - the bands, the centring, the
% percentiles, the levels - and none of them describes what the trace measures.  A twin
% that drops nothing is what a quantity-agnostic wrapper looks like from here.
t = intensityTwinOf(stepById(reg,'vasomotion'));
% THE GLOB AND THE BRANCH ARE SET BY HAND, and the helper cannot derive either.  From
% '*_BFI_r.mat' it yields the branch-agnostic '*_I_r.mat', which would offer this step
% on a cardiac cycle and on a bolus; and it clears branch 'contrast' to '' because there
% is no contrast branch here.  Vasomotion lives on the ANGIOGRAM alone: '_c_I' is one
% averaged beat on a phase axis with no clock a wavelet can run on, and '_b_I' is
% twenty-five to thirty seconds, which is one to seven cycles of the slowest band edge
% and inside the cone of influence for every one of them.
t.inGlob = '*_a_I_r.mat';
t.branch = 'angiogram';
% requires NAMES ITS ENTRY STEP, unlike 'dynamicSegmentationI' two blocks up.  There the
% cascade through segmentationI reaches the right producer by itself; here it would not
% necessarily - segmentationI is satisfied by ANY of the three entry steps, and a chain
% that satisfied it through the cardiac or the bolus one would leave this step with no
% '_a_I' to read.  The speckle original says requires={'BFI'}, which is the product
% token this branch does not have [D10]; the honest prerequisite here is the step that
% writes the angiogram plus the step that puts traces on it.
t.requires = {'segmentationI','intensity'};
reg(end+1) = t;

% ---- 34. pulsatilityI (fluorescence pulsatility - its OWN wrapper) -----------
s = base();
s.id='pulsatilityI'; s.label='Pulsatility'; s.wrapper=@runIntensityPulsatility;
% WRITTEN OUT BY HAND AND NOT AN intensityTwin, and it is the only step in this file
% that keeps a twin's ID while running a different wrapper [D13, author 07-Aug-2026].
% The ID stays because a step declares ONE inGlob, which is the whole reason a twin
% needs its own id, and that reason is unchanged; changing it would ripple into
% wbTypePresets, wbRefBranch and wbStateEngine for nothing.
%
% WHY IT IS NOT A TWIN.  runPulsatility's own header: "a field-name PREFIX encoding the
% physical quantity (ps = pulsatile flow, pd = pulsatile diameter)".  Pointed at a
% '_c_I' product it would write psPI meaning a pulsatile INTENSITY - one column name
% carrying two physical quantities, in a table that is pooled across recordings and
% exported to one spreadsheet.  runIntensityPulsatility writes pv*, for pulsatile
% VOLUME: fluorescence intensity is proportional to the labelled plasma in the light
% path, so that is what was measured.
%
% AND IT WRITES NO pd* COLUMN.  The per-frame diameter runDynamicSegmentation produces
% is quantised at a quarter of a pixel, and the cardiac width change on the reference
% recording is 1.6 to 18 times smaller than one such step.  Where a cardiac width change
% IS measured on this branch is 'motionEnhancement', sub-pixel and against a matched
% control.  The two must never be pooled, exactly as an nd* column and an ns* one.
%
% IT DOES NOT MAGNIFY ANYTHING (master plan section 7, Q1).  One motion step, not two:
% a trace step has nowhere to keep a matched control, and two implementations of the
% magnifier would drift apart in silence.  The Cardiac preset ticks BOTH steps.
s.inGlob='*_c_I_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runIntensityPulsatility';
% requires NAMES THE CARDIAC ENTRY STEP EXPLICITLY, and here that is not optional:
% segmentationI's own requiresAny puts 'intensity' FIRST, so a one-click chain that
% reached this step through segmentationI alone would pull in the angiogram and leave it
% with no '_c_I' to read.
s.requires={'segmentationI','intensityCycle'}; s.produces={'pulsatility'};
s.artifacts={'_rep_pulsatility.jpg'};   % no legacy tail: nothing on this branch
                                        % predates the unified reporting rename, and
                                        % the speckle pulsatility step writes no page
                                        % at all, so no column can claim another's
s.modalities={'EPFL'}; s.branch='cardiac';
s.settingGroups={ 'Harmonic model',{'nHarm'};
                  'Analysis levels',{'segPulsReturn','ppxPulsReturn'} };
s.basicFields={'segPulsReturn','ppxPulsReturn'};
s.sharedKeys={'libraryFolder'};
% THE MODEL IS OFF BY DEFAULT, and that is measured rather than tasteful.  A fluorescence
% beat at these frame rates carries about ten independent samples, which supports the
% fundamental and four harmonics; a five-harmonic model has eleven free parameters, so
% it would interpolate the beat rather than describe it and its R2 would read about 1
% whatever the data said.  There is no per-pixel fit at all for the same reason.
s.presets=struct('default',struct('nHarm',5, ...
    'segPulsReturn',{{'markers'}},'ppxPulsReturn',{{'markers'}}));
s.labels=struct('nHarm','Harmonics in the model', ...
                'segPulsReturn','What to keep per vessel', ...
                'ppxPulsReturn','What to picture per pixel');
s.tips=struct('nHarm',['how many harmonics the optional model uses.  About ten points ' ...
        'of a fluorescence beat are actually measured, so five is the most a model can ' ...
        'be told apart from the beat itself'], ...
    'segPulsReturn',['which numbers to keep per vessel.  The plain measurements are ' ...
        'the default and need no model; adding the model fits a sum of sinusoids to ' ...
        'every beat, which on this kind of recording follows the beat rather than ' ...
        'describing it'], ...
    'ppxPulsReturn',['whether to also make a picture of each number, one value per ' ...
        'pixel.  Leave it empty for the vessels alone, which is what a table is read ' ...
        'from']);
s.note=['How much the dye in each vessel rises and falls over one heartbeat, as a ' ...
    'fraction of its own average - a volume, not a speed and not a width.  How much ' ...
    'each vessel''s width changed over the same beat is Wall motion''s answer and is ' ...
    'reported in its own columns; the two are different instruments and their numbers ' ...
    'do not belong in one column.'];
reg(end+1)=s;

% ---- 35. vesselTypesI (the fluorescence twin of the vessel-typing step) ------
% A GENUINE TWIN.  setVesselTypes paints artery / vein / uncertain onto segments and
% names them; nothing about painting is a fact about how the picture was made.  What IS
% modality-specific is the automatic guess underneath the painting, and that is resolved
% per PRODUCT by getVascularCues rather than per step here - so the registry entry is a
% twin and the science sits one layer down.
t = intensityTwinOf(stepById(reg,'vesselTypes'));
% THE REFERENCE BRANCH IS THE BOLUS, and this is the field the helper cannot derive: it
% carries 'cardiac' over verbatim, which is right on the speckle branch (runPulsatility
% writes ps* on the internal cycle) and wrong here.  On the fluorescence branch the
% column that separates an artery from a vein is the TRACER arrival runCTTH writes on the
% bolus product; a cardiac intensity product would drive the guess from pv* instead,
% which works but is confined to one beat.  The glob stays branch-agnostic '*_I_r.mat'
% and branchScope stays 'all', so a recording's angiogram and cardiac cycle are painted
% or inherited in the same call - refBranch only says which of them opens the editor.
t.refBranch = 'bolus';
% requires NAMES ONLY WHAT THE STEP CANNOT RUN WITHOUT, which is the shape the speckle
% original uses and the reason it does not name pulsatility: painting types by hand is a
% legitimate way to use this step, and a product with no timing on it still opens, still
% paints and still saves.  What changed at S10 is that it SAYS SO - setVesselTypes warns
% by name, writes the reason into the slider panel and prints it on the report page -
% instead of quietly labelling every vessel "Uncertain".  The note below is the same
% sentence before the box is ticked.  'BFI' is gone from the list because there is no BFI
% twin on this branch [D10]; the honest prerequisite is the segmentation alone.
t.requires = {'segmentationI'};
t.note = ['Which vessels are arteries and which are veins.  The automatic guess is ' ...
    'driven by when the tracer reaches each vessel, which Transit time measures.  ' ...
    'Without it the step still opens and the types can be painted by hand, and it ' ...
    'says so rather than reporting every vessel as undecided.'];
reg(end+1) = t;

% ---- 36. vascularTreeI (the fluorescence twin of the vascular-tree step) -----
% A GENUINE TWIN OF THE WRAPPER, ON A DIFFERENT PRODUCT.  setVascularTree orders the
% segments by a monotone flow potential and searches for parent/daughter links; the
% search knows nothing about what was measured, and defaultFlowParams resolves the
% columns per product.
t = intensityTwinOf(stepById(reg,'vascularTree'));
% THE GLOB AND BOTH BRANCH FIELDS ARE SET BY HAND, and this is the decision S10 took.
% From '*_c_BFI_r.mat' the helper yields '*_c_I_r.mat' - the cardiac intensity product,
% where pv* lives - and the master plan puts this step on the BOLUS instead.  The bolus
% wins, measured rather than assumed: a tracer's arrival spreads over seconds across the
% field (Delay p5-p95 0.11-0.77 s, Mtt 0.22-1.39 s over 5215 segments) while a pulse
% phase is confined inside one ~0.1 s beat AND WRAPS, so two vessels a beat apart in true
% arrival are indistinguishable on it.  The bolus is also the product on which this
% preparation's arterial-to-capillary-to-venous ordering has actually been measured.
% claude-docs/intensity-branch/10-vessel-types-tree.md is the argument in full.
t.inGlob    = '*_b_I_r.mat';
t.branch    = 'bolus';
t.refBranch = 'bolus';
% AND THE PARTNERS THAT INHERIT IT.  The helper carries the preset over verbatim, and
% {'t','s'} names two stages that do not exist on this branch - the copy would have found
% no such file and moved on in silence, leaving the angiogram and the cardiac cycle
% without a hierarchy and saying nothing.  setVascularTree resolves this per product when
% the field is unset; the preset says it out loud so the panel shows which products will
% inherit before the step is run.
t.presets.default.propagatePartners = {'a','c'};
t.tips.propagatePartners = 'after the injection recording is derived, copy the hierarchy to these partners';
% requires MIRRORS THE SPECKLE ORIGINAL ELEMENT FOR ELEMENT - {'vesselTypes','pulsatility'}
% becomes {'vesselTypesI','ctth'} - because the transit-time step is what occupies the
% timing slot on this branch.  Both are hard: getVascularTree REFUSES by name when no
% arrival parameter has data, so unlike the step above this one cannot degrade politely.
t.requires = {'vesselTypesI','ctth'};
t.note = ['Which vessel feeds which, from the arterial inlet to the venous outlet.  ' ...
    'The direction comes from when the tracer reached each vessel, so this runs on the ' ...
    'injection recording and the result is copied onto the others.  It stops rather ' ...
    'than guessing if the transit times are not there.'];
reg(end+1) = t;

validateRegistry(reg);          % a registry edit is checked before anyone sees it

% ---- optional modality filter, and the pruning that must follow it ----------
if nargin>=1 && ~isempty(modality), reg = filterTo(reg, modality); end
end

% =====================================================================
function t = intensityTwinOf(s, alsoDrop)
%intensityTwinOf  The FLUORESCENCE twin of a speckle step, derived in ONE place.
%
%   A registry step declares one inGlob, so a step that must serve both '_K'/'_BFI'
%   and '_I' needs two IDS driving the SAME WRAPPER.  Written out by hand twenty times
%   those two would drift; written here once they cannot.  The label is deliberately
%   NOT changed - both say 'Segmentation' - because the modality filter means a user
%   only ever sees one of them.
%
%   WHAT IT CARRIES OVER
%     wrapper, arity, branchScope, fanOut, fileOrder, interactive, needsRaw,
%     artifacts, note, produces, and every setting that is not dropped below.
%     gatingField is carried VERBATIM, and that is the point: the settings field is
%     named after the WRITER, the writer is the same function, so a recording
%     processed on either branch records the same field and wbStateEngine needs to
%     know nothing about which branch it was.
%
%   WHAT IT REWRITES
%     id        gains a trailing 'I'  ('segmentation' -> 'segmentationI')
%     modalities becomes {'EPFL'}
%     inGlob    the PRODUCT TOKEN becomes I and the stage flag is kept:
%                 '*_K_r.mat'      -> '*_I_r.mat'      (any branch)
%                 '*_BFI_r.mat'    -> '*_I_r.mat'      - there is no BFI twin, because
%                                     fluorescence intensity is already proportional
%                                     to plasma volume and there is nothing to invert
%                 '*_c_BFI_r.mat'  -> '*_c_I_r.mat'
%               A glob naming '_t_' or '_s_' is REFUSED: those two stages exist only
%               on the speckle side, so there is no honest rewrite of them and the
%               caller has to say which intensity branch it meant.
%
%   WHAT IT REFUSES TO GUESS, and why each one is cleared rather than carried
%     requires / requiresAny / conflictsWith  - these name STEP IDS, and a twin's
%       producers are different steps.  Carried over they would name LSCI ids that do
%       not survive the EPFL filter, pruneReferences would strike them silently, and
%       the twin would end up with no prerequisite at all while its source still LOOKS
%       correct.  Cleared, the registry entry visibly says what has not been decided.
%     branch / refBranch when they say 'contrast' - there is no contrast branch on the
%       fluorescence side.  'cardiac' and '' are carried, because those are the same
%       branch on both.
%     legacyArtifacts - the pre-rename tails are a speckle-era fact; no fluorescence
%       run ever wrote one, so carrying them could only ever list another step's file.
%
%   AND ONE THING IT THROWS ON.  A step declaring an outSuffix is a RAW PRODUCER or a
%   step that writes a new triplet, and none of the planned twins is: an entry step is
%   written out by hand because its output name, its branch and its settings are all
%   its own.  Deriving one here would produce a name nothing reads, which is the
%   quietest failure in this file.
%
%   SETTINGS DROPPED, and why the list is written out rather than pattern-matched.
%   The contrast-only settings are named one by one below because a name rule cannot
%   do it: 'trustLimitsK' must go and 'trustLimitsI' must stay, and 'sMaxKK' is a
%   curvature bound that ends in K and has nothing to do with contrast.  A field is
%   removed from settingGroups, basicFields, sharedKeys, every preset, tips, enums and
%   labels together - anything less leaves a tooltip for a row that is not there.
%   Pass alsoDrop to remove more; the twin GAINS nothing here, because anything an
%   intensity step needs that its speckle source does not is a decision the session
%   that writes that step makes, not one this helper can derive.
%
%   Syntax:
%      twin = wbStepRegistry('intensityTwin', step)
%      twin = wbStepRegistry('intensityTwin', step, {'extraField',...})
%
%   Example:
%      t = wbStepRegistry('intensityTwin', stepById(reg,'segmentation'));
%      t.requiresAny = {'intensity','intensityCycle','intensityBolus'};

CONTRAST_ONLY = {'trustLimitsK','maskLimitsK','contrastType','contrastKernel', ...
                 'contrastKernelS','contrastKernelT','contrastKernelPreproc'};

if nargin<2, alsoDrop = {}; end
if ~isstruct(s) || ~isscalar(s) || ~isfield(s,'id')
    error('wbStepRegistry:twinNotAStep', ...
        'intensityTwin takes ONE fully-formed step struct.');
end
if ~isempty(s.outSuffix)
    error('wbStepRegistry:twinWritesAProduct', ...
        ['Step "%s" writes its own triplet (%s), so it is an entry step and its twin ' ...
         'has to be written out by hand - its output name, branch and settings are ' ...
         'all its own.'], s.id, strjoin(reshape(s.outSuffix,1,[]),', '));
end

t = s;
t.id         = [char(s.id) 'I'];
t.modalities = {'EPFL'};
t.inGlob     = intensityGlob(s.inGlob, s.id);

% the three fields that name STEP IDS, and the two that name a speckle-only branch
t.requires      = {};
t.requiresAny   = {};
t.conflictsWith = {};
if strcmp(t.branch,'contrast'),    t.branch    = ''; end
if strcmp(t.refBranch,'contrast'), t.refBranch = ''; end
t.legacyArtifacts = {};

t = dropSettings(t, [reshape(CONTRAST_ONLY,1,[]), reshape(alsoDrop,1,[])]);
end

% =====================================================================
function st = stepById(reg, id)
%stepById  The one step of reg with this id.  The twins are derived from the steps
%   already in the array rather than from a second copy of their definitions, so a
%   later edit to setRegions or segmentation reaches its twin by construction.
st = reg(strcmp({reg.id}, id));
if isempty(st)
    error('wbStepRegistry:noSuchStep','No step "%s" to derive a twin from.', id);
end
st = st(1);
end

% =====================================================================
function g = intensityGlob(g, id)
%intensityGlob  The same glob with the product token swapped for I, or a refusal.
g = char(g);
if isempty(g), return; end
if ~endsWith(g,'.mat')
    error('wbStepRegistry:twinRawGlob', ...
        'Step "%s" reads a raw recording (%s); an entry step is not twinned.', id, g);
end
if contains(g,'_t_') || contains(g,'_s_')
    error('wbStepRegistry:twinStageOnlySpeckle', ...
        ['Step "%s" reads %s, and the temporal/spatial stages exist only on the ' ...
         'speckle side.  Set the twin''s inGlob yourself to the intensity branch ' ...
         'it means.'], id, g);
end
g = strrep(g,'_BFI_','_I_');
g = strrep(g,'_K_','_I_');
end

% =====================================================================
function st = dropSettings(st, names)
%dropSettings  Remove settings from EVERY place one step mentions them.  A field left
%   in tips or in a preset after leaving settingGroups is a tooltip and a default for a
%   row nobody can see, which is how a dead setting survives a redesign.
names = names(~cellfun(@isempty,names));
if isempty(names), return; end

% settingGroups is Nx2 {label,{fields}}; a group emptied by the removal goes with it
if ~isempty(st.settingGroups)
    keep = true(size(st.settingGroups,1),1);
    for k = 1:size(st.settingGroups,1)
        f = reshape(st.settingGroups{k,2},1,[]);
        f = f(~ismember(f,names));
        st.settingGroups{k,2} = f;
        keep(k) = ~isempty(f);
    end
    st.settingGroups = st.settingGroups(keep,:);
end

for f = {'basicFields','sharedKeys'}
    l = reshape(st.(f{1}),1,[]);
    st.(f{1}) = l(~ismember(l,names));
end

% every NAMED preset bundle, not just 'default'
pn = fieldnames(st.presets);
for i = 1:numel(pn)
    b = st.presets.(pn{i});
    for k = 1:numel(names)
        if isfield(b,names{k}), b = rmfield(b,names{k}); end
    end
    st.presets.(pn{i}) = b;
end

for f = {'tips','enums','labels'}
    b = st.(f{1});
    for k = 1:numel(names)
        if isfield(b,names{k}), b = rmfield(b,names{k}); end
    end
    st.(f{1}) = b;
end
end

% =====================================================================
function reg = filterTo(reg, modality)
%filterTo  Keep the steps ANY of these modalities exposes, then prune.  A char is
%   one modality, a cellstr the union over several - in registry order either way.
if isempty(modality) || isempty(reg), return; end
want = modality; if ~iscell(want), want = {char(want)}; end
keep = arrayfun(@(st) any(ismember(want, st.modalities)), reg);
reg  = pruneReferences(reg(keep));
end

% =====================================================================
function reg = pruneReferences(reg)
%pruneReferences  Strike every id the modality filter removed out of the steps
%   that survived it, so no surviving step names one that is not there (see the
%   note above).  Nothing is renamed or re-ordered - only names of absent steps go.
if isempty(reg), return; end
have = {reg.id};
for k = 1:numel(reg)
    for f = {'requires','requiresAny','conflictsWith'}
        l = reg(k).(f{1});
        if isempty(l), continue; end
        reg(k).(f{1}) = reshape(l(ismember(l, have)),1,[]);
    end
end
end

% =====================================================================
function validateRegistry(reg)
%validateRegistry  The two things a step array may not do, checked where a
%   registry edit is made rather than where a user clicks.
%
%   A CONFLICT MUST BE SYMMETRIC: ticking A unticks B only if A names B, so a
%   one-sided declaration would make the rule depend on which box was clicked.
%   A CONFLICT MAY NOT ALSO BE A PREREQUISITE: the tick cascade would pull the
%   step in and the conflict rule would push it straight back out.
if isempty(reg) || ~isstruct(reg) || ~isfield(reg,'conflictsWith'), return; end
ids = {reg.id};
for k = 1:numel(reg)
    cw = reg(k).conflictsWith;
    if isempty(cw), continue; end
    cw = reshape(cw,1,[]);
    for i = 1:numel(cw)
        j = find(strcmp(cw{i}, ids),1);
        if isempty(j)
            error('wbStepRegistry:unknownConflict', ...
                'Step "%s" conflicts with "%s", which is not a step.', ids{k}, cw{i});
        end
        if ~any(strcmp(ids{k}, reg(j).conflictsWith))
            error('wbStepRegistry:asymmetricConflict', ...
                'Step "%s" conflicts with "%s" but "%s" does not say so; declare it on both sides.', ...
                ids{k}, cw{i}, cw{i});
        end
    end
    clash = intersect(cw, wbPrereqs('all', reg(k)), 'stable');
    if ~isempty(clash)
        error('wbStepRegistry:conflictIsPrerequisite', ...
            'Step "%s" both requires and conflicts with %s.', ids{k}, strjoin(clash,', '));
    end
end
end

% =====================================================================
function s = base()
%base  A fully-formed step struct with every field defaulted, so appending
%   elements to the registry array never triggers a field-mismatch error.
s = struct( ...
    'id','', 'label','', 'wrapper',[], 'arity','perFile', ...
    'inGlob','', 'outSuffix',{{}}, 'outKind','inplace', 'outTransform',[], ...
    'gatingField','', 'requires',{{}}, 'requiresAny',{{}}, 'conflictsWith',{{}}, ...
    'produces',{{}}, ...
    'interactive',false, 'needsRaw',false, 'artifacts',{{}}, 'legacyArtifacts',{{}}, ...
    'settingGroups',{{}}, 'basicFields',{{}}, 'sharedKeys',{{}}, ...
    'presets',struct('default',struct()), 'tips',struct(), 'enums',struct(), ...
    'labels',struct(), 'note','', ...
    'modalities',{{'LSCI'}}, 'branch','', 'branchScope','one', 'fanOut','flat', ...
    'fileOrder','independent', 'refBranch','');
end

% =====================================================================
function t = parforTip(what)
%parforTip  One sentence, same shape for every parallel switch: what it parallelises,
%   and what turning it off actually costs.  The trade is time against worker
%   processes, and that is the only thing the user has to decide.
t = sprintf('on: run %s on several CPU workers.  off: one item at a time in this MATLAB - slower, but no worker processes and no extra memory.', what);
end
