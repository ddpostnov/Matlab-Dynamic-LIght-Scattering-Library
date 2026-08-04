%wbStepRegistry - Declarative spec of the Processing-Workbench steps.
%
%   Returns the ordered array of pipeline-step specifications that drives every
%   part of the workbench: the matrix columns, the on-disk gating, the settings
%   panels, execution, and the report links.  Adding a modality or a step later
%   is a data edit here, not a code rewrite.  The specs are transcribed FROM THE
%   CODE (the Wrappers/ functions and the Launchers/ %Example-of-s blocks) as of
%   branch main - see claude-docs/processing-workbench/01-pipeline-map.md.
%
%   Steps are listed in dependency order, speckle first and then the myograph:
%     contrast · internalCycle · externalCycle · setRegions · segmentation ·
%     dynamicSegmentation · guided · registration · BFI · vasomotion ·
%     pulsatility · vesselTypes · vascularTree ·
%     myoVideo · labChart · myoPresetIntervals · myoCrop · myoDiameter ·
%     myoIntervals · myoPropagation · myoVasomotion
%
% Syntax:
%    reg = wbStepRegistry()
%    reg = wbStepRegistry(modality)      % filter to steps exposing that modality
%    reg = wbStepRegistry(modalities)    % cellstr: the UNION over several, in
%                                        %   registry order (a mixed working set)
%    reg = wbStepRegistry('filter', reg, modality)   % the same filter + prune,
%                                        %   applied to a step array you hold
%    wbStepRegistry('validate', reg)     % throw on a malformed step array
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
%       branch       'contrast' | 'cardiac' | '' (which file branch this step's
%                    column belongs to; wbFileModel also derives 'epoch'/'bolus')
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
%   * THREE MODALITY FAMILIES ARE HERE: speckle (LSCI), the PRESSURE MYOGRAPH
%     (PMYO, a video) and the WIRE MYOGRAPH (WMYO, a LabChart '.adicht' recording).
%     The .cxd family (bolus / intensity / CTTH) is still deferred; the schema
%     supports it.
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
%   * REAL gating fields (differ from the function name, 01 A1):
%     runBFI->calculateBFI, runExternalCycle->externalCycle,
%     runCTTH->ctthCalculation.  Encoded below.
%   * BRANCH: 'contrast' is the temporal OR spatial contrast side (t|s) - the user
%     picks per s.contrastType and the analysis (segmentation/BFI/vasomotion/...)
%     runs on either; 'cardiac' is the internal-cycle side (c); '' is branch-
%     agnostic.  vasomotion is 'contrast' (t|s only, NOT c or e).  State stays at
%     the recording level here; per-file column filtering by branch/stage is a
%     Phase-3 concern (this field is the hint it uses).
%   * STAGE flag is a SINGLE token: the external cycle is '_e_K' (NOT '_e_t_K')
%     and the internal cycle '_c_K'.  A cycle's contrast base (t|s) is recorded in
%     the SETTINGS, not the name, since a project is not expected to mix bases -
%     so the suffix stays simple, showing only what the next step needs
%     (t_K / s_K / c_K / e_K).  NOTE runExternalCycle currently WRITES the legacy
%     '_t_e_K' (strrep _K_d->_e_K_d on a _t_K_d input); wbFileModel still parses it
%     to the right identity, but reconciling the wrapper to emit '_e_K' is a
%     recommended separate fix.
%   * 'guided' maps to runGuidedContrast (intensity variant deferred).
%     setRegions is modelled perFile (only registration/vesselTypes are perAnimal).
%   * EXCEL EXPORT IS NOT A STEP (author, 2026-07-28).  It is neither configured in
%     the Constructor nor run by the Processor: it is a STANDALONE tool that reads a
%     finished session (spec §7/D9).  Everything here writes data products; the
%     export writes a report, and mixing the two put a workbook column in a
%     processing matrix that nobody wanted to tick.
%   * BRANCHSCOPE - one raw recording yields SEVERAL co-registered products ('_t_K',
%     '_c_K', '_e_K', ...), and the launchers are explicit about how many of them a
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
%     RESOLVED (Phase 6, 2026-07-29): registration and vesselTypes are 'all'.
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
%     externalCycle is perFile / flat / branchScope 'one' / no refBranch - the exact
%     profile of an independent step - yet it indexes a per-file stimulus list, so
%     every rule that guesses gets it wrong (author, 2026-07-31).
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
%     prefers '_c' too.  '' = no preference.  It is the reference twin of
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
% Last revision: 01-August-2026

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
s.id='internalCycle'; s.label='Internal cycle'; s.wrapper=@runInternalCycle;
s.inGlob='*.rls'; s.outSuffix={'_c_K_d','_c_K_r','_c_K_s'}; s.outKind='new';
s.outTransform=struct('from','.rls','to','_c_K_r.mat');
s.gatingField='runInternalCycle'; s.requires={}; s.produces={'cardiac'};
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

% ---- 3. externalCycle (NVC) -------------------------------------------------
s = base();
s.id='externalCycle'; s.label='External cycle'; s.wrapper=@runExternalCycle;
% external/epoch cycle of the contrast side (t or s): the stage flag BECOMES e,
% i.e. the single suffix _e_K (the contrast base is kept in the settings, not the
% name - see wbFileModel rationale)
s.inGlob='*_K_r.mat'; s.outSuffix={'_e_K_d','_e_K_r','_e_K_s'}; s.outKind='new';
s.outTransform=struct('from','_t_K_r.mat','to','_e_K_r.mat');   % replaces the t|s flag with e
s.gatingField='externalCycle';                       % REAL field (differs from fn name)
s.requires={'contrast'}; s.produces={'epochAvg'};
s.interactive=@(ss) isfield(ss,'enablelRejectionModification') && isequal(ss.enablelRejectionModification,1);
s.artifacts={'_rep_epochs.jpg','_rep_epoch-average.jpg'};
s.legacyArtifacts={'_ec.jpg','_ec2.jpg'}; s.branch='contrast';
s.fileOrder='ordered';                               % the stimulus timing is a PER-FILE list
                                                     % indexed by position in fNames

s.settingGroups={ 'Stimulation',{'stimStartType','stimOffset','epochsN','epochDurationSec', ...
                     'epochBaselineSec','epochStimStartSec','stimDurationSec','epochFinaleSec'};
                  'Masking',{'maskType','enablelRejectionModification'};
                  'Rejection',{'rejectBlCoef','rejectEpochCoef','rejectFinCoef','rejectPeakCoef', ...
                     'rejectBlSimCoef','rejectSimCoef','rejectTimeLoss','rejectFirstEpoch'} };
s.basicFields={'stimStartType','stimOffset','epochsN','epochDurationSec','epochBaselineSec','epochStimStartSec','stimDurationSec','epochFinaleSec','maskType','enablelRejectionModification'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('stimStartType',{{'offset','manual'}},'maskType',{{'basic','cMask','selection'}});
s.presets=struct('default',struct('stimStartType','offset','stimOffset',0,'epochsN',20, ...
    'epochDurationSec',30,'epochBaselineSec',[0 10],'epochStimStartSec',10,'stimDurationSec',5, ...
    'epochFinaleSec',[-5 0], ...
    'maskType','cMask','enablelRejectionModification',1,'rejectBlCoef',1,'rejectEpochCoef',1, ...
    'rejectFinCoef',1,'rejectPeakCoef',1,'rejectBlSimCoef',1,'rejectSimCoef',1, ...
    'rejectTimeLoss',0.5,'rejectFirstEpoch',1));
s.tips=struct('maskType','''cMask'' reads the segmentation mask; ''basic'' the contrast mask', ...
    'enablelRejectionModification','1 = interactive epoch-rejection editing', ...
    'stimStartType','''offset'' (fixed delay) or ''manual'' (timestamp list)', ...
    'stimDurationSec',['how long the stimulus lasts, seconds.  This step does not use ' ...
       'it - it is recorded here so the response fit reads the protocol instead of ' ...
       'being told it again']);
reg(end+1)=s;

% ---- 4. setRegions ----------------------------------------------------------
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

% ---- 5. segmentation --------------------------------------------------------
s = base();
s.id='segmentation'; s.label='Segmentation'; s.wrapper=@runSegmentation;
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runSegmentation'; s.requires={}; s.requiresAny={'contrast','internalCycle'};
s.produces={'segmentation'};
s.artifacts={'_rep_categories.jpg','_rep_segments.jpg'};
s.legacyArtifacts={'_cm.jpg','_vs.jpg'};              % '_vs' came from showSegmentsPreview
s.branch=''; s.branchScope='copy';                % computed on contrast side, copied to cardiac
s.settingGroups={ 'Contrast',{'trustLimitsK'};
                  'Categorization',{'lSizeN','sSizeN','sens','sSizeScale','deSens','lThinN','imOpen','iEdge','eEdge'};
                  'Labelling & traces',{'sStat','sMinL','prchNSize','correctNodes','simR','difR'};
                  'Parallel',{'parforSegmentationLabels'} };
                  % no 'Copy to siblings' panel: branchScope 'copy' means wbExecutor
                  % derives s.fNamesCopyTo from the recording's own branch products
s.basicFields={'lSizeN','sSizeN','sens','sMinL'};
% parforSegmentationLabels is SHARED with dynamicSegmentation: both steps drive the
% same core (getSegmentationLabels), so the machine-level choice is one tick, not two.
s.sharedKeys={'trustLimitsK','prchNSize','sMinL','correctNodes','simR','difR','sStat', ...
    'parforSegmentationLabels','libraryFolder'};
s.enums=struct('sStat',{{'median','mean'}});
s.presets=struct('default',struct('trustLimitsK',[0.001 0.99],'lSizeN',141,'sSizeN',9, ...
    'sens',0.1,'sSizeScale',1,'deSens',1,'lThinN',2,'imOpen',0,'iEdge',2,'eEdge',2, ...
    'categories',{{'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}}, ...
    'sStat','median','sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3,'difR',0.4, ...
    'parforSegmentationLabels',true));
s.labels=struct('parforSegmentationLabels','Use parfor: segment labels');
s.tips=struct('lSizeN','odd, ~2x the largest vessel', ...
    'sSizeN','odd, ~2x the small-vessel diameter', ...
    'sens','segmentation sensitivity (raise to catch faint vessels)', ...
    'sMinL','minimum segment length', ...
    'prchNSize','parenchymal pixel neighbourhood', ...
    'parforSegmentationLabels',parforTip('the per-segment label growing'));
reg(end+1)=s;

% ---- 6. dynamicSegmentation -------------------------------------------------
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

% ---- 7. guided --------------------------------------------------------------
s = base();
s.id='guided'; s.label='Guided'; s.wrapper=@runGuidedContrast;
s.inGlob='*_K_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runGuidedContrast'; s.requires={'segmentation'}; s.produces={'guidedTraces'};
s.needsRaw=true; s.branch='contrast';
s.settingGroups={};                                   % uses sMap + raw; no tunable params
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct());
reg(end+1)=s;

% ---- 8. registration (per animal) --------------------------------------------
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

% ---- 9. BFI -----------------------------------------------------------------
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

% ---- 10. vasomotion ---------------------------------------------------------
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

% ---- 11. pulsatility --------------------------------------------------------
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

% ---- 12. fitNVC (stimulus-locked response fit) -------------------------------
s = base();
s.id='fitNVC'; s.label='NVC fit'; s.wrapper=@runFitNVC;
s.inGlob='*_e_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runFitNVC';
% REQUIRES NAMES STEPS, NOT CAPABILITIES.  wbPrereqs matches these against step
% IDS, so the epoch product's producer is 'externalCycle' - the token it PRODUCES
% ('epochAvg') is a separate vocabulary nothing in the dependency graph reads.
s.requires={'BFI','externalCycle'}; s.produces={'nvcFit'};
% THE ROW IS 'contrast', NOT 'epoch', AND THAT IS NOT A TYPO.  There are two
% different senses of "branch" in this library and only one of them is a row.
% wbFileModel>branchOf classifies a FILE by its stage flag, and an '_e' product is
% indeed on the 'epoch' branch there.  A ROW, though, exists only where a RAW
% PRODUCER declares one (wbTypeSelection>branchesOf) - contrast, cardiac, myograph -
% and externalCycle is a DERIVED consumer of '*_K_r.mat', so it creates no row of
% its own.  Declaring 'epoch' here made wbTypeSelection>offers false on every row
% that exists and the step unreachable in the workbench.  The epoch stage is part of
% the contrast PIPELINE, which is exactly what wbStateEngine>samePipeline already
% says - "a stage the pipeline grew into later (the external cycle's '_e')" - so the
% NVC chain is ticked on the contrast row, beside the externalCycle that starts it.
s.branch='contrast';
s.artifacts={'_rep_nvcfit.jpg'};
s.settingGroups={ 'Model',{'nvcModel','nvcDip','stimDurationSec'};
                  'Analysis levels',{'segNvcReturn','ppxNvcReturn'};
                  'Fitting',{'nvcStarts'};
                  'Parallel',{'parforNvcSegments','parforNvcPixels'} };
s.basicFields={'stimDurationSec','nvcModel'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('nvcModel',{{'secondOrder','doubleGamma'}});
s.presets=struct('default',struct('nvcModel','secondOrder','nvcDip',false, ...
    'stimDurationSec',5,'segNvcReturn',{{'markers','model','reconstruction'}}, ...
    'ppxNvcReturn',{{'markers'}},'nvcStarts',16, ...
    'parforNvcSegments',true,'parforNvcPixels',true));
s.labels=struct('parforNvcSegments','Use parfor: segments', ...
                'parforNvcPixels','Use parfor: per pixel');
s.tips=struct( ...
    'nvcModel',['''secondOrder'' fits the flow response of a damped autoregulatory ' ...
       'system - its undershoot is mechanistic and its damping compares across ' ...
       'protocols.  ''doubleGamma'' only describes the shape'], ...
    'nvcDip','also fit a fast constriction before the rise, and report both fits side by side', ...
    'stimDurationSec',['how long the stimulus lasts, seconds.  Taken from the ' ...
       'external cycle step unless it is set here'], ...
    'segNvcReturn','per-segment levels: markers / model / reconstruction', ...
    'ppxNvcReturn',['per-pixel maps ([] = off).  Adding ''model'' fits every pixel ' ...
       'of the field and takes minutes rather than seconds'], ...
    'nvcStarts',['how many start points the fit is tried from.  Fewer is faster and ' ...
       'more likely to settle on a wrong answer'], ...
    'parforNvcSegments',parforTip('the per-segment fit'), ...
    'parforNvcPixels',parforTip('the per-pixel analysis'));
reg(end+1)=s;

% ---- 13. fitVasoreactivity (drug / gas challenge response fit) ---------------
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

% ---- 14. vesselTypes (per animal) --------------------------------------------
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

% ---- 15. vascularTree -------------------------------------------------------
s = base();
s.id='vascularTree'; s.label='Vascular tree'; s.wrapper=@setVascularTree;
s.inGlob='*_c_BFI_r.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setVascularTree'; s.requires={'vesselTypes','pulsatility'}; s.produces={'hierarchy'};
s.interactive=@(ss) ~(isfield(ss,'autoOnly') && isscalar(ss.autoOnly) && isequal(ss.autoOnly,true));
s.branch='cardiac'; s.refBranch='cardiac';   % the hierarchy is derived on the _c side
s.settingGroups={ 'Hierarchy derivation',{'autoOnly','phiWeights','useHarmonicPhase','propagatePartners'} };
s.basicFields={'autoOnly'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('autoOnly',false,'phiWeights',[1 1 1],'useHarmonicPhase',false, ...
    'propagatePartners',{{'t','s'}}));
s.tips=struct('autoOnly','true = derive & save without opening the tree editor', ...
    'phiWeights','relative weight of [foot peak -PI] in the flow potential', ...
    'propagatePartners','after a _c file is derived, copy the hierarchy to these partners');
reg(end+1)=s;

% ---- 16. myoVideo (pressure myograph, entry step) ----------------------------
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

% ---- 17. labChart (wire myograph, entry step) --------------------------------
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

% ---- 18. myoPresetIntervals (before anything is measured) --------------------
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

% ---- 19. myoCrop -------------------------------------------------------------
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

% ---- 20. myoDiameter ---------------------------------------------------------
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

% ---- 21. myoIntervals --------------------------------------------------------
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

% ---- 22. myoPropagation ------------------------------------------------------
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

% ---- 23. myoVasomotion -------------------------------------------------------
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

validateRegistry(reg);          % a registry edit is checked before anyone sees it

% ---- optional modality filter, and the pruning that must follow it ----------
if nargin>=1 && ~isempty(modality), reg = filterTo(reg, modality); end
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
