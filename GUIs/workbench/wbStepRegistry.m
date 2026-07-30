%wbStepRegistry - Declarative spec of the v1 (LSCI) Processing-Workbench steps.
%
%   Returns the ordered array of pipeline-step specifications that drives every
%   part of the workbench: the matrix columns, the on-disk gating, the settings
%   panels, execution, and the report links.  Adding a modality or a step later
%   is a data edit here, not a code rewrite.  The specs are transcribed FROM THE
%   CODE (the Wrappers/ functions and the Launchers/ %Example-of-s blocks) as of
%   branch main - see claude-docs/processing-workbench/01-pipeline-map.md.
%
%   Steps are listed in dependency order:
%     contrast · internalCycle · externalCycle · setRegions · splitRegions ·
%     segmentation · dynamicSegmentation · guided · registration · BFI ·
%     vasomotion · pulsatility · vesselTypes · vascularTree
%
% Syntax:
%    reg = wbStepRegistry()
%    reg = wbStepRegistry(modality)      % filter to steps exposing that modality
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
%       outTransform strrep rule on the _d name      (struct from/to, or [])
%       gatingField  settings.<field> written        ('runContrastFromRLS')
%       requires     upstream step ids, ALL needed  ({} | {'setRegions'} ...)
%       requiresAny  upstream step ids, ANY one is enough - the branch-agnostic
%                    middle of the pipeline (see below); wbPrereqs owns the rule
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
%       modalities   which modalities expose it       ({'LSCI'})
%       branch       'contrast' | 'cardiac' | '' (which file branch this step's
%                    column belongs to; wbFileModel also derives 'epoch'/'bolus')
%       branchScope  'one' | 'all' | 'copy'  - HOW MANY of a recording's branch
%                    products the step consumes (see below)
%       fanOut       'flat' | 'animal' - the SHAPE of the fNames a call gets when
%                    several recordings are batched into it (see below)
%       refBranch    'contrast' | 'cardiac' | '' - which BRANCH of the animal's
%                    REFERENCE RECORDING this step prefers (see below)
%
% Notes (for the author):
%   * The v1 registry is LSCI only.  The .cxd (bolus/intensity/CTTH) and .avi
%     (myograph) steps are deferred to the next steps; the schema supports them.
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
%     step touches: '*_t_K_d.mat' means one branch, '*_K_d.mat' means all of them.
%     A workbench ROW is the recording, not the file, so each step declares the
%     fan-out its launcher cell implies (wbExecutor>buildFNames does the resolving):
%       'one'  - a single file, disambiguated by branch/desiredStage.  The default,
%                and correct for a step that is meaningful on ONE branch only
%                (vasomotion on _t, pulsatility/vascularTree on _c, the entry steps
%                that read the raw recording).
%       'all'  - every branch product of the recording, as an Nx1 fNames column, the
%                way the launcher passes '*_K_d.mat' (splitRegions, BFI,
%                dynamicSegmentation).  The wrapper's own per-file loop then
%                covers the branches, one workbench cell for the lot.
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
%     Their launcher cells are branch-wide ('Roi*_K_d.mat' and 'Roi*_BFI_d.mat',
%     both of which match _t AND _c), and leaving them at 'one' did not just under-
%     cover - it silently picked the WRONG file, because resolveStepInputs falls
%     back to the first dir() match when the step declares no branch, and
%     '_c_K_d.mat' sorts before '_t_K_d.mat'.  A two-pipeline recording was
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
%     getFileNamesList(rootFolder,'*_t_K_d.mat','[A-Z]+\d+') - whose grouping
%     pattern is the animal token.  Left flat, every recording would arrive as its
%     own row and the drawn region could never persist from one recording of an
%     animal to the next, which is the whole point of the step (author, 2026-07-29).
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
%   * REQUIRES vs REQUIRESANY.  '*_K_d.mat' steps (setRegions, segmentation, BFI,
%     registration) read ANY branch product of a recording, so their real
%     prerequisite is "some entry step has run", not "contrast has run".  A purely
%     pulsatile protocol starts at internalCycle and never computes a contrast
%     product - listing 'contrast' as a hard requirement made it tick a step it
%     does not run.  They therefore declare requiresAny={'contrast','internalCycle'};
%     the FIRST entry is the default producer a one-click chain pulls in.  Hard
%     chains (splitRegions after setRegions, pulsatility after BFI+internalCycle)
%     stay in 'requires'.  wbPrereqs is the single definition of "satisfied".
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
%
% See also: wbFileModel, wbStateEngine, wbSettingsModel, wbInvalidate,
%           wbPrereqs, wbRefBranch, wbTypeSelection, wbTypePresets
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function reg = wbStepRegistry(modality)

reg = base();  reg(1) = [];   % empty 1x0 struct array with the right fields

% ---- 1. contrast (temporal) --------------------------------------------------
s = base();
s.id='contrast'; s.label='Contrast'; s.wrapper=@runContrastFromRLS;
% s.contrastType chooses the flag: temporal -> _t_K, spatial -> _s_K (same branch)
s.inGlob='*.rls'; s.outSuffix={'_t_K_d','_t_K_r','_t_K_s'}; s.outKind='new';
s.outTransform=struct('from','.rls','to','_t_K_d.mat');   % spatial: '_s_K_d'
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
    'contrastKernel','25 for temporal, 5-7 for spatial', ...
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
s.outTransform=struct('from','.rls','to','_c_K_d.mat');
s.gatingField='runInternalCycle'; s.requires={}; s.produces={'cardiac'};
s.artifacts={'_rep_cycle-detect.jpg','_rep_cycle-average.jpg'};
s.legacyArtifacts={'_ic1.jpg','_ic2.jpg'}; s.branch='cardiac';
s.settingGroups={ 'Contrast calculation',{'trustLimitsK','trustLimitsI','contrastKernelS'};
                  'Frequency band',{'maxFrqIni','minFrqIni','rangeFrq'};
                  'Exclusion criteria',{'excludeFirstNCycles','coeffsSTD','coeffsRel','coeffsAbs'};
                  'Cycle calculation',{'method','decimationSpace','framesToAverage', ...
                     'contrastKernelT','contrastKernelPreproc','interpFactor','smoothCoef1','minPromCoef'} };
s.basicFields={'method','maxFrqIni','minFrqIni','excludeFirstNCycles'};
s.sharedKeys={'trustLimitsK','trustLimitsI','libraryFolder'};
s.enums=struct('method',{{'sLSCIMM','tLSCIMM','ltLSCIMM','sLSCIMMM'}});
s.presets=struct('default',struct('trustLimitsK',[0.01 0.3],'trustLimitsI',[5 250], ...
    'contrastKernelS',5,'maxFrqIni',20,'minFrqIni',1,'rangeFrq',1,'excludeFirstNCycles',0, ...
    'coeffsSTD',[3 2 2 2 2 3 3 2 2],'coeffsRel',[0.5 0.1],'coeffsAbs',2,'method','sLSCIMM', ...
    'decimationSpace',4,'framesToAverage',1,'contrastKernelT',25,'contrastKernelPreproc',5, ...
    'interpFactor',10,'smoothCoef1',1/3,'minPromCoef',1/4));
s.tips=struct('method','sLSCIMM recommended; ltLSCIMM for high-quality data', ...
    'maxFrqIni','initial max frequency of interest, Hz', ...
    'minFrqIni','initial min frequency of interest, Hz', ...
    'coeffsSTD','pulse-rejection coefficients relative to feature std');
reg(end+1)=s;

% ---- 3. externalCycle (NVC) -------------------------------------------------
s = base();
s.id='externalCycle'; s.label='External cycle'; s.wrapper=@runExternalCycle;
% external/epoch cycle of the contrast side (t or s): the stage flag BECOMES e,
% i.e. the single suffix _e_K (the contrast base is kept in the settings, not the
% name - see wbFileModel rationale)
s.inGlob='*_K_d.mat'; s.outSuffix={'_e_K_d','_e_K_r','_e_K_s'}; s.outKind='new';
s.outTransform=struct('from','_t_K_d.mat','to','_e_K_d.mat');   % replaces the t|s flag with e
s.gatingField='externalCycle';                       % REAL field (differs from fn name)
s.requires={'contrast'}; s.produces={'epochAvg'};
s.interactive=@(ss) isfield(ss,'enablelRejectionModification') && isequal(ss.enablelRejectionModification,1);
s.artifacts={'_rep_epochs.jpg','_rep_epoch-average.jpg'};
s.legacyArtifacts={'_ec.jpg','_ec2.jpg'}; s.branch='contrast';
s.settingGroups={ 'Stimulation',{'stimStartType','stimOffset','epochsN','epochDurationSec', ...
                     'epochBaselineSec','epochStimStartSec','epochFinaleSec'};
                  'Masking',{'maskType','enablelRejectionModification'};
                  'Rejection',{'rejectBlCoef','rejectEpochCoef','rejectFinCoef','rejectPeakCoef', ...
                     'rejectBlSimCoef','rejectSimCoef','rejectTimeLoss','rejectFirstEpoch'} };
s.basicFields={'stimStartType','stimOffset','epochsN','epochDurationSec','epochBaselineSec','epochStimStartSec','epochFinaleSec','maskType','enablelRejectionModification'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('stimStartType',{{'offset','manual'}},'maskType',{{'basic','cMask','selection'}});
s.presets=struct('default',struct('stimStartType','offset','stimOffset',0,'epochsN',20, ...
    'epochDurationSec',30,'epochBaselineSec',[0 10],'epochStimStartSec',10,'epochFinaleSec',[-5 0], ...
    'maskType','cMask','enablelRejectionModification',1,'rejectBlCoef',1,'rejectEpochCoef',1, ...
    'rejectFinCoef',1,'rejectPeakCoef',1,'rejectBlSimCoef',1,'rejectSimCoef',1, ...
    'rejectTimeLoss',0.5,'rejectFirstEpoch',1));
s.tips=struct('maskType','''cMask'' reads the segmentation mask; ''basic'' the contrast mask', ...
    'enablelRejectionModification','1 = interactive epoch-rejection editing', ...
    'stimStartType','''offset'' (fixed delay) or ''manual'' (timestamp list)');
reg(end+1)=s;

% ---- 4. setRegions ----------------------------------------------------------
s = base();
s.id='setRegions'; s.label='Regions'; s.wrapper=@setRegions;
s.inGlob='*_K_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setRegions'; s.requires={}; s.requiresAny={'contrast','internalCycle'};
s.produces={'regionsMask'};
s.artifacts={'_rep_regions.jpg'};                     % NEW in Session 3: the drawn
                                                      % regions, drawn or copied - no
                                                      % legacy tail, it wrote none
s.interactive=true;                                   % always opens the ROI editor
s.branch=''; s.branchScope='copy';                    % draw on _t, inherit onto _c/_e
s.fanOut='animal';                                    % rows = animals: ROIs carry along a row
s.settingGroups={};                                   % fully interactive, no numeric params
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct());
reg(end+1)=s;

% ---- 5. splitRegions --------------------------------------------------------
s = base();
s.id='splitRegions'; s.label='Split regions'; s.wrapper=@splitRegions;
s.inGlob='*_K_d.mat'; s.outSuffix={}; s.outKind='prefix';
s.outTransform=struct('from','','to','Roi');          % RoiN_ prefix, N per region
s.gatingField='splitRegions'; s.requires={'setRegions'}; s.produces={'regionCrops'};
s.branch=''; s.branchScope='all';                     % crop every branch of the recording
s.settingGroups={ 'Options',{'deleteOriginal'} };
s.basicFields={'deleteOriginal'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('deleteOriginal',false));
s.tips=struct('deleteOriginal','true only if you will not re-define regions');
reg(end+1)=s;

% ---- 6. segmentation --------------------------------------------------------
s = base();
s.id='segmentation'; s.label='Segmentation'; s.wrapper=@runSegmentation;
s.inGlob='*_K_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runSegmentation'; s.requires={}; s.requiresAny={'contrast','internalCycle'};
s.produces={'segmentation'};
s.artifacts={'_rep_categories.jpg','_rep_segments.jpg'};
s.legacyArtifacts={'_cm.jpg','_vs.jpg'};              % '_vs' came from showSegmentsPreview
s.branch=''; s.branchScope='copy';                % computed on contrast side, copied to cardiac
s.settingGroups={ 'Contrast',{'trustLimitsK'};
                  'Categorization',{'lSizeN','sSizeN','sens','sSizeScale','deSens','lThinN','imOpen','iEdge','eEdge'};
                  'Labelling & traces',{'sStat','sMinL','prchNSize','correctNodes','simR','difR'} };
                  % no 'Copy to siblings' panel: branchScope 'copy' means wbExecutor
                  % derives s.fNamesCopyTo from the recording's own branch products
s.basicFields={'lSizeN','sSizeN','sens','sMinL'};
s.sharedKeys={'trustLimitsK','prchNSize','sMinL','correctNodes','simR','difR','sStat','libraryFolder'};
s.enums=struct('sStat',{{'median','mean'}});
s.presets=struct('default',struct('trustLimitsK',[0.001 0.99],'lSizeN',141,'sSizeN',9, ...
    'sens',0.1,'sSizeScale',1,'deSens',1,'lThinN',2,'imOpen',0,'iEdge',2,'eEdge',2, ...
    'categories',{{'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}}, ...
    'sStat','median','sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3,'difR',0.4));
s.tips=struct('lSizeN','odd, ~2x the largest vessel', ...
    'sSizeN','odd, ~2x the small-vessel diameter', ...
    'sens','segmentation sensitivity (raise to catch faint vessels)', ...
    'sMinL','minimum segment length', ...
    'prchNSize','parenchymal pixel neighbourhood');
reg(end+1)=s;

% ---- 7. dynamicSegmentation -------------------------------------------------
s = base();
s.id='dynamicSegmentation'; s.label='Dynamic segmentation'; s.wrapper=@runDynamicSegmentation;
s.inGlob='*_K_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runDynamicSegmentation'; s.requires={'segmentation'}; s.produces={'dynamicSeg'};
s.artifacts={'_rep_segments.jpg'};                    % re-emitted, overwriting the static one
s.legacyArtifacts={'_vs.jpg'}; s.branch=''; s.branchScope='all';
s.settingGroups={ 'Labelling (match segmentation)',{'sMinL','prchNSize','correctNodes','simR','difR'};
                  'Dynamic segmentation',{'sMinP2R2','sMaxLBI','sMaxCLR','sMaxKK','iniNSize','sMaxP2D'};
                  'Quality & interpolation',{'gSizeN','minOverlapMask','minOverlapSelf','pInterpF'} };
s.basicFields={'sMinL','minOverlapMask'};
s.sharedKeys={'sMinL','prchNSize','correctNodes','simR','difR','libraryFolder'};
s.presets=struct('default',struct('sMinL',15,'prchNSize',50,'correctNodes',true,'simR',0.3, ...
    'difR',0.4,'sMinP2R2',0.95,'sMaxLBI',(1/7)/15,'sMaxCLR',1.3,'sMaxKK',0.3,'iniNSize',7, ...
    'sMaxP2D',3,'gSizeN',3,'minOverlapMask',0.6,'minOverlapSelf',0.2,'pInterpF',4));
s.tips=struct('sMinP2R2','min accepted R^2 of the 3-degree polynomial fit', ...
    'sMaxCLR','max chord-length ratio (1 straight, 2 coiled)', ...
    'minOverlapMask','min overlap of centre line with the per-frame mask');
reg(end+1)=s;

% ---- 8. guided --------------------------------------------------------------
s = base();
s.id='guided'; s.label='Guided'; s.wrapper=@runGuidedContrast;
s.inGlob='*_K_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runGuidedContrast'; s.requires={'segmentation'}; s.produces={'guidedTraces'};
s.needsRaw=true; s.branch='contrast';
s.settingGroups={};                                   % uses sMap + raw; no tunable params
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct());
reg(end+1)=s;

% ---- 9. registration (per animal) --------------------------------------------
s = base();
s.id='registration'; s.label='Registration'; s.wrapper=@runRegistration;
s.arity='perAnimal';                                   % 2-D fNames row; col 1 = template
s.inGlob='*_K_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runRegistration';                      % code writes runRegistration (header drift A4)
s.requires={}; s.requiresAny={'contrast','internalCycle'}; s.produces={'registered'};
s.interactive=@(ss) ~(isfield(ss,'silent') && isscalar(ss.silent) && isequal(ss.silent,true));
s.artifacts={'_rep_registration.jpg'};                % was a .png; now a report page
s.legacyArtifacts={'_registration.png'};
s.branch=''; s.refBranch='contrast';                  % template = the reference's _t|_s
s.branchScope='all';                                   % EVERY product of every member (launcher: 'Roi*_K_d.mat')
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

% ---- 10. BFI ----------------------------------------------------------------
s = base();
s.id='BFI'; s.label='BFI'; s.wrapper=@runBFI;
s.inGlob='*_K_d.mat'; s.outSuffix={'_BFI_d','_BFI_r','_BFI_s'}; s.outKind='new';
s.outTransform=struct('from','_K_d.mat','to','_BFI_d.mat');   % branch-preserving strrep
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

% ---- 11. vasomotion ---------------------------------------------------------
s = base();
s.id='vasomotion'; s.label='Vasomotion'; s.wrapper=@runVasomotion;
s.inGlob='*_BFI_d.mat'; s.outSuffix={}; s.outKind='inplace';   % contrast-side BFI (_t_BFI or _s_BFI)
s.gatingField='runVasomotion'; s.requires={'BFI'}; s.produces={'vasomotion'};
s.branch='contrast';
s.settingGroups={ 'Bands',{'vFR','cFR','wFR','wVPO'};
                  'Normalisation',{'normalisation','normsize','tgtFS'};
                  'Peaks & percentiles',{'pcts','otsuMaxN','otsuElbow','nPeakProm'};
                  'Signals & levels',{'vsmSignals','segVsmReturn','ppxVsmReturn'} };
s.basicFields={'vFR','cFR','tgtFS','vsmSignals'};
s.sharedKeys={'libraryFolder'};
s.enums=struct('normalisation',{{'mean','median','mmean','mmedian'}});
s.presets=struct('default',struct('vFR',[0.05 0.25],'cFR',[0.4 0.6],'wFR',[0.01 1],'wVPO',10, ...
    'normalisation','median','normsize',101,'tgtFS',1,'pcts',0:10:100,'otsuMaxN',5, ...
    'otsuElbow',0.05,'nPeakProm',0.10, ...
    'vsmSignals',{{'sData','dvsData','dvsDiameter','gsData'}}, ...
    'segVsmReturn',{{'bands','moments','series','clustering','spectrum'}},'ppxVsmReturn',[]));
s.tips=struct('vFR','vasomotion frequency band [lo hi], Hz', ...
    'cFR','control (cardiac) frequency band [lo hi], Hz', ...
    'nPeakProm','VB peak-count prominence as a fraction of the band range', ...
    'segVsmReturn','which per-segment levels to store (set of tokens)', ...
    'ppxVsmReturn','per-pixel analysis ([] = off; e.g. {''bands''} to enable)');
reg(end+1)=s;

% ---- 12. pulsatility --------------------------------------------------------
s = base();
s.id='pulsatility'; s.label='Pulsatility'; s.wrapper=@runPulsatility;
s.inGlob='*_c_BFI_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='runPulsatility'; s.requires={'BFI','internalCycle'}; s.produces={'pulsatility'};
s.branch='cardiac';
s.settingGroups={ 'Harmonic model',{'nHarm'};
                  'Analysis levels',{'segPulsReturn','ppxPulsReturn'} };
s.basicFields={'nHarm'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('nHarm',5, ...
    'segPulsReturn',{{'markers','model','reconstruction'}},'ppxPulsReturn',{{'markers'}}));
s.tips=struct('nHarm','number of harmonics in the sinusoidal cardiac model', ...
    'segPulsReturn','per-segment levels: markers / model / reconstruction', ...
    'ppxPulsReturn','per-pixel maps ([] = off; non-empty enables marker maps)');
reg(end+1)=s;

% ---- 13. vesselTypes (per animal) --------------------------------------------
s = base();
s.id='vesselTypes'; s.label='Vessel types'; s.wrapper=@setVesselTypes;
s.arity='perAnimal';                                   % 2-D fNames row; col 1 = reference
s.inGlob='*_BFI_d.mat'; s.outSuffix={}; s.outKind='inplace';
s.gatingField='setVesselTypes'; s.requires={'BFI','segmentation'}; s.produces={'vesselTypes'};
s.artifacts={'_rep_vesseltypes.jpg'};                 % NEW in Session 3: the artery/vein
                                                      % map, painted or inherited - no
                                                      % legacy tail, it wrote none
s.interactive=true;                                   % paint GUI (per-file skipped under useReference)
s.branch=''; s.refBranch='cardiac';                   % paint target = the reference's _c
s.branchScope='all';                                  % EVERY product of every member (launcher: 'Roi*_BFI_d.mat')
s.settingGroups={ 'Reference',{'useReference'} };
s.basicFields={'useReference'};
s.sharedKeys={'libraryFolder'};
s.presets=struct('default',struct('useReference',false));
s.tips=struct('useReference','true = paint the first (reference) file only, inherit to the rest');
reg(end+1)=s;

% ---- 14. vascularTree -------------------------------------------------------
s = base();
s.id='vascularTree'; s.label='Vascular tree'; s.wrapper=@setVascularTree;
s.inGlob='*_c_BFI_d.mat'; s.outSuffix={}; s.outKind='inplace';
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

% ---- optional modality filter ----------------------------------------------
if nargin>=1 && ~isempty(modality)
    keep = arrayfun(@(st) any(strcmp(modality,st.modalities)), reg);
    reg = reg(keep);
end
end

% =====================================================================
function s = base()
%base  A fully-formed step struct with every field defaulted, so appending
%   elements to the registry array never triggers a field-mismatch error.
s = struct( ...
    'id','', 'label','', 'wrapper',[], 'arity','perFile', ...
    'inGlob','', 'outSuffix',{{}}, 'outKind','inplace', 'outTransform',[], ...
    'gatingField','', 'requires',{{}}, 'requiresAny',{{}}, 'produces',{{}}, ...
    'interactive',false, 'needsRaw',false, 'artifacts',{{}}, 'legacyArtifacts',{{}}, ...
    'settingGroups',{{}}, 'basicFields',{{}}, 'sharedKeys',{{}}, ...
    'presets',struct('default',struct()), 'tips',struct(), 'enums',struct(), ...
    'modalities',{{'LSCI'}}, 'branch','', 'branchScope','one', 'fanOut','flat', ...
    'refBranch','');
end
