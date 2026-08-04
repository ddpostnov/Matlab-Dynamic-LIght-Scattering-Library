%runFitVasoreactivity  Drug / gas-challenge vasoreactivity fit into results.vasoreactivity
%
%   runFitVasoreactivity(s,fNames) loads every contrast-branch *_t_BFI_r.mat (or
%   *_s_BFI_r.mat) whole-recording product in fNames (runContrastFromRLS ->
%   runSegmentation -> runBFI: results.sData/dvsData/dvsDiameter [nT x nSeg] on
%   results.time, plus the flow cube source.data [Y x X x nT]) and reduces every trace
%   to model-free response markers and a fitted vasoreactivity response.
%
%   THIS STEP MODELS PHARMACOKINETICS DRIVING A QUASI-STATIC VASCULAR GAIN, not
%   vascular dynamics.  A rate constant from here is a drug rate, not a vessel time
%   constant, and nothing it reports is comparable with anything runFitNVC reports -
%   the two timescales differ by two to three orders of magnitude.  fitVasoreactivity's
%   header carries the full argument; read it before reading a number out of this tree.
%
%   The per-trace science lives in the shared core fitVasoreactivity; this producer is
%   the multi-file GLUE: it resolves the per-file injection time, builds the fit layout
%   once per recording, runs the core over every segment (parfor) and every masked
%   pixel, and assembles the results into the RESULTS.vasoreactivity tree (mirroring
%   RESULTS.nvc, RESULTS.pulsatility and RESULTS.vasomotion).
%
%   ONE CURVE PER RECORDING, ANCHORED PER FILE.  A drug challenge is not an epoch
%   average: it is one continuous trace with one injection, so s.injectionSec is a
%   PER-FILE LIST indexed by position in fNames (a scalar means the same time for every
%   file).  That is why the registry declares fileOrder 'ordered' - cut the list and
%   the files up differently and every fit silently anchors to the wrong time.
%
%   OUTPUT TREE  (RESULTS.vasoreactivity)
%     Shared root axes - the DECIMATED base every scalar and every fitted curve below
%     lives on, plus the geometry, so a downstream plot is self-describing without
%     re-deriving anything:
%       .time       [nT x 1]  decimated recording clock, seconds, single
%       .fs         achieved sampling rate after decimation, Hz, single
%       .injection  the injection time on that clock, s, single
%       .baseline   [1 x 2]   the baseline window on that clock, s, single
%       .map        [nT x 1]  the measured arterial-pressure covariate resampled onto
%                             that clock, single - ONLY when s.mapTrace was supplied
%     One sub-tree per analysed signal, all floats SINGLE, with a field-name PREFIX
%     encoding the physical quantity (rs = reactivity flow, rd = reactivity diameter):
%       .sData        { scalars(rs*), timeVectors.fData }   (-> results.sMetrics)
%       .dvsData      { scalars(rs*), timeVectors.fData }   (-> results.dvsMetrics)
%       .dvsDiameter  { scalars(rd*), timeVectors.fData }   (-> results.dvsMetrics)
%       .ppx          { scalars(rs*) [Y x X] maps, timeVectors.fData [Y x X x nT] }
%     scalars.* (per segment [nSeg x 1], per pixel [Y x X]) carry the core's scalar
%     names with the prefix attached - rsPeak, rsTMax, rsCVR, rsSteal, rdPeak, ... -
%     and the core owns which names exist at which level (see fitVasoreactivity).  An
%     INVALID trace (non-finite) is NaN in EVERY field.
%
%   THE MEASURED PRESSURE COVARIATE IS STORED, NOT REGRESSED.  Drugs change mean
%   arterial pressure, heart rate and PaCO2, and cerebral autoregulation then acts to
%   OPPOSE the resulting flow change - a confound that is negligible in neurovascular
%   coupling and potentially dominant here.  s.mapTrace / s.mapTime (per-file, on their
%   own clock) are resampled onto the decimated base and kept beside the fit and on the
%   report page, so the confound is VISIBLE in the product.  Nothing subtracts it: with
%   'bateman' there is no input channel to put it in, and a covariate quietly removed
%   is worse than one plainly displayed.
%
%   LEVELS  (selector cell arrays; a subset of {'markers','model','reconstruction'})
%     s.segVrcReturn  per-segment levels; default {'markers','model','reconstruction'}.
%     s.ppxVrcReturn  per-pixel levels + GATE: NON-EMPTY runs the per-pixel path,
%                     [] skips it; DEFAULT [] - OFF.  This differs from runFitNVC,
%                     whose epoch cube is 40 s long; a drug-challenge cube is
%                     45 minutes, so enabling it is a deliberate act.
%     THE LEVELS GATE THE TREE, NOT THE METRICS TABLE.  The model-free markers cost no
%     fit, so they are always computed and always duplicated into the metrics tables;
%     dropping 'markers' drops them from the tree only.  Same rule as runFitNVC and
%     runPulsatility, and it keeps a metrics column from appearing and disappearing
%     with a display setting.
%
%   METRICS-TABLE DUPLICATION (a small key set, same names as the tree, single)
%     results.sMetrics    (from sData)      : rsPeak rsTMax rsTWash rsCVR rsSteal rsR2
%     results.dvsMetrics  (from dvsData)    : the same six rs*
%                         (from dvsDiameter): the six rd*
%     CVR and R2 come from the fit, so they appear only when 'model' was computed.
%     rsSteal IS IN THE TABLE ON PURPOSE: a fit with the steal flag raised is a fit not
%     to trust, and a trust flag that only lives in the tree is one nobody reads.
%     BFI / std(BFI) are runBFI's columns and are LEFT UNTOUCHED.
%
%   INPUTS
%     s        parameter struct with fields
%                • injectionSec  REQUIRED.  When the drug was given, seconds from the
%                                recording start.  Scalar, or one entry per file.
%                • baselineSec   REQUIRED.  [t1 t2] pre-injection baseline window,
%                                seconds; must end before every file's injection.
%                • vrcModel      'bateman' (default; the only model in v1)
%                • tgtFS         rate the traces are decimated to before fitting, Hz
%                                (default 1; [] keeps the recording's own rate)
%                • vrcAucSec     AUC horizon after the injection, s (default 2700)
%                • vrcStealK     steal threshold in baseline SDs (default 2)
%                • vrcStealFrac  steal threshold floor, as a fraction of the peak
%                                (default 0.15); the larger of the two is used
%                • vrcStealSec   how long it must stay below, s (default 60)
%                • vrcWashFrac   washout level as a fraction of the peak (default 0.1)
%                • vrcStarts     multi-start start points (default 16)
%                • segVrcReturn  per-segment level cell (default all three)
%                • ppxVrcReturn  per-pixel level cell / gate (default [] = off)
%                • mapTrace / mapTime  (optional) measured arterial pressure and its
%                                own clock, per file, stored beside the fit
%                • parforVrcSegments logical, default true
%                • parforVrcPixels   logical, default true - both are WORKER BOUNDS,
%                                not branches: false runs the identical loop body
%                                serially in the client and starts no pool.
%     fNames   cell array of *_t_BFI_r.mat (or *_s_BFI_r.mat) paths, IN THE ORDER
%              s.injectionSec is written in.
%                • Optional workbench hooks in s (no-op when absent):
%                  s.stageFcn(stage,detail), s.cancelFcn()->tf.  Cancel is checked
%                  between files (never inside a parfor).
%
%   OUTPUT FILES (side-effects) - NON-DESTRUCTIVE: the flow cube is preserved
%       *_t_BFI_d.mat   SOURCE   - NEVER re-saved (source.data read, never modified)
%       *_t_BFI_r.mat   RESULTS  - RESULTS.vasoreactivity.* + the rs*/rd* metrics columns
%       *_t_BFI_s.mat   SETTINGS - field settings.runFitVasoreactivity added
%       *_rep_vasoreactivityfit.jpg - the mean trace with the fit, the injection mark,
%                                   the baseline window, the residuals and the
%                                   parameters, with the steal flag called out
%
%   EXAMPLE
%     s.injectionSec=[600 615 590];      % one per recording, in the fNames order
%     s.baselineSec=[60 540];
%     s.tgtFS=1;
%     s.segVrcReturn={'markers','model','reconstruction'};
%     fNames = getFileNamesList(resultsFolder,'*_t_BFI_r.mat');
%     runFitVasoreactivity(s, fNames(:));
%
%   DEPENDS ON
%     Core/Vasoreactivity/fitVasoreactivity (shared response core), MATLAB
%     Optimization Toolbox (lsqcurvefit); core LSCI library utilities.
%
% See also: fitVasoreactivity, runFitNVC, runVasomotion, runPulsatility, getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - When the drug was given
% s.injectionSec=600;     % seconds from the start of the recording.  One number per
%                         % recording, in the same order as the file list - a single
%                         % number is used for all of them.
% s.baselineSec=[60,540]; % the quiet stretch before the injection the response is
%                         % measured against.  It has to end before the injection.
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Response model
% s.vrcModel='bateman';   % a rising and a falling phase on a drifting baseline; the
%                         % reported amplitude is the peak change in the trace's units.
% s.tgtFS=1;              % Hz.  The traces are averaged down to this before fitting -
%                         % the response takes minutes, so 1 Hz is plenty and 45
%                         % minutes of raw frames would not fit in memory.
% %ADJUSTED IF NECESSARY - What counts as a response
% s.vrcAucSec=2700;       % how long after the injection the area is integrated over, s
% s.vrcStealK=2;          % how far below baseline, in baseline noise levels, counts
% s.vrcStealFrac=0.15;    % ...or this fraction of the response, whichever is larger
% s.vrcStealSec=60;       % ...and for how long, before it is flagged as a steal
% s.vrcWashFrac=0.1;      % washout time = when the response is back to this fraction
%                         % of its peak
% %ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
% s.segVrcReturn={'markers','model','reconstruction'};  % markers (model-free scalars),
%                         % model (fitted parameters - runs the fit), reconstruction
%                         % (the fitted curve). Default = all three.
% %ADJUSTED IF NECESSARY - Per-pixel maps/cubes (GATE + selector)
% s.ppxVrcReturn=[];      % [] = off (the default).  NON-EMPTY fits every pixel of the
%                         % field over the whole recording - expect many minutes.
% %ADJUSTED IF NECESSARY - Fitting
% s.vrcStarts=16;         % how many start points the fit is tried from.  Fewer is
%                         % faster and more likely to settle on a wrong answer.
% %ADJUSTED IF NECESSARY - Parallel execution
% s.parforVrcSegments=true;
% s.parforVrcPixels=true; % false fits one at a time in this MATLAB and starts no
%                         % parallel pool - slower, but it costs no worker processes.

%------------- BEGIN CODE --------------
function runFitVasoreactivity(s,fNames)

%THE INPUT IS A WHOLE RECORDING ON THE CONTRAST BRANCH.  An epoch average ('_e_') and a
%cardiac average ('_c_') both collapse the recording clock, so there is no time axis
%left to place an injection on - the same argument runExternalCycle makes for refusing
%the cardiac product.  The guard is written positively, naming the two stages that DO
%carry a recording clock.
if ~all( cellfun(@(x) isempty(x) || contains(x,'_t_BFI_r.mat') || ...
        contains(x,'_s_BFI_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries are not a whole-recording BFI product.  ' ...
        'This step places a drug injection on the RECORDING clock, so it needs the ' ...
        'whole trace: an epoch average ("_e_") and a cardiac average ("_c_") carry ' ...
        'no recording clock at all.  List the right files with ' ...
        'getFileNamesList(resultsFolder,''*_t_BFI_r.mat'').']);
end

%resolve defaults so they are recorded in the saved settings (the core defaults the
%same values internally).  segVrcReturn defaults when absent OR empty; ppxVrcReturn
%defaults ONLY when the field is ABSENT - an explicit [] must stay empty, exactly like
%runPulsatility's s.ppxPulsReturn.  ITS DEFAULT IS [] EITHER WAY here: a 45-minute
%per-pixel cube is not something to switch on by omission.
if ~isfield(s,'vrcModel') || isempty(s.vrcModel),      s.vrcModel='bateman';   end
if ~isfield(s,'tgtFS'),                                s.tgtFS=1;              end
if ~isfield(s,'vrcAucSec')   || isempty(s.vrcAucSec),  s.vrcAucSec=2700;       end
if ~isfield(s,'vrcStealK')   || isempty(s.vrcStealK),  s.vrcStealK=2;          end
if ~isfield(s,'vrcStealFrac')|| isempty(s.vrcStealFrac),s.vrcStealFrac=0.15;   end
if ~isfield(s,'vrcStealSec') || isempty(s.vrcStealSec),s.vrcStealSec=60;       end
if ~isfield(s,'vrcWashFrac') || isempty(s.vrcWashFrac),s.vrcWashFrac=0.1;      end
if ~isfield(s,'vrcStarts')   || isempty(s.vrcStarts),  s.vrcStarts=16;         end
if ~isfield(s,'segVrcReturn') || isempty(s.segVrcReturn)
    s.segVrcReturn={'markers','model','reconstruction'};
end
if ~isfield(s,'ppxVrcReturn'), s.ppxVrcReturn=[]; end
%PARALLELISM IS OPTIONAL.  Each switch is a BOUND on its parfor (Inf workers or 0),
%never a branch: parfor(...,0) runs the identical loop body serially IN THE CLIENT and
%starts no pool at all.
if ~isfield(s,'parforVrcSegments') || isempty(s.parforVrcSegments)
    s.parforVrcSegments=true;
end
if ~isfield(s,'parforVrcPixels') || isempty(s.parforVrcPixels)
    s.parforVrcPixels=true;
end

%the two REQUIRED settings are checked ONCE, before any file is opened, so a mistyped
%injection list is a message rather than a partly-processed data set
requireInjection(s,numel(fNames));

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  Cancel is only
% checked between files (a hook inside a parfor would broadcast oddly).
rep=reportOpen(s,'Vasoreactivity fit',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        %SOURCE is read only for the per-pixel path (and never written back).
        if ~isempty(s.ppxVrcReturn)
            load(getProductPath(s.fName,'d'),'source')
        end
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        if ~isfield(results,'time') || isempty(results.time)
            error('runFitVasoreactivity:noTime', ...
                '%s carries no results.time, so there is no recording clock to fit on.', ...
                s.fName);
        end

        %THIS file's injection time (see the header on why the list is per file)
        sf=s;
        sf.injectionSec=perFile(s.injectionSec,fidx);

        %runFitVasoreactivity owns RESULTS.vasoreactivity entirely; rebuild it from
        %scratch so no stale sub-branch survives a re-run with different levels.
        if isfield(results,'vasoreactivity')
            results=rmfield(results,'vasoreactivity');
        end

        %TWO LAYOUTS, ONE TIME BASE.  layoutCalc is what the core is actually called
        %with - the requested levels PLUS 'markers', because the markers cost no fit
        %and the metrics tables always carry them.  layoutTree is the user's own
        %selection, and its scalarNames (a subsequence of layoutCalc's) is exactly what
        %enters the saved tree.  Neither the loop nor the assembly below has to know a
        %single field name.
        sCalc=sf; sCalc.segVrcReturn=withMarkers(sf.segVrcReturn);
        layoutCalc=fitVasoreactivity(results.time,sCalc);
        layoutTree=fitVasoreactivity(results.time,sf);
        names=layoutCalc.scalarNames;
        nT=layoutCalc.nT;

        %shared root axes (single; the tree's OWN copy - results.time stays double).
        %The DECIMATED base is the one every scalar and every fitted curve below lives
        %on, so it is the one the tree publishes.
        results.vasoreactivity.time     =single(layoutCalc.time(:));
        results.vasoreactivity.fs       =single(layoutCalc.fs);
        results.vasoreactivity.injection=single(layoutCalc.tInj);
        results.vasoreactivity.baseline =single(layoutCalc.t1+layoutCalc.bl);
        mapTrace=resampleCovariate(sf,fidx,layoutCalc);
        if ~isempty(mapTrace)
            results.vasoreactivity.map=single(mapTrace);
        end

        % =============================================================
        % Per-segment analysis: sData/dvsData (rs prefix), dvsDiameter (rd prefix).
        % =============================================================
        nwSeg=0; if s.parforVrcSegments, nwSeg=Inf; end   %worker bound, not a branch
        nRecSeg=1; if layoutTree.want.reconstruction, nRecSeg=nT; end

        sigNames={'sData','dvsData','dvsDiameter'};
        for kSig=1:numel(sigNames)
            sigName=sigNames{kSig};
            if ~isfield(results,sigName), continue; end
            if strcmp(sigName,'dvsDiameter'), prefix='rd'; else, prefix='rs'; end

            sigMat=results.(sigName);
            if isempty(sigMat), continue; end
            nSeg=size(sigMat,2);

            %ONE ROW PER SEGMENT, ONE COLUMN PER SCALAR.  The core emits some thirty
            %scalars, so named accumulators would be unreadable; the row/column form
            %also lets the tree and the metrics duplication below be driven by the name
            %list rather than typed out twice.  BOTH slices are assigned on EVERY
            %iteration, so nothing about the parfor depends on which branch the core
            %took.  The traces go in at the RECORDING's rate - the core decimates.
            M=nan(nSeg,numel(names),'single');
            F=nan(nRecSeg,nSeg,'single');
            %the workbench hooks are function handles bound to a GUI; they are
            %transport rather than parameters, the core never reads them, and
            %broadcasting them to workers is the one thing in s that can fail to
            %serialise.  reportSettings strips exactly those two, here as at the save.
            sPar=reportSettings(sCalc);
            parfor (i=1:nSeg, nwSeg)
                [M(i,:),F(:,i)]=vrcTraceRow(sigMat(:,i),layoutCalc,sPar,names,nRecSeg);
            end

            %assemble the sub-tree from the USER's level selection
            T=struct();
            for k=1:numel(layoutTree.scalarNames)
                nm=layoutTree.scalarNames{k};
                T.scalars.([prefix nm])=M(:,strcmp(names,nm));
            end
            if layoutTree.want.reconstruction
                T.timeVectors.fData=F;
            end
            results.vasoreactivity.(sigName)=T;

            %duplicate the key set into the metrics tables (DATA-MODEL section 11).
            %Row order matches the 1:nSeg loop (segment invariant).  BFI/std(BFI) are
            %runBFI's columns and are intentionally NOT written here.
            switch sigName
                case 'sData',                    tbl='sMetrics';
                case {'dvsData','dvsDiameter'},  tbl='dvsMetrics';
            end
            dup=intersect({'Peak','TMax','TWash','CVR','Steal','R2'},names,'stable');
            for k=1:numel(dup)
                results.(tbl).([prefix dup{k}])=M(:,strcmp(names,dup{k}));
            end
        end

        % =============================================================
        % Per-pixel analysis (gated by non-empty s.ppxVrcReturn; rs prefix, since
        % source.data is a flow cube).  Every masked pixel goes through the same core
        % the segments did, so a per-pixel map and a per-segment scalar cannot drift.
        % =============================================================
        if ~isempty(s.ppxVrcReturn)
            if ~isfield(results,'sMap')
                error('runFitVasoreactivity:noSMap', ...
                    ['results.sMap not found; runSegmentation must be run before ' ...
                     'per-pixel vasoreactivity.']);
            end
            sPix=sf; sPix.segVrcReturn=s.ppxVrcReturn;
            sPixCalc=sPix; sPixCalc.segVrcReturn=withMarkers(s.ppxVrcReturn);
            layoutPixCalc=fitVasoreactivity(results.time,sPixCalc);
            layoutPixTree=fitVasoreactivity(results.time,sPix);
            namesPix=layoutPixCalc.scalarNames;

            sz=size(source.data); Y=sz(1); X=sz(2); npx=Y*X;
            Dpix=reshape(source.data,npx,sz(3));     %[npx x nT] pixel timecourses
            sMapLin=results.sMap(:);

            %NaN-preallocated, so a background pixel - which no iteration writes -
            %stays NaN in every map.  An unfitted pixel has no response, and a zero
            %there would be a measurement that was never made.
            Mp=nan(npx,numel(namesPix),'single');
            nRecPix=1; if layoutPixTree.want.reconstruction, nRecPix=layoutPixCalc.nT; end
            Fp=nan(npx,nRecPix,'single');

            nwPix=0; if s.parforVrcPixels, nwPix=Inf; end   %worker bound, not a branch
            sParPix=reportSettings(sPixCalc);               %see the note above
            parfor (p=1:npx, nwPix)
                if sMapLin(p)==0, continue; end              %background stays NaN
                [Mp(p,:),Fp(p,:)]=vrcTraceRow(Dpix(p,:).',layoutPixCalc,sParPix, ...
                    namesPix,nRecPix);
            end

            ppx=struct();
            for k=1:numel(layoutPixTree.scalarNames)
                nm=layoutPixTree.scalarNames{k};
                ppx.scalars.(['rs' nm])=reshape(Mp(:,strcmp(namesPix,nm)),Y,X);
            end
            if layoutPixTree.want.reconstruction
                ppx.timeVectors.fData=reshape(Fp,Y,X,nRecPix);
            end
            results.vasoreactivity.ppx=ppx;
        end

        % =============================================================
        % Report page: the recording's mean trace, its fit, and the parameters.
        % =============================================================
        drawVrcPage(rep,results,sf,mapTrace);

        settings.runFitVasoreactivity=reportSettings(sf);
        reportWriting(rep);
        %NON-DESTRUCTIVE: SOURCE (_d) is never re-saved - only RESULTS and SETTINGS.
        save(s.fName,'results','-v7.3','-nocompression');
        save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);

end

% =====================================================================
function [row,fdata]=vrcTraceRow(y,layout,s,names,nRec)
%vrcTraceRow  One trace through the core, flattened to a fixed-width row.
%   The row carries the layout's scalars IN ORDER, so the caller indexes by name once
%   (when it assembles the tree) instead of once per segment.  nRec is the decimated
%   length when the reconstruction was asked for and 1 otherwise, which keeps the
%   parfor's second sliced output the same shape on every iteration without allocating
%   a full cube nobody wants.
m=fitVasoreactivity(y,layout,s);
row=nan(1,numel(names),'single');
for k=1:numel(names)
    row(k)=m.(names{k});
end
fdata=nan(nRec,1,'single');
if isfield(m,'fData') && numel(m.fData)==nRec
    fdata=m.fData;
end
end

% =====================================================================
function lv=withMarkers(lv)
%withMarkers  The level set the core is CALLED with: the user's plus 'markers'.
%   The model-free markers cost no fit and are always duplicated into the metrics
%   tables, so they are always computed; the user's own selection still decides what
%   reaches the saved tree (see the header).
if ischar(lv)||isstring(lv), lv=cellstr(lv); end
if isempty(lv), lv={}; end
lv=reshape(lv,1,[]);
if ~ismember('markers',lv), lv=[lv {'markers'}]; end
end

% =====================================================================
function requireInjection(s,nFiles)
%requireInjection  The two settings with no upstream owner, checked before any file is
%   opened.  NOTHING IN THE LIBRARY WRITES AN INJECTION TIME - there is no producing
%   step to inherit it from the way runFitNVC inherits its stimulus geometry from the
%   external cycle - so it lives here, and a missing or mis-sized list has to be a
%   named message rather than a subscript error inside the loop.
if ~isfield(s,'injectionSec') || isempty(s.injectionSec)
    error('runFitVasoreactivity:noInjection', ...
        ['s.injectionSec is required - when the drug was given, in seconds from the ' ...
         'start of each recording.  Give one number, or one per file in the order ' ...
         'fNames lists them.']);
end
v=double(s.injectionSec(:));
if ~all(isfinite(v))
    error('runFitVasoreactivity:badInjection', ...
        's.injectionSec must be finite - it is a time in seconds, one per recording.');
end
if ~(isscalar(v) || numel(v)==nFiles)
    error('runFitVasoreactivity:injectionCount', ...
        ['s.injectionSec has %d entries but there are %d files.  Give one number for ' ...
         'all of them, or exactly one per file in the order fNames lists them.'], ...
        numel(v),nFiles);
end
if ~isfield(s,'baselineSec') || numel(s.baselineSec)~=2
    error('runFitVasoreactivity:noBaseline', ...
        ['s.baselineSec is required - the [start end] of the quiet stretch before ' ...
         'the injection, in seconds.  It must end before the earliest injection.']);
end
end

% =====================================================================
function v=perFile(list,fidx)
%perFile  One entry of a per-file setting; a scalar covers every file.
list=double(list(:));
if isscalar(list), v=list; else, v=list(fidx); end
end

% =====================================================================
function m=resampleCovariate(s,fidx,L)
%resampleCovariate  The measured pressure trace on the fit's own time base, or [].
%   STORED, NEVER SUBTRACTED (see the header).  It arrives on its own clock - a
%   physiological recorder samples nothing like an LSCI camera - so it is interpolated
%   rather than decimated, and it is NOT extrapolated: outside its own span it is NaN,
%   which says "not measured here" instead of inventing a pressure.
m=[];
if ~isfield(s,'mapTrace') || isempty(s.mapTrace), return, end
tr=perFileVector(s.mapTrace,fidx);
if isempty(tr), return, end
if ~isfield(s,'mapTime') || isempty(s.mapTime)
    error('runFitVasoreactivity:noMapTime', ...
        ['s.mapTrace was given without s.mapTime, so there is no clock to place the ' ...
         'pressure on.  Give the seconds each pressure sample was taken at.']);
end
tt=perFileVector(s.mapTime,fidx);
if numel(tt)~=numel(tr)
    error('runFitVasoreactivity:mapLength', ...
        ['s.mapTime has %d samples and s.mapTrace has %d; they describe the same ' ...
         'measurement and must be the same length.'],numel(tt),numel(tr));
end
m=interp1(double(tt(:)),double(tr(:)),L.time,'linear',NaN);
end

% =====================================================================
function v=perFileVector(item,fidx)
%perFileVector  One file's vector from a cell of vectors, or the same vector for all.
if iscell(item)
    if isscalar(item), v=item{1}; else, v=item{fidx}; end
else
    v=item;
end
v=double(v(:));
end

% =====================================================================
function drawVrcPage(rep,results,s,mapTrace)
%drawVrcPage  The one report page: the recording's mean trace with its fit, the
%   injection and baseline geometry, the residuals, and the fitted parameters.
%
%   THE PAGE FITS THE MEAN TRACE, ONCE.  Averaging the per-segment fits would give a
%   curve that is not a fit of anything and whose parameters could not be tabulated
%   beside it; fitting the mean gives one trace, one fit, one residual and one
%   parameter set that all describe the same thing - the recording as a whole.
%
%   AND IT CANNOT KILL THE RUN.  reportSave already swallows a failed export; this
%   swallows a failed DRAWING, because the results are computed by the time the page is
%   attempted and losing a whole recording's analysis to a plotting fault would be
%   absurd.  The figure is deleted on every path.
fh=[];
try
    [y,what]=meanTrace(results);
    if isempty(y), return, end

    sPage=s; sPage.segVrcReturn={'markers','model','reconstruction'};
    LP=fitVasoreactivity(results.time,sPage);
    m =fitVasoreactivity(y,LP,sPage);
    t =double(LP.time);
    %the measured curve is decimated by the SAME function the core used, so the page
    %cannot draw the fit against a slightly different rule than it was fitted to
    yd=blockDecimate(double(y),LP.avgN);
    fd=double(m.fData);

    fh=reportFigure(rep,'vasoreactivityfit');
    tl=tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');

    % ---- the trace and the fit ----
    ax=nexttile(tl,1,[1 2]);
    hold(ax,'on')
    xregion(ax,LP.t1+LP.bl(1),LP.t1+LP.bl(2), ...
        'FaceColor',[0.45 0.45 0.45],'FaceAlpha',0.10);
    plot(ax,t,yd,'-','Color',[0.45 0.45 0.45],'LineWidth',1)
    plot(ax,t,fd,'-','Color',[0.85 0.20 0.15],'LineWidth',1.5)
    xline(ax,LP.tInj,'-','Color',[0.30 0.55 0.85],'LineWidth',1.5);
    hold(ax,'off')
    xlabel(ax,'Time, s'); ylabel(ax,what)
    legend(ax,{'baseline','measured','fitted','injection'}, ...
        'Location','best','Box','off')
    axis(ax,'tight'); grid(ax,'on')
    title(ax,'Vasoreactivity response')

    % ---- the residuals, with the pressure covariate when there is one ----
    %   THE COVARIATE GOES ON THE RESIDUAL PANEL, not beside the trace, because that is
    %   where it is diagnostic: a residual that tracks the pressure is the systemic
    %   confound the header warns about, visible as a shape rather than as a caveat.
    ax=nexttile(tl,4,[1 2]);
    hasMap=~isempty(mapTrace) && any(isfinite(mapTrace));
    if hasMap
        yyaxis(ax,'right')
        plot(ax,t,double(mapTrace),'-','Color',[0.55 0.30 0.65])
        ylabel(ax,'Arterial pressure')
        yyaxis(ax,'left')
    end
    plot(ax,t,yd-fd,'-','Color',[0.20 0.35 0.60])
    yline(ax,0,'-','Color',[0.6 0.6 0.6]);
    xlabel(ax,'Time, s'); ylabel(ax,'Measured - fitted')
    axis(ax,'tight'); grid(ax,'on')
    if hasMap
        title(ax,'Residuals, and the pressure they were not corrected for')
    else
        title(ax,'Residuals')
    end

    % ---- the parameters ----
    ax=nexttile(tl,3,[2 1]);
    axis(ax,'off')
    text(ax,0.02,0.98,pageLines(m,LP),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName','Courier New','FontSize',10,'Interpreter','none')
    title(ax,'Fitted response')

    reportSave(rep,fh,'vasoreactivityfit');   % deletes the figure on every path of its own
    fh=[];
catch
    % a page is a by-product and must never kill a run - reportSave says the same
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function [y,what]=meanTrace(results)
%meanTrace  The recording's average trace, and what it is a trace of.
y=[]; what='';
for c={'sData','Flow, BFI'; 'dvsData','Flow, BFI'; 'dvsDiameter','Diameter, px'}'
    nm=c{1}; lbl=c{2};
    if isfield(results,nm) && ~isempty(results.(nm))
        v=mean(double(results.(nm)),2,'omitnan');
        if all(isfinite(v)), y=v; what=lbl; return, end
    end
end
end

% =====================================================================
function txt=pageLines(m,L)
%pageLines  The parameter block, written for a reader rather than for a parser.
%   Minutes, not seconds, for every timing: a t_max of 780 s is a number to divide, and
%   13 min is a number to read.
lines={ ...
    sprintf('peak change  %.4g',m.Gain), ...
    sprintf('relative     %+.1f %%',100*m.CVR), ...
    sprintf('  +/- 95%%    %.2f %%',100*1.96*m.CVRSE), ...
    sprintf('onset        %.1f min',m.Onset/60), ...
    sprintf('time to peak %.1f min',m.TMaxModel/60), ...
    sprintf('washout      %.1f min',m.TWash/60), ...
    sprintf('rise rate    %.3g /s',m.K), ...
    '', ...
    sprintf('rate, rising %.4g 1/s',m.RateFast), ...
    sprintf('rate, decay  %.4g 1/s',m.RateSlow), ...
    sprintf('drift        %.3g per s',m.Drift), ...
    '', ...
    sprintf('R2           %.4f',m.R2), ...
    sprintf('RMSE         %.4g',m.RMSE), ...
    sprintf('correlated   rho %.2f, %d independent samples',m.Rho,round(m.NEff)), ...
    sprintf('starts agree %d of %d',m.StartsAgree,L.nStarts), ...
    sprintf('samples      %d at %.3g Hz',L.nT,L.fs)};
if m.Steal>0
    % THE ONE THING ON THIS PAGE THAT CHANGES WHAT THE NUMBERS MEAN, so it is not a
    % line in the block - it is called out.
    lines=[lines,{'', ...
        '*** STEAL DETECTED ***', ...
        sprintf('%.1f min below baseline',m.StealSec/60), ...
        'this model cannot describe a', ...
        'sub-baseline dip, so the peak', ...
        'and washout times above are', ...
        'not to be trusted'}];
end
txt=strjoin(lines,newline);
end
