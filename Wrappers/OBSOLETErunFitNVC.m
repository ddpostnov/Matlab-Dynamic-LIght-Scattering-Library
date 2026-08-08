%runFitNVC  Stimulus-locked NVC response fit into the results.nvc tree
%
%   runFitNVC(s,fNames) loads every *_e_BFI_r.mat epoch product in fNames (produced
%   by runExternalCycle -> runSegmentation -> runBFI: results.sData/dvsData/
%   dvsDiameter [nT x nSeg] on results.time, plus the averaged epoch cube
%   source.data [Y x X x nT]) and reduces every stimulus-locked trace to model-free
%   response markers and a fitted neurovascular-coupling response.
%
%   The per-trace science lives in the shared core fitNVC; this producer is the
%   multi-file GLUE: it reads the stimulus geometry off the recording, builds the
%   fit layout once per time base, runs the core over every segment (parfor) and
%   every masked pixel, and assembles the results into the RESULTS.nvc tree
%   (mirroring RESULTS.pulsatility and RESULTS.vasomotion).
%
%   THE STIMULUS GEOMETRY IS READ OFF THE RECORDING, NOT RETYPED.  When the stimulus
%   starts, how long it lasts and where the baseline window sits were all decided at
%   the external-cycle step and are recorded in that recording's own
%   settings.externalCycle.  This step reads them from there, so a protocol is
%   described once.  An explicit s.epochStimStartSec / s.stimDurationSec /
%   s.epochBaselineSec overrides, for a re-analysis that deliberately differs.
%
%   OUTPUT TREE  (RESULTS.nvc)
%     Shared root axes:
%       .time      [nT x 1]  epoch time base (= results.time), single
%       .stimulus  [nT x 1]  the stimulus boxcar, single - storing it makes every
%                            downstream plot self-describing without re-deriving it
%     One sub-tree per analysed signal, all floats SINGLE, with a field-name PREFIX
%     encoding the physical quantity (ns = NVC flow, nd = NVC diameter):
%       .sData        { scalars(ns*), timeVectors.fData }   (-> results.sMetrics)
%       .dvsData      { scalars(ns*), timeVectors.fData }   (-> results.dvsMetrics)
%       .dvsDiameter  { scalars(nd*), timeVectors.fData }   (-> results.dvsMetrics)
%       .ppx          { scalars(ns*) [Y x X] maps, timeVectors.fData [Y x X x nT] }
%     scalars.* (per segment [nSeg x 1], per pixel [Y x X]) carry the core's scalar
%     names with the prefix attached - nsPeak, nsTTP, nsZeta, nsR2, ndPeak, ... - and
%     the core owns which names exist at which level (see fitNVC).  An INVALID trace
%     (non-finite) is NaN in EVERY field.  timeVectors.fData is the fitted curve.
%
%   LEVELS  (selector cell arrays; a subset of {'markers','model','reconstruction'})
%     s.segNvcReturn  per-segment levels; default {'markers','model','reconstruction'}.
%     s.ppxNvcReturn  per-pixel levels + GATE: NON-EMPTY runs the per-pixel path,
%                     [] skips it; default (absent) {'markers'} (maps ON, fit OFF).
%     THE LEVELS GATE THE TREE, NOT THE METRICS TABLE.  The model-free markers cost
%     no fit, so they are always computed and always duplicated into the metrics
%     tables; dropping 'markers' drops them from the tree only.  This is the same
%     rule runPulsatility follows and it keeps a metrics column from appearing and
%     disappearing with a display setting.
%
%   METRICS-TABLE DUPLICATION (a small key set, same names as the tree, single)
%     results.sMetrics    (from sData)      : nsPeak nsTTP nsFWHM nsUndershoot
%                                             nsZeta nsR2
%     results.dvsMetrics  (from dvsData)    : the same six ns*
%                         (from dvsDiameter): the six nd*
%     Zeta and R2 come from the fit, so they appear only when 'model' was computed;
%     under s.nvcModel='doubleGamma' there is no Zeta at all and the other five are
%     written.  BFI / std(BFI) are runBFI's columns and are LEFT UNTOUCHED.
%
%   ONE PATH FOR THE MARKERS, SEGMENTS AND PIXELS ALIKE.  Unlike runPulsatility -
%   which keeps a second, vectorised per-pixel marker path to stay bit-identical to
%   its pre-refactor output - every marker here is computed by the core, so a
%   per-pixel map and a per-segment scalar cannot drift apart.  There is no legacy
%   golden to reproduce, so there is no reason for a second definition.
%
%   INPUTS
%     s        parameter struct with fields
%                • nvcModel        'secondOrder' (default) | 'doubleGamma'
%                • nvcDip          false (default) | true - also fit the optional
%                                  pre-dip model and store both AICc values
%                • stimDurationSec how long the stimulus lasts, seconds.  REQUIRED,
%                                  here or in the recording's settings.externalCycle
%                • segNvcReturn    per-segment level cell (default all three)
%                • ppxNvcReturn    per-pixel level cell / gate (default {'markers'};
%                                  [] = per-pixel analysis off)
%                • nvcStarts       multi-start start points (default 16)
%                • nvcWeights      (optional) [nT x 1] per-timepoint weights
%                • parforNvcSegments logical, default true
%                • parforNvcPixels   logical, default true - both are WORKER BOUNDS,
%                                  not branches: false runs the identical loop body
%                                  serially in the client and starts no pool.
%     fNames   cell array of *_e_BFI_r.mat paths.
%                • Optional workbench hooks in s (no-op when absent):
%                  s.stageFcn(stage,detail), s.cancelFcn()->tf.  Cancel is checked
%                  between files (never inside a parfor).
%
%   OUTPUT FILES (side-effects) - NON-DESTRUCTIVE: the averaged epoch is preserved
%       *_e_BFI_d.mat   SOURCE   - NEVER re-saved (source.data read, never modified)
%       *_e_BFI_r.mat   RESULTS  - RESULTS.nvc.* + the ns*/nd* metrics columns
%       *_e_BFI_s.mat   SETTINGS - field settings.runFitNVC added
%       *_rep_nvcfit.jpg          - the mean trace with the fit, the residuals and
%                                   the fitted parameters
%
%   EXAMPLE
%     s.stimDurationSec=5;
%     s.nvcModel='secondOrder';
%     s.segNvcReturn={'markers','model','reconstruction'};
%     s.ppxNvcReturn={'markers'};
%     fNames = getFileNamesList(resultsFolder,'*_e_BFI_r.mat');
%     runFitNVC(s, fNames(:));
%
%   DEPENDS ON
%     Core/NVC/fitNVC (shared NVC response core), MATLAB Optimization Toolbox
%     (lsqcurvefit); core LSCI library utilities.
%
% See also: fitNVC, runExternalCycle, runPulsatility, runVasomotion, getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 04-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Response model
% s.nvcModel='secondOrder'; % 'secondOrder' fits the flow response of a damped
%                         % autoregulatory system - its undershoot is mechanistic and
%                         % its damping ratio is comparable between protocols.
%                         % 'doubleGamma' is a descriptive alternative.
% s.stimDurationSec=5;    % how long the stimulus lasts, seconds.  Read from the
%                         % external-cycle step when it was set there.
% s.nvcDip=false;         % true also fits a fast pre-constriction and reports both
%                         % models side by side - nothing is selected for you.
% %ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
% s.segNvcReturn={'markers','model','reconstruction'};  % markers (model-free
%                         % scalars), model (fitted parameters - runs the fit),
%                         % reconstruction (the fitted curve). Default = all three.
% %ADJUSTED IF NECESSARY - Per-pixel maps/cubes (GATE + selector)
% s.ppxNvcReturn={'markers'};  % NON-EMPTY = per-pixel maps ON (RESULTS.nvc.ppx),
%                         % [] = off.  Adding 'model' fits every masked pixel, which
%                         % is many thousands of fits - expect minutes, not seconds.
% %ADJUSTED IF NECESSARY - Fitting
% s.nvcStarts=16;         % how many start points the fit is tried from.  Fewer is
%                         % faster and more likely to settle on a wrong answer.
% %ADJUSTED IF NECESSARY - Parallel execution
% s.parforNvcSegments=true;
% s.parforNvcPixels=true;  % false fits one at a time in this MATLAB and starts no
%                         % parallel pool - slower, but it costs no worker processes.

%------------- BEGIN CODE --------------
function runFitNVC(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_e_BFI_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_e_BFI_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(resultsFolder,''*_e_BFI_r.mat'').']);
end

%resolve defaults so they are recorded in the saved settings (the core defaults the
%same values internally).  segNvcReturn defaults when absent OR empty; ppxNvcReturn
%defaults ONLY when the field is ABSENT - an explicit [] must stay empty (per-pixel
%analysis off), exactly like runPulsatility's s.ppxPulsReturn.
if ~isfield(s,'nvcModel') || isempty(s.nvcModel)
    s.nvcModel='secondOrder';
end
if ~isfield(s,'nvcDip') || isempty(s.nvcDip)
    s.nvcDip=false;
end
if ~isfield(s,'nvcStarts') || isempty(s.nvcStarts)
    s.nvcStarts=16;
end
if ~isfield(s,'segNvcReturn') || isempty(s.segNvcReturn)
    s.segNvcReturn={'markers','model','reconstruction'};
end
if ~isfield(s,'ppxNvcReturn')
    s.ppxNvcReturn={'markers'};
end
%PARALLELISM IS OPTIONAL.  Each switch is a BOUND on its parfor (Inf workers or 0),
%never a branch: parfor(...,0) runs the identical loop body serially IN THE CLIENT
%and starts no pool at all.
if ~isfield(s,'parforNvcSegments') || isempty(s.parforNvcSegments)
    s.parforNvcSegments=true;
end
if ~isfield(s,'parforNvcPixels') || isempty(s.parforNvcPixels)
    s.parforNvcPixels=true;
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  Cancel is only
% checked between files (a hook inside a parfor would broadcast oddly).
rep=reportOpen(s,'NVC fit',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        %SOURCE is read only for the per-pixel path (and never written back).
        if ~isempty(s.ppxNvcReturn)
            load(getProductPath(s.fName,'d'),'source')
        end
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        if ~isfield(results,'time') || isempty(results.time)
            error('runFitNVC:noTime', ...
                '%s carries no results.time, so there is no epoch to fit.',s.fName);
        end

        %the stimulus geometry of THIS recording, s overriding what the external
        %cycle recorded (see the header)
        sf=applyGeometry(s,settings);

        %runFitNVC owns RESULTS.nvc entirely; rebuild it from scratch so no stale
        %sub-branch survives a re-run with different levels.
        if isfield(results,'nvc')
            results=rmfield(results,'nvc');
        end

        %TWO LAYOUTS, ONE TIME BASE.  layoutCalc is what the core is actually called
        %with - the requested levels PLUS 'markers', because the markers cost no fit
        %and the metrics tables always carry them.  layoutTree is the user's own
        %selection, and its scalarNames (a subsequence of layoutCalc's) is exactly
        %what enters the saved tree.  Neither the loop nor the assembly below has to
        %know a single field name.
        sCalc=sf; sCalc.segNvcReturn=withMarkers(sf.segNvcReturn);
        layoutCalc=fitNVC(results.time,sCalc);
        layoutTree=fitNVC(results.time,sf);
        names=layoutCalc.scalarNames;
        nT=layoutCalc.nT;

        %shared root axes (single; the tree's OWN copy - results.time stays double)
        results.nvc.time    =single(results.time(:));
        results.nvc.stimulus=single(layoutCalc.u);

        % =============================================================
        % Per-segment analysis: sData/dvsData (ns prefix), dvsDiameter (nd prefix).
        % =============================================================
        nwSeg=0; if s.parforNvcSegments, nwSeg=Inf; end   %worker bound, not a branch
        nRecSeg=1; if layoutTree.want.reconstruction, nRecSeg=nT; end

        sigNames={'sData','dvsData','dvsDiameter'};
        for kSig=1:numel(sigNames)
            sigName=sigNames{kSig};
            if ~isfield(results,sigName), continue; end
            if strcmp(sigName,'dvsDiameter'), prefix='nd'; else, prefix='ns'; end

            sigMat=results.(sigName);
            if isempty(sigMat), continue; end
            nSeg=size(sigMat,2);

            %ONE ROW PER SEGMENT, ONE COLUMN PER SCALAR.  The core emits some twenty
            %to thirty scalars, so the twenty named accumulators runPulsatility uses
            %would be unreadable here; the row/column form also lets the tree and the
            %metrics duplication below be driven by the name list rather than typed
            %out twice.  BOTH slices are assigned on EVERY iteration, so nothing
            %about the parfor depends on which branch the core took.
            M=nan(nSeg,numel(names),'single');
            F=nan(nRecSeg,nSeg,'single');
            %the workbench hooks are function handles bound to a GUI; they are
            %transport rather than parameters, the core never reads them, and
            %broadcasting them to workers is the one thing in s that can fail to
            %serialise.  reportSettings strips exactly those two, here as at the save.
            sPar=reportSettings(sCalc);
            parfor (i=1:nSeg, nwSeg)
                [M(i,:),F(:,i)]=nvcTraceRow(sigMat(:,i),layoutCalc,sPar,names,nRecSeg);
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
            results.nvc.(sigName)=T;

            %duplicate the key set into the metrics tables (DATA-MODEL section 11).
            %Row order matches the 1:nSeg loop (segment invariant).  BFI/std(BFI) are
            %runBFI's columns and are intentionally NOT written here.
            switch sigName
                case 'sData',                    tbl='sMetrics';
                case {'dvsData','dvsDiameter'},  tbl='dvsMetrics';
            end
            dup=intersect({'Peak','TTP','FWHM','Undershoot','Zeta','R2'},names,'stable');
            for k=1:numel(dup)
                results.(tbl).([prefix dup{k}])=M(:,strcmp(names,dup{k}));
            end
        end

        % =============================================================
        % Per-pixel analysis (gated by non-empty s.ppxNvcReturn; ns prefix, since
        % source.data is a flow cube).  Every masked pixel goes through the same core
        % the segments did - see the header on why there is no second marker path.
        % =============================================================
        if ~isempty(s.ppxNvcReturn)
            if ~isfield(results,'sMap')
                error('runFitNVC:noSMap', ...
                    'results.sMap not found; runSegmentation must be run before per-pixel NVC.');
            end
            sPix=sf; sPix.segNvcReturn=s.ppxNvcReturn;
            sPixCalc=sPix; sPixCalc.segNvcReturn=withMarkers(s.ppxNvcReturn);
            layoutPixCalc=fitNVC(results.time,sPixCalc);
            layoutPixTree=fitNVC(results.time,sPix);
            namesPix=layoutPixCalc.scalarNames;

            sz=size(source.data); Y=sz(1); X=sz(2); nTp=sz(3); npx=Y*X;
            Dpix=reshape(source.data,npx,nTp);       %[npx x nT] pixel timecourses
            sMapLin=results.sMap(:);

            %NaN-preallocated, so a background pixel - which no iteration writes -
            %stays NaN in every map.  An unfitted pixel has no response, and a zero
            %there would be a measurement that was never made.
            Mp=nan(npx,numel(namesPix),'single');
            nRecPix=1; if layoutPixTree.want.reconstruction, nRecPix=nTp; end
            Fp=nan(npx,nRecPix,'single');

            nwPix=0; if s.parforNvcPixels, nwPix=Inf; end   %worker bound, not a branch
            sParPix=reportSettings(sPixCalc);               %see the note above
            parfor (p=1:npx, nwPix)
                if sMapLin(p)==0, continue; end              %background stays NaN
                [Mp(p,:),Fp(p,:)]=nvcTraceRow(Dpix(p,:).',layoutPixCalc,sParPix, ...
                    namesPix,nRecPix);
            end

            ppx=struct();
            for k=1:numel(layoutPixTree.scalarNames)
                nm=layoutPixTree.scalarNames{k};
                ppx.scalars.(['ns' nm])=reshape(Mp(:,strcmp(namesPix,nm)),Y,X);
            end
            if layoutPixTree.want.reconstruction
                ppx.timeVectors.fData=reshape(Fp,Y,X,nTp);
            end
            results.nvc.ppx=ppx;
        end

        % =============================================================
        % Report page: the recording's mean trace, its fit, and the parameters.
        % =============================================================
        drawNvcPage(rep,results,sf,nT);

        settings.runFitNVC=reportSettings(sf);
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
function [row,fdata]=nvcTraceRow(y,layout,s,names,nRec)
%nvcTraceRow  One trace through the core, flattened to a fixed-width row.
%   The row carries the layout's scalars IN ORDER, so the caller indexes by name
%   once (when it assembles the tree) instead of once per segment.  nRec is nT when
%   the reconstruction was asked for and 1 otherwise, which keeps the parfor's second
%   sliced output the same shape on every iteration without allocating a full cube
%   nobody wants.
m=fitNVC(y,layout,s);
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
function s=applyGeometry(s,settings)
%applyGeometry  The stimulus geometry of THIS recording (see the header).
%   The external-cycle step decided when the stimulus starts, how long it lasts and
%   where the baseline window is, and recorded all three in settings.externalCycle.
%   They are taken from there unless s names them, so a protocol is described once
%   and a re-analysis can still deviate on purpose.
ec=struct();
if isfield(settings,'externalCycle') && isstruct(settings.externalCycle)
    ec=settings.externalCycle;
end

s.epochStimStartSec=preferred(s,ec,'epochStimStartSec');
s.epochBaselineSec =preferred(s,ec,'epochBaselineSec');
s.stimDurationSec  =preferred(s,ec,'stimDurationSec');

if isempty(s.epochStimStartSec) || isempty(s.epochBaselineSec)
    error('runFitNVC:noStimGeometry', ...
        ['This product records no stimulus geometry, so the response cannot be ' ...
         'placed on the epoch.  Fit an epoch product written by the external-cycle ' ...
         'step, or set s.epochStimStartSec and s.epochBaselineSec explicitly.']);
end
if isempty(s.stimDurationSec)
    error('runFitNVC:noStimDuration', ...
        ['How long the stimulus lasts is not recorded for this recording and was ' ...
         'not given.  Set s.stimDurationSec here, or set it at the external-cycle ' ...
         'step so every later step inherits it.']);
end
end

% =====================================================================
function v=preferred(s,ec,name)
%preferred  s wins, the recording's own settings are the default, [] means neither.
v=[];
if isfield(ec,name) && ~isempty(ec.(name)), v=ec.(name); end
if isfield(s,name)  && ~isempty(s.(name)),  v=s.(name);  end
end

% =====================================================================
function drawNvcPage(rep,results,s,nT)
%drawNvcPage  The one report page: the recording's mean trace with its fit, the
%   residuals, and the fitted parameters.
%
%   THE PAGE FITS THE MEAN TRACE, ONCE.  Averaging the per-segment fits would give a
%   curve that is not a fit of anything and whose parameters could not be tabulated
%   beside it; fitting the mean gives one trace, one fit, one residual and one
%   parameter set that all describe the same thing - the recording as a whole.  It
%   is one extra call to the core and it runs whatever levels the user chose, so a
%   markers-only run still gets a page with a fit on it.
%
%   AND IT CANNOT KILL THE RUN.  reportSave already swallows a failed export; this
%   swallows a failed DRAWING, because the results are computed by the time the page
%   is attempted and losing a whole recording's analysis to a plotting fault would be
%   absurd.  The figure is deleted on every path - reportSave does it when it gets
%   that far, the tail below when it does not.
fh=[];
try
    [y,what]=meanTrace(results);
    if isempty(y), return, end

    sPage=s; sPage.segNvcReturn={'markers','model','reconstruction'};
    L=fitNVC(results.time,sPage);
    m=fitNVC(y,L,sPage);

    fh=reportFigure(rep,'nvcfit');
    tl=tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');

    % ---- the trace and the fit ----
    ax=nexttile(tl,1,[1 2]);
    hold(ax,'on')
    xregion(ax,L.tStim,L.tStim+L.D,'FaceColor',[0.30 0.55 0.85],'FaceAlpha',0.12);
    plot(ax,results.time,y,'-','Color',[0.45 0.45 0.45],'LineWidth',1)
    plot(ax,results.time,double(m.fData),'-','Color',[0.85 0.20 0.15],'LineWidth',1.5)
    hold(ax,'off')
    xlabel(ax,'Time, s'); ylabel(ax,what)
    legend(ax,{'stimulus','measured','fitted'},'Location','best','Box','off')
    axis(ax,'tight'); grid(ax,'on')
    title(ax,'Averaged response')

    % ---- the residuals ----
    ax=nexttile(tl,4,[1 2]);
    plot(ax,results.time,y-double(m.fData),'-','Color',[0.20 0.35 0.60])
    yline(ax,0,'-','Color',[0.6 0.6 0.6]);
    xlabel(ax,'Time, s'); ylabel(ax,'Measured - fitted')
    axis(ax,'tight'); grid(ax,'on')
    title(ax,'Residuals')

    % ---- the parameters ----
    ax=nexttile(tl,3,[2 1]);
    axis(ax,'off')
    text(ax,0.02,0.98,pageLines(m,L,nT),'Units','normalized', ...
        'VerticalAlignment','top','HorizontalAlignment','left', ...
        'FontName','Courier New','FontSize',10,'Interpreter','none')
    title(ax,'Fitted response')

    reportSave(rep,fh,'nvcfit');    % deletes the figure on every path of its own
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
function txt=pageLines(m,L,nT)
%pageLines  The parameter block, written for a reader rather than for a parser.
lines={sprintf('model        %s',L.model)};
if strcmp(L.model,'secondOrder')
    lines=[lines,{ ...
        sprintf('gain         %.4g',m.Gain), ...
        sprintf('onset        %.3g s',m.Onset), ...
        sprintf('tau signal   %.3g s',m.TauS), ...
        sprintf('tau feedback %.3g s',m.TauF), ...
        sprintf('damping      %.3f',m.Zeta), ...
        sprintf('ring period  %.3g s',m.RingPeriod)}];
else
    lines=[lines,{ ...
        sprintf('gain         %.4g',m.Gain), ...
        sprintf('onset        %.3g s',m.Onset), ...
        sprintf('peak shape   %.3g / %.3g',m.A1,m.Beta1), ...
        sprintf('dip shape    %.3g / %.3g',m.A2,m.Beta2), ...
        sprintf('dip fraction %.3f',m.CRatio)}];
end
lines=[lines,{'', ...
    sprintf('peak         %.4g',m.Peak), ...
    sprintf('time to peak %.3g s',m.TTP), ...
    sprintf('width        %.3g s',m.FWHM), ...
    sprintf('undershoot   %.3g',m.URatio), ...
    '', ...
    sprintf('R2           %.4f',m.R2), ...
    sprintf('RMSE         %.4g',m.RMSE), ...
    sprintf('starts agree %d of %d',m.StartsAgree,L.nStarts), ...
    sprintf('samples      %d',nT)}];
if L.dip
    lines=[lines,{'', ...
        sprintf('with pre-dip R2 %.4f',m.DipR2), ...
        sprintf('AICc  %.1f vs %.1f',m.AICc,m.AICcDip)}];
end
txt=strjoin(lines,newline);
end
