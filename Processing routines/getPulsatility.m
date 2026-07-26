%getPulsatility  Harmonic pulsatility analysis into the results.pulsatility tree
%
%   getPulsatility(s,fNames) loads every *_BFI_d.mat cardiac-cycle triplet in
%   fNames (produced by getInternalCycle/getExternalCycle -> getSegmentation ->
%   getBFI: results.sData/dvsData/dvsDiameter [nT x nSeg] on results.time, plus
%   the source cube source.data [Y x X x nT]) and reduces every averaged cycle to
%   classical pulsatility markers and, optionally, an nHarm-harmonic sine fit
%
%        y(x) = SUM_{n=1..nHarm} a_n*sin(2*pi*n*x + b_n) + c
%
%   The per-cycle "science" lives in the shared core getPulsatilityMetrics; this
%   producer is the multi-file GLUE: it builds the fit layout once, runs the core
%   over every segment (parfor) and every masked pixel, and assembles the results
%   into the RESULTS.pulsatility tree (mirroring RESULTS.vasomotion).
%
%   OUTPUT TREE  (RESULTS.pulsatility)
%     Shared root axes:
%       .time       [nT x 1]     one-cycle time base (= results.time), single
%       .harmonics  [nHarm x 1]  [1;2;...;nHarm] column index for scalars.hAmp/hPhase
%     One sub-tree per analysed signal, all floats SINGLE, with a field-name PREFIX
%     encoding the physical quantity (ps = pulsatile flow, pd = pulsatile diameter):
%       .sData        { scalars(ps*), timeVectors.fData }   (-> results.sMetrics)
%       .dvsData      { scalars(ps*), timeVectors.fData }   (-> results.dvsMetrics)
%       .dvsDiameter  { scalars(pd*), timeVectors.fData }   (-> results.dvsMetrics)
%       .ppx          { scalars(ps*) [Y x X] maps, timeVectors.fData [Y x X x nT] }
%     scalars.* (per segment [nSeg x 1], per pixel [Y x X]):
%       <p>Min <p>Max <p>TimeMin <p>TimeMax <p>Mean <p>Std <p>PI <p>RI <p>SymRatio
%       (markers, always computed) and hAmp/hPhase [.. x nHarm] + <p>R2 (the fit;
%       hAmp/hPhase are bare - the sub-tree key disambiguates flow vs diameter).
%       PI=(max-min)/mean (Gosling), RI=(max-min)/max (Pourcelot), SymRatio=
%       descend/ascend cycle asymmetry.  timeVectors.fData [nT x nSeg] is the model
%       reconstruction.  An INVALID cycle (non-finite / all-NaN) is NaN in EVERY field.
%
%   LEVELS  (selector cell arrays; a subset of {'markers','model','reconstruction'})
%     s.segPulsReturn  per-segment levels; default {'markers','model','reconstruction'}.
%     s.ppxPulsReturn  per-pixel levels + GATE: NON-EMPTY runs the per-pixel path,
%                      [] skips it; default (absent) {'markers'} (maps ON, fit OFF).
%     'markers' stores the model-free scalars (also duplicated into the metrics
%     tables regardless); 'model'/'reconstruction' run the fit (coefficients /
%     reconstruction cube).  s.nHarm (default 5) sets the number of harmonics.
%
%   METRICS-TABLE DUPLICATION (a small key set, same names as the tree, single)
%     results.sMetrics    (from sData)      : psPI psRI psMean psTimeMin psTimeMax psSymRatio
%     results.dvsMetrics  (from dvsData)    : the six ps*
%                         (from dvsDiameter): the six pd*
%     BFI / std(BFI) are getBFI's columns and are LEFT UNTOUCHED (not written here).
%
%   INPUTS
%     s        parameter struct with fields
%                • nHarm          number of harmonics in the sine model (default 5)
%                • segPulsReturn  per-segment level cell (default all three)
%                • ppxPulsReturn  per-pixel level cell / gate (default {'markers'};
%                                 [] = per-pixel analysis off)
%                • fitOptions     (optional) fitoptions override for the model
%     fNames   cell array of *_BFI_d.mat paths.
%
%   OUTPUT FILES (side-effects) - NON-DESTRUCTIVE: the raw cycle is preserved
%       *_BFI_d.mat   SOURCE   - NEVER re-saved (source.data read, never modified)
%       *_BFI_r.mat   RESULTS  - RESULTS.pulsatility.* + the ps*/pd* metrics columns
%       *_BFI_s.mat   SETTINGS - field settings.getPulsatility added
%
%   EXAMPLE
%     s.nHarm=5;
%     s.segPulsReturn={'markers','model','reconstruction'};
%     s.ppxPulsReturn={'markers'};
%     files = dir(fullfile(dataRoot,'*_BFI_d.mat'));
%     getPulsatility(s, fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     Basic functions/getPulsatilityMetrics (shared harmonic pulsatility core),
%     MATLAB Curve Fitting Toolbox (fit); core LSCI library utilities.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 25-July-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Harmonic model
% s.nHarm=5;              % number of harmonics in y=SUM a_n*sin(2*pi*n*x+b_n)+c
% %ADJUSTED IF NECESSARY - Which per-segment analysis levels to compute/store
% s.segPulsReturn={'markers','model','reconstruction'};  % markers (scalars),
%                         % model (hAmp/hPhase/R2 - runs the fit), reconstruction
%                         % (timeVectors.fData - runs the fit). Default = all three.
% %ADJUSTED IF NECESSARY - Per-pixel maps/cubes (GATE + selector)
% s.ppxPulsReturn={'markers'};  % NON-EMPTY = per-pixel maps ON (RESULTS.pulsatility.ppx),
%                         % [] = off. {'markers'} gives the light [Y x X] marker maps
%                         % (reproduces the old always-computed imgPI/imgRI); add
%                         % 'model'/'reconstruction' to fit every masked pixel (large
%                         % cubes at full field resolution).

function getPulsatility(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_BFI_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_BFI_d.mat".');
end

%resolve defaults so they are recorded in the saved settings (the core defaults
%the same values internally).  nHarm and segPulsReturn default when absent OR empty;
%ppxPulsReturn defaults ONLY when the field is ABSENT - an explicit [] must stay
%empty (per-pixel analysis off), exactly like getVasomotion's s.ppxVsmReturn.
if ~isfield(s,'nHarm') || isempty(s.nHarm)
    s.nHarm=5;
end
if ~isfield(s,'segPulsReturn') || isempty(s.segPulsReturn)
    s.segPulsReturn={'markers','model','reconstruction'};
end
if ~isfield(s,'ppxPulsReturn')
    s.ppxPulsReturn={'markers'};
end

for fidx=1:1:numel(fNames)
    if ~isempty(fNames{fidx})
        disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        s.fName=fNames{fidx};
        clearvars results source settings
        %SOURCE is read only for the per-pixel path (and never written back).
        if ~isempty(s.ppxPulsReturn)
            load(s.fName,'source')
        end
        load(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings');
        load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

        %getPulsatility owns RESULTS.pulsatility entirely; rebuild it from scratch so
        %no stale sub-branch survives a re-run with different levels.
        if isfield(results,'pulsatility')
            results=rmfield(results,'pulsatility');
        end

        %SETUP once for the results.time base (compiling the fittype is the one
        %non-trivial per-time-base cost).  All three per-segment signals share it.
        layout=getPulsatilityMetrics(results.time,s);
        nHarm=layout.nHarm; want=layout.want; nT=numel(results.time);

        %shared root axes (single; the tree's OWN copy - results.time stays double)
        results.pulsatility.time=single(results.time(:));
        results.pulsatility.harmonics=single((1:nHarm)');

        % =============================================================
        % Per-segment analysis: sData/dvsData (ps prefix), dvsDiameter (pd prefix).
        % =============================================================
        sigNames={'sData','dvsData','dvsDiameter'};
        for kSig=1:numel(sigNames)
            sigName=sigNames{kSig};
            if ~isfield(results,sigName), continue; end
            if strcmp(sigName,'dvsDiameter'), prefix='pd'; else, prefix='ps'; end

            %pull the signal matrix into a plain local (one dynamic-field access rather
            %than one per iteration; the marker source cube is chosen inside the core).
            sigMat=results.(sigName);
            if isempty(sigMat), continue; end
            nSeg=size(sigMat,2);

            %NaN-preallocate every accumulator; an invalid segment is written by no
            %branch and therefore stays NaN across every field (invalid -> NaN).
            [mMin,mMax,mTimeMin,mTimeMax,mMean,mStd,mPI,mRI,mSym,mR2]=deal(nan(nSeg,1,'single'));
            mHAmp=nan(nSeg,nHarm,'single'); mHPhase=nan(nSeg,nHarm,'single');
            mFData=nan(nT,nSeg,'single');

            %The per-segment harmonic fit runs SERIALLY (a plain for), matching the
            %pre-refactor code so the coefficients / reconstruction stay bit-identical
            %at single precision: the Trust-Region NLS is BLAS-threading sensitive, so a
            %parfor worker converges ~1 single ULP away from the serial client and would
            %fail the session-1 golden.  The per-segment fit is cheap (nSeg is small); the
            %expensive per-pixel fit below keeps its parfor (as it did pre-refactor).
            for i=1:nSeg
                m=getPulsatilityMetrics(sigMat(:,i),layout,s);
                if m.valid   %markers/fit written only for a valid cycle; invalid stays NaN
                    mMin(i)=m.min; mMax(i)=m.max; mTimeMin(i)=m.timeMin; mTimeMax(i)=m.timeMax;
                    mMean(i)=m.mean; mStd(i)=m.std; mPI(i)=m.PI; mRI(i)=m.RI; mSym(i)=m.symRatio;
                    if want.model
                        mHAmp(i,:)=m.hAmp; mHPhase(i,:)=m.hPhase; mR2(i)=m.R2;
                    end
                    if want.reconstruction
                        mFData(:,i)=m.fData;
                    end
                end
            end

            %assemble the sub-tree; each level is gated by `want` (markers scalars,
            %model scalars, reconstruction timeVectors).
            T=struct();
            if want.markers
                T.scalars.([prefix 'Min'])     =mMin;
                T.scalars.([prefix 'Max'])     =mMax;
                T.scalars.([prefix 'TimeMin']) =mTimeMin;
                T.scalars.([prefix 'TimeMax']) =mTimeMax;
                T.scalars.([prefix 'Mean'])    =mMean;
                T.scalars.([prefix 'Std'])     =mStd;
                T.scalars.([prefix 'PI'])      =mPI;
                T.scalars.([prefix 'RI'])      =mRI;
                T.scalars.([prefix 'SymRatio'])=mSym;
            end
            if want.model
                T.scalars.hAmp        =mHAmp;
                T.scalars.hPhase      =mHPhase;
                T.scalars.([prefix 'R2'])=mR2;
            end
            if want.reconstruction
                T.timeVectors.fData=mFData;
            end
            results.pulsatility.(sigName)=T;

            %duplicate the key marker set into the metrics tables (DATA-MODEL section 11).
            %Markers are always computed by the core, so this runs regardless of the
            %`markers` level.  Row order matches the 1:nSeg loop (segment invariant).
            %BFI/std(BFI) are getBFI's columns and are intentionally NOT written here.
            switch sigName
                case 'sData'
                    results.sMetrics.psPI=mPI;             results.sMetrics.psRI=mRI;
                    results.sMetrics.psMean=mMean;         results.sMetrics.psTimeMin=mTimeMin;
                    results.sMetrics.psTimeMax=mTimeMax;   results.sMetrics.psSymRatio=mSym;
                case 'dvsData'
                    results.dvsMetrics.psPI=mPI;           results.dvsMetrics.psRI=mRI;
                    results.dvsMetrics.psMean=mMean;       results.dvsMetrics.psTimeMin=mTimeMin;
                    results.dvsMetrics.psTimeMax=mTimeMax; results.dvsMetrics.psSymRatio=mSym;
                case 'dvsDiameter'
                    results.dvsMetrics.pdPI=mPI;           results.dvsMetrics.pdRI=mRI;
                    results.dvsMetrics.pdMean=mMean;       results.dvsMetrics.pdTimeMin=mTimeMin;
                    results.dvsMetrics.pdTimeMax=mTimeMax; results.dvsMetrics.pdSymRatio=mSym;
            end
        end

        % =============================================================
        % Per-pixel analysis (gated by non-empty s.ppxPulsReturn; ps prefix, since
        % source.data is a flow cube).  Markers stay VECTORISED over the whole cube
        % (today's fast path); the fit is a flat pixel parfor when model/reconstruction
        % is requested.  This replaces the old imgPI/imgRI/extendedMetrics.img*/fCoeffs.
        % =============================================================
        if ~isempty(s.ppxPulsReturn)
            if ~isfield(results,'sMap')
                error('results.sMap not found; getSegmentation must be run before per-pixel pulsatility.');
            end
            sPix=s; sPix.segPulsReturn=s.ppxPulsReturn;
            layoutPix=getPulsatilityMetrics(results.time,sPix);
            wantPix=layoutPix.want;
            fitRanPix=wantPix.model || wantPix.reconstruction;

            sz=size(source.data); Y=sz(1); X=sz(2); nTp=sz(3); npx=Y*X;
            T=results.time(end);                 %cycle duration for the wrap/symmetry
            invPix=all(isnan(source.data),3);    %[Y x X] all-NaN cycles (invalid -> NaN)

            if fitRanPix
                %fit every masked pixel (skip sMap==0 background).  The reconstruction
                %cube doubles as the marker source, so always evaluate it here even if
                %only 'model' was requested (a fit-then-eval; coefficients are unchanged).
                sFit=sPix; lv=sFit.segPulsReturn;
                if ischar(lv)||isstring(lv), lv=cellstr(lv); end
                if ~ismember('reconstruction',lv), lv=[lv(:)' {'reconstruction'}]; end
                sFit.segPulsReturn=lv;
                layoutFit=getPulsatilityMetrics(results.time,sFit);
                wFitModel=wantPix.model;

                Dpix=reshape(source.data,npx,nTp);   %[npx x nT] pixel timecourses (plain local)
                sMapLin=results.sMap(:);
                %model maps ZERO-init so background (unfitted) stays 0, reproducing the
                %old results.fCoeffs (zeros); an invalid masked pixel is explicitly NaN.
                hAmpAcc=zeros(npx,nHarm,'single'); hPhaseAcc=zeros(npx,nHarm,'single');
                r2Acc=zeros(npx,1,'single');
                %reconstruction cube seeded from the RAW cube so background stays raw
                %(the old fit overwrote only masked pixels); masked pixels overwritten below.
                fDataAcc=single(Dpix);

                tic
                parfor p=1:npx
                    if sMapLin(p)==0, continue; end          %background: model 0, fData raw
                    mp=getPulsatilityMetrics(Dpix(p,:).',layoutFit,sFit);
                    if mp.valid
                        fDataAcc(p,:)=mp.fData.';
                        if wFitModel
                            hAmpAcc(p,:)=mp.hAmp; hPhaseAcc(p,:)=mp.hPhase; r2Acc(p)=mp.R2;
                        end
                    else                                     %invalid masked pixel -> NaN
                        fDataAcc(p,:)=NaN;
                        hAmpAcc(p,:)=NaN; hPhaseAcc(p,:)=NaN; r2Acc(p)=NaN;
                    end
                end
                fprintf('Per-pixel pulsatility fit: %d masked pixels in %.2fs\n',nnz(sMapLin>0),toc);

                fDataCube=reshape(fDataAcc,Y,X,nTp);
                W=fDataCube;                                 %markers on the fitted cube
            else
                W=source.data;                               %markers on the raw cube
            end

            %vectorised model-free markers over the whole cube (identical expressions to
            %the pre-refactor imgPI/imgRI/extendedMetrics; times index source.time).
            [maxMap,maxIdx]=max(W,[],3);
            [minMap,minIdx]=min(W,[],3);
            tMaxMap=source.time(maxIdx);
            tMinMap=source.time(minIdx);
            wrap=tMinMap>tMaxMap;
            tMinMap(wrap)=tMinMap(wrap)-T;                   %wrap-shift so tMin<=tMax
            meanMap=mean(W,3,'omitnan');
            stdMap =std(W,0,3,'omitnan');
            piMap  =(maxMap-minMap)./meanMap;
            riMap  =(maxMap-minMap)./maxMap;
            ascMap =tMaxMap-tMinMap;
            symMap =(T-ascMap)./ascMap;

            %invalid pixels get time(1)/Inf from the vectorised max/min - override every
            %marker to NaN (invalid -> NaN, matching the per-segment rule).
            minMap(invPix)=NaN; maxMap(invPix)=NaN;
            tMinMap(invPix)=NaN; tMaxMap(invPix)=NaN;
            meanMap(invPix)=NaN; stdMap(invPix)=NaN;
            piMap(invPix)=NaN; riMap(invPix)=NaN; symMap(invPix)=NaN;

            %assemble the ppx sub-tree (all single; marker times cast from double).
            ppx=struct();
            ppx.scalars.psMin     =single(minMap);
            ppx.scalars.psMax     =single(maxMap);
            ppx.scalars.psTimeMin =single(tMinMap);
            ppx.scalars.psTimeMax =single(tMaxMap);
            ppx.scalars.psMean    =single(meanMap);
            ppx.scalars.psStd     =single(stdMap);
            ppx.scalars.psPI      =single(piMap);
            ppx.scalars.psRI      =single(riMap);
            ppx.scalars.psSymRatio=single(symMap);
            if fitRanPix && wantPix.model
                ppx.scalars.hAmp  =reshape(hAmpAcc,Y,X,nHarm);
                ppx.scalars.hPhase=reshape(hPhaseAcc,Y,X,nHarm);
                ppx.scalars.psR2  =reshape(r2Acc,Y,X);
            end
            if fitRanPix && wantPix.reconstruction
                ppx.timeVectors.fData=fDataCube;
            end
            results.pulsatility.ppx=ppx;
        end

        settings.getPulsatility=s;
        disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        %NON-DESTRUCTIVE: SOURCE (_d) is never re-saved - only RESULTS and SETTINGS.
        save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
        save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');
    end
end

end
