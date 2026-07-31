%runRegistration  Rigid/affine registration of several *_K_d.mat datasets
%
%   runRegistration(s,fNames) aligns every contrast cube in the cell
%   array fNames to the first file in the list.  For each follower it
%   offers three candidate transforms—identity, intensity-based
%   imregtform (tform1), and correlation-based imregcorr (tform2)—plus an
%   optional manual landmark stage.  The chosen affine2d object is applied
%   to every X × Y × T image in SOURCE and RESULTS.  Masks are combined
%   across files to produce hard/soft consensus maps, and (if
%   s.matchSegmentation==true) vessel/parenchyma segment IDs are reconciled
%   with the reference segmentation.
%
%   INPUTS
%     s        parameter structure
%                • optimizer, metric   objects from imregconfig("monomodal")
%                • tFormType           {'affine','similarity',…}
%                • matchSegmentation   logical, true = reconcile segment IDs
%                • prchNSize           parenchyma seed-grid spacing (pixels)
%                • silent              logical, true = auto-select the best
%                                      transform with no user interaction
%                • forceMethod         (silent) 'intensity' | 'correlation' to
%                                      override the smallest-Delta choice
%                • rotationLimit       degrees; reject an intensity/correlation
%                                      transform rotating > this ([]=no limit)
%     fNames   cell array of *_K_d.mat paths.  First file is the template.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf (before work + between files).
%
%   OUTPUT (side-effects)
%       For every file k
%         *_K_d.mat   SOURCE   – warped
%         *_K_r.mat   RESULTS  – warped + consensus masks
%         *_K_s.mat   SETTINGS – field settings.registration added
%       plus the segments report pages already present in the workflow
%       In silent mode also *_rep_registration.jpg - overlay of the warped image on
%       the reference, annotated with the chosen method and the Delta values
%
%   EXAMPLE
%     [opt,met]     = imregconfig("monomodal");
%     opt.MaximumIterations = 500;
%     p.optimizer   = opt;
%     p.metric      = met;
%     p.tFormType   = 'affine';
%     p.matchSegmentation = true;
%     D = dir(fullfile(dataRoot,'*_K_d.mat'));
%     runRegistration(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     registerToReference (candidate generation + transform selection),
%     enhanceForRegistration (the vessel-emphasis registration proxy), LSCI
%     segmentation utilities, and MATLAB's Image Processing Toolbox (imwarp,
%     imref2d, imshowpair, bwdist, etc.).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 23-July-2026


% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
% [s.optimizer,s.metric] = imregconfig("monomodal");
% s.optimizer.MaximumIterations=500;
% s.tFormType='affine';
% s.matchSegmentation=true;
% s.prchNSize=30; % Parenchymal pixels neighbourhoud.
% s.silent=true;
% s.forceMethod='correlation'; % 'intensity' or 'correlation' to force in silent mode
% s.rotationLimit=45; % degrees; reject registrations rotating > 45 ([] = none)

function runRegistration(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  Registration is
% caller-looped per group and can block on manual landmarks, so cancel is checked before
% any work and between files in the final save loop.
rep=reportOpen(s,'Registration',fNames);
if reportCancelled(rep), return; end

% ---- silent (non-interactive) mode ----------------------------------------
% When s.silent==true the best registration is chosen automatically instead of
% asking the user.  Files that are the same acquisition but a different LSCI
% product (differ only in the _c_K_ / _t_K_ / _s_K_ token) reuse the same
% transform, and a *_rep_registration overlay report is written for each file.
silentMode = isfield(s,'silent') && isscalar(s.silent) && s.silent==true;
normNames  = cell(size(fNames));
for i=1:numel(fNames)
    if ~isempty(fNames{i})
        [~,nm,ex]    = fileparts(fNames{i});
        normNames{i} = regexprep([nm ex],'_[cts]_K_','_X_K_');
    end
end

% s.forceMethod (silent mode, Delta-decided case only) may be 'intensity' or
% 'correlation' to always use that method instead of the smallest-Delta choice.
% Absent, empty or unrecognised -> Delta decides as usual.
forceMethod = '';
if isfield(s,'forceMethod') && (ischar(s.forceMethod) || isstring(s.forceMethod))
    forceMethod = lower(char(s.forceMethod));
end

% s.rotationLimit (degrees): if the intensity- or correlation-based transform
% rotates by more than this it is treated as an obvious failure (e.g. a ~180
% deg flip) - dropped from the choice in silent mode, flagged in interactive
% mode.  Empty [] or absent -> no limit is applied.
rotationLimit = [];
if isfield(s,'rotationLimit') && isnumeric(s.rotationLimit) && isscalar(s.rotationLimit)
    rotationLimit = abs(s.rotationLimit);
end

tforms=cell(size(fNames));
views=cell(size(fNames));      % dead: never assigned; retained so saved s.view stays [] (unchanged)
commonMask=cell(size(fNames));

% ---- processing order -----------------------------------------------------
% Files are normally visited in input order.  In silent mode the _t_K_d files
% are registered first (right after the reference) so that, for each
% acquisition, the _t product drives the registration and its _c_K_d / _s_K_d
% siblings reuse that transform.  The first (reference) file always stays first.
procOrder = (1:size(fNames,1))';
if silentMode && numel(procOrder)>1
    rest = procOrder(2:end);
    isT  = cellfun(@(f) ~isempty(f) && contains(f,'_t_K_d'), fNames(rest));
    procOrder = [procOrder(1); rest(isT); rest(~isT)];
end

% ---- registration proxies + masks -----------------------------------------
% Each file's vessel-emphasised registration proxy (dark vessels -> bright,
% background removed, within-mask histogram-equalised) is built by the shared
% primitive enhanceForRegistration - the recipe used to be inline here and is
% provably identical to it (Tests/Registration/testEnhanceForRegistration).  The
% transforms estimated from these proxies are applied to the pristine SOURCE/
% RESULTS data, which is re-loaded further down.  commonMask keeps the category
% values (not the binary support) for the consensus stage that follows.
proxy=cell(size(fNames));
for oidx=1:numel(procOrder)
    fidx=procOrder(oidx);
    if ~isempty(fNames{fidx})
        clearvars results source
        load(fNames{fidx},'source');
        load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

        img=mean(source.data,3);
        mask=ones(size(img));
        commonMask{fidx}=ones(size(img));
        if isfield(results,'cMask')
            mask=results.cMask>0;
            commonMask{fidx}=results.cMask;
        elseif isfield(results,'mask')
            mask=results.mask>0;
            commonMask{fidx}=results.mask;
        end
        proxy{fidx}=enhanceForRegistration(img,'mask',mask,'equalizeInMask',true);
    end
end
imgRefIni=proxy{1};                          % enhanced reference proxy

% ---- acquisition grouping: representatives vs reused siblings --------------
% Files that are the same acquisition but a different LSCI product differ only in
% the _c_K_ / _t_K_ / _s_K_ token (normalised in normNames).  In silent mode the
% first file of each acquisition (in procOrder, so the _t product leads) is its
% representative and is registered; its siblings reuse the representative's
% transform.  matchIdx reproduces the old "reused from file N" report exactly -
% the smallest already-registered input index of the same acquisition.  In GUI
% mode there is no reuse: every non-empty file is its own representative.
repOf=zeros(size(fNames));                   % representative fidx for each file
matchIdx=zeros(size(fNames));                % reuse-source file# for the report (0=none)
assigned=false(size(fNames));                % set true as procOrder is walked
for oidx=1:numel(procOrder)
    fidx=procOrder(oidx);
    if isempty(fNames{fidx}), continue; end
    m=0;
    if silentMode
        for j=1:numel(fNames)
            if assigned(j) && strcmp(normNames{fidx},normNames{j})
                m=j; break;                  % first already-registered file of this acquisition
            end
        end
    end
    if m>0
        repOf(fidx)=repOf(m);                % chain to the true representative
        matchIdx(fidx)=m;
    else
        repOf(fidx)=fidx;                    % representative (the reference included)
    end
    assigned(fidx)=true;
end

% ---- delegate candidate generation + selection to registerToReference ------
% enhance=false: the proxies above are already registration-ready and are reused
% for the report overlays, so the engine must not re-enhance.  The engine forms
% the identity / intensity (imregtform) / correlation (imregcorr) candidates and
% applies the same silent auto-pick / forceMethod / rotationLimit logic (or the
% GUI 4-panel chooser with the manual fallback), returning modern *tform2d
% objects + per-image diagnostics.
repList=find(~cellfun(@isempty,fNames(:)) & repOf(:)==(1:numel(fNames))');
repList=[1; repList(repList~=1)];            % reference first
sEng=struct();
if silentMode, sEng.mode='silent'; else, sEng.mode='gui'; end
sEng.enhance=false;
sEng.allowManual=true;
switch forceMethod
    case 'intensity',   sEng.method='intensity';
    case 'correlation', sEng.method='correlation';
    otherwise,          sEng.method='auto';  % '' / unrecognised -> smallest Delta
end
sEng.rotationLimit=rotationLimit;
if isfield(s,'tFormType'), sEng.tFormType=s.tFormType; end
if isfield(s,'optimizer'), sEng.optimizer=s.optimizer; end
if isfield(s,'metric'),    sEng.metric=s.metric;       end
[repTforms,repDiag]=registerToReference(proxy(repList),sEng);

repPos=zeros(size(fNames));                  % fidx -> its row in repList / repDiag
for m=1:numel(repList)
    tforms{repList(m)}=repTforms{m};
    repPos(repList(m))=m;
end
for fidx=1:numel(fNames)                     % siblings inherit the representative's transform
    if ~isempty(fNames{fidx}) && repOf(fidx)~=fidx
        tforms{fidx}=tforms{repOf(fidx)};
    end
end

% ---- per-file registration report overlays (silent mode only) --------------
% The file's warped proxy over the reference proxy, titled with the chosen method
% and the Delta(s) - same content, now on the shared canvas.  Representatives report the
% engine's choice (from diag); reused siblings report "reused from file N".
if silentMode
    nameLong={'original','intensity registration','correlation registration'};
    for oidx=1:numel(procOrder)
        fidx=procOrder(oidx);
        if isempty(fNames{fidx}) || fidx==1, continue; end
        if matchIdx(fidx)>0
            warped=imwarp(proxy{fidx},tforms{fidx},"OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
            if matchIdx(fidx)==1
                methodStr='original (same acquisition as reference)';
            else
                methodStr=sprintf('reused from file %d (same acquisition)',matchIdx(fidx));
            end
            deltaStr=sprintf('\\Delta applied=%d  (intensity/correlation not evaluated)', ...
                round(sum(abs(imgRefIni(:)-warped(:)))));
        else
            d=repDiag(repPos(fidx));
            warped=d.warped;
            mi=find(strcmp(d.method,{'original','intensity','correlation'}),1);
            methodStr=sprintf('%s (%s)',nameLong{mi},d.reason);
            r2=''; if d.rejected(2), r2=' (rot>lim)'; end
            r3=''; if d.rejected(3), r3=' (rot>lim)'; end
            deltaStr=sprintf('\\Delta original=%d, intensity=%d%s, correlation=%d%s', ...
                round(d.delta(1)),round(d.delta(2)),r2,round(d.delta(3)),r3);
        end

        % reporting overlay: warped image on the reference + chosen method/Delta.
        % This loop runs before the per-file save loop that arms the banner, so the
        % page is named from ITS file explicitly.
        fRep=reportFigure(rep,'registration','single');
        tlR=tiledlayout(fRep,1,1,'TileSpacing','compact','Padding','compact');
        axR=nexttile(tlR);
        imshowpair(rot90(imgRefIni , size(imgRefIni,2) > size(imgRefIni,1)), ...
            rot90(warped , size(warped,2) > size(warped,1)),'Parent',axR)
        axis(axR,'image')
        title(axR,{['Silent registration: ',methodStr],deltaStr})
        reportSave(rep,fRep,'registration',fNames{fidx});
    end
end

for fidx=1:1:size(fNames,1)
    if ~isempty(fNames{fidx})
        if ~(fidx==1)
            commonMask{fidx}=imwarp(commonMask{fidx},tforms{fidx},"nearest","OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
        end
    end

end


commonMask=cat(3,commonMask{:});
commonMask=round(commonMask);
commonMaskSoft=commonMask;
commonMaskSoft(commonMaskSoft(:)==4)=5;
commonMaskSoftest=commonMask;
commonMaskSoftest(commonMaskSoftest(:)==2)=1;
commonMaskSoftest(commonMaskSoftest(:)==4)=5;
commonMask = (max(commonMask,[],3) == min(commonMask,[],3)) & mean(commonMask,3)>0;
commonMaskSoft = (max(commonMaskSoft,[],3) == min(commonMaskSoft,[],3)) & mean(commonMaskSoft,3)>0;
commonMaskSoftest = (max(commonMaskSoftest,[],3) == min(commonMaskSoftest,[],3)) & mean(commonMaskSoftest,3)>0;

for fidx=1:1:size(fNames,1)
    if reportCancelled(rep), break; end          % cooperative cancel between files

    if ~isempty(fNames{fidx})
        reportFile(rep,fidx,fNames{fidx});
        load(fNames{fidx},'source');
        load(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings');
        load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');
        img=mean(source.data,3);


        if fidx==1
            if s.matchSegmentation
                results.sMetrics.('RegID')=NaN(height(results.sMetrics),1);
                results.sMetrics.('RegOverlap')=NaN(height(results.sMetrics),1);

                if isfield(results,"dvsMetrics")
                    results.dvsMetrics.('RegID')=NaN(height(results.dvsMetrics),1);
                    results.dvsMetrics.('RegOverlap')=NaN(height(results.dvsMetrics),1);
                end

                idxs=results.sMetrics.category==1;
                results.sMap(ismember(results.sMap, results.sMetrics.idx(idxs))) = 0;
                results.sData(:,idxs)=[];
                results.sMetrics(idxs,:)=[];

                sMap=results.sMap;
                [M,N]       = size(results.cMask);                  % image dimensions
                step        = double(s.prchNSize);                  % nominal cell width
                R           = step/2;                               % half-offset
                rows        = 1 : R*sqrt(3) : M;                    % √3·R vertical pitch
                cols        = 1 : step         : N;
                [C,Rr]      = meshgrid(cols,rows);                  % full grid of centres
                C(2:2:end,:)= C(2:2:end,:) + R;                     % shift odd rows by R
                C           = round(C);  Rr = round(Rr);            % integer coordinates
                inFrame     =  Rr>=1 & Rr<=M & C>=1 & C<=N;         % guard ①
                idxAll      = sub2ind([M N], Rr(inFrame), C(inFrame));
                idxSeeds    = idxAll(results.cMask(idxAll) >= 0);    % guard ②: inside mask
                seed        = false(M,N);  seed(idxSeeds) = true;   % binary seed image
                [~,lbl]     = bwdist(seed,'euclidean');             % nearest-seed label
                valid       = results.cMask>=1;     % area to overwrite
                lbl(~valid) = 0;
                nz                  = lbl>0;
                [~,~,lbl(nz)]       = unique(lbl(nz));              % consecutive IDs
                lbl                 = int32(lbl) + max(sMap(:));   % avoid clashes
                sMap(valid & results.cMask==1)  = lbl(valid & results.cMask==1);                   % updated label map
                results.sMap=sMap;

                data=reshape(source.data,[],size(source.data,3));
                for i=unique(results.sMap(results.cMask==1))'
                    area=sum(results.sMap(:)==i);
                    c=unique(results.cMask(results.sMap==i));
                    results.sMetrics(i,:)={i,c,NaN,NaN,NaN,NaN,area,NaN,NaN,NaN};
                    results.sData(:,i)=mean(data(results.sMap(:)==i,:),1,'omitnan');
                end
                clearvars data

                results.sMetrics.('RegID')=results.sMetrics.('idx');
                results.sMetrics.('RegOverlap')=results.sMetrics.('area');

                if isfield(results,"dvsMetrics")
                    results.dvsMetrics.('RegID')=results.dvsMetrics.('idx');
                    results.dvsMetrics.('RegOverlap')=results.sMetrics.area(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
                end

                cMaskRef=results.cMask;
                sMapRef=results.sMap;
                sMapRef2=sMapRef;
                for i=unique(sMapRef(cMaskRef==5))'
                    sMapRef2(sMapRef>=i & sMapRef<=i+2)=i;
                end

            end
        else
            fn = fieldnames(results);
            for k=1:numel(fn)
                if size(results.(fn{k}),1)==size(img,1) & size(results.(fn{k}),2)==size(img,2)
                    tmp=zeros(size(imgRefIni,1),size(imgRefIni,2),size(results.(fn{k}),3),class(results.(fn{k})));
                    for i=1:1:size(results.(fn{k}),3)
                        if ~(isa(tmp,'double') || isa(tmp,'single'))
                            tmp(:,:,i)=imwarp(results.(fn{k})(:,:,i),tforms{fidx},"nearest","OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
                        else
                            tmp(:,:,i)=imwarp(results.(fn{k})(:,:,i),tforms{fidx},"OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
                        end
                    end
                    results.(fn{k})=tmp;
                end
            end

            fn = fieldnames(source);
            for k=1:numel(fn)
                if size(source.(fn{k}),1)==size(img,1) & size(source.(fn{k}),2)==size(img,2)
                    tmp=zeros(size(imgRefIni,1),size(imgRefIni,2),size(source.(fn{k}),3),class(source.(fn{k})));
                    for i=1:1:size(source.(fn{k}),3)
                        if ~(isa(tmp,'double') || isa(tmp,'single'))
                            tmp(:,:,i)=imwarp(source.(fn{k})(:,:,i),tforms{fidx},"nearest","OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
                        else
                            tmp(:,:,i)=imwarp(source.(fn{k})(:,:,i),tforms{fidx},"OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
                        end
                    end
                    source.(fn{k})=tmp;
                end
            end


            if s.matchSegmentation
                results.sMetrics.('RegID')=NaN(height(results.sMetrics),1);
                results.sMetrics.('RegOverlap')=NaN(height(results.sMetrics),1);

                if isfield(results,"dvsMetrics")
                    results.dvsMetrics.('RegID')=NaN(height(results.dvsMetrics),1);
                    results.dvsMetrics.('RegOverlap')=NaN(height(results.dvsMetrics),1);
                end

                idxs=results.sMetrics.idx(results.sMetrics.category==5);
                for i=1:1:numel(idxs)
                    lbl=idxs(i);
                    BW=results.sMap>=lbl & results.sMap<=lbl+2;
                    m=mode(sMapRef2(BW(:)==1 & cMaskRef(:)>=3));

                    if m>0
                        results.sMetrics.RegID(results.sMetrics.idx==lbl)=m;
                        results.sMetrics.RegID(results.sMetrics.idx==lbl+1)=m+1;
                        results.sMetrics.RegID(results.sMetrics.idx==lbl+2)=m+2;
                        results.sMetrics.RegOverlap(results.sMetrics.idx==lbl)=sum(sMapRef(:)==m & results.sMap(:)==lbl);
                        results.sMetrics.RegOverlap(results.sMetrics.idx==lbl+1)=sum(sMapRef(:)==(m+1) & results.sMap(:)==(lbl+1));
                        results.sMetrics.RegOverlap(results.sMetrics.idx==lbl+2)=sum(sum(sMapRef(:)==(m+2) & results.sMap(:)==(lbl+2)));

                        if isfield(results,"dvsMetrics")
                            results.dvsMetrics.RegID(results.dvsMetrics.idx==lbl)=m;
                            results.dvsMetrics.RegOverlap(results.dvsMetrics.idx==lbl)=sum(sMapRef(:)==m & results.sMap(:)==lbl);
                        end
                    end
                end

                idxs=results.sMetrics.category==1;
                results.sMap(ismember(results.sMap, results.sMetrics.idx(idxs))) = 0;
                results.sData(:,idxs)=[];
                results.sMetrics(idxs,:)=[];

                sMap=results.sMap;
                [M,N]       = size(results.cMask);                  % image dimensions
                step        = double(s.prchNSize);                  % nominal cell width
                R           = step/2;                               % half-offset
                rows        = 1 : R*sqrt(3) : M;                    % √3·R vertical pitch
                cols        = 1 : step         : N;
                [C,Rr]      = meshgrid(cols,rows);                  % full grid of centres
                C(2:2:end,:)= C(2:2:end,:) + R;                     % shift odd rows by R
                C           = round(C);  Rr = round(Rr);            % integer coordinates
                inFrame     =  Rr>=1 & Rr<=M & C>=1 & C<=N;         % guard ①
                idxAll      = sub2ind([M N], Rr(inFrame), C(inFrame));
                idxSeeds    = idxAll(results.cMask(idxAll) >= 0);    % guard ②: inside mask
                seed        = false(M,N);  seed(idxSeeds) = true;   % binary seed image
                [~,lbl]     = bwdist(seed,'euclidean');             % nearest-seed label
                valid       = results.cMask>=1;     % area to overwrite
                lbl(~valid) = 0;
                nz                  = lbl>0;
                [~,~,lbl(nz)]       = unique(lbl(nz));              % consecutive IDs
                lbl                 = int32(lbl) + max(sMap(:));   % avoid clashes
                sMap(valid & results.cMask==1)  = lbl(valid & results.cMask==1);                   % updated label map
                results.sMap=sMap;

                data=reshape(source.data,[],size(source.data,3));
                for i=unique(results.sMap(results.cMask==2 | results.cMask==1))'
                    area=sum(results.sMap(:)==i);
                    c=unique(results.cMask(results.sMap==i));
                    [m,ma]=mode(sMapRef(results.sMap(:)==i & cMaskRef(:)==c));
                    results.sMetrics(i,:)={i,c,NaN,NaN,NaN,NaN,area,NaN,m,ma};
                    results.sData(:,i)=mean(data(results.sMap(:)==i,:),1,'omitnan');
                end
                clearvars data
            end
        end

        results.commonMask=cat(3,commonMask,commonMaskSoft,commonMaskSoftest) ;
        s.tForm=tforms{fidx};
        s.view=views{fidx};
        s.imgRefIni=imgRefIni;
        settings.runRegistration=reportSettings(s);
        %Save the data
        reportWriting(rep);
        save(fNames{fidx},'source','-v7.3');
        save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
        save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');
        reportSaved(rep);
    end
end
reportClose(rep);

end