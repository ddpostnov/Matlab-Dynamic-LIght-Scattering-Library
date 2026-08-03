%setVesselTypes  Semi-automatic vessel/ROI labelling for *_BFI_r.mat datasets
%
%   setVesselTypes(s,fNames) opens an interactive GUI for each BFI dataset in
%   fNames, displays the BFI image together with a provisional artery/vein
%   guess (derived from psPI heuristics), and lets the user paint or
%   relabel regions.  The routine then:
%       • writes the chosen type ("Artery", "Vein", "Uncertain"),
%         free-text label, and confidence score into RESULTS.sMetrics
%         (and dvsMetrics, if present)
%       • stores the signed confidence map in RESULTS.mapType (NaN outside
%         mask, −maxGuess … +maxGuess in vessels)
%       • updates *_BFI_r.mat and *_BFI_s.mat in place
%
%   When s.useReference is true the first file in fNames is treated as the
%   reference: subsequent files inherit its type/label definitions via
%   REGID matching, skipping the GUI.
%
%   INPUTS
%     s        parameter struct
%                • useReference   true / false
%                • refFName       path to reference file (if used)
%                • prchNSize      grid spacing for parenchyma fill-in
%     fNames   cell array of *_BFI_r.mat paths.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail) and
%     s.cancelFcn()->tf (between files).
%
%   SIDE-EFFECTS
%       *_BFI_r.mat   RESULTS  – fields mapType, .type, .label, etc.
%       *_BFI_s.mat   SETTINGS – sub-field settings.setVesselTypes added
%       *_rep_vesseltypes.jpg  the record of the labelling: the artery/vein map over
%                     the BFI image with the NAMED vessels outlined and annotated.
%                     Written for EVERY file - under s.useReference only the first
%                     opens the GUI, and confirming that the inherited labelling
%                     landed correctly on the rest is what the page is for.
%
%   EXAMPLE
%     p.useReference = false;
%     D = dir(fullfile(dataRoot,'*_BFI_r.mat'));
%     setVesselTypes(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     enhanceForDisplay (background-subtracted preview), Image Processing Toolbox
%     (bwskel, visboundaries, etc.) and core LSCI utility functions.
%
% See also: setRegions, runSegmentation, runPulsatility, runDynamicSegmentation
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
% %If no reference is used
% % s.useReference=false;
% % s.refFName='';
% %If reference is used without proper file naming
% % s.useReference=true; %ASSUMES PRE-REGISTERED FILES
% % s.refFName=REFERENCE FILE NAME;
% %If reference is used in an automated setting
% s.useReference=true; %Assumes PRE-registered files
% %Settings tabs with a reference has to be done ROI by ROI if split ROIS are used. It also has
% %to be done ANIMAL by ANIMAL if multiple animals are compared. It can be
% %convienintly done in a loop as below, but REQUIRES proper file naming.
% for idxA=5 %animals index list
%     for idxR=1:2 %ROIs index list
%         files      = dir(fullfile(rootFolder,'**',sprintf('Roi%d*LH%03d*c_BFI_r.mat', idxR, idxA))); %<---ALWAYS REFER TO "_K_r.mat" files, but you may use regexp to define specific "_K_r.mat" files of interest
%         fNames     = fullfile({files.folder}', {files.name}');
%         s.refFName=fNames{1};
%         setVesselTypes(s,fNames)
%     end
% end



function setVesselTypes(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_BFI_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_BFI_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_c_BFI_r.mat'').']);
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  This step is
% interactive (paint GUI), so only stageFcn (file boundaries) and cancelFcn (between
% files) are wired - progress is not threaded through the GUI.
rep=reportOpen(s,'Vessel types',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');
        if (fidx==1 &&  s.useReference ) || ~s.useReference
            idxs=results.sMetrics{:,'idx'};
            ctgs=results.sMetrics{:,'category'};
            bfi=results.sMetrics{:,'BFI'};
            guess=zeros(size(idxs));
            cMask=results.cMask;
            % maxGuess (the magnitude of a hand-set label, kept beyond the auto-guess
            % range) is computed from the number of guess sliders once they are built.
            N=numel(idxs);
            % Valid category-5 rows the auto-guess may write to (need a 2-endpoint
            % span idxs(i):idxs(i)+1 that stays inside the guess vector).
            m5 = ~isnan(idxs) & idxs>=1 & (idxs+1)<=N & (ctgs==5);
            ii = find(m5);                  % shared by all guess sliders
            vn = results.sMetrics.Properties.VariableNames;   % available feature columns

            % ----- guess sliders 1-11 (computed once) -------------------------------
            % The automatic guess is produced by range sliders (see recomputeGuess);
            % each carries a per-segment value, a full range and initial thumbs, and
            % is sampled at the segment's first endpoint idxs(i):
            %   sld1     psPI gradient (this endpoint - next)
            %   sld2     psPI magnitude
            %   sld3     BFI / diameter ratio
            %   sld4     std(BFI) / BFI ratio
            %   sld5-11  single-column "magnitude" metrics (built by metricSlider):
            %            the time-of-* metrics are inverted.
            % A slider whose source column is absent is inert and not shown. Sliders
            % 1-2 start ON when psPI exists; the rest start OFF. Default thumbs:
            % sld1 at 0, sld2 at the PI medians around thrsh, sld3 at the data
            % quartiles, and sld4-11 at the metric's parenchymal (category-1) median.
            sld1 = struct('seg',nan(N,1),'rng',[-1 1],'lo',0,'hi',0);
            sld2 = struct('seg',nan(N,1),'rng',[ 0 1],'lo',1/3,'hi',2/3);
            sld3 = struct('seg',nan(N,1),'rng',[ 0 1],'lo',1/3,'hi',2/3);
            sld4 = struct('seg',nan(N,1),'rng',[ 0 1],'lo',1/3,'hi',2/3);

            % sliders 1 & 2: psPI gradient and magnitude
            hasPI = any(strcmp('psPI',vn));
            if hasPI
                var=results.sMetrics{:,"psPI"};
                thrsh=median(var(bfi<prctile(bfi(ctgs==5),10) & ctgs==5),'omitnan');
                dv = nan(N,1);  mv = nan(N,1);
                dv(ii) = var(idxs(ii)) - var(idxs(ii)+1);   % this endpoint - next
                mv(ii) = var(idxs(ii));                      % psPI magnitude here
                sld1.seg = dv;  sld2.seg = mv;

                % slider 1: symmetric range spanning the endpoint psPI differences
                % (the category-5 vs category-3 difference along a segment), thumbs at 0
                R1 = max(abs(dv),[],'omitnan');
                if ~isfinite(R1) || R1==0, R1 = 1; end
                sld1.rng = [-R1 R1];  sld1.lo = 0;  sld1.hi = 0;

                % slider 2: range over category-5 psPI; default thumbs at the
                % medians of category-5 psPI below / above thrsh
                lo2 = median(var(ctgs==5 & var<thrsh),'omitnan');
                hi2 = median(var(ctgs==5 & var>thrsh),'omitnan');
                mvLo = min(mv,[],'omitnan');  mvHi = max(mv,[],'omitnan');
                if ~isfinite(mvLo) || ~isfinite(mvHi) || mvLo==mvHi, mvLo=0; mvHi=1; end
                if ~isfinite(lo2), lo2 = mvLo + (mvHi-mvLo)/3;   end
                if ~isfinite(hi2), hi2 = mvLo + 2*(mvHi-mvLo)/3; end
                pad = 0.02*(mvHi-mvLo);
                sld2.rng = [min(mvLo,lo2)-pad, max(mvHi,hi2)+pad];
                sld2.lo  = min(max(lo2,sld2.rng(1)),sld2.rng(2));
                sld2.hi  = min(max(hi2,sld2.rng(1)),sld2.rng(2));
            end

            % slider 3: BFI / diameter ratio (enable-able only if both columns exist)
            has3 = any(strcmp('BFI',vn)) && any(strcmp('diameter',vn));
            if has3
                diam = results.sMetrics{:,'diameter'};
                r3 = nan(N,1);
                r3(ii) = bfi(idxs(ii)) ./ diam(idxs(ii));
                r3(~isfinite(r3)) = NaN;                 % guard divide-by-zero / missing
                sld3.seg = r3;
                [sld3.rng, sld3.lo, sld3.hi] = autoRangeThumbs(r3, m5);
            end

            % slider 4: std(BFI) / BFI ratio; default thumbs at the parenchymal median
            has4 = any(strcmp('BFI',vn)) && any(strcmp('std(BFI)',vn));
            if has4
                sbfi = results.sMetrics{:,"std(BFI)"};
                r4 = nan(N,1);
                r4(ii) = sbfi(idxs(ii)) ./ bfi(idxs(ii));
                r4(~isfinite(r4)) = NaN;
                rP = sbfi(ctgs==1) ./ bfi(ctgs==1);          % ratio over parenchyma
                med4 = median(rP(isfinite(rP)));
                sld4.seg = r4;
                [sld4.rng, sld4.lo, sld4.hi] = thumbsAtValue(r4, m5, med4);
            end

            % sliders 5-11: single-column magnitude sliders (parenchymal-median
            % default). The time-of-* metrics are inverted.
            % An absent column yields an inert slider (.has == false), so it won't show.
            sld5  = metricSlider('psTimeMin',   true);
            sld6  = metricSlider('psTimeMax',   true);
            sld7  = metricSlider('t0B',         true);
            sld8  = metricSlider('t5B',         true);
            sld9  = metricSlider('tUB',         true);
            sld10 = metricSlider('t90',         true);
            sld11 = metricSlider('tMB',         true);

            % A slider is only built/shown when the variable it relies on is present;
            % absent ones are omitted entirely (not just disabled). Sliders 1-2 share
            % psPI; the rest each rely on their own column.
            present = [hasPI hasPI has3 has4 ...
                sld5.has sld6.has sld7.has sld8.has sld9.has sld10.has sld11.has];
            nVis = nnz(present);

            map=[0;guess];
            map=map(results.sMap+1);
            map(cMask<4)=NaN;
            cmap=vesselTypeColormap();   % blue (vein) -> green (0) -> red (artery)
            img=sqrt(results.imgBFI);
            img=mat2gray(img,double(prctile(img(cMask(:)>0),[5,99])));
            fSize=floor((min(size(img))./20))*2+1;
            img=enhanceForDisplay(img,fSize).*(cMask>0);
            rois=repmat("",1,numel(guess));


            f  = figure('Name','Guided labeling','Color','w','WindowState','maximized', ...
                'CloseRequestFcn',@(~,~)uiresume(gcbf));
            TL = tiledlayout(f,1,5,'TileSpacing','none','Padding','compact');
            axL = nexttile(TL,[1 2]);
            imagesc(axL,img), axis(axL,'image','off'), colormap(axL,'parula');
            clim(axL,prctile(img(cMask(:)>0),[1,99]));
            title('BFI');
            axR = nexttile(TL,[1 2]);
            hMap = imagesc(axR,map,'AlphaData',~isnan(map));
            axis(axR,'image','off'), setRightClim(axR,f), colormap(axR,cmap);
            title('Labels');
            axCtl = nexttile(TL); axis(axCtl,'off');
            axis tight
            panelWH  = [150 300];
            sliderWH = [150 max(1,nVis)*63];       % one ~63px band per visible slider
            panelGap = 12;                         % gap between slider and control panels
            p    = uipanel(f,'Units','pixels','BorderType','none','BackgroundColor','w');
            pSld = uipanel(f,'Units','pixels','BorderType','none','BackgroundColor','w');
            repositionPanels(f,axCtl,pSld,sliderWH,p,panelWH,panelGap)
            f.SizeChangedFcn = @(~,~)repositionPanels(f,axCtl,pSld,sliderWH,p,panelWH,panelGap);
            uicontrol(p,'Style','text','String','Left panel image', ...
                'Units','normalized','Position',[.05 .90 .90 .08], ...
                'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
            pop = uicontrol(p,'Style','popup','Units','normalized', ...
                'Position',[.05 .81 .90 .08]);
            uicontrol(p,'Style','text','String','Labelling controls', ...
                'Units','normalized','Position',[.05 .7 .90 .06], ...
                'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
            btnY = [.60 .50 .40 .30 .50 .40];
            btnX = [.05 .05 .05 .05 .55 .55];
            lbl   = {'Artery','Vein','Uncertain','Clear','ROI Vessel','ROI Parench'};
            val = [1 -1 0 9 7 8];   % action codes: +1 artery, -1 vein, 0 uncertain, 9 clear, 7/8 ROI
            btn   = gobjects(1,numel(lbl));
            for k = 1:numel(lbl)
                btn(k) = uicontrol(p,'String',lbl{k},'Units','normalized', ...
                    'Position',[btnX(k) btnY(k) .45 .08], ...
                    'Callback',@(src,~)pickCat(src,val(k)));
            end
            txt=uicontrol(p,'Style','edit','String','Label?', ...
                'Units','normalized','Position',[.55 .6 .45 .08], ...
                'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
            uicontrol(p,'String','Remove ROI','Units','normalized', ...
                'Position',[.55 .30 .45 .08], ...
                'Callback',@(src,~)removeROI(src));
            hoverTxt = uicontrol(p,'Style','text','String','Label: (hover a segment)', ...
                'Units','normalized','Position',[.05 .20 .90 .08], ...
                'HorizontalAlignment','left','BackgroundColor',p.BackgroundColor);
            finishBtn = uicontrol(p,'String','Finish','Units','normalized', ...
                'Position',[.05 .00 .90 .08], ...
                'Callback',@(~,~)uiresume(gcbf));
            setappdata(f,'btnHandles',btn);
            setappdata(f,'txtHandle',txt);
            setappdata(f,'popHandle',pop);     % needed for left-panel redraw / de-labelling
            setappdata(f,'labelOrder',strings(1,0));   % vessel labels, in order of first use
            setappdata(f,'hoverTxtHandle',hoverTxt);    % label inspector (updated on hover)

            % ---------- build selector list (2-D fields == size(map)) ---------------
            names={}; data={};  chk=@(v) isnumeric(v)&&ismatrix(v)&&isequal(size(v),size(map));
            for fn = fieldnames(results).',                  v=results.(fn{1});
                if chk(v), names{end+1}=fn{1}; data{end+1}=v; end, end
            if isfield(results,'extendedMetrics')
                for fn = fieldnames(results.extendedMetrics).', v=results.extendedMetrics.(fn{1});
                    if chk(v), names{end+1}=['extendedMetrics.' fn{1}]; data{end+1}=v; end, end
            end
            % relocated per-pixel pulsatility maps (results.pulsatility.ppx.scalars.*):
            % the [Y x X] marker maps (psPI/psRI/psTimeMin/... - the hAmp/hPhase cubes
            % are 3-D and skipped by chk) that used to live as top-level imgPI/imgRI.
            if isfield(results,'pulsatility') && isfield(results.pulsatility,'ppx') ...
                    && isfield(results.pulsatility.ppx,'scalars')
                for fn = fieldnames(results.pulsatility.ppx.scalars).', v=results.pulsatility.ppx.scalars.(fn{1});
                    if chk(v), names{end+1}=['pulsatility.ppx.scalars.' fn{1}]; data{end+1}=v; end, end
            end
            idx0 = find(strcmp(names,'imgBFI'),1);   % try to start at "imgBFI"
            if isempty(idx0), idx0 = 1; end
            pop.UserData = struct('ax',axL,'data',{data},'names',{names});
            set(pop,'String',names,'Value',idx0,'Callback',@updateImg);   % ← set start

            % ---------- range sliders (only those whose variable is present) -------
            % Each slider: a variable-name label with an On/Off switch beside it, a
            % two-thumb track (right thumb constrained >= left, splitting the range
            % into 3 segments), and a Reset button. All shown sliders feed the guess;
            % 1-2 start ON (psPI present), the rest start OFF. Sliders whose source
            % column is missing are not built at all.
            sldNames   = {'psPI: this - next','psPI: magnitude', ...
                'BFI / diameter','std(BFI) / BFI', ...
                '-psTimeMin','-psTimeMax', ...
                '-t0B','-t5B','-tUB','-t90','-tMB'};
            sldRange   = {sld1.rng,sld2.rng,sld3.rng,sld4.rng,sld5.rng,sld6.rng,sld7.rng,sld8.rng,sld9.rng,sld10.rng,sld11.rng};
            sldLo      = {sld1.lo, sld2.lo, sld3.lo, sld4.lo, sld5.lo, sld6.lo, sld7.lo, sld8.lo, sld9.lo, sld10.lo, sld11.lo};
            sldHi      = {sld1.hi, sld2.hi, sld3.hi, sld4.hi, sld5.hi, sld6.hi, sld7.hi, sld8.hi, sld9.hi, sld10.hi, sld11.hi};
            sldSeg     = {sld1.seg,sld2.seg,sld3.seg,sld4.seg,sld5.seg,sld6.seg,sld7.seg,sld8.seg,sld9.seg,sld10.seg,sld11.seg};
            sldGuess   = true(1,11);                                          % all feed the guess
            sldEnabled = [hasPI hasPI false false false false false false false false false]; % on by default?
            sliders = struct([]);                  % only present sliders get built
            vi = 0;                                % visible-row counter
            for k = 1:numel(sldNames)
                if ~present(k), continue, end      % variable missing -> not visualised
                vi = vi + 1;
                bandH = 1/nVis;
                yBand = 1 - vi*bandH;              % top row first, by visible position
                rs = addRangeSlider(f,pSld,[0 yBand 1 bandH],sldNames{k}, ...
                    sldRange{k},sldLo{k},sldHi{k},sldEnabled(k),false);
                rs.guesses = sldGuess(k);          % whether it contributes to the guess
                rs.segVals = sldSeg{k};            % per-segment variable (or [])
                if vi==1, sliders = rs; else, sliders(vi) = rs; end
            end
            setappdata(f,'sliders',sliders);
            setappdata(f,'activeSlider',[]);
            setappdata(f,'hoverSlider',[]);          % range slider currently under the cursor

            % A hand-set label is stored at +/-maxGuess so it is always one step
            % beyond the largest possible auto-guess (= number of guess sliders shown).
            if isempty(sliders), maxGuess = 1; else, maxGuess = nnz([sliders.guesses]) + 1; end

            % ---------- share data for painting routine -----------------------------
            guidata(f,struct('map',map,'cMask',cMask,'guess',guess,'maxGuess',maxGuess,'rois',rois,'sMap',results.sMap,'hMap',hMap,'axL',axL,'axR',axR, ...
                'userSet',false(size(guess)),'idxs',idxs,'m5',m5));
            setappdata(f,'paintVal',[]);
            set(f,'WindowButtonDownFcn',@startPaint,'WindowButtonUpFcn',@stopPaint, ...
                'WindowButtonMotionFcn',@onHover, ...     % hover: label inspector + slider focus
                'WindowKeyPressFcn',@onSliderKey, ...     % arrows nudge the hovered slider
                'WindowScrollWheelFcn',@onSliderScroll);  % scroll nudges the hovered slider

            recomputeGuess(f);                       % initial guess from the sliders

            uiwait(f);                               % pause until Finish pressed
            guess = guidata(f).guess;
            rois = guidata(f).rois;
            delete(f);    % retrieve & close

            map=[0;guess];
            map=map(results.sMap+1);
            map(cMask<3)=NaN;
            results.mapType=map;


            type = repmat("Uncertain",size(guess));  % default = 0
            type(guess < 0) = "Vein";                % negative values
            type(guess > 0) = "Artery";
            idxs=results.sMetrics.idx;

            results.sMetrics.type=type;
            results.sMetrics.label=rois';
            results.sMetrics.typeConfidence=guess;
            if isfield(results,"dvsMetrics")
                results.dvsMetrics.label=rois(ismember(results.sMetrics.idx, results.dvsMetrics.idx))';
                results.dvsMetrics.type=type(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
                results.dvsMetrics.typeConfidence=guess(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
            end
        else
            results.mapType=map;
            results.sMetrics.type = strings(height(results.sMetrics),1);
            results.sMetrics.label = strings(height(results.sMetrics),1);
            results.sMetrics.typeConfidence = zeros(height(results.sMetrics),1);
            if isfield(results,"dvsMetrics")
                results.dvsMetrics.type = strings(height(results.dvsMetrics),1);
                results.dvsMetrics.label = strings(height(results.dvsMetrics),1);
                results.dvsMetrics.typeConfidence = zeros(height(results.dvsMetrics),1);
            end

            if ismember('RegID', results.sMetrics.Properties.VariableNames)
                [isHit, loc]   = ismember(results.sMetrics.RegID, idxs);
                results.sMetrics.type(isHit) = type(loc(isHit));
                results.sMetrics.typeConfidence(isHit) = guess(loc(isHit));
                results.sMetrics.label(isHit) = rois(loc(isHit))';

                if isfield(results,"dvsMetrics")
                    [isHit, loc]   = ismember(results.dvsMetrics.RegID, idxs);
                    results.dvsMetrics.type(isHit) = type(loc(isHit));
                    results.dvsMetrics.typeConfidence(isHit) = guess(loc(isHit));
                    results.dvsMetrics.label(isHit) = rois(loc(isHit))';
                end
            else
                % no registration ID: the segmentation map is assumed identical
                results.sMetrics.type=type;
                results.sMetrics.label=rois';
                results.sMetrics.typeConfidence=guess;
                if isfield(results,"dvsMetrics")
                    results.dvsMetrics.label=rois(ismember(results.sMetrics.idx, results.dvsMetrics.idx))';
                    results.dvsMetrics.type=type(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
                    results.dvsMetrics.typeConfidence=guess(ismember(results.sMetrics.idx, results.dvsMetrics.idx));
                end

            end
        end

        % The record of the labelling, for EVERY file.  Under s.useReference only
        % the first file opens the GUI and the rest inherit the labels - and
        % checking that an inherited labelling landed correctly on file 7 is
        % precisely what this page is for.
        writeVesselTypesReport(rep,s.fName,results);

        settings.setVesselTypes=reportSettings(s);
        reportWriting(rep);
        save(fNames{fidx},'results','-v7.3','-nocompression');
        save(getProductPath(fNames{fidx},'s'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);


%% =========================  LOCAL FUNCTIONS  ========================= %%
    function startPaint(fig,~)
        [idx,thumb] = hitSlider(fig);            % did we grab a range slider?
        if ~isempty(idx)
            sld = getappdata(fig,'sliders');
            if ~sld(idx).enabled, return, end    % disabled slider -> ignore
            setappdata(fig,'activeSlider',struct('idx',idx,'thumb',thumb));
            set(fig,'WindowButtonMotionFcn',@dragSlider);
            dragSlider(fig);                     % jump the thumb to the click
            return
        end
        set(fig,'WindowButtonMotionFcn',@doPaint);
        doPaint(fig);
    end
    function stopPaint(fig,~)
        setappdata(fig,'activeSlider',[]);           % end any slider drag
        set(fig,'WindowButtonMotionFcn',@onHover);   % resume hover inspection
    end
    function [idx,thumb] = hitSlider(fig)
        % Index of the range slider under the cursor ([] if none) and which
        % thumb (1 = low, 2 = high) should move. Picks the nearer thumb; ties
        % break by the side of the value that was clicked.
        idx = []; thumb = [];
        sld = getappdata(fig,'sliders');
        if isempty(sld), return, end
        ax = ancestor(hittest(fig),'axes');
        if isempty(ax), return, end
        for i = 1:numel(sld)
            if ax==sld(i).ax, idx = i; break, end
        end
        if isempty(idx), return, end
        cp = get(sld(idx).ax,'CurrentPoint'); x = cp(1,1);
        if abs(x-sld(idx).lo) < abs(x-sld(idx).hi)
            thumb = 1;
        elseif abs(x-sld(idx).hi) < abs(x-sld(idx).lo)
            thumb = 2;
        elseif x <= sld(idx).lo
            thumb = 1;
        else
            thumb = 2;
        end
    end
    function dragSlider(fig,~)
        ad = getappdata(fig,'activeSlider');
        if isempty(ad), return, end
        sld = getappdata(fig,'sliders');
        S = sld(ad.idx);
        if ~S.enabled, return, end
        cp = get(S.ax,'CurrentPoint'); x = cp(1,1);
        x = min(max(x,S.rng(1)),S.rng(2));       % clamp to range
        if ad.thumb==1
            S.lo = min(x,S.hi);                  % left thumb stays <= right
        else
            S.hi = max(x,S.lo);                  % right thumb stays >= left
        end
        drawSlider(S);
        sld(ad.idx) = S;
        setappdata(fig,'sliders',sld);
        if S.guesses, recomputeGuess(fig); end   % live-update the guess + right image
        drawnow limitrate
    end
    function onHover(fig,~)
        % Two jobs on every mouse move (both cheap): (1) record which range
        % slider the cursor is over, so the arrow-key / scroll-wheel handlers
        % know which one to nudge; (2) update the label inspector with the
        % segment currently under the cursor, refreshing only when it changes.
        ax = ancestor(hittest(fig),'axes');

        % (1) which range slider is the cursor over? ([] if none)
        sld = getappdata(fig,'sliders');
        hov = [];
        if ~isempty(ax) && ~isempty(sld)
            for i = 1:numel(sld)
                if ax==sld(i).ax, hov = i; break, end
            end
        end
        setappdata(fig,'hoverSlider',hov);

        % (2) label inspector for the two image panels
        gui = guidata(fig);
        if ~isstruct(gui), return, end
        label = "";
        if ~isempty(ax) && (ax==gui.axL || ax==gui.axR)
            cp = get(ax,'CurrentPoint');
            x  = round(cp(1,1));  y = round(cp(1,2));
            if x>=1 && y>=1 && x<=size(gui.sMap,2) && y<=size(gui.sMap,1)
                lbl = gui.sMap(y,x);
                if lbl>=1 && lbl<=numel(gui.rois)
                    label = gui.rois(lbl);
                end
            end
        end
        if strlength(label)==0
            target = 'Label: (hover a segment)';
        else
            target = char("Label: " + label);
        end
        h = getappdata(fig,'hoverTxtHandle');
        if ~strcmp(get(h,'String'),target)
            set(h,'String',target);
            drawnow limitrate
        end
    end
    function doPaint(fig,~)
        val = getappdata(fig,'paintVal');
        if isempty(val); return; end

        gui = guidata(fig);
        if ~isstruct(gui); return; end

        % Which axes is the mouse currently over?
        axUnder = ancestor(hittest(fig),'axes');
        if isempty(axUnder) || ~(axUnder==gui.axL || axUnder==gui.axR)
            return                                   % clicked outside the images
        end

        % Pixel → matrix indices
        cp = get(axUnder,'CurrentPoint');
        x  = round(cp(1,1));
        y  = round(cp(1,2));
        if x<1 || y<1 || x>size(gui.map,2) || y>size(gui.map,1)
            return
        end

        lbl = gui.sMap(y,x);
        if lbl==0, return, end

        txt = getappdata(fig,'txtHandle');
        txt = get(txt,'String');
        isEmptyLabel = isempty(strtrim(txt));   % empty field ⇒ de-label / deselect
        if val==7 && lbl>0 && gui.cMask(y,x)==5
            if isEmptyLabel
                if all(strlength(gui.rois(lbl:lbl+1))==0), return, end  % nothing here
                gui.rois(lbl:lbl+1) = "";                              % drop label
                guidata(fig,gui);
                redrawLeft(fig);                                       % remove outline
                drawnow limitrate
            else
                gui.rois(lbl:lbl+1)=txt;
                mask = ismember(gui.sMap,lbl);
                c = vesselColor(fig,txt);
                hold(gui.axL,'on')
                visboundaries(gui.axL,mask,'Color',c,'LineWidth',2);
                hold(gui.axL,'off')
                hold(gui.axR,'on')
                set(visboundaries(gui.axR,mask,'Color',c,'LineWidth',1),'Tag','roiOutline');
                hold(gui.axR,'off')
                guidata(fig,gui);
                drawnow limitrate
            end
        elseif val==8 && lbl>0 && gui.cMask(y,x)==1
            if isEmptyLabel
                if strlength(gui.rois(lbl))==0, return, end            % nothing here
                gui.rois(lbl) = "";                                    % drop label
                guidata(fig,gui);
                redrawLeft(fig);                                       % remove outline
                drawnow limitrate
            else
                gui.rois(lbl)=txt;
                mask = (gui.sMap==lbl);
                c = parenchymaColor(fig,gui,txt);
                hold(gui.axL,'on')
                visboundaries(gui.axL,mask,'Color',c,'LineWidth',2);
                hold(gui.axL,'off')
                hold(gui.axR,'on')
                set(visboundaries(gui.axR,mask,'Color',c,'LineWidth',1),'Tag','roiOutline');
                hold(gui.axR,'off')
                guidata(fig,gui);
                drawnow limitrate
            end
        elseif val<7 && lbl>0 && gui.cMask(y,x)==5
            g = val * gui.maxGuess;              % +1 artery -> +maxGuess, -1 vein -> -maxGuess, 0 uncertain
            gui.guess(lbl:lbl+1)   = g;
            gui.userSet(lbl:lbl+1) = true;       % hand-labelled: shield from auto-guess
            gui.map(gui.sMap==lbl)   = g;
            gui.map(gui.sMap==lbl+1) = g;
            set(gui.hMap,'CData',gui.map);
            setRightClim(gui.axR, fig);
            guidata(fig,gui);
            drawnow limitrate
        elseif val==9 && lbl>0 && gui.cMask(y,x)==5
            gui.userSet(lbl:lbl+1) = false;      % release: hand off control back to the sliders
            guidata(fig,gui);
            recomputeGuess(fig);                 % segment falls back to the current slider guess
            drawnow limitrate
        else
            return;
        end
    end

    function updateImg(src,~)
        redrawLeft(ancestor(src,'figure'));           % redraw image + ROI outlines
    end

    function redrawLeft(fig)
        % Redraws the left-panel image (current popup selection) and its ROI
        % outlines, and re-syncs the same outlines onto the right panel. Called
        % when the displayed image is switched and when an ROI is de-labelled,
        % so clearing gui.rois(j) removes that outline from BOTH panels.
        gui = guidata(fig);
        pop = getappdata(fig,'popHandle');
        u   = pop.UserData;  sel = pop.Value;
        cm  = u.data{find(strcmp(u.names,'cMask'),1)};
        img = u.data{sel};
        imagesc(gui.axL,img), axis(gui.axL,'image','off')
        clim(gui.axL,prctile(img(cm>0),[1,99]));
        title(gui.axL,u.names{sel},'Interpreter','none')

        drawOutlines(gui.axL,fig);                    % left panel (imagesc cleared old ones)
        delete(findobj(gui.axR,'Tag','roiOutline'));  % right panel: drop old outlines but
        drawOutlines(gui.axR,fig);                    %   keep the label image (hMap) itself
    end

    function drawOutlines(ax, fig)
        % Draws the outline of every currently-labelled ROI onto AX and tags
        % each 'roiOutline' so it can be cleared again. Vessels (and any
        % parenchyma sharing a vessel's label) use a per-label colour; other
        % parenchyma is green. Shared by both panels.
        gui = guidata(fig);
        for j = 1:numel(gui.rois)
            if strlength(gui.rois(j))==0, continue, end   % only labelled ROIs

            hold(ax,'on')
            if any(gui.cMask(gui.sMap==j)==5)
                h = visboundaries(ax,ismember(gui.sMap,j),'Color',vesselColor(fig,gui.rois(j)),'LineWidth',2);
                set(h,'Tag','roiOutline');
            elseif  any(gui.cMask(gui.sMap==j)==1)
                h = visboundaries(ax,gui.sMap==j,'Color',parenchymaColor(fig,gui,gui.rois(j)),'LineWidth',2);
                set(h,'Tag','roiOutline');
            end
            hold(ax,'off')
        end
    end
    function c = vesselColor(fig, label)
        % Returns the outline colour for a vessel ROI label. All segments that
        % share a label get the same colour; each new label takes the next
        % palette entry (cycling). The palette deliberately avoids the parula
        % (blue/green/yellow) and blue-green-red colormaps so outlines stay
        % distinct; visboundaries' two-tone rendering keeps them visible on any
        % background.
        palette = [1.00 0.00 1.00;   % magenta
            1.00 0.50 0.00;   % orange
            1.00 0.20 0.60;   % rose / pink
            0.60 0.20 1.00;   % violet
            0.90 0.40 0.10;   % burnt orange
            1.00 0.55 0.85];  % light pink
        key  = strtrim(string(label));
        seen = getappdata(fig,'labelOrder');
        idx  = find(seen==key,1);
        if isempty(idx)
            seen(end+1) = key;                  % first time this label is used
            setappdata(fig,'labelOrder',seen);
            idx = numel(seen);
        end
        c = palette(mod(idx-1,size(palette,1))+1,:);
    end
    function c = parenchymaColor(fig, gui, label)
        % Parenchyma shares the vessel colour when the same label is also used
        % on a vessel segment (cMask==5); otherwise it stays green. A label used
        % only on parenchyma never consumes a vessel-palette slot, because
        % vesselColor is called only when a matching vessel exists.
        v   = unique(gui.sMap(gui.cMask==5));   % segment ids touching vessel pixels
        v   = v(v>=1 & v<=numel(gui.rois));
        key = strtrim(string(label));
        if any(strtrim(gui.rois(v))==key)
            c = vesselColor(fig,label);         % match → same colour as that vessel
        else
            c = [0 1 0];                        % green
        end
    end
    function setRightClim(ax, fig)
        % Symmetric colour limits for the signed label map, centred on green (0).
        % The half-range is (number of ENABLED guess sliders) + 1, so the colour
        % scale grows with the active sliders: the strongest auto-guess (= that
        % count) sits just inside the extreme, and a hand-set label (+/-maxGuess,
        % always beyond) saturates at the pure red / blue end.
        sld = getappdata(fig,'sliders');
        if isempty(sld)
            L = 1;
        else
            L = nnz([sld.guesses] & [sld.enabled]) + 1;
        end
        clim(ax,[-L L]);
    end
    function removeROI(src)
        % Clears every ROI whose stored label equals the current label field,
        % both visually (outline) and in the stored gui.rois variable. The
        % cleared labels propagate to results.sMetrics/dvsMetrics at save time.
        fig    = ancestor(src,'figure');
        target = strtrim(get(getappdata(fig,'txtHandle'),'String'));
        if isempty(target), return, end            % empty field ⇒ nothing to match
        gui = guidata(fig);
        hit = strcmp(strtrim(gui.rois), target);    % ROIs sharing this label
        if ~any(hit), return, end
        gui.rois(hit) = "";                         % remove stored label
        guidata(fig,gui);
        redrawLeft(fig);                            % remove their outlines
        drawnow
    end
    function pickCat(src,val)
        fig = ancestor(src,'figure');
        setappdata(fig,'paintVal',val)
        bh  = getappdata(fig,'btnHandles');
        set(bh,'Enable','on'), set(src,'Enable','off')
    end
    function repositionPanels(figH, tileAx, topPanel, topWH, botPanel, botWH, gap)
        % Stacks the slider panel (top) above the control panel (bottom) and
        % centres the pair within the control tile. Replaces the single-panel
        % 'reposition' so both panels track figure resizes together.
        figPix  = getpixelposition(figH);
        tilePos = get(tileAx,'Position');
        tilePix = [tilePos(1:2).*figPix(3:4)  tilePos(3:4).*figPix(3:4)];
        W    = max(topWH(1), botWH(1));
        Htot = topWH(2) + gap + botWH(2);
        x = tilePix(1) + max(0,(tilePix(3)-W)/2);
        y = tilePix(2) + max(0,(tilePix(4)-Htot)/2);
        set(botPanel,'Position',[x+(W-botWH(1))/2, y,              botWH]);
        set(topPanel,'Position',[x+(W-topWH(1))/2, y+botWH(2)+gap, topWH]);
    end
    function S = addRangeSlider(fig, parent, band, name, rng, lo, hi, enabled, lockSwitch)
        % One labelled, dual-thumb range slider with an enable switch.
        %   band       = [x y w h] region (parent-normalized) for the whole unit
        %   rng        = [min max] full range of the variable
        %   lo, hi     = initial thumb positions ([] => split the range in thirds;
        %                right thumb is kept at or above the left thumb)
        %   enabled    = logical, whether the slider starts active
        %   lockSwitch = logical, if true the On/Off switch cannot be toggled
        %                (used when the driving variable, e.g. psPI, is absent)
        % The two thumbs split the range into three segments. Handles + state are
        % returned in S (stored in the 'sliders' appdata).
        bx=band(1); by=band(2); bw=band(3); bh=band(4);
        bg = get(parent,'BackgroundColor');
        labelPos  = [bx+0.04*bw, by+0.56*bh, 0.58*bw, 0.36*bh];  % top-left
        switchPos = [bx+0.64*bw, by+0.56*bh, 0.32*bw, 0.36*bh];  % top-right, by the label
        axPos     = [bx+0.04*bw, by+0.12*bh, 0.60*bw, 0.40*bh];  % track, bottom-left
        resetPos  = [bx+0.66*bw, by+0.14*bh, 0.30*bw, 0.34*bh];  % bottom-right (old switch slot)

        span = rng(2)-rng(1);
        if nargin<6 || isempty(lo), lo = rng(1) + span/3;   end
        if nargin<7 || isempty(hi), hi = rng(1) + 2*span/3; end
        if nargin<8 || isempty(enabled),    enabled = true;     end
        if nargin<9 || isempty(lockSwitch), lockSwitch = false; end
        lo = min(max(lo,rng(1)),rng(2));
        hi = min(max(hi,rng(1)),rng(2));
        if hi < lo, hi = lo; end

        S.name    = name;
        S.rng     = rng;
        S.lo      = lo;
        S.hi      = hi;
        S.lo0     = lo;                       % default thumb positions (for Reset)
        S.hi0     = hi;
        S.enabled = logical(enabled);

        % variable-name label, sitting just above the slider
        S.hLabel = uicontrol(parent,'Style','text','String',name, ...
            'Units','normalized','Position',labelPos, ...
            'HorizontalAlignment','left','FontWeight','bold','BackgroundColor',bg);

        % axes hosting the track + thumbs
        S.ax = axes('Parent',parent,'Units','normalized','Position',axPos);
        hold(S.ax,'on');
        set(S.ax,'XLim',rng,'YLim',[0 1],'XTick',[],'YTick',[], ...
            'Color','none','XColor','none','YColor','none','Box','off');
        disableDefaultInteractivity(S.ax);    % no built-in pan/zoom on drag
        S.ax.Toolbar = [];                    % no hover toolbar

        % transparent patch so a click anywhere on the axes grabs the slider
        S.hHit = patch(S.ax,'XData',[rng(1) rng(2) rng(2) rng(1)],'YData',[0 0 1 1], ...
            'FaceColor',[1 1 1],'FaceAlpha',0,'EdgeColor','none','PickableParts','all');
        uistack(S.hHit,'bottom');

        % three track segments (left / middle / right) + the two thumbs
        yMid = 0.5;
        S.hTrackL = line(S.ax,[rng(1) S.lo],[yMid yMid],'LineWidth',6,'Color',[0.80 0.80 0.84]);
        S.hTrackM = line(S.ax,[S.lo S.hi],[yMid yMid],'LineWidth',6,'Color',[0.20 0.45 0.95]);
        S.hTrackR = line(S.ax,[S.hi rng(2)],[yMid yMid],'LineWidth',6,'Color',[0.80 0.80 0.84]);
        S.hLo = line(S.ax,S.lo,yMid,'Marker','o','MarkerSize',12,'LineStyle','none', ...
            'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.20 0.30 0.55],'LineWidth',1.5);
        S.hHi = line(S.ax,S.hi,yMid,'Marker','o','MarkerSize',12,'LineStyle','none', ...
            'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.20 0.30 0.55],'LineWidth',1.5);

        % enable/disable switch
        if S.enabled, swStr='On'; swVal=1; else, swStr='Off'; swVal=0; end
        S.hSwitch = uicontrol(parent,'Style','togglebutton','String',swStr,'Value',swVal, ...
            'Units','normalized','Position',switchPos, ...
            'Callback',@(src,~)toggleSlider(fig,src));
        if lockSwitch, set(S.hSwitch,'Enable','off'); end

        % reset-to-default button (where the switch used to be)
        S.hReset = uicontrol(parent,'Style','pushbutton','String','Reset', ...
            'Units','normalized','Position',resetPos, ...
            'Callback',@(src,~)resetSlider(fig,src));

        setSliderEnabledLook(S);              % grey the track if it starts disabled
    end
    function drawSlider(S)
        % Repaint a slider's segments and thumbs from its stored lo/hi values.
        set(S.hTrackL,'XData',[S.rng(1) S.lo]);
        set(S.hTrackM,'XData',[S.lo S.hi]);
        set(S.hTrackR,'XData',[S.hi S.rng(2)]);
        set(S.hLo,'XData',S.lo,'YData',0.5);
        set(S.hHi,'XData',S.hi,'YData',0.5);
    end
    function toggleSlider(fig, src)
        % Switch callback: flips the slider's enabled flag and restyles it.
        sld = getappdata(fig,'sliders');
        idx = find(arrayfun(@(z) isequal(z.hSwitch,src), sld),1);
        if isempty(idx), return, end
        sld(idx).enabled = logical(get(src,'Value'));
        if sld(idx).enabled, set(src,'String','On'); else, set(src,'String','Off'); end
        setappdata(fig,'sliders',sld);
        setSliderEnabledLook(sld(idx));
        if sld(idx).guesses, recomputeGuess(fig); end   % its contribution changed
    end
    function setSliderEnabledLook(S)
        % Full colour when enabled, greyed when disabled.
        if S.enabled
            cMid=[0.20 0.45 0.95]; cSide=[0.80 0.80 0.84];
            cEdge=[0.20 0.30 0.55]; cFace=[1 1 1];
        else
            cMid=[0.78 0.78 0.82]; cSide=[0.88 0.88 0.90];
            cEdge=[0.72 0.72 0.76]; cFace=[0.95 0.95 0.96];
        end
        set(S.hTrackM,'Color',cMid);
        set([S.hTrackL S.hTrackR],'Color',cSide);
        set([S.hLo S.hHi],'MarkerEdgeColor',cEdge,'MarkerFaceColor',cFace);
    end
    function nudgeSlider(fig, idx, delta, useHi)
        % Move one thumb of slider idx by delta, keeping the right thumb at or
        % above the left and both inside the range, then redraw and (for guess
        % sliders) refresh the guess. Shared by the key and scroll handlers.
        sld = getappdata(fig,'sliders');
        if isempty(idx) || idx<1 || idx>numel(sld), return, end
        S = sld(idx);
        if ~S.enabled, return, end
        if useHi
            S.hi = min(max(S.hi + delta, S.lo), S.rng(2));   % upper thumb
        else
            S.lo = min(max(S.lo + delta, S.rng(1)), S.hi);   % lower thumb
        end
        drawSlider(S);
        sld(idx) = S;
        setappdata(fig,'sliders',sld);
        if S.guesses, recomputeGuess(fig); end
    end
    function onSliderKey(fig, evt)
        % Left/right arrows nudge the LOWER thumb of the slider under the
        % cursor; hold Shift to nudge the UPPER thumb instead.
        idx = getappdata(fig,'hoverSlider');
        if isempty(idx), return, end
        switch evt.Key
            case 'leftarrow',  sgn = -1;
            case 'rightarrow', sgn =  1;
            otherwise, return
        end
        sld = getappdata(fig,'sliders');
        if idx<1 || idx>numel(sld), return, end
        step  = (sld(idx).rng(2)-sld(idx).rng(1))/100;
        useHi = any(strcmp(evt.Modifier,'shift'));
        nudgeSlider(fig, idx, sgn*step, useHi);
    end
    function onSliderScroll(fig, evt)
        % Scroll wheel nudges the LOWER thumb of the slider under the cursor
        % (wheel up = increase); hold Shift to nudge the UPPER thumb instead.
        idx = getappdata(fig,'hoverSlider');
        if isempty(idx), return, end
        sgn = -sign(evt.VerticalScrollCount);    % wheel up (count<0) => increase
        if sgn==0, return, end
        sld = getappdata(fig,'sliders');
        if idx<1 || idx>numel(sld), return, end
        step  = (sld(idx).rng(2)-sld(idx).rng(1))/100;
        useHi = any(strcmp(get(fig,'CurrentModifier'),'shift'));
        nudgeSlider(fig, idx, sgn*step, useHi);
    end
    function resetSlider(fig, src)
        % Reset button: restore a slider's thumbs to their default positions.
        sld = getappdata(fig,'sliders');
        idx = find(arrayfun(@(z) isequal(z.hReset,src), sld),1);
        if isempty(idx), return, end
        sld(idx).lo = sld(idx).lo0;
        sld(idx).hi = sld(idx).hi0;
        setappdata(fig,'sliders',sld);
        drawSlider(sld(idx));
        if sld(idx).guesses, recomputeGuess(fig); end
    end
    function recomputeGuess(fig)
        % Re-derives the automatic part of the guess from the enabled,
        % guess-driving range sliders, then refreshes the right-hand label image.
        %   value above a slider's HIGH thumb -> +1 (artery)
        %   value below a slider's LOW  thumb -> -1 (vein)
        %   value between the two thumbs       ->  0 (unknown)
        % Contributions from every enabled guess slider are summed; a disabled
        % slider contributes nothing. Segments the user labelled by hand
        % (gui.userSet) keep their value and are never overwritten.
        gui = guidata(fig);
        if ~isstruct(gui) || ~isfield(gui,'guess'), return, end
        sld = getappdata(fig,'sliders');

        votes = zeros(numel(gui.guess),1);       % per-row +1 / -1 / 0 tally
        for si = 1:numel(sld)
            if sld(si).guesses && sld(si).enabled && ~isempty(sld(si).segVals)
                vv = sld(si).segVals(:);
                votes = votes + double(vv > sld(si).hi) - double(vv < sld(si).lo);
            end
        end

        % spread each segment's vote over its two endpoints idxs(i):idxs(i)+1
        selRows = gui.m5;                          % valid category-5 rows
        autoG = accumarray(gui.idxs(selRows),   votes(selRows), [numel(gui.guess) 1]) + ...
            accumarray(gui.idxs(selRows)+1, votes(selRows), [numel(gui.guess) 1]);

        gNew = gui.guess;
        gNew(~gui.userSet) = autoG(~gui.userSet); % auto only where not hand-labelled
        gui.guess = gNew;

        lutG = [0; gNew];                          % rebuild the displayed label map
        mpG  = lutG(gui.sMap+1);
        mpG(gui.cMask<4) = NaN;
        gui.map = mpG;
        set(gui.hMap,'CData',mpG);
        setRightClim(gui.axR, fig);
        guidata(fig,gui);
    end
    function [rng, lo, hi] = autoRangeThumbs(vals, vmask)
        % Pick a range and default thumbs for a guess slider from the spread of
        % its per-segment values (over the valid rows in vmask): range = observed
        % min/max (padded), thumbs = medians of the lower / upper halves (roughly
        % the quartiles). Falls back to [0 1] split in thirds without enough data.
        rr = vals(vmask);
        rr = rr(isfinite(rr));
        if isempty(rr)
            rng = [0 1]; lo = 1/3; hi = 2/3; return
        end
        rmin = min(rr);  rmax = max(rr);
        if rmin==rmax, rmin = rmin-0.5; rmax = rmax+0.5; end
        rmed = median(rr);
        lo = median(rr(rr<rmed));            % ~ lower quartile
        hi = median(rr(rr>rmed));            % ~ upper quartile
        if ~isfinite(lo), lo = rmin + (rmax-rmin)/3;   end
        if ~isfinite(hi), hi = rmin + 2*(rmax-rmin)/3; end
        pd  = 0.02*(rmax-rmin);
        rng = [rmin-pd, rmax+pd];
        lo  = min(max(lo,rng(1)),rng(2));
        hi  = min(max(hi,rng(1)),rng(2));
        if hi<lo, hi=lo; end
    end
    function [rng, lo, hi] = thumbsAtValue(vals, vmask, defVal)
        % Range spans the per-segment values over vmask (extended to include
        % defVal), and BOTH thumbs start at defVal -- e.g. the parenchymal median,
        % giving a single threshold there by default. Falls back to the range
        % midpoint / [0 1] when defVal or the data is unavailable.
        rr = vals(vmask);
        rr = rr(isfinite(rr));
        haveData = ~isempty(rr);
        haveDef  = isscalar(defVal) && isfinite(defVal);
        if ~haveData && ~haveDef
            rng = [0 1]; lo = 0.5; hi = 0.5; return
        end
        if haveData, rmin = min(rr); rmax = max(rr); else, rmin = defVal; rmax = defVal; end
        if haveDef,  rmin = min(rmin,defVal); rmax = max(rmax,defVal); end
        if rmin==rmax, rmin = rmin-0.5; rmax = rmax+0.5; end
        pd  = 0.02*(rmax-rmin);
        rng = [rmin-pd, rmax+pd];
        if haveDef, vc = defVal; else, vc = (rmin+rmax)/2; end
        vc  = min(max(vc,rng(1)),rng(2));
        lo  = vc;  hi = vc;
    end
    function S = metricSlider(col, invert)
        % Build a single-column "magnitude" guess slider for the table column
        % named col: per-segment value is that metric (negated when invert is
        % true) sampled at the segment's first endpoint idxs(i), with both thumbs
        % defaulting to its median over parenchymal (category-1) segments. Returns
        % an inert slider with S.has == false when the column is absent, so the
        % caller can skip building/showing it. Uses the shared setup variables
        % vn / results / idxs / ii / ctgs / m5 / N.
        S = struct('seg',nan(N,1),'rng',[0 1],'lo',1/3,'hi',2/3,'has',false);
        if ~any(strcmp(col,vn)), return, end
        xv = results.sMetrics{:,col};
        if invert, xv = -xv; end
        sv = nan(N,1);  sv(ii) = xv(idxs(ii));
        sv(~isfinite(sv)) = NaN;
        mPar = xv(ctgs==1);  medv = median(mPar(isfinite(mPar)));
        S.seg = sv;  S.has = true;
        [S.rng, S.lo, S.hi] = thumbsAtValue(sv, m5, medv);
    end

end

% =====================================================================
% Local functions, OUTSIDE the main function so they do not share its workspace -
% everything above this line is a nested callback that does.
% =====================================================================

function cmap=vesselTypeColormap()
%vesselTypeColormap  The one artery/vein colour convention: blue for vein, green
%   for undecided, red for artery.  Used by the labelling GUI and by the report
%   page, so the page a reviewer looks at is the map the operator worked on.
cmap=zeros(11,3);
cmap(1:5,3)=(1:5)./10+0.5;
cmap(6,2)=1;
cmap(7:11,1)=(5:-1:1)./10+0.5;
end

% =====================================================================
function writeVesselTypesReport(rep,fName,results)
%writeVesselTypesReport  Save <name>_rep_vesseltypes.jpg: the artery/vein map over
%   the BFI image, with the NAMED vessels outlined and annotated.
%
%   Only the non-empty labels count as "labelled vessels" - the ones somebody typed
%   a name for.  Everything else is left to the colour map, which already says
%   artery or vein.  The page is a record, so it is written even when nothing was
%   named: an unlabelled map is still the map that was accepted.
%
%   A FAILURE HERE COSTS A LINE, NOT THE RUN.  This page is new: before it existed
%   setVesselTypes saved its labelling and moved on, and it must still do that.
if ~isfield(results,'mapType') || ~isfield(results,'sMap') || ~isfield(results,'imgBFI')
    return
end

cMask=results.cMask;
img=sqrt(results.imgBFI);
img=mat2gray(img,double(prctile(img(cMask(:)>0),[5,99])));
fSize=floor((min(size(img))./20))*2+1;
img=enhanceForDisplay(img,fSize).*(cMask>0);

map=results.mapType;
L=max(abs(map(:)),[],'omitnan');
if ~isfinite(L) || L==0, L=1; end

% the named vessels: non-empty labels, mapped back onto sMap through their idx
labelled=false(0,1); idxs=[];
if istable(results.sMetrics) && ismember('label',results.sMetrics.Properties.VariableNames)
    idxs=results.sMetrics.idx;
    labelled=strlength(strtrim(string(results.sMetrics.label)))>0 & ~isnan(idxs);
end

f=reportFigure(rep,'vesseltypes','single');
try
    t=tiledlayout(f,1,1,'TileSpacing','compact','Padding','compact');
    ax=nexttile(t);
    imagesc(ax,img); axis(ax,'image','off'); colormap(ax,'gray')
    setClim(ax,prctile(img(cMask(:)>0),[1,99]));
    hold(ax,'on')
    % The type map on top, transparent outside the vessels, on ITS own colour scale.
    % IT GOES IN THE TILE, not on the figure at a copy of the tile's position: a
    % second axes pinned to the same tile follows it when the layout re-flows, and
    % the layout DOES re-flow - reportSave puts the recording's name across the top
    % of every page.  Copied coordinates would leave the map behind, out of register
    % with the image under it and with its caption on top of the page title.
    axMap=axes(t,'Color','none');
    axMap.Layout.Tile=ax.Layout.Tile;
    imagesc(axMap,map,'AlphaData',~isnan(map));
    axis(axMap,'image','off'); colormap(axMap,vesselTypeColormap()); clim(axMap,[-L L]);
    hold(axMap,'on')
    nNamed=0;
    for i=1:1:numel(labelled)
        if ~labelled(i), continue; end
        bw=(results.sMap==idxs(i));
        if ~any(bw(:)), continue; end
        nNamed=nNamed+1;
        visboundaries(axMap,bw,'Color','w','LineWidth',1);
        [y,x]=find(bw);
        text(axMap,mean(x),mean(y),char(strtrim(string(results.sMetrics.label(i)))), ...
            'Color','w','FontWeight','bold','HorizontalAlignment','center', ...
            'VerticalAlignment','middle','Interpreter','none');
    end
    hold(axMap,'off'); hold(ax,'off')
    title(axMap,sprintf('Vessel types (blue vein / red artery), %d named',nNamed))
catch
    delete(f);              % no page, no line: a report is a by-product
    return
end
reportSave(rep,f,'vesseltypes',fName);
end

% =====================================================================
function setClim(ax,lims)
%setClim  clim() with the degenerate case handled: a flat image gives equal
%   percentiles, and clim rejects a non-increasing pair.
lims=double(lims);
if any(~isfinite(lims)) || lims(2)<=lims(1), return; end
clim(ax,lims);
end