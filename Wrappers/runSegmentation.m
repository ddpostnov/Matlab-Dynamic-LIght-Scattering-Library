%runSegmentation  Categorize + label vessels/parenchyma from LSCI mean images
%
%   runSegmentation(s,fNames) is the fully automatic static-segmentation step: for
%   every *_K_r.mat / *_I_r.mat file it builds the five-level categorical mask, turns
%   it into an indexed vessel/parenchyma label map, extracts per-structure scalar
%   metrics and mean/median time-series, saves the updated *_r.mat / *_s.mat, and
%   writes two report pages (_rep_categories, _rep_segments).  It takes NO
%   user input: the interactive region selection now lives in setRegions (which
%   writes results.regionsMask), and the automatic dynamic-segmentation loop lives
%   in runDynamicSegmentation.  runSegmentation replaces the retired runCategories
%   (its categorization is now the getPixelCategories core) and the static half of
%   the old runSegmentation.
%
%   PER-FILE PIPELINE
%     read results.regionsMask (from setRegions; whole window if absent)
%       -> getPixelCategories  (edge size + enhancement + trust mask + 5-level cMask)
%       -> _rep_categories page
%       -> getSegmentationLabels (centerlines -> sMap/pMap; merges inner walls)
%       -> one grouping pass over sMap -> sMetrics (geometry) + sData (traces)
%       -> _rep_segments page -> save
%
%   WHICH WAY ROUND THE VESSELS ARE IS RESOLVED IN ONE PLACE, getVesselPolarity, and
%   never by a file-name test here: a contrast product's vessels are dark as a matter
%   of physics and an intensity product's are bright by this library's stated
%   assumption.  It matters that the two are asked as one question, because inverting
%   a picture that should not have been inverted segments the tissue BETWEEN the
%   vessels and nothing about the run looks wrong afterwards.
%
%   INPUTS
%     s        parameter structure.
%              Categorization (read by getPixelCategories): trustLimitsK, iniSizeN,
%                diffusionSchedule, lSizeN, sSizeN, sens, deSens, sSizeScale,
%                lThinN, imOpen, iEdge, eEdge.
%              Labelling (read by getSegmentationLabels): correctNodes, simR, difR,
%                sMinL, prchNSize.
%              Traces: sStat ('mean' or, default, 'median').
%              s.fNamesCopyTo (optional, default {}): assign the computed segmentation
%                to co-registered sibling files - see below.
%     fNames   FLAT (order-independent) cell array of *_K_r.mat / *_I_r.mat paths;
%              iterate element by element (grouping was setRegions' job).  Each file
%              must have matching *_d.mat / *_s.mat siblings.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   s.fNamesCopyTo - assign the segmentation to siblings (replaces assignCategories)
%     Mirrors fNames with one extra dimension: s.fNamesCopyTo(i,:) lists the sibling
%     files that inherit the segmentation computed for fNames{i} (0..K targets per
%     source; use a COLUMN fNames so row i lines up with source i, e.g.
%     s.fNamesCopyTo = regexprep(fNames,'_t_K_r.mat$','_c_K_r.mat')).  For each target
%     the SHARED SPATIAL products (cMask, regionsMask, mask, sMap, pMap, sMetrics) are
%     copied verbatim and the target's OWN sData is RE-EXTRACTED from its OWN cube using
%     the copied sMap (siblings are different temporal data on the same FOV, so masks/
%     labels are shared but traces are per-file).  Assumes source and target are
%     co-registered / same FOV (the old assignCategories assumption).
%
%   OUTPUT SIDE-EFFECTS (per file, and per copy target)
%     <name>_r.mat   results.{cMask,regionsMask,mask,sMap,pMap,sMetrics,sData}
%     <name>_s.mat   settings.runSegmentation = s (now carrying edgeSize and sStat)
%     <name>_rep_categories.jpg  categories report page
%     <name>_rep_segments.jpg    segments report page
%
%   EXAMPLE
%     fNames = getFileNamesList(root,'*_t_K_r.mat');           % flat column list
%     setRegions(s,fNames);                                    % define regions first
%     s.fNamesCopyTo = regexprep(fNames,'_t_K_r.mat$','_c_K_r.mat');
%     runSegmentation(s,fNames);                               % segment _t, copy onto _c
%     runDynamicSegmentation(s,fNames);                        % optional heavy loop
%
%   DEPENDS ON
%     getVesselPolarity, getPixelCategories, getSegmentationLabels,
%     enhanceForDisplay, showSegmentsPreview (Core); plus MATLAB's Image
%     Processing Toolbox.
%
% See also: setRegions, runDynamicSegmentation, getVesselPolarity,
%           getPixelCategories, getSegmentationLabels, splitRegions
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026


%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CATEGORIZATION
% s.trustLimitsK=[0.001,0.5]; % min (fastest flows) and max (slowest flows) expected contrast
% s.lSizeN=61; % Odd, ~2x the largest vessel
% s.sSizeN=15; % Odd, ~2x the small-vessel diameter
% s.sens=0.3;  % segmentation sensitivity (raise if missing vessels, lower to cut noise)
% s.deSens=1;  % integer >=1; small-vessel sensitivity
% s.sSizeScale=1; % small-object background/unrecognised assignment scaler
% s.lThinN=2;  % large-vessel thinning
% s.imOpen=2;  % small-vessel thinning
% s.iniSizeN=7;% Odd number equal or larger than the spatial contrast kernel
% s.iEdge=2; s.eEdge=2; % internal / external wall edge widths
% s.diffusionSchedule='multiscale'; % 'multiscale' (nine extra coarse-to-fine passes,
%              % written for speckle contrast's noise) or 'single'.  Default when the
%              % field is absent is 'multiscale', which is what every launcher here
%              % was already doing.
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
% s.correctNodes=true; s.simR=0.3; s.difR=0.4; % branch-node correction
% s.sMinL=15;      % minimum segment length
% s.prchNSize=30;  % parenchymal cell neighbourhood
% s.sStat='median';% per-segment trace statistic ('median' default, or 'mean')
% s.fNamesCopyTo={}; % optional: assign the segmentation to co-registered siblings

function runSegmentation(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_K_r.mat')|| contains(x,'_I_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_K_r.mat" or "_I_r.mat".  ' ...
        'Every step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_t_K_r.mat'').']);
end
if ~isfield(s,'fNamesCopyTo'), s.fNamesCopyTo={}; end
% Resolved here, not only in the core, so the choice is RECORDED in the saved
% settings like every other tunable.  getSegmentationLabels reads it as a worker
% bound on its per-label parfor; true (the default) is today's behaviour.
if ~isfield(s,'parforSegmentationLabels') || isempty(s.parforSegmentationLabels)
    s.parforSegmentationLabels=true;
end
% The same, for the diffusion schedule that used to be chosen by the file's product
% token.  'multiscale' is what every launcher in this repository was already getting,
% so an s written before the field existed keeps its numbers exactly; the workbench
% supplies it from the step's own preset, which is where the two branches differ.
if ~isfield(s,'diffusionSchedule') || isempty(s.diffusionSchedule)
    s.diffusionSchedule='multiscale';
end
if ~any(strcmp(s.diffusionSchedule,{'multiscale','single'}))
    error('runSegmentation:diffusionSchedule', ...
        ['The vessel-smoothing schedule is either "multiscale" or "single"; ' ...
         '"%s" is neither.'], char(s.diffusionSchedule));
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Segmentation',fNames);

for fidx=1:1:numel(fNames)
     if reportCancelled(rep), break; end        % cooperative cancel between files
     if ~isempty(fNames{fidx})
    s.fName=fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results source settings
    load(getProductPath(s.fName,'d'),'source')
    load(getProductPath(s.fName,'s'),'settings');
    load(s.fName,'results');

    % --- which way round the vessels are, and which mean image holds them ---
    [polarity,product]=getVesselPolarity(s.fName);
    if strcmp(product,'K')
        if isfield(results,'imgK')
            imgIni=results.imgK;
        else
            imgIni=mean(source.data,3,'omitmissing');
        end
    else
        imgIni=results.imgI;
    end

    % --- region mask from setRegions; a MISSING regionsMask (setRegions skipped, or no
    %     ROI drawn) means the whole window: an all-ones mask the size of the image ---
    if isfield(results,'regionsMask')
        regionsMask=results.regionsMask;
    else
        regionsMask=ones(size(imgIni));
    end

    % --- categorize (automatic core: edge size + enhancement + trust mask + cMask) ---
    % An intensity product carries no results.mask of its own - neither entry step
    % writes one, deliberately (see 03-regions-segmentation.md) - so existingMask is
    % empty there and the trust mask is derived open, restricted only by the regions.
    if isfield(results,'mask'), existingMask=results.mask; else, existingMask=[]; end
    [cMask,s.edgeSize,maskOut,imgVis]=getPixelCategories(imgIni,regionsMask,existingMask,product,polarity,s);
    results.regionsMask=regionsMask;
    results.mask=(maskOut==(regionsMask>0));
    results.cMask=cMask;                          % store the 5-level (un-merged) mask
    writeCategoriesPreview(rep,s.fName,imgVis,cMask);

    % --- indexed label maps (mask algebra; merges inner walls into outer) ---
    edgeSize=s.edgeSize;
    [sMap,pMap,sLines,cMask,~,dMask,~]=getSegmentationLabels(cMask,edgeSize,s); % cMask now merged
    results.pMap=pMap;

    %Build segments table
    % ONE pass over the labelled pixels gives every label its pixel list and its area;
    % the metrics and the traces then index that grouping instead of re-scanning the
    % whole image once per label (5215 labels x 1.69e6 pixels on a 1300 x 1300 field).
    nSeg=double(max(sMap(:)));
    grp=labelPixelGroups(sMap,nSeg);
    sMetrics=segmentMetrics(sLines,cMask,dMask,grp,nSeg);

    dataSize=size(source.data);
    source.data=reshape(source.data,[],dataSize(3));
    sData=groupedTraces(source.data,grp,nSeg,dataSize(3),s.sStat);
    source.data=reshape(source.data,dataSize);

    results.sMap=sMap;
    results.sMetrics=sMetrics;
    results.sData=sData;

    % --- segments preview: the wrapper owns the canvas, the core only draws ---
    fh=reportFigure(rep,'segments');
    showSegmentsPreview(fh,source.data,cMask,sMap,polarity);
    reportSave(rep,fh,'segments',s.fName);

    %Save the data
    settings.runSegmentation=reportSettings(s);
    reportWriting(rep);
    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
    save(s.fName,'results','-v7.3','-nocompression');
    reportSaved(rep);

    % --- assign the segmentation to co-registered siblings (s.fNamesCopyTo) ---
    shared=struct('cMask',results.cMask,'regionsMask',regionsMask,'mask',results.mask, ...
                  'sMap',sMap,'pMap',pMap,'sMetrics',sMetrics);
    tgts=copyTargets(s,fidx);
    for t=1:1:numel(tgts)
        if ~isempty(tgts{t})
            copySegmentationOnto(s,tgts{t},shared,rep,fidx);
        end
    end
     end
end
reportClose(rep);
end

% =====================================================================
function tgts=copyTargets(s,fidx)
%copyTargets  Row fidx of s.fNamesCopyTo (the copy targets for source fidx), or {}.
tgts={};
if isfield(s,'fNamesCopyTo') && ~isempty(s.fNamesCopyTo) && fidx<=size(s.fNamesCopyTo,1)
    tgts=s.fNamesCopyTo(fidx,:);
end
end

% =====================================================================
function copySegmentationOnto(s,targetName,shared,rep,fidx)
%copySegmentationOnto  Copy the shared spatial products onto a co-registered sibling
%   and RE-EXTRACT the target's own sData from its own cube (replaces assignCategories).
%   The sibling is a RECORDING OF ITS OWN, so it gets its own three lines - it is
%   not a sub-stage of the file that was just segmented.
reportFile(rep,fidx,targetName);
sT=s; sT.fName=targetName;
clearvars results source settings
load(getProductPath(targetName,'d'),'source')
load(getProductPath(targetName,'s'),'settings');
load(targetName,'results');

% The TARGET's own polarity and product token - a copy can cross from an intensity
% file to a contrast one, and each has to read its own mean image the right way up.
[polarity,product]=getVesselPolarity(targetName);
if strcmp(product,'K')
    if isfield(results,'imgK'), imgIni=results.imgK; else, imgIni=mean(source.data,3,'omitmissing'); end
else
    imgIni=results.imgI;
end

% shared spatial products (geometry is identical on a co-registered FOV)
results.cMask=shared.cMask;
results.regionsMask=shared.regionsMask;
results.mask=shared.mask;
results.sMap=shared.sMap;
results.pMap=shared.pMap;
results.sMetrics=shared.sMetrics;

% the target's OWN traces from the target's OWN cube, using the copied sMap
results.sData=extractTraces(source.data,shared.sMap,sT.sStat);

% previews from the target's own image with the copied overlays
imgVis=categoryPreviewBackground(imgIni,polarity);
writeCategoriesPreview(rep,targetName,imgVis,shared.cMask);
fh=reportFigure(rep,'segments');
showSegmentsPreview(fh,source.data,shared.cMask,shared.sMap,polarity);
reportSave(rep,fh,'segments',targetName);

settings.runSegmentation=reportSettings(sT);  % carry edgeSize / sStat onto the sibling
reportWriting(rep);
save(getProductPath(targetName,'s'),'settings','-v7.3','-nocompression');
save(targetName,'results','-v7.3','-nocompression');
reportSaved(rep);
end

% =====================================================================
function sData=extractTraces(dataCube,sMap,sStat)
%extractTraces  Per-label mean/median trace over a cube - the same pass runSegmentation
%   runs on its own file, for a copy target whose sMap arrives from the source.
dataSize=size(dataCube);
n=double(max(sMap(:)));
grp=labelPixelGroups(sMap,n);
dataCube=reshape(dataCube,[],dataSize(3));
sData=groupedTraces(dataCube,grp,n,dataSize(3),sStat);
end

% =====================================================================
function g=labelPixelGroups(map,n)
%labelPixelGroups  ONE pass over a label map -> each label's pixel list and area.
%   g.idx holds every labelled linear index sorted by label, g.first(i):g.last(i) is
%   label i's slice of it, and g.area(i) is its pixel count.  This replaces the
%   sum(map(:)==i) / map(:)==i that used to run once per label: on a 1300 x 1300 field
%   with 5215 labels those were 8.8e9 element compares for information one pass holds.
%
%   sort is STABLE, so within a label the indices stay ASCENDING - the exact order, and
%   therefore the exact floating-point summation order, a logical map(:)==i would give.
lin=find(map>0);
lin=lin(:);
[lab,ord]=sort(double(map(lin)));
g.idx  = lin(ord);
g.area = accumarray(lab,1,[n 1]);
g.last = cumsum(g.area);
g.first= g.last-g.area+1;
end

% =====================================================================
function sData=groupedTraces(dataFlat,g,n,nT,sStat)
%groupedTraces  Per-label mean/median trace, one gather per label from the grouping.
%   dataFlat is the cube reshaped to [nPixels x nT]; empty labels get a NaN column.
sData=zeros(nT,n);
useMean=strcmp(sStat,'mean');
for i=1:1:n
    if g.area(i)>0
        rows=g.idx(g.first(i):g.last(i));
        if useMean
            sData(:,i)=mean(dataFlat(rows,:),1,'omitnan');
        else %DEFAULT
            sData(:,i)=median(dataFlat(rows,:),1,'omitnan');
        end
    else
        sData(:,i)=NaN;
    end
end
end

% =====================================================================
function sMetrics=segmentMetrics(sLines,cMask,dMask,g,n)
%segmentMetrics  Per-segment scalars: length, chord-length ratio, diameter, area.
%   Every number is what the per-label whole-image loop produced; only the cost moved.
%     - area and category come from the labelPixelGroups pass, not a full-image compare;
%     - bwdist(~dMask) is computed ONCE and indexed, not recomputed per segment;
%     - bwmorph / bwdistgeodesic run inside the segment's own bounding box, padded by
%       one pixel.  Geodesic distance is a property of the mask and every pixel of the
%       segment is inside its own bounding box, so the distances are unchanged; the pad
%       makes bwmorph see the neighbourhood it would see on the full frame.  A crop is
%       a RECTANGLE, so column-major order inside it matches the full frame's order
%       restricted to it - the endpoint enumeration and the diameter's mean/std run
%       over the same values in the same sequence.
varTypes = ["double","double","double","double","double","double","double","double"];
varNames = ["idx","category","length","CLR","diameter","std(diameter)","area","nearestVesIdx"];
sMetrics=table('Size',[n,8],'VariableTypes',varTypes,'VariableNames',varNames);

dtWall=bwdist(~dMask);                        % ONCE, not once per segment
[H,W]=size(cMask);

% --- every centre-line's bounding box, in one pass over the skeleton pixels ---
bbY1=zeros(n,1); bbY2=zeros(n,1); bbX1=zeros(n,1); bbX2=zeros(n,1);
idxL=find(sLines>0);
if ~isempty(idxL)
    labL=double(sLines(idxL));
    [ry,rx]=ind2sub([H W],idxL(:));
    bbY1=accumarray(labL,ry,[n 1],@min,0);  bbY2=accumarray(labL,ry,[n 1],@max,0);
    bbX1=accumarray(labL,rx,[n 1],@min,0);  bbX2=accumarray(labL,rx,[n 1],@max,0);
end

for i=1:1:n
    area=g.area(i);
    if area>0
        cv=cMask(g.idx(g.first(i):g.last(i)));
        if ~all(cv==cv(1))
            error("Mix of categories detected in the region");
        end
        c=double(cv(1));

        if c==5
            % The segment's own bounding box, padded by one.  A label with no
            % centre-line pixels leaves the box at zero, which degenerates to a 1x1
            % empty crop and reaches the same "Centerline extraction failed" error
            % the whole-image loop reached.
            r1=max(1,bbY1(i)-1); r2=min(H,bbY2(i)+1);
            c1=max(1,bbX1(i)-1); c2=min(W,bbX2(i)+1);

            %Measure segment length
            tmp=(sLines(r1:r2,c1:c2)==i);
            ends=bwmorph(tmp,'endpoints');
            [y,x]  = find(ends);
            if numel(y)<2
                error(['Centerline extraction failed in segment ',num2str(i)]);
            elseif numel(y)>2
                maxDistance=0;
                for pointFirst = 1:numel(y)
                    seedMask = false(size(tmp));
                    seedMask(y(pointFirst), x(pointFirst)) = true;

                    % Calculate distance along the curve using quasi-euclidean metric
                    geodesicMap = bwdistgeodesic(tmp, seedMask, 'quasi-euclidean');

                    for pointSecond = (pointFirst + 1):numel(y)
                        pathLength = geodesicMap(y(pointSecond), x(pointSecond));

                        % Exclude infinite values occurring if points reside on disconnected segments
                        if pathLength > maxDistance && ~isinf(pathLength)
                            maxDistance = pathLength;
                            pointIndexA = pointFirst;
                            pointIndexB = pointSecond;
                        end
                    end
                end

                x = [x(pointIndexA), x(pointIndexB)];
                y = [y(pointIndexA), y(pointIndexB)];
            end
            d        = bwdistgeodesic(tmp, x(1), y(1), 'quasi-euclidean');
            l= max(d(~isnan(d)));

            %Measure CLR
            clr=l./hypot(x(2)-x(1), y(2)-y(1));

            %Get segment diameter
            d=dtWall(r1:r2,c1:c2).*tmp;
            d=[mean(d(d(:)>0)),std(d(d(:)>0))]*2;

            sMetrics(i,:)={i,c,l,clr,d(1),d(2),area,i};
        elseif c==3
            sMetrics(i,:)={i,c,NaN,NaN,NaN,NaN,area,i-1};
        else
            sMetrics(i,:)={i,c,NaN,NaN,NaN,NaN,area,NaN};

        end
    else
        sMetrics(i,:)={NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
    end
end
end

% =====================================================================
function imgVis=categoryPreviewBackground(imgIni,polarity)
%categoryPreviewBackground  Enhanced background for the _rep_categories page.
%   Mirrors getPixelCategories' internal imgVis prep so a copied file's page
%   matches a directly-segmented one (keep in sync with getPixelCategories).
img=imgIni;
if strcmp(polarity,'dark')
    img(img(:)>prctile(img(:),99))=prctile(img(:),99);
    img(img(:)<prctile(img(:),1))=prctile(img(:),1);
    img=imcomplement(img);
end
fSize=floor((min(size(img))./20))*2+1;
img(isnan(img))=0;
imgVis=enhanceForDisplay(img,fSize,min(15,fSize));
end

% =====================================================================
function writeCategoriesPreview(rep,fName,imgVis,cMask)
%writeCategoriesPreview  Save <name>_rep_categories.jpg: enhanced image + category
%   boundaries.  fName is passed explicitly so a copy target gets its own page.
f=reportFigure(rep,'categories');
try
    t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(t);
    imagesc(ax,rot90(imgVis, size(imgVis,2) > size(imgVis,1)))
    hold(ax,'on')
    visboundaries(ax,rot90(cMask>4, size(cMask,2) > size(cMask,1)),'Color','m')
    visboundaries(ax,rot90(cMask>2, size(cMask,2) > size(cMask,1)),'Color','w')
    hold(ax,'off')
    clim(ax,prctile(imgVis(:),[10,99]))
    axis(ax,'image')
    ax=nexttile(t);
    imagesc(ax,rot90(cMask, size(cMask,2) > size(cMask,1)))
    axis(ax,'image')
catch ME
    delete(f); rethrow(ME);      % reportSave never runs, so the figure goes here
end
reportSave(rep,f,'categories',fName);
end
