%runSegmentation  Categorize + label vessels/parenchyma from LSCI mean images
%
%   runSegmentation(s,fNames) is the fully automatic static-segmentation step: for
%   every *_K_d.mat / *_I_d.mat file it builds the five-level categorical mask, turns
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
%       -> single cube pass -> sMetrics (geometry) + sData (traces)
%       -> _rep_segments page -> save
%
%   INPUTS
%     s        parameter structure.
%              Categorization (read by getPixelCategories): trustLimitsK, iniSizeN,
%                lSizeN, sSizeN, sens, deSens, sSizeScale, lThinN, imOpen, iEdge, eEdge.
%              Labelling (read by getSegmentationLabels): correctNodes, simR, difR,
%                sMinL, prchNSize.
%              Traces: sStat ('mean' or, default, 'median').
%              s.fNamesCopyTo (optional, default {}): assign the computed segmentation
%                to co-registered sibling files - see below.
%     fNames   FLAT (order-independent) cell array of *_K_d.mat / *_I_d.mat paths;
%              iterate element by element (grouping was setRegions' job).  Each file
%              must have matching *_s.mat / *_r.mat siblings.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   s.fNamesCopyTo - assign the segmentation to siblings (replaces assignCategories)
%     Mirrors fNames with one extra dimension: s.fNamesCopyTo(i,:) lists the sibling
%     files that inherit the segmentation computed for fNames{i} (0..K targets per
%     source; use a COLUMN fNames so row i lines up with source i, e.g.
%     s.fNamesCopyTo = regexprep(fNames,'_t_K_d.mat$','_c_K_d.mat')).  For each target
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
%     fNames = getFileNamesList(root,'*_t_K_d.mat');           % flat column list
%     setRegions(s,fNames);                                    % define regions first
%     s.fNamesCopyTo = regexprep(fNames,'_t_K_d.mat$','_c_K_d.mat');
%     runSegmentation(s,fNames);                               % segment _t, copy onto _c
%     runDynamicSegmentation(s,fNames);                        % optional heavy loop
%
%   DEPENDS ON
%     getPixelCategories, getSegmentationLabels, enhanceForDisplay,
%     showSegmentsPreview (Core/LSCI); plus MATLAB's Image Processing Toolbox.
%
% See also: setRegions, runDynamicSegmentation, getPixelCategories,
%           getSegmentationLabels, splitRegions
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026


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
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING & TRACES
% s.correctNodes=true; s.simR=0.3; s.difR=0.4; % branch-node correction
% s.sMinL=15;      % minimum segment length
% s.prchNSize=30;  % parenchymal cell neighbourhood
% s.sStat='median';% per-segment trace statistic ('median' default, or 'mean')
% s.fNamesCopyTo={}; % optional: assign the segmentation to co-registered siblings

function runSegmentation(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_K_d.mat')|| contains(x,'_I_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat" or "_I_d.mat".');
end
if ~isfield(s,'fNamesCopyTo'), s.fNamesCopyTo={}; end
% Resolved here, not only in the core, so the choice is RECORDED in the saved
% settings like every other tunable.  getSegmentationLabels reads it as a worker
% bound on its per-label parfor; true (the default) is today's behaviour.
if ~isfield(s,'parforSegmentationLabels') || isempty(s.parforSegmentationLabels)
    s.parforSegmentationLabels=true;
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
    load(s.fName,'source')
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');

    % --- mean image + modality ---
    isK=contains(s.fName,'_K_d.mat');
    if isK
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
    if isfield(results,'mask'), existingMask=results.mask; else, existingMask=[]; end
    [cMask,s.edgeSize,maskOut,imgVis]=getPixelCategories(imgIni,regionsMask,existingMask,isK,s);
    results.regionsMask=regionsMask;
    results.mask=(maskOut==(regionsMask>0));
    results.cMask=cMask;                          % store the 5-level (un-merged) mask
    writeCategoriesPreview(rep,s.fName,imgVis,cMask);

    % --- indexed label maps (mask algebra; merges inner walls into outer) ---
    edgeSize=s.edgeSize;
    [sMap,pMap,sLines,cMask,~,dMask,~]=getSegmentationLabels(cMask,edgeSize,s); % cMask now merged
    results.pMap=pMap;

    %Build segments table
    varTypes = ["double","double","double","double","double","double","double","double"];
    varNames = ["idx","category","length","CLR","diameter","std(diameter)","area","nearestVesIdx"];
    sMetrics=table('Size',[max(sMap(:)),8],'VariableTypes',varTypes,'VariableNames',varNames);
    sData=zeros(size(source.data,3),max(sMap(:)));
    dataSize=size(source.data);
    source.data=reshape(source.data,[],dataSize(3));
    for i=1:1:max(sMap(:))
        area=sum(sMap(:)==i);
        if area>0
            if strcmp(s.sStat,'mean')
            sData(:,i)=mean(source.data(sMap(:)==i,:),1,'omitnan');
            else %DEFAULT
            sData(:,i)=median(source.data(sMap(:)==i,:),1,'omitnan');
            end
            c=unique(cMask(sMap==i));
            if numel(c)~=1
                error("Mix of categories detected in the region");
            end

            if c==5
                %Measure segment length
                tmp=(sLines==i);
                ends=bwmorph(sLines==i,'endpoints');
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
                d        = bwdistgeodesic(sLines==i, x(1), y(1), 'quasi-euclidean');
                l= max(d(~isnan(d)));

                %Measure CLR
                clr=l./hypot(x(2)-x(1), y(2)-y(1));

                %Get segment diameter
                d=bwdist(~dMask).*(sLines==i);
                d=[mean(d(d(:)>0)),std(d(d(:)>0))]*2;

                sMetrics(i,:)={i,c,l,clr,d(1),d(2),area,i};
            elseif c==3
                sMetrics(i,:)={i,c,NaN,NaN,NaN,NaN,area,i-1};
            else
                sMetrics(i,:)={i,c,NaN,NaN,NaN,NaN,area,NaN};

            end
        else
            sData(:,i)=NaN;
            sMetrics(i,:)={NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
        end
    end
    source.data=reshape(source.data,dataSize);

    results.sMap=sMap;
    results.sMetrics=sMetrics;
    results.sData=sData;

    % --- segments preview: the wrapper owns the canvas, the core only draws ---
    fh=reportFigure(rep,'segments');
    showSegmentsPreview(fh,source.data,cMask,sMap,isK);
    reportSave(rep,fh,'segments',s.fName);

    %Save the data
    settings.runSegmentation=reportSettings(s);
    reportWriting(rep);
    save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3','-nocompression');
    save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3','-nocompression');
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
load(targetName,'source')
load(strrep(targetName,'_d.mat','_s.mat'),'settings');
load(strrep(targetName,'_d.mat','_r.mat'),'results');

isK=contains(targetName,'_K_d.mat');
if isK
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
imgVis=categoryPreviewBackground(imgIni,isK);
writeCategoriesPreview(rep,targetName,imgVis,shared.cMask);
fh=reportFigure(rep,'segments');
showSegmentsPreview(fh,source.data,shared.cMask,shared.sMap,isK);
reportSave(rep,fh,'segments',targetName);

settings.runSegmentation=reportSettings(sT);  % carry edgeSize / sStat onto the sibling
reportWriting(rep);
save(strrep(targetName,'_d.mat','_s.mat'),'settings','-v7.3','-nocompression');
save(strrep(targetName,'_d.mat','_r.mat'),'results','-v7.3','-nocompression');
reportSaved(rep);
end

% =====================================================================
function sData=extractTraces(dataCube,sMap,sStat)
%extractTraces  Per-label mean/median trace over a cube - matches runSegmentation's
%   own sData pass exactly (mean/median over sMap(:)==i; NaN column for empty labels).
dataSize=size(dataCube);
n=double(max(sMap(:)));
sData=zeros(dataSize(3),n);
dataCube=reshape(dataCube,[],dataSize(3));
for i=1:1:n
    if sum(sMap(:)==i)>0
        if strcmp(sStat,'mean')
            sData(:,i)=mean(dataCube(sMap(:)==i,:),1,'omitnan');
        else %DEFAULT
            sData(:,i)=median(dataCube(sMap(:)==i,:),1,'omitnan');
        end
    else
        sData(:,i)=NaN;
    end
end
end

% =====================================================================
function imgVis=categoryPreviewBackground(imgIni,isK)
%categoryPreviewBackground  Enhanced background for the _rep_categories page.
%   Mirrors getPixelCategories' internal imgVis prep so a copied file's page
%   matches a directly-segmented one (keep in sync with getPixelCategories).
if isK
    img=imgIni;
    img(img(:)>prctile(img(:),99))=prctile(img(:),99);
    img(img(:)<prctile(img(:),1))=prctile(img(:),1);
    img=imcomplement(img);
else
    img=imgIni;
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
