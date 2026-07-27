%runCategories  Pre-process LSCI *_K_d.mat files and write categorical masks
%
%   runCategories(s,fNames) iterates over every file name in the cell array
%   fNames, builds a five-level categorical mask (background, parenchyma,
%   unsegmented pixels, vessel walls, lumen), updates the companion
%   *_s.mat (settings) and *_r.mat (results) files, and saves a JPEG
%   preview of the segmentation overlay.
%
%   INPUTS
%     s        – parameter structure created by the LSCI processing library.
%                Fields read here include: iniSizeN, trustLimitsK, regionsN,
%                lSizeN, sSizeN, sSizeScale, sens, deSens, lThinN, imOpen,
%                iEdge, eEdge.
%     fNames   – cell array of char vectors or strings containing full paths
%                to *_K_d.mat files.  Each file must have matching *_s.mat
%                and *_r.mat siblings in the same directory.
%
%   SIDE-EFFECTS (per file)
%     * <name>_s.mat   updated ‘settings’ structure (settings.runCategories,
%                      including the computed edgeSize)
%     * <name>_r.mat   updated ‘results’ structure with .cMask, .regionsMask,
%                      and .mask fields added/overwritten
%     * <name>_cm.jpg  segmentation preview (300 dpi)
%
%   EXAMPLE
%     p = defaultParams();                              % user helper
%     D = dir(fullfile(dataRoot,'**','*_K_d.mat'));
%     runCategories(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     getPixelCategories (categorization core; computes edgeSize and calls
%     getEdgeSizeSLSCI internally), enhanceForDisplay (ROI-editor preview).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

% %Example of s structure parametrisation
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
% s.trustLimitsK=[0.001,0.5]; %minimum (first value, fastest flows) and maximum (second value, slowest flows) expected contrast. Usually [0.01,0.3], but can be e.g. [0.01,0.5] for stroke
% s.regionsN=2; %Numer of regions for manual selection. 0 if using entire window.
% s.lSizeN=61; % Odd, approximately 2 times larger than the largest vessel
% s.sSizeN=15; % Odd, approximately 2 times larger than small vessels diameter
% s.sens=0.3; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
% s.deSens=1; % integer, 1 or higher. 1 - maximally sensitive to smallest vessels but also noisy larger values reduce small vessel sensitivity
% s.sSizeScale=1; % scaler for small objects assignment to background or to unregognized regions
% % %ADJUSTED IF NECESSARY - SEGMENTATION ADJUSTEMNTS
% s.lThinN=2; % Large vessels thinning (appears as internal edges)
% s.imOpen=2; % Small vessels thinning (appears as internal edges)
% s.iniSizeN=7; % Odd number equal or larger than the spatial contrast kernel
% % %DO NOT CHANGE - META DATA
% s.categories={'background','parenchyma','unsegmented','externalWall','internalWall','lumen'}; %CATEGORIES

function runCategories(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat')|| contains(s,'_I_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat" or "_I_d.mat".');
end

for fidx=1:1:numel(fNames)
     if ~isempty(fNames{fidx})
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');

    % --- mean image + modality ---
    isK=contains(s.fName,'_K_d.mat');
    if isK
        if isfield(results,'imgK')
            imgIni=results.imgK;
        else
            load(s.fName,'source');
            imgIni=mean(source.data,3,'omitmissing');
            clearvars source
        end
    else
        imgIni=results.imgI;
    end

    % --- interactive region selection (draws on the display-enhanced image) ---
    regionsMask=zeros(size(imgIni));
    if s.regionsN>0
        % Display enhancement for the ROI editor - mirrors getPixelCategories' own
        % prep so the user draws on the same image the core categorizes.
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
        img=enhanceForDisplay(img,fSize,min(15,fSize));

        f=figure(1);
        f.WindowState='maximized';
        tiledlayout(1,2,"TileSpacing",'compact','Padding','compact');
        nexttile
        imagesc(rot90(imgIni , size(imgIni,2) > size(imgIni,1)))
        clim(prctile(imgIni(:),[1,99]))
        axis image
        title('Original image')
        nexttile
        imagesc(rot90(img, size(img,2) > size(img,1)))
        clim(prctile(img(:),[10,99]))
        axis image
        title('Select regions for segmentation')
        for i=1:1:s.regionsN
            if fidx==1
                [BW,xi,yi] = roipoly;
                x{i}=xi;
                y{i}=yi;
                h = drawpolygon(gca, ...         % recreate an editable polygon
                    'Position',[x{i} y{i}], ...
                    'Color', 'white',...
                    'FaceAlpha',0);
                h.InteractionsAllowed = 'none';
            else
                hold on
                h = drawpolygon(gca, ...         % recreate an editable polygon
                    'Position',[x{i} y{i}], ...
                    'Color', 'white',...
                    'StripeColor', 'blue',...
                    'FaceAlpha',0);
                hold off
                wait(h);
                BW=createMask(h);
                x{i}=h.Position(:,1);
                y{i}=h.Position(:,2);
            end
            regionsMask=regionsMask+BW.*i;
        end
        close(f);
    else
        if isfield(results,'regionsMask')
            regionsMask=results.regionsMask;
        else
            regionsMask=regionsMask+1;
        end
    end

    % --- categorize (automatic core) ---
    if isfield(results,'mask'), existingMask=results.mask; else, existingMask=[]; end
    [cMask,s.edgeSize,maskOut,imgVis]=getPixelCategories(imgIni,regionsMask,existingMask,isK,s);

    results.regionsMask=regionsMask;
    results.mask=(maskOut==(regionsMask>0));

    f=figure(1);
    f.WindowState='maximized';
    tiledlayout(1,2,"TileSpacing",'compact','Padding','compact');
    nexttile
    imagesc(rot90(imgVis , size(imgVis,2) > size(imgVis,1)))
    hold on
    visboundaries(rot90(cMask>4, size(cMask,2) > size(cMask,1)),'Color','m')
    visboundaries(rot90(cMask>2, size(cMask,2) > size(cMask,1)),'Color','w')
    hold off
    clim(prctile(imgVis(:),[10,99]))
    axis image
    nexttile
    imagesc(rot90(cMask, size(cMask,2) > size(cMask,1)))
    axis image
    print(f,strrep(s.fName,'.mat','_cm.jpg'), '-djpeg', '-r300');
    results.cMask=cMask;
    settings.runCategories=s;
    %Save the data
    disp(['Saving the results. Elapsed time ',num2str(round(toc)),'s']);
    save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
    save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
    disp('Saving complete');
     end
end
end
