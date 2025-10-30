%getCategories  Pre-process LSCI *_K_d.mat files and write categorical masks
%
%   getCategories(s,fNames) iterates over every file name in the cell array
%   fNames, builds a five-level categorical mask (background, parenchyma,
%   unsegmented pixels, vessel walls, lumen), updates the companion
%   *_s.mat (settings) and *_r.mat (results) files, and saves a JPEG
%   preview of the segmentation overlay.
%
%   INPUTS
%     s        – parameter structure created by the LSCI processing library.
%                Fields read here include: iniSizeN, minK, maxK, regionsN,
%                lSizeN, sSizeN, sThinN, sens, lThinN.
%     fNames   – cell array of char vectors or strings containing full paths
%                to *_K_d.mat files.  Each file must have matching *_s.mat
%                and *_r.mat siblings in the same directory.
%
%   SIDE-EFFECTS (per file)
%     * <name>_s.mat   updated ‘settings’ structure
%     * <name>_r.mat   updated ‘results’ structure with .cMask, .regionsMask,
%                      and .mask fields added/overwritten
%     * <name>_cm.jpg  segmentation preview (300 dpi)
%
%   EXAMPLE
%     p = defaultParams();                              % user helper
%     D = dir(fullfile(dataRoot,'**','*_K_d.mat'));
%     getCategories(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     getEdgeSizeSLSCI, plus other core functions in the LSCI toolbox.
%
%   --------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: This header was generated with ChatGPT and may contain minor
%   inconsistencies—please verify before release.
%   --------------------------------------------------------------------

% %Example of s structure parametrisation
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
% s.maxK=0.3; % Maximum valid contrast - helps with initial masking
% s.minK=0.0001; % Minimum valid contrast
% s.regionsN=2; %Numer of regions for manual selection. 0 if using entire window.
% s.lSizeN=61; % Odd, approximately 2 times larger than the largest vessel
% s.sSizeN=15; % Odd, approximately 2 times larger than small vessels diameter
% s.sens=0.3; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
% % %ADJUSTED IF NECESSARY - SEGMENTATION ADJUSTEMNTS
% s.lThinN=2; % Large vessels thinning (appears as internal edges)
% s.sThinN=2; % Small vessels thinning (appears as internal edges)
% s.iniSizeN=7; % Odd number equal or larger than the spatial contrast kernel
% % %DO NOT CHANGE - META DATA
% s.categories={'background','parenchyma','unsegmented','externalWall','internalWall','lumen'}; %CATEGORIES

function getCategories(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

for fidx=1:1:numel(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    load(s.fName,'source')
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');

    imgIni=mean(source.data,3);
    s.edgeSize=getEdgeSizeSLSCI(imgIni,0.8);
    img=imgIni;
    img(img(:)>prctile(img(:),99))=prctile(img(:),99);
    img(img(:)<prctile(img(:),1))=prctile(img(:),1);
    img=imcomplement(img);
    fSize=floor((min(size(img))./20))*2+1;
    img(isnan(img))=0;
    img=img-imopen(medfilt2(img,[fSize,fSize],"symmetric"),strel('disk',fSize));
    imgVis=img;

    mask=ones(size(img));
    if isfield(results,'mask') && size(results.mask,1)==size(img,1) && size(results.mask,2)==size(img,2)
        mask=double(results.mask>0);
    end
    mask=mask.*ordfilt2(imgIni,1,ones(s.iniSizeN),'symmetric')>=s.minK & ordfilt2(imgIni,s.iniSizeN.*s.iniSizeN,ones(s.iniSizeN),'symmetric')<=s.maxK;

    regionsMask=zeros(size(mask));
    if s.regionsN>0
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
    mask=(regionsMask>0).*mask;
    results.regionsMask=regionsMask;
    results.mask=(mask==(regionsMask>0));

    tmp=zeros(size(img));
    tmp(s.edgeSize+1:end-s.edgeSize,s.edgeSize+1:end-s.edgeSize)=mask(s.edgeSize+1:end-s.edgeSize,s.edgeSize+1:end-s.edgeSize);
    mask=(tmp==1);
    maskIni=mask;

    for i=1:1:max(regionsMask(:))
        img(regionsMask(:)==i)=mat2gray(img(regionsMask(:)==i),double(prctile(img(regionsMask(:)==i),[5,99])));
    end
    img=img.*mask+(1-mask);
    img(isnan(img))=0;


    img=padarray(img,[s.lSizeN,s.lSizeN],'symmetric');
    mask=padarray(mask,[s.lSizeN,s.lSizeN],0);
    tmp=zeros(size(img));
    for i=s.sSizeN:2:s.lSizeN
        tmp=tmp+imbinarize(img,adaptthresh(img,s.sens,'NeighborhoodSize',i)).*mask;
    end
    maskV=tmp>0;
    maskV=maskV+(conv2(maskV,[1,1,1;1,0,1;1,1,1],'same')>4); %connect pixels with more than 4 neighbours
    maskV=bwareaopen(maskV,s.sSizeN.*s.sSizeN); % remove tiny specs
    maskV   = imclose(maskV,strel('disk',s.sThinN));
    maskV   = imopen(maskV,strel('disk',s.sThinN));
    maskV=bwareaopen(maskV,s.lSizeN.*s.sSizeN); % remove large specs
    maskVEE=imdilate(maskV,strel("disk",s.sThinN)).*mask;
    maskVIE=maskV;
    maskV=bwmorph(maskV,'thin',s.lThinN);

    cMask=int32(mask);
    cMask(tmp(:)>0)=2;
    cMask(maskVEE(:)==1)=3;
    cMask(maskVIE(:)==1)=4;
    cMask(maskV(:)==1)=5;

    img=img(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN);
    cMask=cMask(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN);

    f=figure(1);
    f.WindowState='maximized';
    tiledlayout(1,2,"TileSpacing",'compact','Padding','compact');
    nexttile
    imagesc(rot90(imgVis , size(imgVis,2) > size(imgVis,1)))
    hold on
    visboundaries(rot90(cMask>3, size(cMask,2) > size(cMask,1)))
    hold off
    clim(prctile(imgVis(:),[10,99]))
    axis image
    nexttile
    imagesc(rot90(cMask, size(cMask,2) > size(cMask,1)))
    axis image
    print(f,strrep(s.fName,'.mat','_cm.jpg'), '-djpeg', '-r300');
    results.cMask=cMask;
    settings.categoricalMask=s;
    %Save the data
    disp(['Saving the results. Elapsed time ',num2str(round(toc)),'s']);
    save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
    save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
    disp('Saving complete');
end
end