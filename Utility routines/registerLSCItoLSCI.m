%registerLSCItoLSCI  Rigid/affine registration of several *_K_d.mat datasets
%
%   registerLSCItoLSCI(s,fNames) aligns every contrast cube in the cell
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
%     fNames   cell array of *_K_d.mat paths.  First file is the template.
%
%   OUTPUT (side-effects)
%       For every file k
%         *_K_d.mat   SOURCE   – warped
%         *_K_r.mat   RESULTS  – warped + consensus masks
%         *_K_s.mat   SETTINGS – field settings.registration added
%       plus *_vs.jpg previews already present in the workflow
%
%   EXAMPLE
%     [opt,met]     = imregconfig("monomodal");
%     opt.MaximumIterations = 500;
%     p.optimizer   = opt;
%     p.metric      = met;
%     p.tFormType   = 'affine';
%     p.matchSegmentation = true;
%     D = dir(fullfile(dataRoot,'*_K_d.mat'));
%     registerLSCItoLSCI(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     manualByPointRegistration, LSCI segmentation utilities, MATLAB’s
%     Image Processing Toolbox (imregtform, imregcorr, etc.).
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: Header generated with ChatGPT; minor inconsistencies may remain.
%   ----------------------------------------------------------------------


% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
% [s.optimizer,s.metric] = imregconfig("monomodal");
% s.optimizer.MaximumIterations=500;
% s.tFormType='affine';
% s.matchSegmentation=true;
% s.prchNSize=30; % Parenchymal pixels neighbourhoud.

function registerLSCItoLSCI(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

tforms=cell(size(fNames));
views=cell(size(fNames));
commonMask=cell(size(fNames));

for fidx=1:1:size(fNames,1)
    clearvars results source settings
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
    img(img(:)>prctile(img(mask(:)==1),99))=prctile(img(mask(:)==1),99);
    img(img(:)<prctile(img(mask(:)==1),1))=prctile(img(mask(:)==1),1);
    img=imcomplement(img);
    fSize=floor((min(size(img))./20))*2+1;
    img(isnan(img))=0;
    img=img-imopen(medfilt2(img,[fSize,fSize],"symmetric"),strel('disk',fSize));
    img=mat2gray(img).*mask;
    img(mask(:)==1)=histeq(img(mask(:)==1));
    imgRef=img;

    if fidx==1
        imgRefIni=imgRef;
        tforms{fidx}=affine2d(eye(3));
    else
        tform1 = imregtform(imgaussfilt(img,fSize./5),imgaussfilt(imgRef,fSize./5),s.tFormType,s.optimizer,s.metric);
        tmp1=imwarp(img,tform1,"OutputView",imref2d(size(imgRefIni)),'FillValues', 0);
        tform2 =  imregcorr(img,imgRefIni,'similarity');
        tmp=imwarp(img,tform2,"OutputView",imref2d(size(imgRefIni)),'FillValues', 0);

        f=figure(1);
        f.WindowState='Maximized';
        tiledlayout(1,4,"TileSpacing","compact","Padding","compact");
        nexttile
        imshowpair(rot90(imgRefIni , size(imgRefIni,2) > size(imgRefIni,1)),rot90(img , size(img,2) > size(img,1)))
        axis image
        if size(imgRefIni)==size(img)
            title(['Original, \Delta=',num2str(round(sum(abs(imgRefIni(:)-img(:)))))])
        else
            title('Original, size mismatch')
        end
        nexttile
        imshowpair(rot90(imgRefIni , size(imgRefIni,2) > size(imgRefIni,1)),rot90(tmp1 , size(tmp1,2) > size(tmp1,1)))
        axis image
        title(['Intensity registration, \Delta=',num2str(round(sum(abs(imgRefIni(:)-tmp1(:)))))]);
        nexttile
        imshowpair(rot90(imgRefIni , size(imgRefIni,2) > size(imgRefIni,1)),rot90(tmp , size(tmp,2) > size(tmp,1)))
        axis image
        title(['Correlation registration, \Delta=',num2str(round(sum(abs(imgRefIni(:)-tmp(:)))))]);

        ax=nexttile;
        axis(ax,'off');
        p = uipanel(f,'Units','normalized', ...
            'Position', ax.Position, ...
            'BorderType','none');

        vals = {'Use original','Use intensity registration','Use correlation registration','Start manual'};
        for k = 1:4
            label = vals{k};                                      % capture once
            uicontrol(p,'Style','pushbutton','String',label, ...
                'Units','normalized','Position',[0.1 1-0.2*k 0.8 0.15], ...
                'Callback',@(src,evt)buttonDone(f,label));
        end
        uiwait(f);                          % wait for a click
        choice = getappdata(f,'choice');    % 1, 2, or 3
        close(f);

        if find(strcmp(vals,choice))==1
            tforms{fidx}=affine2d(eye(3));
        elseif find(strcmp(vals,choice))==2
            tforms{fidx}=tform1;
        elseif find(strcmp(vals,choice))==3
            tforms{fidx}=tform2;
        else
            tforms{fidx}=manualByPointRegistration(imgRefIni,img,'sideBySide');
        end
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

    if ~isempty(fNames{fidx})
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
        settings.registration=s;
        %Save the data
        disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        save(fNames{fidx},'source','-v7.3');
        save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
        save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');
    end
end
disp('Saving complete');

    function buttonDone(f,sel)
        setappdata(f,'choice',sel);           % store selection
        uiresume(f);                          % release uiwait
    end

end