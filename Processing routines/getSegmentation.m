%getSegmentation  Segment vessels & parenchyma from categorical LSCI masks
%
%   getSegmentation(s,fNames) loads each *_K_d.mat dataset, converts the
%   pre-computed categorical mask (results.cMask) into a detailed vessel/
%   parenchyma label map (results.sMap), extracts scalar metrics and mean
%   time-series for every labelled structure, and—optionally—performs an
%   automated dynamic-segmentation (centre-line tracking) to obtain diameter
%   and flow traces (results.dvs* fields).  Updated *_r.mat and *_s.mat
%   files overwrite the originals, and a JPEG preview of the segmentation
%   is saved alongside the data.
%
%   INPUTS
%     s        parameter structure  
%              ├─ attmemptDS           true = attempt dynamic segmentation  
%              ├─ sMinL, prchNSize     geometric thresholds  
%              ├─ sMinP2R2‥sMaxP2D     dynamic-segment quality limits  
%              ├─ minOverlapMask/Self  acceptance criteria  
%              └─ pInterpF, …          interpolation & kernel sizes
%     fNames   cell array of *_K_d.mat file paths (contrast cubes).
%
%   OUTPUT SIDE-EFFECTS
%     <run>_K_r.mat   RESULTS with fields  
%                       • sMap,  sMetrics,  sData  
%                       • dvsMap, dvsMetrics, dvsDiameter, dvsData  (if DS)  
%     <run>_K_s.mat   SETTINGS   (field settings.vesselsSegmentation added)
%     <run>_vs.jpg    Preview image: raw contrast + vessel labels
%
%   EXAMPLE
%     p = defaultSegmentationParams();
%     files = dir(fullfile(dataRoot,'*_K_d.mat'));
%     getSegmentation(p, fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     Functions from the LSCI processing library: readRLS, getFFT, etc.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: This header was generated with ChatGPT and may contain minor
%   inconsistencies—please verify before release.
%   ----------------------------------------------------------------------


%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - BASIC PARAMETERS
% s.attmemptDS=true; %attempt to perform automated dynamic segmentation or not
% s.sMinL=15; % Minimum length for segments
% s.prchNSize=30; % Parenchymal pixels neighbourhoud.
%
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
% s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
% s.sMaxLBI=(1/5)./s.sMinL; %Max local bending (0 to pi per pixel)
% s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
% s.sMaxDK=0.2; %Max accepted std/mean for the initial diameter estimation
% s.sMaxKK=0.3; %Max accepted std/mean for the initial contrast estimation
% s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
% s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate
%
% %ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
% s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
% s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
% s.pInterpF=10; % leave as is

function getSegmentation(s,fNames)

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

%temporary fix - diameter is estimated including inner walls, rest of
%parameters estimated for lumen and merged walls
    cMask=results.cMask;
    dMask=cMask>3;
    cMask(cMask==4)=3;


    edgeSize=settings.categoricalMask.edgeSize;
    sLines=bwskel(cMask==5);
    tmp=zeros(size(sLines));
    tmp(edgeSize*2+1:end-edgeSize*2,edgeSize*2+1:end-edgeSize*2)=sLines(edgeSize*2+1:end-edgeSize*2,edgeSize*2+1:end-edgeSize*2);
    sLines=(tmp==1);
    nodes=logical((conv2(sLines,[1,1,1;1,0,1;1,1,1],'same')>2).*sLines);
    sLines=sLines-nodes;
    sLines=bwareaopen(sLines,max(2,floor(s.sMinL./2)),8);
    sLines=int32(bwlabel(sLines)); %labeled segment center lines

    labels=nonzeros(unique(sLines));
    distStack=inf([size(cMask) numel(labels)],'single');
    for k = 1:numel(labels)
        distStack(:,:,k) = bwdistgeodesic(cMask>2, sLines == labels(k), 'quasi-euclidean');
    end
    [~,tmp] = min(distStack,[],3);
    vsMap     = zeros(size(cMask),'like',sLines);
    vsMap(cMask>2) = labels(tmp(cMask>2)); %labeled segments
    vsMap(vsMap>0)=(vsMap(vsMap>0)-1).*2+1;
    sLines(sLines>0)=(sLines(sLines>0)-1).*2+1;

    sMap=vsMap;
    sMap(cMask==3)=sMap(cMask==3)+1;
    idxs=int32(bwlabel(cMask==2))+max(sMap(:));
    sMap(cMask==2)=idxs(cMask==2);


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
    sMap(valid & cMask==1)         = lbl(valid & cMask==1);                   % updated label map




    %Build segments table
    varTypes = ["double","double","double","double","double","double","double","double"];
    varNames = ["idx","category","length","CLR","diameter","std(diameter)","area","nearestVesIdx"];
    sMetrics=table('Size',[max(sMap(:)),8],'VariableTypes',varTypes,'VariableNames',varNames);
    sData=zeros(size(source.data,3),max(sMap(:)));
    data=reshape(source.data,[],size(source.data,3));
    for i=1:1:max(sMap(:))
        area=sum(sMap(:)==i);
        if area>0
            sData(:,i)=mean(data(sMap(:)==i,:),1,'omitnan');
            c=unique(cMask(sMap==i));
            if numel(c)~=1
                error("Mix of categories detected in the region");
            end

            if c==5
                %Measure segment length
                ends=bwmorph(sLines==i,'endpoints');
                [y,x]  = find(ends);
                if numel(y)~=2
                    error('Incorrect number of ends found')
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

    nodesD=bwdist(~dMask).*nodes;
    [tmp, idxs] = bwdist(nodes);
    rNode = nodesD(nodes==1);
    r=zeros(numel(nodes),1);
    r(nodes==1)=rNode;
    rP = r(idxs);
    nodesD = (tmp <= rP);
    sLines(nodesD)=0;
    sLines=sLines.*int32(bwareaopen(sLines>0,s.sMinL,8));


    tmp=ones(size(cMask));
    tmp(cMask==0)=nan;
    data=source.data.*tmp;
    data=imcomplement(data);
    d2C=bwdist(~dMask).*(dMask);
    d2MY=islocalmax(mean(data,3,'omitnan'),1).*(cMask==5);
    [~,d2MY]=bwdist(d2MY);
    d2MX=islocalmax(mean(data,3,'omitnan'),2).*(cMask==5);
    [~,d2MX]=bwdist(d2MX);
    dvsDiameter=zeros(size(source.data,3),max(sLines(:)));
    varTypes = ["int32","single","single","single","single","single"];
    varNames = ["idx","length","CLR","R2","overlapMask","overlapSelf"];
    dvsMetrics=table('Size',[max(sLines(:)),6],'VariableTypes',varTypes,'VariableNames',varNames);
    dvsData=zeros(size(source.data,3),max(sLines(:)));
    dvsMap=zeros(size(cMask),'int32');


    counter=1;
    if s.attmemptDS
        for lineIdx=unique(sLines(sLines(:)>0))'
            disp(['Checking segment ',num2str(lineIdx),' out of possible ',num2str(max(sLines(:)))])
            [y,x]=find(sLines==lineIdx);
            if (max(x)-min(x))>=(max(y)-min(y))
                sLines(d2MY(sub2ind(size(sLines),y,x)))=lineIdx;
                [y,x]=find(sLines==lineIdx);
                [p,S,mu] = polyfit(x,y,3);
                xx=min(x):1/s.pInterpF:max(x);
                yy=polyval(p,xx,S,mu);
            else
                sLines(d2MX(sub2ind(size(sLines),y,x)))=lineIdx;
                [y,x]=find(sLines==lineIdx);
                [p,S,mu] = polyfit(y,x,3);
                yy=min(y):1/s.pInterpF:max(y);
                xx=polyval(p,yy,S,mu);
            end
            kappa=abs(6.*p(1).*yy+2.*p(2))./((1+(3.*p(1).*(yy.^2)+2.*p(2).*yy+p(3)).^2).^(3./2));
            sLines(sLines==lineIdx)=0;
            dd=min(hypot(x-xx,y-yy),[],1);
            idxs=dd<=s.sMaxP2D & kappa<=s.sMaxLBI & xx>=1 & yy>=1 & xx<=size(sLines,2) & yy<=size(sLines,1);
            xx=xx(idxs);
            yy=yy(idxs);

            %Leaving only the longest part of the segment
            tmp=zeros(size(sLines),'logical');
            tmp(sub2ind(size(sLines),round(yy),round(xx)))=1;
            tmp=bwareafilt(tmp,1,"largest",8);
            [y,x]=find(tmp);
            idxs=sum(round(yy)==y & round(xx)==x,1)>0;
            xx=xx(idxs);
            yy=yy(idxs);
            if numel(xx)>1
                sL=sum(hypot(diff(xx),diff(yy)));
                sD=ceil(2*[mean(d2C(sub2ind(size(sLines),round(yy),round(xx)))),std(d2C(sub2ind(size(sLines),round(yy),round(xx))))]);
                sCLR=sL/ hypot(xx(end)-xx(1), yy(end)-yy(1));
                sR2=S.rsquared;
                sLines(sub2ind(size(sLines),round(yy),round(xx)))=lineIdx;
                if (max(xx)-min(xx))>=(max(yy)-min(yy))
                    limX=([floor(min(xx)),ceil(max(xx))]);
                    limY=([floor(min(yy))- sD(1)-2*sD(2),ceil(max(yy))+ sD(1)+2*sD(2)]);
                else
                    limY=([floor(min(yy)),ceil(max(yy))]);
                    limX=round([floor(min(xx))- sD(1)-2*sD(2),ceil(max(xx))+ sD(1)+2*sD(2)]);
                end

                if sL>=s.sMinL && sR2>=s.sMinP2R2 && sCLR<=s.sMaxCLR && sD(2)/sD(1)<=s.sMaxKK && limX(1)>edgeSize && limY(1)>edgeSize && limX(2)<size(cMask,2)-edgeSize && limY(2)<size(cMask,1)-edgeSize
                    disp('Fitting a segment')

                    sD=sD.*s.pInterpF;
                    xx = round( (xx - limX(1)) * s.pInterpF ) + 1;
                    yy = round( (yy - limY(1)) * s.pInterpF ) + 1;
                    v     = pca([xx(:) yy(:)]);
                    theta = atan2d(v(2,1), v(1,1));

                    maskROI=single(~((vsMap(limY(1):limY(2),limX(1):limX(2))~=lineIdx & vsMap(limY(1):limY(2),limX(1):limX(2))~=0) | cMask(limY(1):limY(2),limX(1):limX(2))==2));
                    dataROI=source.data(limY(1):limY(2),limX(1):limX(2),:);
                    compVal=max(dataROI(:));
                    maskROI(maskROI==0)=NaN;
                    dataROI=compVal-dataROI;
                    dataROI=dataROI.*maskROI;
                    dataROI=imresize3(dataROI,[size(dataROI,1)*s.pInterpF,size(dataROI,2)*s.pInterpF,size(dataROI,3)]);


                    if (max(xx)-min(xx))>=(max(yy)-min(yy))
                        dataProfile=nan(numel(xx),sum(sD)*2+1,size(dataROI,3));
                        if (min(yy)-sum(sD))>0 && (max(yy)+sum(sD))<size(dataROI,1)
                            for i=1:1:numel(xx)
                                dataProfile(i,:,:)=dataROI(yy(i)-sum(sD):yy(i)+sum(sD),xx(i),:);
                            end
                        else
                            disp('Segmentation failed - out of bounds')
                            continue;
                        end
                    else
                        dataProfile=nan(numel(yy),sum(sD)*2+1,size(dataROI,3));
                        if (min(xx)-sum(sD))>0 && (max(xx)+sum(sD))<size(dataROI,2)
                            for i=1:1:numel(yy)
                                dataProfile(i,:,:)=dataROI(yy(i),xx(i)-sum(sD):xx(i)+sum(sD),:);
                            end

                        else
                            disp('Segmentation failed - out of bounds')
                            continue;
                        end
                    end

                    dataProfile=squeeze(mean(dataProfile,1,'omitnan'));

                    idxIni=zeros(1,size(dataProfile,2));
                    idxL=zeros(1,size(dataProfile,2));
                    idxR=zeros(1,size(dataProfile,2));
                    for t=1:1:size(dataProfile,2)
                        ts=squeeze(dataProfile(:,t));
                        idxsFrgrd=zeros(1,length(ts));
                        [~,idxIni(t)]=max(ts);
                        idxCur=idxIni(t);
                        idxsFrgrd(idxCur)=1;

                        idxL(t)=idxCur;
                        idxR(t)=idxCur;

                        stdFrgrd=std(ts(idxsFrgrd==1));
                        stdBkgrd=std(ts(idxsFrgrd==0));

                        stdSumIni=sqrt(sum(idxsFrgrd==1).*stdFrgrd.^2+sum(idxsFrgrd==0).*stdBkgrd.^2);
                        stdSumCur=stdSumIni;
                        stdSum2=stdSumCur;
                        stdSumCur2=stdSum2;
                        while stdSumCur2<=stdSum2 && idxL(t)>1 && idxR(t)<length(idxsFrgrd)
                            stdSum2=stdSumCur2;
                            stdSum=stdSumCur;
                            idxCur=idxL(t);
                            while stdSumCur<=stdSum && idxCur>1
                                stdSum=stdSumCur;
                                idxCur=idxCur-1;
                                idxsFrgrd(idxCur)=1;
                                stdFrgrd=std(ts(idxsFrgrd==1));
                                stdBkgrd=std(ts(idxsFrgrd==0));
                                stdSumCur=sqrt(sum(idxsFrgrd==1).*stdFrgrd.^2+sum(idxsFrgrd==0).*stdBkgrd.^2);
                            end
                            idxL(t)=idxCur;

                            stdSum=stdSumCur;
                            idxCur=idxR(t);
                            while stdSumCur<=stdSum && idxCur<length(idxsFrgrd)
                                stdSum=stdSumCur;
                                idxCur=idxCur+1;
                                idxsFrgrd(idxCur)=1;
                                stdFrgrd=std(ts(idxsFrgrd==1));
                                stdBkgrd=std(ts(idxsFrgrd==0));
                                stdSumCur=sqrt(sum(idxsFrgrd==1).*stdFrgrd.^2+sum(idxsFrgrd==0).*stdBkgrd.^2);
                            end
                            idxR(t)=idxCur;
                            stdSumCur2=stdSumCur;
                        end
                    end

                    tmp=zeros(size(dataROI));
                    if (max(xx)-min(xx))>=(max(yy)-min(yy))
                        for i=1:1:numel(xx)
                            for ii=1:1:size(dataROI,3)
                                tmp(idxL(ii)+yy(i)-sum(sD):idxR(ii)+yy(i)-sum(sD),xx(i),ii)=1;
                            end
                        end
                    else
                        for i=1:1:numel(xx)
                            for ii=1:1:size(dataROI,3)
                                tmp(yy(i),idxL(ii)+xx(i)-sum(sD):idxR(ii)+xx(i)-sum(sD),ii)=1;
                            end
                        end
                    end
                    tmp=imresize3(tmp,[numel(limY(1):limY(2)),numel(limX(1):limX(2)),size(dataROI,3) ],'nearest');

                    maskROI=vsMap(limY(1):limY(2),limX(1):limX(2))==lineIdx & cMask(limY(1):limY(2),limX(1):limX(2))>3;
                    test1=(1-sum(abs(maskROI-mean(tmp,3)),'all')./sum(maskROI+mean(tmp,3),'all')./2);
                    test2=(1-sum(abs(mean(tmp,3)-(mean(tmp,3)>0)),'all')./sum((mean(tmp,3)>0),'all'));
                    [~,test3]=bwlabel(mean(tmp,3)>0.9,4);
                    test3=(test3==1);

                    dataROI=source.data(limY(1):limY(2),limX(1):limX(2),:);
                    compVal=max(dataROI(:));
                    dataROI=compVal-dataROI;



                    if test1>=s.minOverlapMask && test2>=s.minOverlapSelf && test3
                        tmp=mean(tmp,3)>0;
                        [y,x]=find(tmp>0);
                        x=x+limX(1)-1;
                        y=y+limY(1)-1;

                        dvsMap(sub2ind(size(cMask), y, x))=lineIdx;
                        dvsMetrics(counter,:)={lineIdx,sL,sCLR,sR2,test1,test2};
                        dvsDiameter(:,counter)=(idxR-idxL)* abs(sind(theta))./s.pInterpF;
                        for i=1:1:size(dataProfile,2)
                            dvsData(i,counter)=compVal-mean(dataProfile(idxL(i):idxR(i),i),1,'omitnan');
                        end
                        counter=counter+1;
                    else
                        disp('Segmentation failed - bad quality')
                    end
                end
            end
        end
        dvsData(:,counter:end)=[];
        dvsDiameter(:,counter:end)=[];
        dvsMetrics(counter:end,:)=[];
    end

    results.sMap=sMap;
    results.sMetrics=sMetrics;
    results.sData=sData;
    if s.attmemptDS
        results.dvsMap=dvsMap;
        results.dvsMetrics=dvsMetrics;
        results.dvsDiameter=dvsDiameter;
        results.dvsData=dvsData;
    end

    img=mean(source.data,3);
    img=mat2gray(img,double(prctile(img(cMask(:)>0),[5,99])));
    img=imcomplement(img);
    fSize=floor((min(size(img))./20))*2+1;
    img=(img-imopen(medfilt2(img,[fSize,fSize],"symmetric"),strel('disk',fSize))).*(cMask>0);

    n=double(max(sMap(:)));
    phi=(sqrt(5)-1)/2;
    cmap=hsv2rgb([mod((0:n-1)'*phi,1) 0.8*ones(n,1) 0.9-0.2*mod((0:n-1)',2)]);

    f=figure(1);
    f.WindowState='maximized';
    t=tiledlayout(1,2,"TileSpacing",'compact','Padding','compact');
    t1=nexttile(t);
    imagesc(img,'Parent', t1)
    axis image
    if s.attmemptDS
        hold on
        visboundaries(dvsMap>0);
        hold off
    end
    t2=nexttile(t);
    sMap=single(sMap);
    sMap(sMap==0)=NaN;
    h=imagesc(sMap,'Parent', t2);
    set(h,'AlphaData',~isnan(sMap));
    axis image
    colormap(t1,parula);
    colormap(t2,cmap)
    t2.Colormap=cmap;
    set(t1,'color',[1 1 1])
    set(t2,'color',[1 1 1])
    set(gcf,'Color','w')
    drawnow
    print(f,strrep(s.fName,'.mat','_vs.jpg'), '-djpeg', '-r300');

    %Save the data
    settings.vesselsSegmentation=s;
    disp(['Saving the results. Elapsed time ',num2str(round(toc)),'s']);
    save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
    save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
    disp('Saving complete');
end
end