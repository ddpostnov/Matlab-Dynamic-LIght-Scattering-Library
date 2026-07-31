%runDynamicSegmentation  Track vessel centre-lines over time -> dynamic diameter/flow
%
%   runDynamicSegmentation(s,fNames) is the optional heavy dynamic-segmentation step,
%   lifted verbatim from the old runSegmentation `attmemptDS` loop.  For every
%   *_K_d.mat / *_I_d.mat file it re-derives the labelled centre-lines from the static
%   segmentation, fits each straight-enough vessel with a 3-degree polynomial, walks a
%   perpendicular intensity profile frame by frame to estimate the lumen edges, keeps
%   only segments that pass the geometry / overlap quality gates, and writes their
%   per-frame diameter and flow traces.
%
%   SELF-CONTAINED (reads only the persisted seam).  It needs no transient state from
%   runSegmentation: it reloads the source cube (_d), reads results.cMask and
%   settings.runSegmentation.edgeSize (_r/_s), and RECOMPUTES the labelling transients
%   (sLines, cMask-merged, vsMap, dMask, nodes) with getSegmentationLabels exactly as
%   the static step did - getSegmentationLabels is deterministic, so these are the same
%   values the old inline block consumed.  Run it AFTER runSegmentation on the same
%   file(s).
%
%   INPUTS
%     s        parameter structure.
%              Labelling (read by getSegmentationLabels, must match the runSegmentation
%                call): correctNodes, simR, difR, sMinL, prchNSize.
%              Dynamic segmentation: sMinL, sMinP2R2, sMaxLBI, sMaxCLR, sMaxKK, sMaxP2D,
%                pInterpF, gSizeN, minOverlapMask, minOverlapSelf.
%     fNames   cell array of *_K_d.mat / *_I_d.mat paths already processed by
%              runSegmentation (each with matching *_s.mat / *_r.mat siblings).
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUT SIDE-EFFECTS (per file)
%     <name>_r.mat   results.{dvsMap,dvsMetrics,dvsDiameter,dvsData}
%     <name>_s.mat   settings.runDynamicSegmentation = s
%     <name>_rep_segments.jpg  segments page, RE-EMITTED with the accepted dynamic
%                    segments (visboundaries(dvsMap>0)) overlaid on the enhanced mean
%                    image - overwrites the static page runSegmentation wrote, so the
%                    final artifact matches the pre-refactor attmemptDS preview.
%
%   EXAMPLE
%     runSegmentation(s,fNames);
%     runDynamicSegmentation(s,fNames);
%
%   DEPENDS ON
%     getSegmentationLabels, showSegmentsPreview (Core/LSCI) + MATLAB's Image
%     Processing Toolbox (bwdist, islocalmax, bwareaopen, bwareafilt, polyfit, pca,
%     imresize3, ...).
%
% See also: runSegmentation, setRegions, getSegmentationLabels
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 27-July-2026


%%Example of s structure parametrisation
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
% s.sMinL=15;       % Minimum length for segments (also gates getSegmentationLabels)
% s.sMinP2R2=0.95;  % Min accepted R2 of the 3-degree polynomial fit
% s.sMaxLBI=(1/5)./s.sMinL; % Max local bending (0 to pi per pixel)
% s.sMaxCLR=1.3;    % Maximum accepted CLR (1 straight, 1.5 slow bend, 2 coil)
% s.sMaxKK=0.3;     % Max accepted std/mean for the initial contrast estimation
% s.sMaxP2D=3;      % Max accepted deviation of the fit from the centre estimate
% s.pInterpF=10;    % interpolation factor (leave as is)
% s.gSizeN=3;       % gap-fill size
% s.minOverlapMask=0.6; % min overlap between the centre line and the segmentation mask
% s.minOverlapSelf=0.2; % min segmented area vs the initial ROI

function runDynamicSegmentation(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_K_d.mat')|| contains(x,'_I_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat" or "_I_d.mat".');
end
% Resolved here, not only in the core, so the choice is RECORDED in the saved
% settings like every other tunable.  getSegmentationLabels reads it as a worker
% bound on its per-label parfor; true (the default) is today's behaviour.
if ~isfield(s,'parforSegmentationLabels') || isempty(s.parforSegmentationLabels)
    s.parforSegmentationLabels=true;
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Dynamic segmentation',fNames);

for fidx=1:1:numel(fNames)
     if reportCancelled(rep), break; end        % cooperative cancel between files
     if ~isempty(fNames{fidx})
    s.fName=fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results source settings
    load(s.fName,'source')
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');

    % --- re-derive the labelling transients from the persisted seam (option a) ---
    edgeSize = settings.runSegmentation.edgeSize;
    [~,~,sLines,cMask,vsMap,dMask,nodes] = getSegmentationLabels(results.cMask,edgeSize,s);

    % --- dynamic segmentation (verbatim from the old runSegmentation attmemptDS loop) ---
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
    img=mean(source.data,3,'omitnan');
    img=img.*tmp;
    img=imcomplement(img);
    d2C=bwdist(~dMask).*(dMask);
    d2MY=islocalmax(img,1).*(cMask==5);
    [~,d2MY]=bwdist(d2MY);
    d2MX=islocalmax(img,2).*(cMask==5);
    [~,d2MX]=bwdist(d2MX);
    dvsDiameter=zeros(size(source.data,3),max(sLines(:)));
    varTypes = ["int32","single","single","single","single","single"];
    varNames = ["idx","length","CLR","R2","overlapMask","overlapSelf"];
    dvsMetrics=table('Size',[max(sLines(:)),6],'VariableTypes',varTypes,'VariableNames',varNames);
    dvsData=zeros(size(source.data,3),max(sLines(:)));
    dvsMap=zeros(size(cMask),'int32');


    segIdxs=unique(sLines(sLines(:)>0))';
    nSeg=numel(segIdxs);
    counter=1;
    for segNo=1:nSeg
        lineIdx=segIdxs(segNo);
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
                sD=sD.*s.pInterpF;
                xx = round( (xx - limX(1)) * s.pInterpF ) + 1;
                yy = round( (yy - limY(1)) * s.pInterpF ) + 1;

                xx=xx(:);
                yy=yy(:);
                stepDist = hypot(diff(xx), diff(yy));

                % Identify gaps larger than adjacent diagonals but smaller than the threshold
                gapIndices = find(stepDist > 1.5 & stepDist < s.gSizeN.*s.pInterpF);

                for k = flip(gapIndices(:)')
                    pointCount = ceil(stepDist(k));
                    % Generate linear coordinate insertions
                    fillX = round(linspace(xx(k), xx(k+1), pointCount + 1)');
                    fillY = round(linspace(yy(k), yy(k+1), pointCount + 1)');
                    % Splice arrays to insert interpolated pixels, excluding existing endpoints
                    xx = [xx(1:k); fillX(2:end-1); xx(k+1:end)];
                    yy = [yy(1:k); fillY(2:end-1); yy(k+1:end)];
                end




                v     = pca([xx(:) yy(:)]);
                theta = atan2d(v(2,1), v(1,1));

                maskROI=single(~((vsMap(limY(1):limY(2),limX(1):limX(2))~=lineIdx & vsMap(limY(1):limY(2),limX(1):limX(2))~=0) | cMask(limY(1):limY(2),limX(1):limX(2))==2));
                dataROI=single(source.data(limY(1):limY(2),limX(1):limX(2),:));
                compVal=max(dataROI(:));
                maskROI(maskROI==0)=NaN;
                if contains(s.fName,'_K_d.mat')
                    dataROI=compVal-dataROI;
                elseif contains(s.fName,'_I_d.mat')
                    dataROI=dataROI-min(dataROI(dataROI(:)>0));
                end
                dataROI=dataROI.*maskROI;
                dataROI=imresize3(dataROI,[size(dataROI,1)*s.pInterpF,size(dataROI,2)*s.pInterpF,size(dataROI,3)]);


                if (max(xx)-min(xx))>=(max(yy)-min(yy))
                    try
                        dataProfile=nan(numel(xx),sum(sD)*2+1,size(dataROI,3));
                    catch
                        continue;               % too large to fit in memory - skipped
                    end

                    if (min(yy)-sum(sD))>0 && (max(yy)+sum(sD))<size(dataROI,1)
                        for i=1:1:numel(xx)
                            dataProfile(i,:,:)=dataROI(yy(i)-sum(sD):yy(i)+sum(sD),xx(i),:);
                        end
                    else
                        continue;               % reaches outside the image - skipped
                    end
                else
                    try
                        dataProfile=nan(numel(yy),sum(sD)*2+1,size(dataROI,3));
                    catch
                        continue;               % too large to fit in memory - skipped
                    end
                    if (min(xx)-sum(sD))>0 && (max(xx)+sum(sD))<size(dataROI,2)
                        for i=1:1:numel(yy)
                            dataProfile(i,:,:)=dataROI(yy(i),xx(i)-sum(sD):xx(i)+sum(sD),:);
                        end

                    else
                        continue;               % reaches outside the image - skipped
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


                if test1>=s.minOverlapMask && test2>=s.minOverlapSelf && test3
                    tmp=mean(tmp,3)>0;
                    [y,x]=find(tmp>0);
                    x=x+limX(1)-1;
                    y=y+limY(1)-1;

                    dvsMap(sub2ind(size(cMask), y, x))=lineIdx;
                    dvsMetrics(counter,:)={lineIdx,sL,sCLR,sR2,test1,test2};
                    dvsDiameter(:,counter)=(idxR-idxL)* abs(sind(theta))./s.pInterpF;
                    for i=1:1:size(dataProfile,2)
                        if contains(s.fName,'_K_d.mat')
                            dvsData(i,counter)=compVal-mean(dataProfile(idxL(i):idxR(i),i),1,'omitnan');
                        elseif contains(s.fName,'_I_d.mat')
                            dvsData(i,counter)=mean(dataProfile(idxL(i):idxR(i),i),1,'omitnan')-min(dataROI(dataROI(:)>0));
                        end
                    end
                    counter=counter+1;
                end                             % else: did not fit well enough - skipped
            end
        end
    end
    dvsData(:,counter:end)=[];
    dvsDiameter(:,counter:end)=[];
    dvsMetrics(counter:end,:)=[];

    results.dvsMap=dvsMap;
    results.dvsMetrics=dvsMetrics;
    results.dvsDiameter=dvsDiameter;
    results.dvsData=dvsData;

    % --- segments preview: re-emit the page with the accepted dynamic segments overlaid ---
    % (overwrites the static one runSegmentation wrote, matching the old attmemptDS
    %  artifact) using the merged cMask (recomputed above) and the persisted results.sMap.
    isK=contains(s.fName,'_K_d.mat');
    fh=reportFigure(rep,'segments');
    showSegmentsPreview(fh,source.data,cMask,results.sMap,isK,dvsMap);
    reportSave(rep,fh,'segments',s.fName);

    %Save the data
    settings.runDynamicSegmentation=reportSettings(s);
    reportWriting(rep);
    save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
    save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
    reportSaved(rep);
     end
end
reportClose(rep);

end
