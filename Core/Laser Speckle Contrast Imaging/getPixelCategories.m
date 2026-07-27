%getPixelCategories - binarize an LSCI mean image into a 5-level category mask
%
%   [cMask,edgeSize,mask,imgVis] = getPixelCategories(imgIni,regionsMask, ...
%                                                     existingMask,isK,s)
%
% DESCRIPTION
%   The automatic categorization core lifted verbatim from runCategories: it turns
%   a mean LSCI image into the five-level categorical mask cMask
%     0 background | 1 parenchyma | 2 unsegmented | 3 external wall |
%     4 internal wall | 5 lumen
%   The pipeline is: getEdgeSizeSLSCI -> per-modality prep -> enhanceForDisplay ->
%   trust mask -> per-region normalisation -> anisotropic diffusion -> multiscale
%   top-hat + adaptive threshold -> morphological cleanup.  Only the automatic math
%   lives here; the file I/O and the interactive region GUI stay in the wrapper,
%   which supplies regionsMask.
%
%   The trust mask is automatic (not interactive): an existing same-size mask
%   (existingMask, typically results.mask from runContrastFromRLS) is reused as-is,
%   otherwise - for contrast (_K) data - it is derived from the ordfilt2 min/max
%   of imgIni against s.trustLimitsK.  The mean image imgIni (NOT the enhanced
%   image) drives that ordfilt2, so it is passed alongside the enhancement.
%
% INPUT
%   imgIni        2-D mean image (results.imgK for contrast, results.imgI for
%                 intensity).  Used both for edge-size / enhancement and, for _K,
%                 for the trust-mask ordfilt2.
%   regionsMask   integer label image (0 = excluded).  1..N select regions to
%                 categorize and normalise; the wrapper builds it (GUI or whole
%                 window).  Its size sets the working frame.
%   existingMask  a mask to reuse verbatim (e.g. results.mask), or [] to derive
%                 one.  Reused only when it matches size(imgIni).
%   isK           true for contrast (_K) data, false for intensity (_I).  Selects
%                 the imcomplement + multi-scale diffusion schedule (_K) versus the
%                 plain path (_I); do NOT sniff the file name inside the core.
%   s             parameter struct.  Reads: iniSizeN, trustLimitsK (K only),
%                 lSizeN, sSizeN, sens, deSens, sSizeScale, lThinN, imOpen,
%                 iEdge, eEdge.
%
% OUTPUT
%   cMask     int32 five-level category mask, size(imgIni).
%   edgeSize  border width (px) getEdgeSizeSLSCI found (0 for _I); the scalar the
%             wrapper stores in settings and the label stage reuses.
%   mask      double 0/1 trust-and-region mask BEFORE the edge trim; the wrapper
%             turns it into results.mask via (mask == (regionsMask>0)).
%   imgVis    the enhanceForDisplay image (pre-normalisation) for the preview.
%
% DEPENDS ON
%   getEdgeSizeSLSCI, enhanceForDisplay, and MATLAB's Image Processing Toolbox
%   (ordfilt2, imdiffuseest, imdiffusefilt, imtophat, adaptthresh, imbinarize,
%   bwareaopen, imclose, imopen, imdilate, bwmorph, padarray, mat2gray).
%
% See also: runSegmentation, getSegmentationLabels, enhanceForDisplay, getEdgeSizeSLSCI
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 26-July-2026

function [cMask,edgeSize,mask,imgVis] = getPixelCategories(imgIni,regionsMask,existingMask,isK,s)

% --- edge size + per-modality preparation of the mean image ---
if isK
    edgeSize=getEdgeSizeSLSCI(imgIni,0.8);
    img=imgIni;
    img(img(:)>prctile(img(:),99))=prctile(img(:),99);
    img(img(:)<prctile(img(:),1))=prctile(img(:),1);
    img=imcomplement(img);
else
    edgeSize=0;
    img=imgIni;
end

fSize=floor((min(size(img))./20))*2+1;
img(isnan(img))=0;
img=enhanceForDisplay(img,fSize,min(15,fSize));
imgVis=img;

% --- trust mask: reuse an existing same-size mask, else contrast-limit (K only) ---
mask=ones(size(img));
if ~isempty(existingMask) && size(existingMask,1)==size(img,1) && size(existingMask,2)==size(img,2)
    mask=double(existingMask>0);
elseif isK
    mask=mask.*ordfilt2(imgIni,1,ones(s.iniSizeN),'symmetric')>=s.trustLimitsK(1) & ordfilt2(imgIni,s.iniSizeN.*s.iniSizeN,ones(s.iniSizeN),'symmetric')<=s.trustLimitsK(2);
end

% --- restrict to the selected regions; keep this mask for results.mask ---
mask=(regionsMask>0).*mask;
maskOut=mask;

tmp=zeros(size(img));
tmp(edgeSize+1:end-edgeSize,edgeSize+1:end-edgeSize)=mask(edgeSize+1:end-edgeSize,edgeSize+1:end-edgeSize);
mask=(tmp==1);

% --- per-region intensity normalisation ---
for i=1:1:max(regionsMask(:))
    img(regionsMask(:)==i)=mat2gray(img(regionsMask(:)==i),double(prctile(img(regionsMask(:)==i),[5,99])));
end
img=img.*mask+(1-mask);
img(isnan(img))=0;

% --- anisotropic diffusion (padded) ---
img=padarray(img,[s.lSizeN,s.lSizeN],'symmetric');
[gradThresh,numIter] = imdiffuseest(img,'ConductionMethod','quadratic');
if isK
    img = imdiffusefilt(img,'ConductionMethod','quadratic', 'GradientThreshold',[gradThresh,(9/10*min(gradThresh)):(-min(gradThresh)/10):(min(gradThresh)/10)],'NumberOfIterations',numIter+9);
else
    img = imdiffusefilt(img,'ConductionMethod','quadratic', 'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
end
pval = medfilt2(img, [s.lSizeN s.lSizeN], 'symmetric');
pval=min(pval);
mask=padarray(mask,[s.lSizeN,s.lSizeN],0);

% --- multiscale top-hat + adaptive threshold -> unsegmented vessel pixels ---
tmp=zeros(size(img));
for i=s.sSizeN:2:s.lSizeN
    tmp2=imtophat(img,strel('disk',i-s.sSizeN+s.deSens));
    tmp2=imbinarize(tmp2,adaptthresh(tmp2,s.sens,'NeighborhoodSize',i));
    tmp2(img<=pval)=0;
    tmp2=bwareaopen(tmp2,round(s.sSizeN.*s.sSizeN.*s.sSizeScale));
    tmp=tmp+tmp2;
end
cMask=int32(mask);
cMask(tmp(:)>0)=2;

% --- morphological cleanup and wall / lumen levels ---
maskV=tmp>0;
maskV=(maskV+(conv2(maskV,[1,1,1;1,0,1;1,1,1],'same')>4))>0; %connect pixels with more than 4 neighbours
maskV   = imclose(maskV,strel('disk',s.lThinN));
tmp(:)=1;
tmp(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN)=maskV(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN);
tmp   = imopen(tmp,strel('disk',s.imOpen));
maskV=bwareaopen(tmp,s.lSizeN.*s.lSizeN).*mask;
tmp(:)=0;
tmp(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN)=maskV(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN);
maskV=bwareaopen(tmp,s.lSizeN.*s.sSizeN);
tmp=maskV(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN);
maskV=padarray(tmp,[s.lSizeN,s.lSizeN],"symmetric");

maskVEE=maskV;
for i=1:1:s.eEdge
    maskVEE=imdilate(maskVEE,strel("disk",1));
end
maskVIE=maskV;
maskV=bwmorph(maskV,'thin',s.iEdge);

cMask(maskVEE(:)==1)=3;
cMask(maskVIE(:)==1)=4;
cMask(maskV(:)==1)=5;
tmp=bwareaopen(cMask==2,s.sSizeN.*s.sSizeN);
cMask(cMask(:)==2 & ~tmp(:))=1;

cMask=cMask.*int32(mask);
cMask=cMask(s.lSizeN+1:end-s.lSizeN,s.lSizeN+1:end-s.lSizeN);

mask=maskOut;
end
