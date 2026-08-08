%getPixelCategories - binarize an LSCI mean image into a 5-level category mask
%
%   [cMask,edgeSize,mask,imgVis] = getPixelCategories(imgIni,regionsMask, ...
%                                                     existingMask,product,polarity,s)
%
% DESCRIPTION
%   The automatic categorization core lifted verbatim from runCategories: it turns
%   a mean LSCI image into the five-level categorical mask cMask
%     0 background | 1 parenchyma | 2 unsegmented | 3 external wall |
%     4 internal wall | 5 lumen
%   The pipeline is: edge size -> polarity prep -> enhanceForDisplay -> trust mask ->
%   per-region normalisation -> anisotropic diffusion -> multiscale top-hat +
%   adaptive threshold -> morphological cleanup.  Only the automatic math lives
%   here; the file I/O and the interactive region GUI stay in the wrapper, which
%   supplies regionsMask.
%
%   FOUR DECISIONS, AND THEY DO NOT MOVE TOGETHER.  Until 2026-08-08 all four were
%   taken by one boolean called isK, which meant "this came from a _K file" - a
%   single name for four unrelated questions that were only ever asked of the same
%   files.  The fluorescence branch broke that coincidence, so each one now says
%   what it follows:
%
%     inversion       follows POLARITY.  In speckle contrast fast flow is dark, so
%                     the image is inverted to look for vessels as bright ridges.
%                     A plasma-labelled fluorescence recording is already the right
%                     way round and inverting it would segment the tissue BETWEEN
%                     the vessels.  Which way round an intensity recording is, is
%                     the operator's answer, not the file name's - see
%                     getVesselPolarity.
%     edge size       follows the PRODUCT TOKEN.  getEdgeSizeSLSCI measures the
%                     border a spatial-contrast kernel leaves round the frame; a
%                     fluorescence frame has no such border whichever way round its
%                     vessels are, so this is 0 for '_I' and stays that way.
%     trust mask      follows the PRODUCT TOKEN.  s.trustLimitsK bounds a CONTRAST,
%                     and there is no contrast on an intensity product to bound -
%                     the intensity twin does not even carry the setting.  An
%                     intensity product's trust mask is left fully open unless the
%                     caller hands one in as existingMask.
%     diffusion       follows s.diffusionSchedule, and is the one of the four that
%                     is a TUNING choice rather than a fact about the file.  The
%                     multi-scale schedule was written for speckle contrast's noise
%                     and may well suit a fluorescence image too; nobody has
%                     measured it, so it is a setting with per-step defaults that
%                     reproduce what each branch did before, and moving it is a
%                     measurement somebody can now make without a code change.
%
%   The trust mask is automatic (not interactive): an existing same-size mask
%   (existingMask, typically results.mask from runContrastFromRLS) is reused as-is,
%   otherwise - for contrast ('K') data - it is derived from the ordfilt2 min/max
%   of imgIni against s.trustLimitsK.  The mean image imgIni (NOT the enhanced
%   image) drives that ordfilt2, so it is passed alongside the enhancement.
%
% INPUT
%   imgIni        2-D mean image (results.imgK for contrast, results.imgI for
%                 intensity).  Used both for edge-size / enhancement and, for 'K',
%                 for the trust-mask ordfilt2.
%   regionsMask   integer label image (0 = excluded).  1..N select regions to
%                 categorize and normalise; the wrapper builds it (GUI or whole
%                 window).  Its size sets the working frame.
%   existingMask  a mask to reuse verbatim (e.g. results.mask), or [] to derive
%                 one.  Reused only when it matches size(imgIni).
%   product       'K' for contrast, 'I' for intensity.  The wrapper reads it off
%                 the file name through getVesselPolarity; do NOT sniff the file
%                 name inside the core.
%   polarity      'bright' or 'dark' - are the vessels brighter or darker than the
%                 tissue.  'dark' clips the extreme percentiles and inverts.
%   s             parameter struct.  Reads: iniSizeN, trustLimitsK ('K' only),
%                 diffusionSchedule, lSizeN, sSizeN, sens, deSens, sSizeScale,
%                 lThinN, imOpen, iEdge, eEdge.
%
% OUTPUT
%   cMask     int32 five-level category mask, size(imgIni).
%   edgeSize  border width (px) getEdgeSizeSLSCI found (0 for 'I'); the scalar the
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
% See also: getVesselPolarity, runSegmentation, getSegmentationLabels,
%           enhanceForDisplay, getEdgeSizeSLSCI
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

function [cMask,edgeSize,mask,imgVis] = getPixelCategories(imgIni,regionsMask,existingMask,product,polarity,s)

isContrast=strcmp(product,'K');

% --- the SLSCI border, which only a contrast image has ---
if isContrast
    edgeSize=getEdgeSizeSLSCI(imgIni,0.8);
else
    edgeSize=0;
end

% --- polarity: everything after this looks for vessels as BRIGHT ridges ---
% The percentile clip belongs to the inversion rather than beside it: one hot pixel
% becomes an extreme low after imcomplement and drags the enhancement with it.
img=imgIni;
if strcmp(polarity,'dark')
    img(img(:)>prctile(img(:),99))=prctile(img(:),99);
    img(img(:)<prctile(img(:),1))=prctile(img(:),1);
    img=imcomplement(img);
end

fSize=floor((min(size(img))./20))*2+1;
img(isnan(img))=0;
img=enhanceForDisplay(img,fSize,min(15,fSize));
imgVis=img;

% --- trust mask: reuse an existing same-size mask, else contrast-limit ('K' only) ---
% s.trustLimitsK bounds a CONTRAST, so it is the product token that decides this and
% not the polarity: a dark-vessel intensity image has no contrast to bound either.
mask=ones(size(img));
if ~isempty(existingMask) && size(existingMask,1)==size(img,1) && size(existingMask,2)==size(img,2)
    mask=double(existingMask>0);
elseif isContrast
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
% A TUNING CHOICE, AND THE ONLY ONE OF THE FOUR THAT IS.  'multiscale' adds nine
% coarse-to-fine passes at falling gradient thresholds; it was written for speckle
% contrast, which is a ratio of two small-kernel statistics and therefore noisy.
% Whether a fluorescence mean image wants it too is a measurement nobody has made,
% so it is a setting and each step's preset reproduces what that branch did before.
if strcmp(s.diffusionSchedule,'multiscale')
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
