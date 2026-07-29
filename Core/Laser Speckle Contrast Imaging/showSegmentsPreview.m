%showSegmentsPreview - save the _rep_segments report page (optionally with dvs overlay)
%
%   showSegmentsPreview(rep,fName,dataCube,cMask,sMap,isK)
%   showSegmentsPreview(rep,fName,dataCube,cMask,sMap,isK,dvsMap)
%
% DESCRIPTION
%   Writes <name>_rep_segments.jpg: a two-tile preview of a segmented recording.
%   The left tile is the enhanced mean image (percentile-scaled, inverted for _K /
%   min-offset for _I, background-subtracted by enhanceForDisplay, then masked to
%   cMask>0); the right tile is the golden-angle colour-mapped segment label map
%   sMap.  This is the segments preview shared by the static runSegmentation and the
%   dynamic runDynamicSegmentation - lifted verbatim from the runSegmentation local
%   it used to be, mirroring how the enhancement idiom became enhanceForDisplay.
%
%   IT TAKES THE REPORTING CONTEXT because the page has to land in the same ledger
%   as every other report image.  This function already drew invisibly, on white,
%   and deleted its figure - it is the pattern Core/Reporting was generalised from -
%   but it wrote the file itself, so the segmentation PDF was silently losing half
%   its pages.  reportFigure/reportSave now own the canvas, the name and the ledger;
%   fName is passed explicitly so a copy target gets its own page (see
%   runSegmentation's copySegmentationOnto).
%
%   When dvsMap is supplied (dynamic segmentation), the accepted-segment boundaries
%   (visboundaries(dvsMap>0)) are overlaid on the left tile, reproducing the
%   pre-refactor dynamic preview.  Omit dvsMap (or pass []) for the static one.
%   runDynamicSegmentation re-emits the page (overwriting the static one that
%   runSegmentation wrote) so the final on-disk artifact matches the old attmemptDS
%   pipeline exactly.
%
% INPUT
%   rep       the reporting context from reportOpen (the caller's).
%   fName     the *_d.mat path; the page is written next to it.
%   dataCube  source data cube; mean over dim 3 builds the background image.
%   cMask     merged category mask; cMask>0 masks the background image.
%   sMap      indexed segment label map (right tile; 0 = background).
%   isK       true for contrast (_K) data (imcomplement), false for intensity (_I).
%
% OPTIONAL
%   dvsMap    accepted dynamic-segment label map; its nonzero boundaries are drawn
%             on the left tile.  Default [] (static preview, no overlay).
%
% OUTPUT
%   none - writes <name>_rep_segments.jpg as a side-effect and records it in rep's
%   ledger, so reportClose can page it into the call's PDF.
%
% DEPENDS ON
%   reportFigure, reportSave, enhanceForDisplay, and MATLAB's Image Processing
%   Toolbox (visboundaries, hsv2rgb, mat2gray, imcomplement).
%
% See also: runSegmentation, runDynamicSegmentation, enhanceForDisplay,
%           reportFigure, reportSave
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

function showSegmentsPreview(rep,fName,dataCube,cMask,sMap,isK,dvsMap)

if nargin<7, dvsMap=[]; end

img=mean(dataCube,3);
img=mat2gray(img,double(prctile(img(cMask(:)>0),[5,99])));
if isK
    img=imcomplement(img);
else
    img=img-min(img(cMask(:)>0 & img(:)>0));
end

fSize=floor((min(size(img))./20))*2+1;
img=enhanceForDisplay(img,fSize).*(cMask>0);

n=double(max(sMap(:)));
phi=(sqrt(5)-1)/2;
cmap=hsv2rgb([mod((0:n-1)'*phi,1) 0.8*ones(n,1) 0.9-0.2*mod((0:n-1)',2)]);

f=reportFigure(rep,'segments');
try
    t=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    t1=nexttile(t);
    imagesc(img,'Parent', t1)
    axis image
    if ~isempty(dvsMap)
        hold(t1,'on')
        visboundaries(t1,dvsMap>0);
        hold(t1,'off')
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
catch ME
    delete(f); rethrow(ME);      % reportSave never runs, so the figure goes here
end
reportSave(rep,f,'segments',fName);
end
