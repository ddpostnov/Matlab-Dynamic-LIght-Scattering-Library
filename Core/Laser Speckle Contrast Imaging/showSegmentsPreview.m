%showSegmentsPreview - draw the segments report page into a figure
%
%   showSegmentsPreview(fh,dataCube,cMask,sMap,polarity)
%   showSegmentsPreview(fh,dataCube,cMask,sMap,polarity,dvsMap)
%
% DESCRIPTION
%   Draws a two-tile preview of a segmented recording into the figure it is given.
%   The left tile is the enhanced mean image (percentile-scaled, then brought to
%   vessels-bright-on-a-zero-floor either way round - inverted when the vessels are
%   dark, floor-subtracted when they are already bright - background-subtracted by
%   enhanceForDisplay, then masked to cMask>0); the right tile is the
%   golden-angle colour-mapped segment label map
%   sMap.  This is the segments preview shared by the static runSegmentation and the
%   dynamic runDynamicSegmentation - lifted verbatim from the runSegmentation local
%   it used to be, mirroring how the enhancement idiom became enhanceForDisplay.
%
%   IT TAKES A FIGURE HANDLE, NOT THE REPORTING CONTEXT.  Only wrappers report: a
%   core that contributes a page is handed a canvas and draws on it, and the
%   wrapper that opened the canvas is the one that names, writes and ledgers the
%   image.  So the caller does
%
%       fh=reportFigure(rep,'segments');
%       showSegmentsPreview(fh,...);
%       reportSave(rep,fh,'segments',fName);
%
%   and this function knows nothing about files, tags or PDFs.  It does not delete
%   the figure either - reportSave does that on every path, including the one
%   where the drawing below throws.
%
%   When dvsMap is supplied (dynamic segmentation), the accepted-segment boundaries
%   (visboundaries(dvsMap>0)) are overlaid on the left tile, reproducing the
%   pre-refactor dynamic preview.  Omit dvsMap (or pass []) for the static one.
%   runDynamicSegmentation re-emits the page (overwriting the static one that
%   runSegmentation wrote) so the final on-disk artifact matches the old attmemptDS
%   pipeline exactly.
%
% INPUT
%   fh        the figure to draw into, from reportFigure.
%   dataCube  source data cube; mean over dim 3 builds the background image.
%   cMask     merged category mask; cMask>0 masks the background image.
%   sMap      indexed segment label map (right tile; 0 = background).
%   polarity  'bright' or 'dark' - are the vessels brighter or darker than the
%             tissue, from getVesselPolarity.  Both arms end with the vessels
%             bright against a zero floor, which is what the colour scale assumes.
%
% OPTIONAL
%   dvsMap    accepted dynamic-segment label map; its nonzero boundaries are drawn
%             on the left tile.  Default [] (static preview, no overlay).
%
% OUTPUT
%   none - the tiles are added to fh, which the caller owns.
%
% DEPENDS ON
%   enhanceForDisplay and MATLAB's Image Processing Toolbox (visboundaries,
%   hsv2rgb, mat2gray, imcomplement).
%
% See also: getVesselPolarity, runSegmentation, runDynamicSegmentation,
%           enhanceForDisplay, reportFigure, reportSave
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

function showSegmentsPreview(fh,dataCube,cMask,sMap,polarity,dvsMap)

if nargin<6, dvsMap=[]; end

img=mean(dataCube,3);
img=mat2gray(img,double(prctile(img(cMask(:)>0),[5,99])));
if strcmp(polarity,'dark')
    img=imcomplement(img);
else
    img=img-min(img(cMask(:)>0 & img(:)>0));
end

fSize=floor((min(size(img))./20))*2+1;
img=enhanceForDisplay(img,fSize).*(cMask>0);

n=double(max(sMap(:)));
phi=(sqrt(5)-1)/2;
cmap=hsv2rgb([mod((0:n-1)'*phi,1) 0.8*ones(n,1) 0.9-0.2*mod((0:n-1)',2)]);

t=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
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
end
