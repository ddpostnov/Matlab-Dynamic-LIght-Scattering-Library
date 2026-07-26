%getEdgeSizeSLSCI - Estimate the edge-artefact border width of a contrast image.
%
%   Spatial-contrast images carry a border of biased values because the
%   contrast kernel runs past the image edge.  getEdgeSizeSLSCI grows a border
%   inward while more than a fraction 'thresh' of the outer pixels are still
%   brighter than the pixels just inside them (the signature of the edge
%   artefact) and returns the width, in pixels, of that biased border so it can
%   be cropped before further analysis.
%
% Syntax:
%    edgeSize = getEdgeSizeSLSCI(imgK, thresh)
%
% Inputs:
%    imgK    - a single spatial-contrast image, 2-D [y x].
%    thresh  - fraction (0-1) of edge pixels that must exceed their inner
%              neighbour for the border row/column to count as still affected.
%
% Outputs:
%    edgeSize - width, in pixels, of the biased border (0 = none detected).
%
% Example:
%    e   = getEdgeSizeSLSCI(getK(raw,'spatial','kernelSize',5,'procType','gpu'), 0.8);
%    img = imgK(e+1:end-e, e+1:end-e);          % crop the biased border
%
% See also: getK, getContrastFromRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function edgeSize = getEdgeSizeSLSCI(imgK,thresh)

edgeSize = 1;
while (sum(imgK(:,edgeSize)>imgK(:,edgeSize+1))/size(imgK,1))>thresh || ...
      (sum(imgK(edgeSize,:)>imgK(edgeSize+1,:))/size(imgK,2))>thresh || ...
      (sum(imgK(:,end-edgeSize+1)>imgK(:,end-edgeSize))/size(imgK,1))>thresh || ...
      (sum(imgK(end-edgeSize+1,:)>imgK(end-edgeSize,:))/size(imgK,2))>thresh
    edgeSize = edgeSize+1;
end
edgeSize = edgeSize-1;
end
%------------- END OF CODE --------------
