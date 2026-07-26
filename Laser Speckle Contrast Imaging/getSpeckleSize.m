%getSpeckleSize - Estimate speckle size from the average spatial autocorrelation.
%
%   Estimates the speckle size as the width at half-prominence of the mean 1-D
%   spatial autocorrelation of an image, averaged over its rows and columns.
%
% Syntax:
%    sSize = getSpeckleSize(img, maxLags, toPlot)
%
% Inputs:
%    img     - a single 2-D image (e.g. one raw speckle frame).
%    maxLags - maximum autocorrelation lag, pixels.
%    toPlot  - logical; when true, plots the image and the autocorrelation peak.
%
% Outputs:
%    sSize - estimated speckle size (full width at half prominence), pixels.
%
% Example:
%    s = getSpeckleSize(rawFrame, 15, true);
%
% Dependencies: Signal Processing Toolbox (xcorr, findpeaks).
% See also: getK, getContrastFromRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function sSize=getSpeckleSize(img,maxLags,toPlot)
counter=1;
for i=maxLags+1:size(img,1)-maxLags
x=squeeze(img(i,:));
x=x-mean(x);
[r(counter,:),lags]=xcorr(x,maxLags,'unbiased');
counter=counter+1;
end
for i=maxLags+1:size(img,2)-maxLags
x=squeeze(img(:,i));
x=x-mean(x);
[r(counter,:),lags]=xcorr(x,maxLags,'unbiased');
counter=counter+1;
end
rMean=mean(r,1);
rMean=rMean/max(rMean);
[~,locs,w,~]=findpeaks(rMean);
sSize=w(locs==maxLags+1);
if toPlot
subplot(1,2,1)
imagesc(img)
title(std(single(img(:)))./mean(single(img(:))))
subplot(1,2,2)
findpeaks(rMean,lags,'Annotate','extents')
title(sSize)
end
end