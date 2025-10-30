% getSpeckelSize - estimates specklesize as width at half prominance of
% the average spatial autocorrelation
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 29 August 2022

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