% getNormalizedG2 - calculates g2 - normalized intensity autocorrelation 
% function drom 3D (X,Y,T) data.
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

function g2=getNormalizedG2(data,lagMax,sampleL,blocksN)

sampleL=sampleL-lagMax-1;
sizeY=size(data,1);
sizeX=size(data,2);
blockSizeY=floor(sizeY/blocksN(1));
blockSizeX=floor(sizeX/blocksN(2));
X=[1:blockSizeX:(sizeX-blockSizeX+1),sizeX];
Y=[1:blockSizeY:(sizeY-blockSizeY+1),sizeY];

g2=gpuArray(zeros(size(data,1),size(data,2),lagMax+1,'single'));
for i=1:1:length(X)-1
    for j=1:1:length(Y)-1
        subdata=gpuArray(single(data(Y(j):Y(j+1),X(i):X(i+1),:)));
        for lag=0:1:lagMax
            g2(Y(j):Y(j+1),X(i):X(i+1),lag+1)=(mean(subdata(:,:,1:sampleL).*subdata(:,:,1+lag:sampleL+lag),3))./(mean(subdata(:,:,1:sampleL),3).^2);
        end
    end
end
g2=gather(g2);
end
