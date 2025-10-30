% getTauc - uses thresholding for simplified evaluation of decorrelation
% time
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

function [tauC, decorThreshold]=getTauc(data,threshold,fPar,lagMax,fps)
decorThreshold=zeros(size(data,1),size(data,2));
if isempty(threshold)
    frame=imgaussfilt(squeeze(data(:,:,end)),fPar);    
    decorThreshold(:)=max(frame(:));
elseif threshold==-1
    frame=squeeze(data(:,:,1))-squeeze(min(data(:,:,:),[],3));
decorThreshold=frame./(2.718*2.718)+squeeze(min(data(:,:,:),[],3));
elseif threshold==-2    
    frame=imgaussfilt(squeeze(data(:,:,1)),fPar)-imgaussfilt(squeeze(data(:,:,end)),fPar);
    decorThreshold=frame./(2.718*2.718)+imgaussfilt(squeeze(data(:,:,end)),fPar);
else
    decorThreshold(:)=threshold;
end

tauC=zeros(size(data,1),size(data,2));
xq=0:1:lagMax;
for x=1:1:size(data,1)
    for y=1:1:size(data,2)
        ts=squeeze(data(x,y,:));
        idx=find(ts<=decorThreshold(x,y),1);
        if ~isempty(idx)
            if idx>1
                yR=ts(idx);
                yL=ts(idx-1);
                cR=(decorThreshold(x,y)-yL)/(yR-yL);
                cL=(yR - decorThreshold(x,y))/(yR-yL);
                tauC(x,y)=(xq(idx-1)*cL+xq(idx)*cR)./fps;
            else
                tauC(x,y)=xq(1)./fps;
            end
        else
            tauC(x,y)=xq(end)./fps;
        end
    end
end

end