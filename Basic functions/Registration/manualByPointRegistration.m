% manualByPointRegistration - performs point by point registration allowing
% user to select any number of anchor points
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

function [tform,imgOut]=manualByPointRegistration(refImg,img,visType)

tform=affine2d(eye(3));
X1=[];
Y1=[];
X2=[];
Y2=[];
imgOut=img;
h=figure(1);
h.WindowState="maximized";
if strcmp(visType,'sideBySide')
    subplot(1,3,1)
    imagesc(refImg)
    axis image
    title('Reference image')
    subplot(1,3,2)
    imagesc(img)
    axis image
    title('Image to register')
    subplot(1,3,3)

    imshowpair(refImg,img)
    title('Overlay');
else
    imshowpair(refImg,img)
end

set(zoom(h),'ActionPostCallback',@(h,evd)adjustColorLimits(h));
set(pan(h),'ActionPostCallback',@(h,evd)adjustColorLimits(h));

while true
    zoom on
    figure(h);
    waitfor(gcf, 'CurrentCharacter', char(13))
    if exist('h') && ishandle(h)
        zoom off
        set(gcf, 'CurrentCharacter', char(12))
        [x,y]=ginput(1);
        x=round(x);
        y=round(y);
        X1=[X1,x];
        Y1=[Y1,y];
        hold on
        plot(X1,Y1,'xm','MarkerSize',20)
        hold off

        zoom on
        waitfor(gcf, 'CurrentCharacter', char(13))
        zoom off
        set(gcf, 'CurrentCharacter', char(12))
        [x,y]=ginput(1);
        x=round(x);
        y=round(y);
        X2=[X2,x];
        Y2=[Y2,y];
        hold on
        plot(X2,Y2,'xg','MarkerSize',20)
        hold off


        if length(X2)>3
            img=imgOut;
            img(isnan(img))=0;
            tform=estgeotform2d([X2',Y2'],[X1',Y1'],"projective");
            img=imwarp(img,tform,"OutputView",imref2d(size(refImg)),'FillValues', 0);
            subplot(1,3,3)
            imshowpair(refImg,img)
            title(['Overlay, \Delta=',num2str(round(sum(abs(refImg(:)-img(:)))))]);
        end
    else
        break;
    end
end
imgOut=img;
end




