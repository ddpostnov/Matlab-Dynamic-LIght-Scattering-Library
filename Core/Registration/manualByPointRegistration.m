%manualByPointRegistration - Interactive landmark-based image registration.
%
%   Lets the user pick any number of matching anchor-point pairs on a reference
%   image and a moving image (zoom/pan, press Enter to place each point), then
%   estimates a projective transform once more than three pairs are available and
%   shows the running overlay.  Returns the transform and the warped image.
%
% Syntax:
%    [tform,imgOut] = manualByPointRegistration(refImg, img, visType)
%
% Inputs:
%    refImg  - reference (fixed) image, 2-D.
%    img     - moving image to register, 2-D.
%    visType - 'sideBySide' for a 3-panel view, otherwise overlay only.
%
% Outputs:
%    tform  - estimated projective transform (moving -> reference).
%    imgOut - the moving image warped into the reference frame.
%
% Usage:
%    Zoom/pan to a landmark, press Enter, then click it on the reference; repeat
%    on the moving image.  Close the figure to finish.
%
% Dependencies: Image Processing Toolbox (estgeotform2d, imwarp, imshowpair).
% See also: translateMask, runRegistration
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
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

while true
    zoom on
    figure(h);
    waitfor(gcf, 'CurrentCharacter', char(13))
    if exist('h','var') && ishandle(h)
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




