%getEllepticalROI - creates an elleptical ROI
%
% Syntax:  ROIs=getEllepticalROI(img)
%
% Inputs:
%    img  - image to select regions from
%
% Outputs:
%    ROIs - structure that contains trapezoid mask and bounding box
%           coordinates for the mask position in the original image
%
% Example:
%    ROIs=getTrapezoidSelection(img)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 17-Jun-2024

%------------- BEGIN CODE --------------

function ROIs=getEllepticalROI(h,img,roisN)
ii=[];
if ~isempty(img)
    h=figure;
    imagesc(img);
    axis image
    colormap gray
    imgSize=size(img);
end
imgSize=size(img);
figure(h);
disp('ROIs selection');
ii=1;
set(gcf, 'CurrentCharacter', char(12))
title('Press enter to add the first ROI')
waitfor(gcf, 'CurrentCharacter', char(13))
while true
    if exist('h') && ishandle(h)
        try
            roi = drawellipse('Color','r');
            
            set(h,'KeyPressFcn',@(H,E)deleteROI(H,E));
            set(gcf, 'CurrentCharacter', char(12))
            title('Press enter to confirm ROI. Press esc to finish')
            waitfor(gcf, 'CurrentCharacter', char(13))
            ROIs(ii).mask=createMask(roi);
            ii=ii+1;
            if ~isempty(roisN) && ii==roisN+1
                delete(h);
            end
        catch
        end
    else
        break;
    end
end
    function deleteROI(H,E)
        switch E.Key           
            case 'escape'
                delete(h);
            otherwise
        end
    end

end