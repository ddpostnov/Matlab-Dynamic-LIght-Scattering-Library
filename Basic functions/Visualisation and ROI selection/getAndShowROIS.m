% getAndShowROIS - allwos user to select ROI, returns the mask. If ROIS are
% passed as an argument - plots them instead. Requires figure handle to be
% passed as an argument
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 12 Feb 2023

function ROIs=getAndShowROIS(h,colors,N,ROIs)
figure(h);
if isempty(ROIs)
    count=1;
    ROIs=[];
    set(gcf, 'CurrentCharacter', char(12))
    title('Press enter to add the first ROI')
    waitfor(gcf, 'CurrentCharacter', char(13))
    while true
        if exist('h') && ishandle(h)
            try
                title('Use mouse to select ROI')
                BW = roipoly;
                s=regionprops(BW,'Orientation','PixelList');
                if ~isempty(N) && length(s.PixelList)>N
                    if abs(s.Orientation)<45
                        s.PixelList=sortrows(s.PixelList,1);
                        for i=N+1:1:length(s.PixelList)
                            BW(s.PixelList(i,2),s.PixelList(i,1))=0;
                        end
                    else
                        s.PixelList=sortrows(s.PixelList,2);
                        for i=N+1:1:length(s.PixelList)
                            BW(s.PixelList(i,2),s.PixelList(i,1))=0;
                        end

                    end
                end
                if isempty(N) || length(s.PixelList)>N
                    hold on
                    if isempty(colors)
                        h1(count)=visboundaries(BW);
                    else
                        h1(count)=visboundaries(BW,'Color',colors(count,:));
                    end
                    hold off
                    ROIs(:,:,count)=BW;
                    set(h,'KeyPressFcn',@(H,E)deleteROI(H,E));
                    set(gcf, 'CurrentCharacter', char(12))
                    title('Press enter or delete to add or remove ROI. Press esc to finish')
                    waitfor(gcf, 'CurrentCharacter', char(13))
                    count=count+1;
                end
            catch
            end
        else
            break;
        end
    end
else
    for i=1:1:size(ROIs,3)
        BW = squeeze(ROIs(:,:,i));
        hold on
        if isempty(colors)
            visboundaries(BW);
        else
            visboundaries(BW,'Color',colors(i,:));
        end
        hold off
    end
end

    function deleteROI(H,E)
        switch E.Key
            case 'delete'
                if count>0
                    delete(h1(count))
                    ROIs=ROIs(:,:,1:count-1);
                    count=count-1;
                end
            case 'escape'
                delete(h);
            otherwise
        end
    end
end