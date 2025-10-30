function translateMask(H,E,offset,name) 
switch E.Key
    case 'downarrow'
        offset.y=offset.y+1;
    case 'uparrow'
        offset.y=offset.y-1;
    case 'rightarrow'
        offset.x=offset.x+1;
    case 'leftarrow'
        offset.x=offset.x-1;
    otherwise
end
assignin('base', name, offset);
end