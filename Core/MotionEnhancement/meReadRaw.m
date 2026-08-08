%meReadRaw - Open one of the working copies written by mePassA.
%
%   THE ONE READER OF THE FLAT FILES.  mePassA writes two of them and every later
%   step reads them, so the sample size, the array shape and the mapping back to
%   the full-resolution frame are read from the sidecar in one place rather than
%   repeated at each call site.  A raw file carries none of that itself: get the
%   shape wrong and it reads as a plausible picture of nothing.
%
%   NOTHING IS LOADED.  A memmapfile is returned, so a step touches only the
%   frames it asks for.  The binned copy of the reference recording is 8 GB and
%   the full-resolution window is another 8; a step that wants two hundred frames
%   of one window should cost two hundred frames of one window.
%
%   THE COPIES DO NOT SHARE A COORDINATE SYSTEM, so geom carries the map.  A
%   position in either copy becomes a position in the original frame through
%   geom.origin and geom.scale, CENTRE TO CENTRE:
%
%       fullRow = geom.origin(1) + (row-1)*geom.scale + (geom.scale-1)/2
%
%   The last term is the half-pixel a binned sample sits at, and leaving it out
%   moves every binned measurement half a pixel against the full-resolution one -
%   which in a block that resolves tenths of a pixel is not a rounding detail.
%   geom.pixelSize is already multiplied out, so a distance measured in either
%   copy is micrometres without any further thought about which copy it came from.
%
% Syntax:
%    [m,geom] = meReadRaw(pass, copy)
%
% Inputs:
%    pass - the sidecar struct from mePassA, or the path to its .mat file.
%    copy - 'bin' for the binned copy of the whole field, 'crop' for the
%           full-resolution copy of one window.
%
% Outputs:
%    m    - a memmapfile whose m.Data.d is [rows columns frames].
%    geom - .size      [rows columns frames]
%           .class     sample type of the copy
%           .origin    [row column] of its first pixel in the original frame
%           .scale     original pixels per pixel of this copy
%           .pixelSize micrometres per pixel of this copy
%           .fs        frames per second
%           .time      seconds, one per frame
%           .lo,.hi    the intensity limits the samples were scaled between, so a
%                      sample can be turned back into counts
%
% Example:
%    pass      = load('rec_ME_passA.mat').pass;
%    [m,geom]  = meReadRaw(pass,'crop');
%    frame     = m.Data.d(:,:,1500);
%
% Dependencies: none.
% See also: mePassA, meSettings, meProbe, memmapfile
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function [m,geom] = meReadRaw(pass, copy)

if ischar(pass) || isstring(pass)
    loaded = load(char(pass),'pass');
    pass   = loaded.pass;
end
copy = lower(char(copy));

switch copy
    case 'bin'
        fPath  = pass.binPath;
        sz     = [pass.binSize pass.sizeT];
        origin = [1 1];
        scale  = pass.binFactor;
    case 'crop'
        fPath  = pass.cropPath;
        sz     = [pass.cropSize pass.sizeT];
        origin = [pass.roi(1,1) pass.roi(2,1)];
        scale  = 1;
    otherwise
        error('meReadRaw:copy','Ask for the binned copy or the cropped one, not %s.',copy);
end

if ~isfile(fPath)
    error('meReadRaw:missing','The working copy is not there: %s',fPath);
end

% A wrong shape reads as a plausible picture of nothing, so the file's own size
% is checked against the shape the sidecar claims before anything is mapped.
nByte  = numel(typecast(cast(0,pass.rawClass),'uint8'));
want   = prod(sz)*nByte;
onDisk = dir(fPath);
if onDisk.bytes ~= want
    error('meReadRaw:size', ...
        'The working copy is %d bytes and the sidecar describes %d. Re-run mePassA.', ...
        onDisk.bytes, want);
end

m = memmapfile(fPath,'Format',{pass.rawClass, sz, 'd'},'Writable',false);

geom           = struct();
geom.size      = sz;
geom.class     = pass.rawClass;
geom.origin    = origin;
geom.scale     = scale;
geom.pixelSize = pass.pixelSize*scale;
geom.fs        = pass.fs;
geom.time      = (0:pass.sizeT-1)'./pass.fs;
geom.lo        = pass.lo;
geom.hi        = pass.hi;
end
%------------- END OF CODE --------------
