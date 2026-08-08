%meGlobal - Take the movement of the whole preparation out of the phase.
%
%   THE DIFFERENTIAL MODE, AND IT IS DONE IN THE PHASE DOMAIN RATHER THAN BY
%   REGISTERING FRAMES.  A rigid translation d of the whole field appears at pyramid
%   level k as a phase that is the SAME EVERYWHERE, roughly omega_k*d.  Subtracting
%   the spatial mean of the phase therefore removes it exactly - and, the part that
%   matters here, RESAMPLES NOTHING.  Registering each frame would interpolate the
%   stack, and interpolation is precisely what destroys a motion of a hundredth of a
%   pixel: the thing we came to measure would be spent paying for the correction.
%
%   THE MEAN IS AMPLITUDE-WEIGHTED because most of the frame is flat tissue where
%   the band is near zero and its phase is the direction of the noise.  An unweighted
%   mean would be dominated by those pixels and would subtract noise from signal.
%
%   ORDER 0 IS A TRANSLATION.  ORDER 1 ADDS A PLANE - scale, rotation and parallax
%   read as a phase that varies linearly across the field, so a weighted plane fit
%   removes them too.  Order 1 also removes any genuine motion that happens to be
%   linear across the frame, which on a field of a few hundred micrometres a vessel
%   can be; it is a setting, not the default, for that reason.
%
%   ORDER -1 LEAVES EVERYTHING IN.  That is the absolute mode, and it is not a
%   fallback: the DIFFERENCE between the absolute and the differential output is a
%   measurement in its own right, because it says how much of the in-band motion was
%   the preparation rather than the vessel.  On this recording that matters - the
%   author sees an artery move with the naked eye and the probe measured the field
%   itself moving 0.078 px between frames, which is larger than any wall motion in
%   it.
%
% Syntax:
%    [pc, ps] = meGlobal(pc, ps, amp, order)
%
% Inputs:
%    pc, ps - the phase, [rows columns frames], as mePhase returns it and meTemporal
%             has filtered it.
%    amp    - the weight, same size, from mePhase.
%    order  - -1 leaves the bulk motion in, 0 removes a rigid translation, 1 also
%             removes scale, rotation and parallax.
%
% Outputs:
%    pc, ps - the same, with the bulk component removed frame by frame.
%
% Example:
%    [pc,ps] = meGlobal(pc, ps, ph.amp{3}, 0);
%
% Dependencies: none beyond core MATLAB.
% See also: mePhase, meTemporal, meShift, meBlur
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function [pc, ps] = meGlobal(pc, ps, amp, order)

if order < 0, return; end

sw = sum(amp,[1 2]);
sw(sw==0) = 1;

if order == 0
    pc = pc - sum(amp.*pc,[1 2])./sw;
    ps = ps - sum(amp.*ps,[1 2])./sw;
    return
end
if order ~= 1
    error('meGlobal:order','The order is -1, 0 or 1, not %g.',order);
end

% A weighted plane, frame by frame.  The coordinates are centred so the 3x3 normal
% matrix stays well conditioned on a large level.
[nr,nc,nT] = size(pc);
[xx,yy] = meshgrid(1:nc, 1:nr);
xx = cast(xx - (nc+1)/2,'like',pc);
yy = cast(yy - (nr+1)/2,'like',pc);

s00 = sw;
s10 = sum(amp.*xx ,[1 2]);
s01 = sum(amp.*yy ,[1 2]);
s20 = sum(amp.*xx.^2,[1 2]);
s11 = sum(amp.*xx.*yy,[1 2]);
s02 = sum(amp.*yy.^2,[1 2]);

bc = cat(4, sum(amp.*pc,[1 2]), sum(amp.*xx.*pc,[1 2]), sum(amp.*yy.*pc,[1 2]));
bs = cat(4, sum(amp.*ps,[1 2]), sum(amp.*xx.*ps,[1 2]), sum(amp.*yy.*ps,[1 2]));

for t = 1:nT
    M = double([s00(t) s10(t) s01(t); s10(t) s20(t) s11(t); s01(t) s11(t) s02(t)]);
    if rcond(M) < 1e-12, continue; end
    a = M\double(squeeze(bc(1,1,t,:)));
    b = M\double(squeeze(bs(1,1,t,:)));
    pc(:,:,t) = pc(:,:,t) - cast(a(1) + a(2).*xx + a(3).*yy,'like',pc);
    ps(:,:,t) = ps(:,:,t) - cast(b(1) + b(2).*xx + b(3).*yy,'like',ps);
end
end
%------------- END OF CODE --------------
