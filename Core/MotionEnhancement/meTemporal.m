%meTemporal - Keep the part of the phase that moves at the frequency of interest.
%
%   ZERO-PHASE, ALWAYS.  THIS IS THE RULE WITH THE LARGEST CONSEQUENCE IN THE BLOCK
%   AND IT IS NOT A PREFERENCE.  A causal filter puts a frequency-dependent lag on
%   whatever passes through it.  The local frequency content is not the same at
%   every point of a vascular field, so the lag is not the same either, and the
%   magnified movie then shows a wave crossing the vasculature THAT DOES NOT EXIST.
%   This library measures pulse-wave propagation - setVascularTree builds the flow
%   graph from pulse-arrival phase - so an invented travelling wave is not an
%   artefact one squints past, it is a wrong scientific result that looks like the
%   right one.  The published pseudocode uses a first-order IIR because it runs
%   online.  We are offline and have no such excuse.
%
%   FOUR FILTERS, ONE PROPERTY.  Every branch below is symmetric in time, either by
%   construction (a real, even frequency mask; a symmetric kernel) or by running the
%   recursion both ways (filtfilt).  Anything added here must be too.
%
%   'ideal'        an FFT band-pass.  The sharpest edge available, which is what a
%                  narrow cardiac band wants, and the default.  The series is
%                  mirrored at both ends before the transform so the wrap-around
%                  does not put a step at the joint and ring across the whole
%                  record.
%   'cyclic'       the same FFT band-pass with NO padding, for a series that is one
%                  whole period - the cine.  There the wrap-around is the physics
%                  and it is the mirror that would be the artefact.
%   'butter'       designfilt plus filtfilt.  Gentler skirts, no ringing.
%   'acceleration' the second derivative of a Gaussian in time, Zhang et al.,
%                  CVPR 2017.  A LARGE motion is locally linear in time and a
%                  second derivative annihilates it, while a small fast one
%                  survives.  Our data have exactly that pathology - vasomotion,
%                  respiration and stage drift underneath the pulsation.
%
%   THE ACCELERATION BRANCH IS NORMALISED TO UNIT GAIN AT THE BAND CENTRE, and that
%   is not a fudge factor on the magnifier.  A second derivative multiplies by
%   -omega^2, so without it the branch would report a gain of its own and the
%   calibration would come out at something other than 1+alpha for a reason that has
%   nothing to do with the magnifier.  The other two branches already pass their
%   band at unity; this makes the three comparable.  It also fixes the sign the
%   second derivative flips.
%
% Syntax:
%    y = meTemporal(x, fs, s)
%
% Inputs:
%    x  - [rows columns frames], filtered along the third dimension.  Any real
%         class; single in and single out.
%    fs - frames per second.
%    s  - settings.  Fields read:
%         .temporalFilter 'ideal' | 'cyclic' | 'butter' | 'acceleration'
%         .passband       [low high] in the units of the series' own axis.
%         .filterOrder    half the IIR order used by 'butter'.
%
% Outputs:
%    y  - the same size and class as x.
%
% Example:
%    s = meSettings;  s.passband = [7.9 14.7];
%    y = meTemporal(ph.cos{3}, 99.96, s);
%
% Dependencies: Signal Processing Toolbox (designfilt, filtfilt); Image Processing
%               Toolbox (imfilter).
% See also: mePhase, meShift, meGlobal, meSettings
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function y = meTemporal(x, fs, s)

nT = size(x,3);
if nT < 4
    error('meTemporal:frames','A temporal filter needs more than %d frames.',nT);
end
band = double(s.passband(:)');
band = [max(band(1), 1e-6), min(band(2), 0.499*fs)];
if band(2) <= band(1)
    error('meTemporal:band','The passband [%g %g] Hz does not fit under %g Hz.', ...
        s.passband(1), s.passband(2), fs/2);
end

switch lower(s.temporalFilter)

    case 'ideal'
        % Mirrored at both ends, endpoints not repeated, by two periods of the
        % lowest frequency that is being kept.  Enough that the brick wall's
        % ringing has decayed before it reaches the real data.
        nPad = min(nT-1, max(64, ceil(2*fs/band(1))));
        xp   = cat(3, x(:,:,nPad+1:-1:2), x, x(:,:,end-1:-1:end-nPad));
        M    = size(xp,3);

        f    = (0:M-1)./M.*fs;
        keep = (f>=band(1) & f<=band(2)) | (f>=fs-band(2) & f<=fs-band(1));
        H    = zeros(1,1,M,'like',x);
        H(1,1,keep) = 1;

        y = real(ifft(fft(xp,[],3).*H,[],3));
        y = y(:,:,nPad+1:nPad+nT);

    case 'cyclic'
        % THE SERIES IS ONE WHOLE PERIOD, SO THE WRAP-AROUND IS THE PHYSICS.  A cine
        % is a cardiac cycle sampled at nBin phase points and its last bin is
        % adjacent to its first; mirroring it, as 'ideal' does, would make it EVEN
        % about both ends, which throws away the odd part of the beat - the systolic
        % upstroke against the slower decay, which is the shape being measured.  So
        % nothing is padded, the axis is phase rather than time, and fs is bins per
        % cycle while the band is in cycles per cycle: [0.5 5.5] keeps the
        % fundamental and four harmonics and drops the static picture.
        f    = (0:nT-1)./nT.*fs;
        keep = (f>=band(1) & f<=band(2)) | (f>=fs-band(2) & f<=fs-band(1));
        H    = zeros(1,1,nT,'like',x);
        H(1,1,keep) = 1;
        y = real(ifft(fft(x,[],3).*H,[],3));

    case 'butter'
        d = designfilt('bandpassiir', ...
            'FilterOrder',        2*s.filterOrder, ...
            'HalfPowerFrequency1',band(1), ...
            'HalfPowerFrequency2',band(2), ...
            'SampleRate',         fs);
        sz = size(x);
        v  = reshape(x, [], nT).';
        v  = filtfilt(d, double(v));
        y  = reshape(cast(v.','like',x), sz);

    case 'acceleration'
        % Peak response of the second derivative of a Gaussian is at
        % omega = sqrt(2)/sigma, so the band centre fixes sigma.
        fc  = mean(band);
        sig = sqrt(2)/(2*pi*fc)*fs;
        L   = max(2, ceil(4*sig));
        t   = (-L:L);
        g   = exp(-t.^2./(2*sig^2));
        h   = (t.^2./sig^4 - 1/sig^2).*g;

        % Unit gain, right sign, at the band centre.  h is real and symmetric so
        % this response is real; dividing by it does both jobs at once.
        h = h ./ real(sum(h.*exp(-1i*2*pi*fc*t./fs)));
        y = imfilter(x, reshape(h,1,1,[]), 'symmetric', 'same');

    otherwise
        error('meTemporal:filter', ...
            'The temporal filter is ideal, cyclic, butter or acceleration, not %s.', ...
            s.temporalFilter);
end
end
%------------- END OF CODE --------------
