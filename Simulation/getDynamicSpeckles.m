%getDynamicSpeckles - Simulate dynamic laser speckle image stacks.
%
% Syntax:
%    [I, info] = getDynamicSpeckles('Name', Value, ...)
%
% Description:
%    Generates a time series of laser speckle intensity images from the
%    coherent superposition of the fields scattered by moving point
%    particles (Zilpelwar et al., Biomed. Opt. Express 13(12):6533, 2022;
%    Postnov et al., Sci. Rep. 9:2542, 2019). Every pixel is the modulus
%    squared of the summed spherical waves, E = sum_n exp(i k r_n)/r_n.
%
%    The speckle-to-pixel ratio is controlled WITHOUT distorting the
%    illumination: the optical field is rendered on a fine grid sampled at
%    'sampPerSpeckle' points per speckle (>= Nyquist) and every output pixel
%    is the area-average of 'sampPerPixel' sub-samples per side. Changing
%    'pixelsPerSpeckle' only changes how the SAME finely sampled field is
%    binned, so the illumination envelope, speckle size and statistics are
%    left intact. The physical speckle size is measured from the rendered
%    field (calibrateSpeckle) so the requested ratio is met accurately.
%
%    Three decorrelation regimes are provided through 'motionType', and each
%    output frame may either be instantaneous ('exposureTime' empty) or the
%    integral of many instantaneous sub-frames over a finite exposure.
%
% Inputs (all optional Name-Value pairs):
%    'pixelsN'          - Output pixels, scalar N or [Ny Nx]. Default 64.
%    'pixelsPerSpeckle' - Speckle-to-pixel ratio s/p (pixels spanned by one
%                         speckle). >1 oversampled/high contrast, <1
%                         undersampled/low contrast. Default 4.
%    'sampPerSpeckle'   - Fine-grid samples across one speckle (Nyquist >=2,
%                         4 is safe). Default 4.
%    'nFrames'          - Number of output frames. Default 100.
%    'motionType'       - 'brownian' (unordered diffusion, slow decorr.),
%                         'ordered' (directed flow) or 'decorrelated'
%                         (independent frame to frame). Default 'brownian'.
%    'tauC'             - Target decorrelation time, us. Default 10.
%    'dT'               - Time step between instantaneous sub-frames, us.
%                         Default 1.
%    'exposureTime'     - Exposure per output frame, us. [] = instantaneous
%                         (one sub-frame); otherwise round(exposureTime/dT)
%                         sub-frames are integrated. Default [].
%    'lambda'           - Wavelength, um. Default 0.785.
%    'particlesN'       - Number of scattering particles. Default 1000.
%    'sizeX','sizeY'    - Scattering volume width, um. Default 1000.
%    'sizeZ'            - Scattering volume depth, um. Default 100.
%    'sensorZ'          - Sensor distance from the volume, um. Default 3500.
%    'addNoise'         - Apply the camera noise model. Default false.
%    'cameraParams'     - Struct of camera settings, missing fields take the
%                         defaults: photonsPerPixel 1000 (mean photons/pixel
%                         at the mean signal - the main SNR knob; set [] to
%                         derive it from flux instead), flux 30
%                         photons/um^2/ms, QE 0.6, darkCurrent 0.014 e/ms,
%                         readNoise 1.4 e, bias 0 e, bitRange 255, meanLevel
%                         0.25 (mean signal as a fraction of bitRange; lower
%                         it to avoid saturating bright speckle).
%    'procType'         - 'gpu' or 'cpu'. Default 'gpu' (falls back to cpu).
%    'calibrateSpeckle' - Measure the speckle size from the field so the
%                         requested ratio is accurate. Default true.
%    'speckleSize'      - Fix the speckle size (um) explicitly, skipping the
%                         estimate and calibration (useful to keep many calls
%                         consistent). [] = auto. Default [].
%    'rngSeed'          - Seed for reproducibility. [] leaves rng untouched.
%                         Default 0.
%    'memoryCoef'       - Fraction of free device memory used per field tile.
%                         Default 0.25.
%    'verbose'          - Print progress. Default false.
%
% Outputs:
%    I    - Speckle intensity stack, single [Ny Nx nFrames]. Mean-normalised
%           to 1 when addNoise is false; integer digital counts in
%           [0 bitRange] when addNoise is true.
%    info - Struct with the achieved parameters (speckleSize um, pixelSize,
%           pixelsPerSpeckle, sampPerPixel, effective sampPerSpeckle, fineN,
%           nSub, exposureTime, geometry, seed, ...).
%
% Examples:
%    % Instantaneous, slowly decorrelating stack (DLSI-style):
%    I = getDynamicSpeckles('motionType','brownian','tauC',10,'nFrames',200);
%
%    % Exposure-integrated frames (LSCI-style):
%    I = getDynamicSpeckles('motionType','brownian','exposureTime',50);
%
%    % Independent frames with a fine field for a pixels/speckle sweep:
%    I = getDynamicSpeckles('motionType','decorrelated', ...
%           'pixelsPerSpeckle',8,'sampPerSpeckle',16,'pixelsN',256);
%
% Dependencies: Parallel Computing Toolbox (GPU, optional); Statistics and
%    Machine Learning Toolbox (poissrnd/normrnd, only when addNoise=true).
% See also: getNormalizedG2, getTauC, getK, Launcher_speckleSimulation
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 20-July-2026

%------------- BEGIN CODE --------------
function [I, info] = getDynamicSpeckles(varargin)

p = inputParser;
p.KeepUnmatched = false;
addParameter(p, 'pixelsN', 64, @(x) isnumeric(x) && any(numel(x)==[1 2]));
addParameter(p, 'pixelsPerSpeckle', 4, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sampPerSpeckle', 4, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p, 'nFrames', 100, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'motionType', 'brownian', @(x) any(validatestring(x, {'brownian','ordered','decorrelated'})));
addParameter(p, 'tauC', 10, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'dT', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'exposureTime', [], @(x) isempty(x) || (isscalar(x) && x>0));
addParameter(p, 'lambda', 0.785, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'particlesN', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'sizeX', 1000, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sizeY', 1000, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sizeZ', 100, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sensorZ', 3500, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'addNoise', false, @(x) islogical(x) || ismember(x,[0 1]));
addParameter(p, 'cameraParams', struct(), @isstruct);
addParameter(p, 'procType', 'gpu', @(x) any(validatestring(x, {'gpu','cpu'})));
addParameter(p, 'calibrateSpeckle', true, @(x) islogical(x) || ismember(x,[0 1]));
addParameter(p, 'speckleSize', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
addParameter(p, 'rngSeed', 0, @(x) isempty(x) || (isscalar(x) && x>=0));
addParameter(p, 'memoryCoef', 0.25, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
addParameter(p, 'verbose', false, @(x) islogical(x) || ismember(x,[0 1]));
parse(p, varargin{:});
r = p.Results;

% Output geometry
if isscalar(r.pixelsN)
    pixelsNy = r.pixelsN; pixelsNx = r.pixelsN;
else
    pixelsNy = r.pixelsN(1); pixelsNx = r.pixelsN(2);
end

% Reproducibility
if ~isempty(r.rngSeed)
    rng(r.rngSeed);
end

% Device selection (gpu requested but unavailable -> cpu)
useGPU = strcmpi(r.procType,'gpu');
if useGPU
    try
        useGPU = gpuDeviceCount > 0;
    catch
        useGPU = false;
    end
end
if useGPU
    toDev = @(x) gpuArray(single(x));
    gd = gpuDevice; memAvail = gd.AvailableMemory;
else
    toDev = @(x) single(x);
    [~, mem] = memory; memAvail = mem.PhysicalMemory.Available;
end

% Camera parameters (defaults merged with user overrides)
cam = mergeStruct(defaultCamera(), r.cameraParams);

% Physical constants
k = 2*pi/r.lambda;                 % wavenumber, 1/um
sensorCx = r.sizeX/2;              % sensor centred on the volume
sensorCy = r.sizeY/2;

% Exposure / integration
if isempty(r.exposureTime)
    nSub = 1; exposureTime = r.dT;
else
    nSub = max(1, round(r.exposureTime/r.dT)); exposureTime = nSub*r.dT;
end

% Brownian step (per axis, um): D = 1/(tauC k^2), step = sqrt(2 D dT)
% -> field autocorrelation g1(tau)=exp(-tau/tauC), so g2=1+beta exp(-2 tau/tauC).
brownianStep = sqrt(2*r.dT/(r.tauC*k^2));

% Initial particle positions [1 x 1 x N]
pX = toDev(rand(1,1,r.particlesN)*r.sizeX);
pY = toDev(rand(1,1,r.particlesN)*r.sizeY);
pZ = toDev(rand(1,1,r.particlesN)*r.sizeZ);

% --- Speckle size ----------------------------------------------------------
speckleEst = r.lambda*r.sensorZ/min(r.sizeX,r.sizeY);   % analytic estimate
if ~isempty(r.speckleSize)
    speckleSize = r.speckleSize;                        % explicit override
elseif r.calibrateSpeckle
    speckleSize = measureSpeckleSize(pX,pY,pZ, speckleEst, k, r.sensorZ, ...
        sensorCx, sensorCy, r.sizeX, r.sizeY, r.sizeZ, toDev, memAvail, r.memoryCoef);
else
    speckleSize = speckleEst;
end
drift = speckleSize/r.tauC*r.dT;       % ordered-flow displacement per sub-frame, um

% --- Fine grid from the requested ratio ------------------------------------
pixelSize   = speckleSize/r.pixelsPerSpeckle;              % um
sampPerPix  = max(1, round(r.sampPerSpeckle/r.pixelsPerSpeckle));
subPixel    = pixelSize/sampPerPix;                        % fine-grid spacing, um
fineNy = pixelsNy*sampPerPix;
fineNx = pixelsNx*sampPerPix;
xFine = toDev(((1:fineNx)-(fineNx+1)/2)*subPixel + sensorCx);   % 1 x fineNx
yFine = toDev(((1:fineNy)-(fineNy+1)/2)*subPixel + sensorCy);   % 1 x fineNy

% Row-tile size to bound field-tensor memory
rowBlock = tileRows(fineNx, r.particlesN, memAvail, r.memoryCoef);

% --- Main loop -------------------------------------------------------------
I = zeros(pixelsNy, pixelsNx, r.nFrames, 'single');
for f = 1:r.nFrames
    acc = toDev(zeros(fineNy, fineNx));
    for s = 1:nSub
        acc = acc + renderFine(pX,pY,pZ, xFine,yFine, k, r.sensorZ, rowBlock);
        [pX,pY,pZ] = stepParticles(pX,pY,pZ, r.motionType, brownianStep, drift, ...
            r.sizeX, r.sizeY, r.sizeZ, toDev);
    end
    I(:,:,f) = gather(binMean(acc, sampPerPix));
    if r.verbose && (mod(f,max(1,round(r.nFrames/10)))==0 || f==r.nFrames)
        fprintf('getDynamicSpeckles: frame %d/%d\n', f, r.nFrames);
    end
end

% --- Normalisation / camera noise -----------------------------------------
I = I/mean(I(:));                    % mean-normalised speckle intensity
if r.addNoise
    I = applyCameraNoise(I, cam, pixelSize, exposureTime);
end

% --- Info ------------------------------------------------------------------
info = struct();
info.speckleSize       = speckleSize;
info.speckleSizeEst    = speckleEst;
info.pixelSize         = pixelSize;
info.pixelsPerSpeckle  = speckleSize/pixelSize;
info.sampPerPixel      = sampPerPix;
info.sampPerSpeckleEff = sampPerPix*r.pixelsPerSpeckle;
info.fineN             = [fineNy fineNx];
info.pixelsN           = [pixelsNy pixelsNx];
info.nSub              = nSub;
info.exposureTime      = exposureTime;
info.dT                = r.dT;
info.tauC              = r.tauC;
info.motionType        = r.motionType;
info.particlesN        = r.particlesN;
info.lambda            = r.lambda;
info.geometry          = struct('sizeX',r.sizeX,'sizeY',r.sizeY,'sizeZ',r.sizeZ,'sensorZ',r.sensorZ);
info.addNoise          = logical(r.addNoise);
info.procType          = ternary(useGPU,'gpu','cpu');
info.rngSeed           = r.rngSeed;
end

%------------- LOCAL FUNCTIONS ---------------------------------------------
function fineI = renderFine(pX,pY,pZ, xFine, yFine, k, sensorZ, rowBlock)
% Modulus-squared of the superposed spherical waves on the fine grid,
% computed in row tiles to bound the [rows x cols x particles] tensor.
fineNy = numel(yFine); fineNx = numel(xFine);
fineI = zeros(fineNy, fineNx, 'like', real(xFine(1)));
for r0 = 1:rowBlock:fineNy
    r1 = min(r0+rowBlock-1, fineNy);
    yB = reshape(yFine(r0:r1), [], 1);                     % rows x 1
    rr = sqrt((pX - xFine).^2 + (pY - yB).^2 + (pZ - sensorZ).^2);  % rows x cols x N
    E  = sum(exp(1i*k*rr)./rr, 3);                         % rows x cols
    fineI(r0:r1,:) = real(E.*conj(E));
end
end

function [pX,pY,pZ] = stepParticles(pX,pY,pZ, motionType, brownianStep, drift, sizeX,sizeY,sizeZ, toDev)
% Advance particle positions by one sub-frame (dT).
N = numel(pX);
switch motionType
    case 'brownian'                       % isotropic diffusion, periodic volume
        pX = mod(pX + brownianStep*randn(1,1,N,'like',pX), sizeX);
        pY = mod(pY + brownianStep*randn(1,1,N,'like',pY), sizeY);
        pZ = mod(pZ + brownianStep*randn(1,1,N,'like',pZ), sizeZ);
    case 'ordered'                        % directed flow along +x, refill at exit
        pX = pX + drift;
        out = pX > sizeX;
        if any(out(:))
            pX(out) = pX(out) - sizeX;
            pY(out) = toDev(rand(1,1,nnz(out)))*sizeY;
            pZ(out) = toDev(rand(1,1,nnz(out)))*sizeZ;
        end
    case 'decorrelated'                   % independent configuration each sub-frame
        pX = toDev(rand(1,1,N))*sizeX;
        pY = toDev(rand(1,1,N))*sizeY;
        pZ = toDev(rand(1,1,N))*sizeZ;
end
end

function out = binMean(fineI, b)
% Area-average b x b fine samples into each output pixel.
if b==1, out = fineI; return; end
[fy, fx] = size(fineI);
out = reshape(mean(reshape(fineI, b, fy/b, fx), 1), fy/b, fx);      % rows
out = reshape(mean(reshape(out.', b, fx/b, fy/b), 1), fx/b, fy/b).'; % cols
end

function speckleSize = measureSpeckleSize(pX,pY,pZ, speckleEst, k, sensorZ, cx, cy, sizeX,sizeY,sizeZ, toDev, memAvail, memoryCoef)
% Speckle size = sqrt of the intensity correlation area, averaged over a few
% independent frames. This correlation-area definition is the one consistent
% with the pixel-integrated contrast K = (s/p)/sqrt(1+(s/p)^2).
nCal = 6; Ncal = 128; dCal = speckleEst/8;
xC = toDev(((1:Ncal)-(Ncal+1)/2)*dCal + cx);
yC = toDev(((1:Ncal)-(Ncal+1)/2)*dCal + cy);
rowBlock = tileRows(Ncal, numel(pX), memAvail, memoryCoef);
ac = zeros(Ncal, Ncal, 'single');
for c = 1:nCal
    fineI = gather(renderFine(pX,pY,pZ, xC,yC, k, sensorZ, rowBlock));
    d = fineI - mean(fineI(:));
    a = fftshift(real(ifft2(abs(fft2(d)).^2)));   % unnormalised autocovariance
    ac = ac + a/max(a(:));                         % normalise to peak 1
    pX = toDev(rand(1,1,numel(pX)))*sizeX;         % fresh speckle realisation
    pY = toDev(rand(1,1,numel(pY)))*sizeY;
    pZ = toDev(rand(1,1,numel(pZ)))*sizeZ;
end
ac = ac/nCal;                                     % normalised autocovariance, peak 1
% Integrate over a central window: this captures the true correlation area
% while excluding the periodic negative pedestal (a full-frame sum is 0 after
% mean subtraction) and the rectified far-field noise a clipped sum would add.
c0 = Ncal/2 + 1;                                  % zero-lag index after fftshift
Rwin = min(floor(Ncal/2)-1, round(4*speckleEst/dCal));
win = ac(c0-Rwin:c0+Rwin, c0-Rwin:c0+Rwin);
corrArea = sum(win(:))*dCal*dCal;                 % um^2 (Goodman speckle area)
speckleSize = sqrt(max(corrArea, dCal*dCal));
end

function rowBlock = tileRows(fineNx, N, memAvail, memoryCoef)
% Rows processed at once so that [rows x fineNx x N] tensors fit the budget.
bytesPerElem = 48;                                % r, exp(complex), ./r temporaries
rowBlock = max(1, floor(memAvail*memoryCoef/(fineNx*N*bytesPerElem)));
end

function I = applyCameraNoise(I, cam, pixelSize, exposureTime)
% Poisson shot noise (with QE) + Poisson dark + Gaussian read, quantised to
% bitRange (Zilpelwar et al. 2022, Eq. 18). The mean signal is placed at
% meanLevel*bitRange (default 0.25, so bright speckle rarely saturates); the
% noise magnitude is set by the mean photons per pixel, given directly
% (cameraParams.photonsPerPixel) or derived from flux, pixel area and exposure.
sz = size(I);
tExpMs = exposureTime/1000;
if isfield(cam,'photonsPerPixel') && ~isempty(cam.photonsPerPixel)
    meanPhotons = cam.photonsPerPixel;
else
    meanPhotons = cam.flux*pixelSize^2*tExpMs;
end
Ephot = double(I)*meanPhotons;                      % expected photons per pixel
elec  = poissrnd(cam.QE*Ephot) ...                  % photoelectron shot noise
      + poissrnd(cam.darkCurrent*tExpMs, sz) ...    % dark current
      + round(normrnd(cam.bias, cam.readNoise, sz));% read noise
meanE = cam.QE*meanPhotons;                         % mean signal electrons
I = double(elec)/meanE*(cam.meanLevel*cam.bitRange);% digital counts
I(I<0) = 0; I(I>cam.bitRange) = cam.bitRange;       % clip to the sensor range
I = single(round(I));
end

function s = defaultCamera()
s = struct('photonsPerPixel',1000,'flux',30,'QE',0.6,'darkCurrent',0.014, ...
    'readNoise',1.4,'bias',0,'bitRange',255,'meanLevel',0.25);
end

function base = mergeStruct(base, over)
f = fieldnames(over);
for i = 1:numel(f), base.(f{i}) = over.(f{i}); end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
