%getFFT - Windowed FFT amplitude and phase spectra, averaged over time windows.
%
%   Splits a signal (per pixel for 3-D input) into non-overlapping windows of
%   fftN samples, computes the single-sided FFT of each window and returns the
%   window-averaged amplitude and phase spectra.  Several per-window
%   normalizations are available.
%
% Syntax:
%    [fftAmp,fftPhase,f] = getFFT(data, fs, fftN, procType)
%    [fftAmp,fftPhase,f] = getFFT(data, fs, fftN, procType, normalizationType)
%
% Inputs:
%    data              - signal, 3-D [y x t] or a 1-D vector.
%    fs                - sampling frequency, Hz.
%    fftN              - window length (samples) and FFT size.
%    procType          - 'cpu' or 'gpu'.
%    normalizationType - (optional) per-window normalization before the FFT:
%                        'subtractMean' (default), 'none', 'samescale' (1-99 pct
%                        rescale then de-mean) or 'normalized' (divide by the
%                        mean, minus 1).
%
% Outputs:
%    fftAmp   - single-sided amplitude spectrum, [y x fftN/2] (or [fftN/2]).
%    fftPhase - phase spectrum (radians), same size.
%    f        - frequency axis, Hz, length fftN/2.
%
% Example:
%    [amp,ph,f] = getFFT(data, 100, 256, 'gpu', 'subtractMean');
%
% Dependencies: Parallel Computing Toolbox (procType 'gpu').
% See also: getPulsatility, getVasomotion
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function [fftAmp,fftPhase,f]=getFFT(data,fs,fftN,procType,varargin)

normalizationType='subtractMean';
if ~isempty(varargin)
    normalizationType = varargin{1};
end

if sum(size(data)>1)==1
    data=reshape(data,1,1,length(data(:)));
end

f = fs.*(0:(fftN/2))/fftN;
f=f(1:end-1);
fftAmp=zeros(size(data,1),size(data,2),fftN/2);
fftPhase=zeros(size(data,1),size(data,2),fftN/2);
for i=1:1:floor(size(data,3)./fftN)    
    subdata=single(data(:,:,(i-1)*fftN+1:i*fftN));
if strcmp(procType,'gpu')
    subdata=gpuArray(subdata);
end
switch normalizationType
    case 'none'
    case 'subtractMean'
        subdata=subdata-mean(subdata,3);
    case 'samescale'
        subdataP1=prctile(subdata,1,3);
    subdataP99=prctile(subdata,99,3);
    subdata=(subdata-subdataP1)./(subdataP99-subdataP1);
    subdata=subdata-mean(subdata,3);
    case 'normalized'        
    subdata=(subdata./mean(subdata,3))-1;
end
        

fftRes = fft(subdata,fftN,3);
fftRes = fftRes(:,:,1:fftN/2);
tmpAmp=abs(fftRes./fftN);
tmpAmp(:,:,2:end)= 2.*tmpAmp(:,:,2:end);
fftAmp = fftAmp+tmpAmp;
fftPhase=fftPhase+angle(fftRes);
end
fftAmp=squeeze(fftAmp./floor(size(data,3)./fftN));
fftPhase=squeeze(fftPhase./floor(size(data,3)./fftN));
if strcmp(procType,'gpu')
    fftAmp=gather(fftAmp);
    fftPhase=gather(fftPhase);
end
end
