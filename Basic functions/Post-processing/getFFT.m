% getFFT - calculates FFT spectrum and phase over averaged over temporal windows
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

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
