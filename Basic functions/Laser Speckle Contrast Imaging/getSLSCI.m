%getSLSCI - Calculates Spatial Laser Speckle Contrast Images (sLSCI).
%
%   This function computes local spatial contrast using a square kernel.
%   It includes memory-safe batch processing, edge-artifact correction
%   (via normalization maps), and optional post-calculation temporal 
%   decimation (averaging).
%
% Syntax:  data = getSLSCI(dataIn, kernelSize, procType)
%          data = getSLSCI(dataIn, kernelSize, procType, nFramesAvg)
%          data = getSLSCI(dataIn, kernelSize, procType, nFramesAvg, memoryCoef)
%
% Inputs:
%    dataIn     - Raw laser speckle data as a 3D [y, x, t] matrix.
%                 Supports uint8, uint16, single, or double.
%    kernelSize - Scalar integer specifying the side length of the square 
%                 kernel (e.g., 5 or 7).
%    procType   - String specifying the processing mode:
%                 'cpu' - Standard processing on the CPU.
%                 'gpu' - Parallel processing on the GPU (requires Parallel Computing Toolbox).
%
% Optional Inputs (Positional):
%    nFramesAvg - (Scalar) Number of contrast frames to average temporally 
%                 after calculation. Default is 1 (no averaging).
%                 Example: Set to 10 to reduce output file size by 10x.
%    memoryCoef - (Scalar, 0.0 to 1.0) Fraction of available RAM/GPU memory
%                 to use for batch allocation. Default is 0.8.
%
% Outputs:
%    data       - Processed contrast data as a 3D [y, x, t_new] matrix 
%                 of type 'single'.
%
% Example 1 (Basic GPU processing):
%    sLSCI = getSLSCI(raw_data, 5, 'gpu');
%
% Example 2 (CPU processing with temporal averaging):
%    % Calculate contrast with 7x7 kernel, then average every 10 frames
%    sLSCI = getSLSCI(raw_data, 7, 'cpu', 10);
%
% Example 3 (Adjusting memory safety):
%    % Use only 50% of available GPU memory
%    sLSCI = getSLSCI(raw_data, 5, 'gpu', 1, 0.5);
%
% See also: getTLSCI.m
% Authors: Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 29-January-2026

%------------- BEGIN CODE --------------

function data=getSLSCI(dataIn,kernelSize,procType,varargin)

memoryCoef=0.8;
nFramesAvg=1; % Default: No temporal decimation
if ~isempty(varargin)
    nFramesAvg = varargin{1};
    if length(varargin) > 1
        memoryCoef = varargin{2};
    end
end

sz = size(dataIn);
outT = floor(sz(3) / nFramesAvg);
outSz = [sz(1), sz(2), outT];
framesToProc=outT*nFramesAvg;

[~, memAvailRAM] = memory;
memAvailRAM = memAvailRAM.PhysicalMemory.Available;
if memAvailRAM*memoryCoef<prod(outSz)*4
    error('Insufficient memory available for keeping the processed file.');
end

kernel=single(ones(kernelSize,kernelSize,1)./(kernelSize^2));
normMap=single(conv2(ones(size(dataIn,1),size(dataIn,2)), gather(kernel), 'same'));

switch procType
    case 'gpu'
        data=zeros(outSz,'single');
        memAvailGPU = gpuDevice;
        memAvailGPU = memAvailGPU.AvailableMemory;
        batchNum=ceil(max((3*numel(dataIn)*4)./(memAvailGPU*memoryCoef),1));
        batchSize=floor(size(dataIn,3)./batchNum);

        % Ensure batchSize is a multiple of nFramesAvg
        batchSize = floor(batchSize/nFramesAvg) * nFramesAvg;
        if batchSize == 0, batchSize = nFramesAvg; end

        kernel=gpuArray(kernel);
        normMap=gpuArray(normMap);
        for i=1:batchSize:framesToProc
            i2=min(i+batchSize-1,framesToProc);
            dataGPU=single(gpuArray(dataIn(:,:,i:i2)));
            m=convn(dataGPU, kernel, 'same')./normMap;
            dataGPU=convn(dataGPU.^2, kernel, 'same')./normMap;
            dataGPU = sqrt(max(dataGPU - m.^2, 0));
            dataGPU=dataGPU./m;

            if nFramesAvg > 1
                [h, w, t] = size(dataGPU);
                dataGPU = reshape(dataGPU, h, w, nFramesAvg, t/nFramesAvg);
                dataGPU = mean(dataGPU, 3);
                dataGPU = reshape(dataGPU, h, w, t/nFramesAvg);
            end
            
            j = (i-1)/nFramesAvg + 1;
            j2 = j + size(dataGPU,3) - 1;
            data(:,:,j:j2)=gather(dataGPU);
        end
    case 'cpu' %convn can be replaced with imboxfilt for a minor (less than 5% boost in performance)
        data=zeros(outSz,'single');
        [~, memAvailRAM] = memory;
        memAvailRAM = memAvailRAM.PhysicalMemory.Available;
        batchNum=ceil(max((4*numel(dataIn)*4)./(memAvailRAM*memoryCoef),1)); % operations below should take aprox 2x batch array size, but, we allocate 3x just for safety
        batchSize=floor(size(dataIn,3)./batchNum);

        % Ensure batchSize is a multiple of nFramesAvg
        batchSize = floor(batchSize/nFramesAvg) * nFramesAvg;
        if batchSize == 0, batchSize = nFramesAvg; end

        for i=1:batchSize:framesToProc
            i2=min(i+batchSize-1,framesToProc);
            chunk = single(dataIn(:,:,i:i2));
            m=convn(chunk, kernel, 'same')./normMap;
            chunk=convn(chunk.^2, kernel, 'same')./normMap;
            chunk = sqrt(max(chunk - m.^2, 0));
            chunk=chunk./m;

            if nFramesAvg > 1
                [h, w, t] = size(chunk);
                chunk = reshape(chunk, h, w, nFramesAvg, t/nFramesAvg);
                chunk = squeeze(mean(chunk, 3));
                chunk = reshape(chunk, h, w, t/nFramesAvg);
            end

            j = (i-1)/nFramesAvg + 1;
            j2 = j + size(chunk,3) - 1;
            data(:,:,j:j2) = chunk;
        end
end
end

%------------- END OF CODE --------------
% Comments: note that there might be up to e-07 order rounding error
% processing types.
%
% Large input data can lead to a memory overflow, particularly
% when uint8 input data is provided. This can be controlled by an additional
% outer loop and/or by converting sLSCI data to scaled integers or by
% allowing downsampling with a specific kernel size. Futhermore, memory overflow
% can be controlled when using 'fastcpu' or 'fastgpu' by specifying a
% smaller batch size.
%
% Finally, note that the 'gpu' and the 'fastgpu' processing types require
% a MATLAB supported GPU and their drivers to be installed.