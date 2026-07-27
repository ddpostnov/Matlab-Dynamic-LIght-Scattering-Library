%runGuidedContrast - Per-segment spatial contrast at full temporal resolution.
%
%   Uses an existing segmentation (results.sMap) as a guide to extract, for
%   every labelled region, the spatial contrast at the ORIGINAL frame rate of
%   the recording.  In the standard pipeline the raw recording is first turned
%   into temporally decimated contrast and then segmented, so results.sData is
%   only available at the decimated frame rate.  runGuidedContrast streams the
%   raw file in batches and, for each frame, computes one contrast value per
%   region as the coefficient of variation (std/mean) over that region's valid
%   pixels.  The result is a set of per-region traces at maximum temporal
%   resolution, stored in results.gsData with the matching time vector
%   results.gsTime (and results.gsType = "contrast").
%
%   "Valid" pixels are those that share the region label in results.sMap AND are
%   kept in results.mask.  The columns of results.gsData are aligned with the
%   rows of results.sMetrics (column i corresponds to region label i); regions
%   with no valid pixels are returned as NaN.
%
%   The raw recording is never loaded in full: frames are read in batches sized
%   to the available memory, so arbitrarily long recordings can be processed.
%
% Syntax:
%    runGuidedContrast(s, fNames)
%    runGuidedContrast(s, fNames, fNamesRaw)
%
% Inputs:
%    s         - parameter struct (uses s.libraryFolder for traceability;
%                optional s.memoryCoef, fraction of free RAM to use, default 0.25).
%    fNames    - cell array of segmented *_K_d.mat (or *_I_d.mat) paths; the
%                companion *_r.mat must already contain results.sMap.
%    fNamesRaw - (optional) cell array of the matching raw recordings (.rls or
%                .cxd), same size as fNames.  If omitted or left empty, the raw
%                file name is derived from each *_d.mat name and expected in the
%                same folder (just like the rest of the pipeline).
%
% Outputs:
%    (none) - updates each *_r.mat with results.gsData [nFrames x nRegions],
%             results.gsTime [nFrames x 1] and results.gsType, and records
%             settings.runGuidedContrast.
%
% Example:
%    s.libraryFolder = libraryFolder;
%    fNames = getFileNamesList(rootFolder,'*_t_K_d.mat');
%    runGuidedContrast(s,fNames(:));
%
% Dependencies: getPointerRLS (.rls); Bio-Formats bfGetReader/bfGetPlane (.cxd).
% See also: runGuidedIntensity, runSegmentation, runContrastFromRLS, getPointerRLS
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 20-July-2026
%------------- BEGIN CODE --------------

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% s.memoryCoef=0.25; %fraction of free RAM used per read batch

function runGuidedContrast(s,fNames,fNamesRaw)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_d.mat".');
end
if nargin<3 || isempty(fNamesRaw)
    fNamesRaw=deriveRawNames(fNames);
end
if ~isfield(s,'memoryCoef') || isempty(s.memoryCoef), s.memoryCoef=0.25; end

for fidx=1:1:numel(fNames)
    if ~isempty(fNames{fidx})
        tic
        s.fName=char(fNames{fidx});
        s.fNameRaw=char(fNamesRaw{fidx});
        disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        clearvars results settings

        load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
        load(strrep(s.fName,'_d.mat','_r.mat'),'results');

        if ~isfield(results,'sMap')
            error('Reference results file is missing sMap. Run runSegmentation first.');
        end
        if ~isfile(s.fNameRaw)
            error(['Raw recording not found: ',s.fNameRaw,newline,...
                'Place it next to the *_d.mat file or pass fNamesRaw explicitly.']);
        end

        % Indicator matrix of valid pixels per region (same sMap index AND kept
        % in the mask).  A sparse [nValidPix x nRegions] matrix lets us reduce a
        % whole batch of frames to per-region sums with a single matrix product.
        sMap=results.sMap;
        if isfield(results,'mask') && ~isempty(results.mask) ...
                && isequal(size(results.mask),size(sMap))
            validMask=(sMap>0) & results.mask;
        else
            validMask=(sMap>0);
        end
        nRegions=double(max(sMap(:)));
        pixIdx=find(validMask);
        labels=double(sMap(pixIdx));
        regInd=sparse(1:numel(pixIdx),labels,1,numel(pixIdx),nRegions); % [P x nRegions]
        countPerRegion=full(sum(regInd,1))';                            % [nRegions x 1]

        % Open the raw recording, verify geometry and size the read batches
        [st,cfg]=openRawStream(s.fNameRaw);
        if cfg.sizeY~=size(sMap,1) || cfg.sizeX~=size(sMap,2)
            closeRawStream(st,cfg);
            error(['Raw frame size (',num2str(cfg.sizeY),'x',num2str(cfg.sizeX),...
                ') does not match sMap (',num2str(size(sMap,1)),'x',num2str(size(sMap,2)),...
                '). Guided extraction needs the original, uncropped recording.']);
        end
        nT=cfg.sizeT;
        batchSize=getBatchSize(numel(pixIdx),cfg.sizeY*cfg.sizeX,nT,s.memoryCoef);

        gsData=nan(nT,nRegions);
        timeStamps=zeros(nT,1);
        done=0;
        while done<nT
            b=min(batchSize,nT-done);
            [chunk,tsB,st]=readRawBatch(st,cfg,b,done);
            vals=double(reshape(chunk,cfg.sizeY*cfg.sizeX,b));
            vals=vals(pixIdx,:);                       % [P x b] valid pixels
            sumX =regInd'*vals;                        % [nRegions x b]
            sumX2=regInd'*(vals.^2);
            meanX=sumX./countPerRegion;
            varX =max(sumX2./countPerRegion-meanX.^2,0);   % population variance (/N), as in getK
            gsData(done+1:done+b,:)=(sqrt(varX)./meanX)';
            timeStamps(done+1:done+b)=tsB;
            done=done+b;
            disp(['   frames ',num2str(done),'/',num2str(nT),', elapsed ',num2str(round(toc)),'s'])
        end
        closeRawStream(st,cfg);

        results.gsData=gsData;
        results.gsTime=toSeconds(timeStamps,cfg);
        results.gsType="contrast";

        % Save the settings and results
        s.rawFrameRate=1./median(diff(results.gsTime));
        settings.runGuidedContrast=s;
        disp(['Saving the results. Elapsed time ',num2str(round(toc)),'s']);
        save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
        save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
        disp('Saving complete');
    end
end
end

%------------- LOCAL FUNCTIONS --------------
function fNamesRaw=deriveRawNames(fNames)
% Map each segmented *_d.mat name to its raw recording in the same folder.
fNamesRaw=cell(size(fNames));
for i=1:1:numel(fNames)
    f=fNames{i};
    if isempty(f), fNamesRaw{i}=''; continue; end
    if contains(f,'_I_d.mat')
        fNamesRaw{i}=regexprep(f,'_I_d\.mat$','.cxd');   % intensity pipeline -> .cxd
    else
        r=regexprep(f,'_[a-z]_K_d\.mat$','.rls');        % contrast pipeline -> .rls
        if strcmp(r,f), r=regexprep(f,'_K_d\.mat$','.rls'); end
        fNamesRaw{i}=r;
    end
end
end

function batchSize=getBatchSize(nValidPix,nPixFrame,nT,memoryCoef)
% Frames per read batch, sized so the working arrays fit in free RAM.
try
    [~,mem]=memory; avail=mem.PhysicalMemory.Available;
catch
    avail=2e9;   % conservative fallback if memory() is unavailable
end
perFrameBytes=nValidPix*8*3 + nPixFrame*2;   % vals + vals.^2 + slack + raw chunk
batchSize=max(1,min(nT,floor(avail*memoryCoef/perFrameBytes)));
end

function [st,cfg]=openRawStream(fNameRaw)
% Open a raw recording for batched reading (.rls or Bio-Formats .cxd).
[~,~,ext]=fileparts(fNameRaw);
cfg.isRLS=strcmpi(ext,'.rls');
if cfg.isRLS
    [st.fp,meta]=getPointerRLS(fNameRaw);
    cfg.sizeY=double(meta.sizeY);
    cfg.sizeX=double(meta.sizeX);
    cfg.sizeT=double(meta.sizeT);
    cfg.dataType=meta.dataType;
else
    st.reader=bfGetReader(fNameRaw);
    omeMeta=st.reader.getMetadataStore();
    cfg.sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
    cfg.sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
    cfg.sizeT=double(omeMeta.getPlaneCount(0));
    cfg.dataType=char(omeMeta.getPixelsType(0).toString());
    cfg.sampling=double(omeMeta.getPixelsTimeIncrement(0).value());  % seconds
end
end

function [chunk,ts,st]=readRawBatch(st,cfg,b,done)
% Read b consecutive frames (and their time stamps) into [sizeY x sizeX x b].
chunk=zeros(cfg.sizeY,cfg.sizeX,b,cfg.dataType);
ts=zeros(b,1);
if cfg.isRLS
    for j=1:1:b
        ts(j)=double(fread(st.fp,1,'*uint64'));
        chunk(:,:,j)=fread(st.fp,[cfg.sizeY,cfg.sizeX],['*',cfg.dataType]);
    end
else
    for j=1:1:b
        chunk(:,:,j)=bfGetPlane(st.reader,done+j);
    end
    ts=((done:done+b-1)').*cfg.sampling;                 % seconds
end
end

function closeRawStream(st,cfg)
if cfg.isRLS, fclose(st.fp); else, st.reader.close(); end
end

function time=toSeconds(timeStamps,cfg)
% Relative time vector in seconds, starting at zero.
if cfg.isRLS
    time=(timeStamps-timeStamps(1))./1000;               % .rls stamps are in ms
else
    time=timeStamps-timeStamps(1);                        % .cxd stamps already in s
end
time=time(:);
end
%------------- END OF CODE --------------
