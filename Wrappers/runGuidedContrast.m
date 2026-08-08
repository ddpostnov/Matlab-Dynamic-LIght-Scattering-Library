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
%   results.gsTime.  That gsData holds CONTRAST is said by
%   settings.runGuidedContrast: the guided step leaves exactly one of
%   runGuidedContrast / runGuidedIntensity in the settings, and that is the
%   library's only record of which quantity the traces are.
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
%    fNames    - cell array of segmented *_K_r.mat (or *_I_r.mat) paths; each must
%                already contain results.sMap.
%    fNamesRaw - (optional) cell array of the matching raw recordings (.rls or
%                .cxd), same size as fNames.  If omitted or left empty, the raw
%                file name is derived from each *_r.mat name and expected in the
%                same folder (just like the rest of the pipeline).
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.
%
% Outputs:
%    (none) - updates each *_r.mat with results.gsData [nFrames x nRegions] and
%             results.gsTime [nFrames x 1], and records settings.runGuidedContrast,
%             removing settings.runGuidedIntensity if an earlier run left one.
%
% Example:
%    s.libraryFolder = libraryFolder;
%    fNames = getFileNamesList(rootFolder,'*_t_K_r.mat');
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

if ~all( cellfun(@(x) isempty(x) || contains(x,'_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_t_K_r.mat'').']);
end
if nargin<3 || isempty(fNamesRaw)
    fNamesRaw=deriveRawNames(fNames,s);
end
if ~isfield(s,'memoryCoef') || isempty(s.memoryCoef), s.memoryCoef=0.25; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Guided',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=char(fNames{fidx});
        s.fNameRaw=char(fNamesRaw{fidx});
        reportFile(rep,fidx,s.fName);
        clearvars results settings

        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        if ~isfield(results,'sMap')
            error('Reference results file is missing sMap. Run runSegmentation first.');
        end
        if ~isfile(s.fNameRaw)
            error(['Raw recording not found: ',s.fNameRaw,newline,...
                'Place it next to the *_r.mat file or pass fNamesRaw explicitly.']);
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
        end
        closeRawStream(st,cfg);

        results.gsData=gsData;
        results.gsTime=toSeconds(timeStamps,cfg);

        % Save the settings and results.  The settings are the only record of which
        % quantity gsData holds, so the sibling wrapper's entry goes when this one is
        % written - re-running the guided step the other way round must leave ONE
        % answer standing, not two.
        s.rawFrameRate=1./median(diff(results.gsTime));
        settings.runGuidedContrast=reportSettings(s);
        if isfield(settings,'runGuidedIntensity')
            settings=rmfield(settings,'runGuidedIntensity');
        end
        reportWriting(rep);
        save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
        save(s.fName,'results','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);
end

%------------- LOCAL FUNCTIONS --------------
function fNamesRaw=deriveRawNames(fNames,s)
% Map each segmented *_r.mat name to its raw recording: back out of the results tree
% first, then swap the extension.  The two are separate steps because only the FOLDER
% ever moved - the stem is the recording's own in both trees.  With no results folder
% set the first step is a no-op and the raw sits beside the product, as it always did.
fNamesRaw=cell(size(fNames));
for i=1:1:numel(fNames)
    f=fNames{i};
    if isempty(f), fNamesRaw{i}=''; continue; end
    f=getResultsPath(f,s,'root');
    if contains(f,'_I_r.mat')
        % EVERY intensity product carries its stage flag - '_a' angiogram, '_c'
        % cardiac, '_b' bolus - so the flag comes off with the product token, exactly
        % as it does for '_K' one line below.  The flagless '_I_r.mat' this used to
        % strip was runIntensity's old name; nothing writes it and nothing reads it,
        % and stripping only '_I_r.mat' turned 'M_a_I_r.mat' into 'M_a.cxd'.
        r=regexprep(f,'_[a-z]_I_r\.mat$','.cxd');        % intensity pipeline -> .cxd
        if strcmp(r,f)
            error('runGuidedContrast:noRawName', ...
                ['Cannot name the recording behind "%s": an intensity product is ' ...
                 '"<stem>_a_I_r.mat", "_c_I_r.mat" or "_b_I_r.mat".'],f);
        end
        fNamesRaw{i}=r;
    else
        r=regexprep(f,'_[a-z]_K_r\.mat$','.rls');        % contrast pipeline -> .rls
        if strcmp(r,f), r=regexprep(f,'_K_r\.mat$','.rls'); end
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
