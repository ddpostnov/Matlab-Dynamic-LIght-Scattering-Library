%runIntensity - Read a .cxd recording plane by plane and write its ANGIOGRAM product.
%
%   THE RECORDING IS READ ONCE, AND A PLANE AT A TIME.  The version this replaces
%   allocated the whole stack before the first frame arrived - 129 GB as single for
%   the 60 GB reference recording, so it could not be run on that recording at all.
%   The loop below holds ONE plane and ONE accumulator the size of a plane, whatever
%   the recording is.  The only large array left is the output cube itself, and that
%   one is sized, compared against free memory and REFUSED BEFORE THE READ rather
%   than discovered twenty minutes into it.
%
%   THE PRODUCT IS THE ANGIOGRAM, '<stem>_a_I_d/_r/_s.mat'.  Stage 'a', product 'I',
%   branch 'angiogram' - the whole-recording intensity stack and its time-mean image.
%   It is one of the three products a .cxd can yield ('_a_I' here, '_c_I' the cardiac
%   cycle, '_b_I' a bolus); which one a recording is for is the protocol's answer, not
%   this step's.  The flagless '_I' this used to write is gone: no code in the library
%   reads that name.
%
%   DECIMATION IS TEMPORAL AND NOTHING ELSE.  s.decimFactor consecutive frames are
%   AVERAGED into one output frame and the blocks do not overlap - the 'sharp' sense of
%   the contrast step's decimMethod - so the output frame rate is the original divided
%   by s.decimFactor.  A trailing partial block is dropped, exactly as the contrast
%   core drops one.  THE FULL SPATIAL RESOLUTION IS ALWAYS KEPT: an angiogram is read
%   for its vessels, and a vessel is a few pixels wide.
%
%   THE TIME OF AN OUTPUT FRAME IS THE TIME OF THE FIRST RAW FRAME IN ITS BLOCK, so
%   time(1) is 0 whatever the decimation, and the acquisition clock is in
%   results.timeStamp.
%
%   WHETHER THE FRAMES ARE KEPT IS A DECISION ABOUT WHAT THE RECORDING IS FOR.
%   s.saveSource stays, and it is TRUE by default.  Off, the step writes '_a_I_r' and
%   '_a_I_s' and no '_a_I_d': the angiogram, the intensity traces and the settings, at
%   a memory cost of one frame however long the recording is.  What that buys is the
%   only way a recording too long to hold as a cube yields an angiogram at all.  What
%   it costs is every later step that reads the frames - velocity, vasomotion, motion
%   enhancement, background removal, every per-pixel path - which then has nothing to
%   read on that recording and cannot be run on it.  What still works is everything
%   derived from the mean image: regions, segmentation, vessel density and topology.
%   Turning it off is not a smaller version of the same product, it is a different one,
%   and the tooltip in the workbench says so.
%
%   IT SAYS HOW FAR IT HAS GOT, AND ONLY WHEN THAT IS WORTH SAYING.  The reporting
%   seam is three lines per recording, which is right for a sixty-file batch and wrong
%   for a single read that outlives the operator's patience.  So this step emits a
%   progress line with an ETA through the same sink - the command window and the
%   workbench log alike - but only once the read has already taken longer than
%   REPORTEVERYSEC below.  A short recording therefore still narrates in exactly three
%   lines, and the twenty-minute one does not go silent.
%
% Syntax:
%    runIntensity(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%               .decimFactor  frames averaged into each saved frame (1 = none).
%               .dataTypeOut  class the frames are stored as, e.g. 'single'.
%                             Empty keeps the recording's own class.
%               .saveSource   logical; write the frames as well as the angiogram.
%    fNames - cell array of .cxd file paths.
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.
%
% Outputs:
%    (none) - writes, through getResultsPath, beside the recording or in the project's
%             results folder:
%               <stem>_a_I_d.mat  SOURCE   source.data [Y x X x T], source.time [T x 1]
%               <stem>_a_I_r.mat  RESULTS  time [T x 1], timeStamp, imgI (the
%                                          angiogram), meanI / minI / maxI [T x 1]
%               <stem>_a_I_s.mat  SETTINGS settings.runIntensity
%               <stem>_rep_intensity.jpg   angiogram + the intensity trace
%
% Example:
%    s.libraryFolder = libraryFolder;
%    s.decimFactor   = 1;
%    s.dataTypeOut   = 'single';
%    s.saveSource    = true;
%    D = dir(fullfile(dataRoot,'*.cxd'));
%    runIntensity(s, fullfile({D.folder}',{D.name}'));
%
% Dependencies: the Bio-Formats MATLAB library (bfGetReader, bfGetPlane).
% See also: runIntensityBolus, runContrastFromRLS, getResultsPath, wbFileModel,
%           wbStepRegistry
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% s.decimFactor=1; %frames averaged into each saved frame. Output framerate = original / s.decimFactor
% s.dataTypeOut='single'; %class the frames are stored as; '' keeps the recording's own
% s.saveSource=true; %false keeps only the angiogram - see the header before using it

%------------- BEGIN CODE --------------
function runIntensity(s,fNames)

REPORTEVERYSEC = 30;    % how often a read that is running long says where it is

if ~all( cellfun(@(f) isempty(f) || contains(f,'.cxd'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".cxd".');
end

% The three settings this step has, defaulted where a hand-written s omits one.  The
% registry preset supplies all three, so this is for a script that names only what it
% wants changed.  saveSource defaults to TRUE: a whole product is the safe direction.
if ~isfield(s,'decimFactor') || isempty(s.decimFactor), s.decimFactor=1;   end
if ~isfield(s,'dataTypeOut'),                           s.dataTypeOut='';  end
if ~isfield(s,'saveSource')  || isempty(s.saveSource),  s.saveSource=true; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Angiogram',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        %set file name to load data
        s.fName=char(fNames{fidx});
        reportFile(rep,fidx,s.fName);
        clearvars results source settings data

        % Everything that can be judged from the settings alone is judged here, BEFORE
        % the recording is opened: a Bio-Formats index of a large .cxd is seventeen
        % minutes, and nobody should pay it to be told about a typo.
        dF=checkDecimFactor(s.decimFactor);
        if ~isempty(s.dataTypeOut), checkOutClass(s.dataTypeOut); end

        reader=bfGetReader(s.fName);
        % The guard is what closes the recording when a refusal below throws.  It is an
        % ANONYMOUS handle capturing the reader by value: an onCleanup holding a nested
        % function handle would pin this workspace and never fire.  Re-assigning it on
        % the next file fires this one, which is why the explicit close is safe to
        % follow it - closeReader swallows a double close.
        guard=onCleanup(@() closeReader(reader));

        omeMeta=reader.getMetadataStore();
        sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
        sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
        sizeTRaw=double(omeMeta.getPlaneCount(0));
        startT=omeMeta.getImageAcquisitionDate(0);
        timeStamp=round(posixtime(datetime(string(startT),'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC')))*1000; %ms
        dataType=matlabClass(char(omeMeta.getPixelsType(0).toString()));
        sampling=double(omeMeta.getPixelsTimeIncrement(0).value()); %s

        %The whole shape of the product, before a single plane has been read.  Doing it
        %here is the point of the rewrite: an allocation this size fails at the end of
        %the read, and the answer the user needs - which decimation would have fitted -
        %is arithmetic that costs nothing and can be given before the wait.
        sizeT=floor(sizeTRaw/dF);
        if sizeT<1
            error('runIntensity:noOutputFrames', ...
                ['s.decimFactor=%d averages more frames than the %d in %s, so there ' ...
                 'would be no output frame; the largest value that yields one is %d.'], ...
                dF,sizeTRaw,shortName(s.fName),sizeTRaw);
        end
        % Asked again on the RESOLVED class, which is the recording's own when the
        % setting was left empty - an OME pixel type MATLAB cannot allocate is then
        % named here rather than surfacing out of zeros().
        outClass=s.dataTypeOut;
        if isempty(outClass), outClass=dataType; end
        outClass=checkOutClass(outClass);
        time=((0:1:(sizeT-1)).*sampling.*dF)';      % [T x 1] COLUMN, library-wide

        %The cube exists only when the frames are kept.  With saveSource off there is
        %nothing here to refuse: the read below holds one plane and one accumulator, so
        %an angiogram of a recording of any length costs the same.
        cubeBytes=sizeY*sizeX*sizeT*bytesPerSample(outClass);
        if s.saveSource
            checkCubeFits(cubeBytes,sizeY,sizeX,sizeTRaw,outClass,dF,s.fName);
            data=zeros(sizeY,sizeX,sizeT,outClass);
        end

        %What the report page and every later step read off the product.  They are
        %accumulated per output frame rather than taken from the cube at the end,
        %because with saveSource off there is no cube to take them from - and because a
        %min over a 129 GB array is a second pass over it.
        imgSum=zeros(sizeY,sizeX);                  % double: exact over 30 000 frames
        meanI=zeros(sizeT,1);
        minI=zeros(sizeT,1);
        maxI=zeros(sizeT,1);
        acc=zeros(sizeY,sizeX);                     % the decimation block, summed in place
        nPix=sizeY*sizeX;

        %read the recording, one plane at a time, and close it
        tRead=tic;  tNext=REPORTEVERYSEC;
        for t=1:1:sizeT
            acc(:)=0;
            for j=1:1:dF
                acc=acc+double(bfGetPlane(reader,(t-1).*dF+j));
            end
            %The frame is CAST FIRST and measured afterwards, so imgI and the traces
            %describe what is in the product rather than what it was rounded from.  With
            %an integer dataTypeOut on a decimated read those two differ, and the one a
            %reader can check against the file is this one.
            frame=cast(acc./dF,outClass);
            if s.saveSource, data(:,:,t)=frame; end
            fd=double(frame);
            imgSum=imgSum+fd;
            meanI(t)=sum(fd(:))./nPix;
            minI(t)=min(fd(:));
            maxI(t)=max(fd(:));

            if toc(tRead)>=tNext
                el=toc(tRead);
                rep.emit(sprintf('  read %d of %d frames, %.0f%%, about %s left', ...
                    t,sizeT,100*t/sizeT,hms(el*(sizeT-t)/t)));
                tNext=el+REPORTEVERYSEC;
            end
        end
        reader.close();

        results.time=time;
        results.timeStamp=timeStamp;
        results.imgI=imgSum./sizeT;                 % the angiogram
        results.meanI=meanI;
        results.minI=minI;
        results.maxI=maxI;
        settings.runIntensity=reportSettings(s);

        fh=reportFigure(rep,'intensity');
        tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl);
        imagesc(ax,results.imgI)
        clim(ax,[prctile(results.imgI(:),1),prctile(results.imgI(:),99)])
        colorbar(ax)
        axis(ax,'image')
        ax=nexttile(tl);
        semilogy(ax,time,meanI);
        hold(ax,'on')
        semilogy(ax,time,minI);
        semilogy(ax,time,maxI);
        hold(ax,'off')
        legend(ax,{'Mean','Min','Max'})
        ylabel(ax,'Intensity')
        xlabel(ax,'Time, s')
        %A log axis cannot show a zero, and a dark pixel gives one; a single output
        %frame gives a degenerate x range.  Both are ordinary recordings, so the limits
        %are made valid rather than left to throw out of a report page.
        xlim(ax,spanOf(time,[0 1]));
        ylim(ax,logSpanOf([minI;meanI;maxI]));
        reportSave(rep,fh,'intensity');     % reportSave titles the page

        % Save the settings and results
        reportWriting(rep);
        % This is the ENTRY step, so this is where a raw path becomes a product path -
        % and the only place that has to know a project may keep its results apart from
        % its recordings.  s.fName is NOT reassigned: it stays the raw recording.  With
        % no results folder set the name comes back verbatim.
        outBase=getResultsPath(s.fName,s);
        if s.saveSource
            source.data=data;
            source.time=time;
            save(strrep(outBase,'.cxd','_a_I_d.mat'),'source','-v7.3','-nocompression');
        end
        save(strrep(outBase,'.cxd','_a_I_r.mat'),'results','-v7.3','-nocompression');
        save(strrep(outBase,'.cxd','_a_I_s.mat'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
        %Both copies of the cube go before the next file allocates its own: a plain
        %re-assignment builds the new one while the old is still resident, which on a
        %batch of large recordings is twice the peak for nothing.
        clear data source
    end
end
reportClose(rep);

end

% =====================================================================
function dF = checkDecimFactor(v)
%checkDecimFactor  A whole number of frames, one or more.  Checked before the
%   recording is opened, because a bad value here is a typo in a settings file and the
%   user should not pay a Bio-Formats index to hear about it.
dF = double(v);
if ~isscalar(dF) || ~isfinite(dF) || dF<1 || mod(dF,1)~=0
    error('runIntensity:decimFactor', ...
        's.decimFactor must be a whole number of frames, 1 or more; it is %s.', ...
        mat2str(v));
end
end

% =====================================================================
function cls = matlabClass(omeType)
%matlabClass  The OME pixel type as a class zeros() will accept.  OME says 'float'
%   where MATLAB says 'single', and nothing else in the vocabulary differs.
cls = char(omeType);
if strcmp(cls,'float'), cls = 'single'; end
end

% =====================================================================
function cls = checkOutClass(cls)
%checkOutClass  The classes a frame may be stored as.  A wrong string here used to
%   surface as an error out of zeros() naming neither the setting nor the file.
cls = char(cls);
ok = {'single','double','uint8','int8','uint16','int16','uint32','int32'};
if ~any(strcmp(cls,ok))
    error('runIntensity:dataTypeOut', ...
        's.dataTypeOut must be one of %s, or empty to keep the recording''s own; it is ''%s''.', ...
        strjoin(ok,', '), cls);
end
end

% =====================================================================
function n = bytesPerSample(cls)
%bytesPerSample  Bytes one sample of a stored frame costs.
switch cls
    case {'uint8','int8'},   n = 1;
    case {'uint16','int16'}, n = 2;
    case {'single','uint32','int32'}, n = 4;
    otherwise, n = 8;                    % double
end
end

% =====================================================================
function checkCubeFits(cubeBytes,sizeY,sizeX,sizeTRaw,outClass,dF,fName)
%checkCubeFits  REFUSE BEFORE THE READ, and say what would have worked.
%
%   Half the free physical memory is the bound runContrastInternalCycle already applies
%   own largest array, and it is the same reasoning: the cube is not the only thing
%   resident while it is written.  The DIFFERENCE from that step is that this one
%   refuses instead of adjusting.  That step raises decimationSpace by itself
%   because a coarser sampling of a spatial average is still that average, so nothing
%   the user asked for changes.  s.decimFactor is not like that - it sets the output
%   frame rate, which is the scientific content of the product - so choosing it is the
%   user's and not this function's.
%
%   The number offered is the decimation that would fit the cube in the same bound,
%   and the second sentence names the other way out, which costs nothing at all.
[~,mem] = memory;
avail   = mem.PhysicalMemory.Available;
if cubeBytes <= avail/2, return; end

frameBytes = sizeY*sizeX*bytesPerSample(outClass);
tFit       = max(1,floor((avail/2)/frameBytes));
dFNeed     = max(dF+1,ceil(sizeTRaw/tFit));
error('runIntensity:cubeTooLarge', ...
    ['The %d frames of %s would be %.1f GB as %s and half the free memory is ' ...
     '%.1f GB.  s.decimFactor=%d would bring it to %.1f GB; s.saveSource=false ' ...
     'keeps the angiogram alone and needs no cube at all.'], ...
    floor(sizeTRaw/dF), shortName(fName), cubeBytes/2^30, outClass, ...
    avail/2/2^30, dFNeed, frameBytes*floor(sizeTRaw/dFNeed)/2^30);
end

% =====================================================================
function lim = spanOf(v,fallback)
%spanOf  [first last] of a vector, widened when the two are the same value.
lim = [v(1) v(end)];
if ~(lim(2)>lim(1)), lim = fallback; end
end

% =====================================================================
function lim = logSpanOf(v)
%logSpanOf  A y range a log axis can take: strictly positive and strictly increasing.
%   Handed EVERY trace on the axis rather than only the minima, because a dark corner
%   pixel makes every frame's minimum zero and the floor then has to come from the
%   traces that are not zero - on a recording already scaled into [0 1] a fixed
%   fall-back of 1 would put the whole plot below the axis.
a = min(v(v>0));
if isempty(a), a = 1; end
b = max(v);
if ~(b>a), b = a*10; end
lim = [a b];
end

% =====================================================================
function s = hms(sec)
%hms  Seconds as minutes and seconds, for a line an operator reads while waiting.
s = sprintf('%d min %02d s', floor(sec/60), round(mod(sec,60)));
end

% =====================================================================
function n = shortName(f)
%shortName  Name and extension - what the user sees on disk.
[~,b,e] = fileparts(char(f));
n = [b e];
end

% =====================================================================
function closeReader(reader)
try
    reader.close();
catch
end
end
%------------- END OF CODE --------------
