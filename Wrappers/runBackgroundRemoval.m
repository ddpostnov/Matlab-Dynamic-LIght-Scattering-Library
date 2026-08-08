%runBackgroundRemoval - Take the glow off the picture of the vessels.
%
%   WHAT THE STEP IS FOR (author, 2026-08-07, [D15]): "Its role is to produce the
%   clearest image of the vessels for the steps that work with vessel structure,
%   topology quantification and wall dynamics analysis.  It is applied to imgI in the
%   Bolus or PIV run, on which segmentation is run, but it is not applied to their
%   source data.  Therefore among the options there should be a choice to apply to
%   imgI only or to source data as well."
%
%   SO IT CLEANS results.imgI AND IT DOES NOT TOUCH THE FRAMES UNLESS IT IS TOLD TO.
%   That is not a compromise between two consumers, it is the observation that they
%   read different members of the same product: runSegmentation on the intensity path
%   reads results.imgI and nothing else, and runCTTH reads results.sData and the cube.
%   Clean the first and the second cannot change - the safety property is a property
%   of the data flow rather than a rule in a document.
%
%   AND IT STORES WHAT IT REMOVED.  results.imgBackground sits beside the cleaned
%   picture, so the operation is EXACTLY REVERSIBLE (imgI_raw = imgI + imgBackground),
%   nothing is destroyed, and a second run at a different radius restores the raw
%   picture first rather than cleaning an already-cleaned one.  It costs a few
%   megabytes, and it is what makes the report page's second panel possible at all.
%
%   s.applyToSource CLEANS EVERY FRAME TOO, AND IT IS FALSE BY DEFAULT.  It exists for
%   wall dynamics, the one named consumer that reads frames rather than the picture -
%   runDynamicSegmentation fits a wall on every frame out of source.data, and so does
%   the motion-enhancement kymograph, and cleaning imgI does nothing for either of
%   them.  It is free for PIV: R1 measured every spatial removal neutral to three
%   significant figures for the lumen texture PIV reads.
%
%   IT IS EXPENSIVE AND IT IS NOT INTERLOCKED [D16].  Cleaning the cube widens uint16
%   to single - +60 GB on the 60 GB reference recording - and R1 measured what happens
%   when a transit time is then read off it: ARRIVAL SHIFTED BY +2171 ms (217 frames),
%   rCBV's segment-to-segment CV 95%, and no calibre is safe (30-50 um vessels have
%   12.5% of segments within tolerance).  IT LOOKS COMPLETELY NORMAL.  The author
%   ruled that the combination is nevertheless the user's to make: "If a user chooses
%   apply to source - they can still run CTTH.  It is their choice.  What is correct
%   is another question."  So nothing is blocked.  What ships instead is the part that
%   costs the user nothing: the choice is RECORDED in
%   settings.runBackgroundRemoval.applyToSource, so a transit time measured on a
%   cleaned cube is auditable afterwards rather than indistinguishable, and runCTTH
%   says so once per file as an ordinary warning.
%
%   IT RUNS ON ANY INTENSITY PRODUCT, on all three branches ('_a_I', '_c_I', '_b_I'),
%   because all three entry steps write results.imgI.  Nothing requires it and it
%   requires nothing but an entry step, so it may be run before or after Regions - the
%   ORDER changes nothing about the answer, since regions are geometry and this step
%   does not move them.  AFTER is the recommendation: the operator draws the regions on
%   the picture they recognise, and the cleaning is then a step whose effect they can
%   see on its own report page.
%
%   A PRODUCT WITH NO FRAMES IS FINE, AND ONLY THE OPTION IS NOT.  runIntensity's
%   s.saveSource=false writes no '_d' member at all; the picture is still there, so
%   the step still does its job, and only s.applyToSource is refused by name.  (S1's
%   contract listed background removal among the steps such a product cannot carry.
%   Under [D15] that is only half true and this is the half that changed.)
%
% Syntax:
%    runBackgroundRemoval(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%               .pixelSize     micrometres per pixel.  REQUIRED, no default - every
%                              length below is in micrometres and this is the only
%                              number that turns them into pixels.
%               .radiusUm      opening radius, micrometres.  Default 200.
%               .medianUm      median-filter window, micrometres.  Default 37.5.
%               .applyToSource logical; clean every frame as well.  Default FALSE.
%    fNames - cell array of *_I_r.mat paths - the RESULTS member of an intensity
%             product.  The SOURCE cube and the SETTINGS are named from it.
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.
%
% Outputs:
%    (none) - writes, per file:
%      <name>_I_r.mat   RESULTS  imgI cleaned in place; imgBackground added
%      <name>_I_d.mat   SOURCE   only when s.applyToSource - every frame cleaned, and
%                                the cube widened to single if it was an integer class
%      <name>_I_s.mat   SETTINGS settings.runBackgroundRemoval
%      <name>_rep_background.jpg   the picture before, the glow, the picture after,
%                                  and the measurement that says which of the two it
%                                  actually removed
%
% Example:
%    s.libraryFolder = libraryFolder;
%    s.pixelSize     = 2.5;      % micrometres per pixel - measure it for your rig
%    s.radiusUm      = 200;
%    s.medianUm      = 37.5;
%    s.applyToSource = false;
%    fNames = getFileNamesList(rootFolder,'*_b_I_r.mat');
%    runBackgroundRemoval(s,fNames);
%
% Dependencies: getBackground, getHaloProfile (Core); the Core Reporting module;
%               MATLAB's Image Processing Toolbox.
% See also: getBackground, getHaloProfile, runSegmentation, runCTTH, runIntensity
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL
% s.pixelSize=2.5;   % micrometres per pixel.  REQUIRED - measured for the rig, not
%                    % guessed: every length below scales with it
% s.radiusUm=200;    % the width of the glow around the vessels.  Measured: 200 um
%                    % leaves 0.2% of the halo where 377.5 um leaves 10.9%.  Set it
%                    % above the widest vessel you want to keep
% s.medianUm=37.5;   % smooths speckle before the glow is measured; keep it below the
%                    % narrowest vessel you care about
% s.applyToSource=false; % clean every frame too, not just the picture - see the
%                    % header before turning it on

%------------- BEGIN CODE --------------
function runBackgroundRemoval(s,fNames)

REPORTEVERYSEC = 30;        % how often a long cube pass says where it has got to
BLOCKBYTES     = 256*2^20;  % the frame block the cube is cleaned in

if ~all( cellfun(@(f) isempty(f) || contains(f,'_I_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_I_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_I_r.mat'').']);
end

% Judged from the settings alone, BEFORE the first product is opened: a missing pixel
% size is a settings mistake and nobody should pay a multi-gigabyte load to hear about
% it.  There is no default, and there deliberately is not one - 2.5 um/px is measured
% for one rig and every micrometre in this step scales with it, so a default would be
% a silent wrong answer on any other magnification rather than a refusal.
if ~isfield(s,'pixelSize') || isempty(s.pixelSize)
    error('runBackgroundRemoval:pixelSizeNotSet', ...
        ['This step works in micrometres, so it needs to know how many micrometres ' ...
         'a pixel is on this recording.  Set it and run again.']);
end
if ~isfield(s,'radiusUm')      || isempty(s.radiusUm),      s.radiusUm=200;       end
if ~isfield(s,'medianUm')      || isempty(s.medianUm),      s.medianUm=37.5;      end
if ~isfield(s,'applyToSource') || isempty(s.applyToSource), s.applyToSource=false; end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Background removal',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end          % cooperative cancel between files
    if ~isempty(fNames{fidx})

    s.fName=fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results source settings
    load(getProductPath(s.fName,'s'),'settings');
    load(s.fName,'results');

    % --- THE PICTURE THIS STEP IS ABOUT, restored if it has already been cleaned ---
    % imgBackground present and imgBackground absent are two CURRENT states of the
    % product - it has been cleaned, or it has not - and one reader covers both by
    % design.  Restoring is what makes a second run at a different radius mean the new
    % radius rather than the new radius on top of the old one.
    imgRaw=double(results.imgI);
    if isfield(results,'imgBackground') && ~isempty(results.imgBackground)
        imgRaw=imgRaw+double(results.imgBackground);
    end

    % --- the cube, only when it is going to be cleaned ---------------------------
    % With s.applyToSource false the frames are not read at all: on the reference
    % recording that is 60 GB not loaded to compute one image.
    if s.applyToSource
        cleanedAlready=isfield(settings,'runBackgroundRemoval') && ...
            isfield(settings.runBackgroundRemoval,'applyToSource') && ...
            isequal(settings.runBackgroundRemoval.applyToSource,true);
        if cleanedAlready
            % The picture keeps its estimate and can be restored; the frames do not,
            % because there is no per-frame background product by decision.  So a
            % second pass over an already-cleaned cube would subtract the glow twice
            % and there would be no way back short of re-reading the recording.
            error('runBackgroundRemoval:cubeAlreadyCleaned', ...
                ['The frames of %s have already had the glow taken off them, and ' ...
                 'unlike the picture they cannot be put back.  Read the recording ' ...
                 'again if a different setting is wanted.'], shortName(s.fName));
        end
        dName=getProductPath(s.fName,'d');
        if ~isfile(dName)
            error('runBackgroundRemoval:noFrames', ...
                ['%s was written without its frames, so there are none to clean.  ' ...
                 'The picture of the vasculature is still cleaned - turn off ' ...
                 'cleaning the recording itself, or read the recording again ' ...
                 'keeping the frames.'], shortName(s.fName));
        end
        load(dName,'source');
    end

    % --- the estimate, and the picture it comes off -------------------------------
    bg=getBackground(imgRaw,s);
    results.imgBackground=bg;               % single, counts - 4 to 7 MB, and the
                                            % operation is exactly reversible with it
    results.imgI=imgRaw-double(bg);

    % --- the frames, in blocks -----------------------------------------------------
    if s.applyToSource
        source.data=cleanCube(source.data,s,rep,BLOCKBYTES,REPORTEVERYSEC);
    end

    % --- the report page ------------------------------------------------------------
    writeBackgroundPage(rep,s,imgRaw,bg,results.imgI);

    settings.runBackgroundRemoval=reportSettings(s);
    reportWriting(rep);
    save(s.fName,'results','-v7.3','-nocompression');
    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
    if s.applyToSource
        save(getProductPath(s.fName,'d'),'source','-v7.3','-nocompression');
    end
    reportSaved(rep);
    end
end
reportClose(rep);
end

% =====================================================================
function data=cleanCube(data,s,rep,blockBytes,everySec)
%cleanCube  Take the glow off every frame, a block at a time.
%
%   NOTHING FULL-SIZE IS ALLOCATED BESIDE THE CUBE.  The obvious form -
%   single(data) - getBackground(data,...) - builds a SECOND cube the size of the
%   first to hold the estimate, which on the reference recording is another 120 GB for
%   an array that is thrown away one frame later.  The estimate is made and subtracted
%   a block at a time instead, so the peak is the cube plus one block.
%
%   THE WIDENING IS THE COST, AND IT IS REFUSED RATHER THAN DISCOVERED.  A cleaned
%   frame has negative values in it, so an integer cube has to become single: +60 GB
%   on the reference recording.  That allocation is checked against free memory before
%   it is attempted, the way runIntensity checks its own.
[nR,nC,nT]=size(data,1,2,3);
if ~isa(data,'single')
    checkCleanCubeFits(nR,nC,nT,class(data));
    data=single(data);
end

blk=max(1,floor(blockBytes/max(nR*nC*4,1)));
tRun=tic;  tNext=everySec;
for k0=1:blk:nT
    idx=k0:min(k0+blk-1,nT);
    data(:,:,idx)=data(:,:,idx)-getBackground(data(:,:,idx),s);
    if toc(tRun)>=tNext
        el=toc(tRun);  done=idx(end);
        rep.emit(sprintf('  cleaned %d of %d frames, %.0f%%, about %s left', ...
            done,nT,100*done/nT,hms(el*(nT-done)/done)));
        tNext=el+everySec;
    end
end
end

% =====================================================================
function checkCleanCubeFits(nR,nC,nT,fromClass)
%checkCleanCubeFits  The single cube the widening needs, against half the free memory.
%   Half is the bound runContrastInternalCycle and runIntensity both apply to their own
%   largest array, and the reasoning is the same: the cube is not the only thing
%   resident while it is written.  The raw cube is already loaded when this is asked,
%   so what is being sized is the ADDITIONAL allocation.
bytes=nR*nC*nT*4;
[~,mem]=memory;
avail=mem.PhysicalMemory.Available;
if bytes<=avail/2, return; end
error('runBackgroundRemoval:cubeTooLarge', ...
    ['A cleaned frame has negative values in it, so the %d frames would have to be ' ...
     'widened from %s to single - %.1f GB, against %.1f GB of free memory to spend.  ' ...
     'Clean the picture only, or work on a shorter recording.'], ...
    nT, fromClass, bytes/2^30, avail/2/2^30);
end

% =====================================================================
function writeBackgroundPage(rep,s,imgRaw,bg,imgClean)
%writeBackgroundPage  <name>_rep_background.jpg - and the fourth tile is the point.
%
%   A BACKGROUND REMOVAL ALWAYS MAKES A PRETTIER PICTURE.  Three panels of before,
%   estimate and after are therefore not evidence of anything: the same three panels
%   look just as convincing when what was removed was the vasculature.  The fourth
%   tile is the measurement that separates the two - the light around the vessels
%   against distance from them, before and after, with what is left of it in the
%   title, beside what is left of the vessels themselves.
%
%   EACH IMAGE ON ITS OWN COLOUR LIMITS.  The estimate is smooth and vessel-free and
%   on the picture's limits it is a flat grey rectangle; the cleaned picture on the
%   raw picture's limits is nearly black.  On its own limits the cleaned picture is
%   the case FOR the step: it is where a small vessel and a large one are legible at
%   the same time.
fh=reportFigure(rep,'background');
try
    tl=tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

    ax=nexttile(tl); showImage(ax,imgRaw);            title(ax,'Before');
    ax=nexttile(tl); showImage(ax,double(bg));        title(ax,'The glow that was removed');
    ax=nexttile(tl); showImage(ax,imgClean);          title(ax,'After');

    ax=nexttile(tl);
    [p0,ref]=getHaloProfile(imgRaw,s.pixelSize);
    p1=getHaloProfile(imgClean,ref);
    if p0.ok && p1.ok && p0.grad~=0
        plot(ax,p0.d,p0.value-p0.floor,'LineWidth',1.2); hold(ax,'on')
        plot(ax,p1.d,p1.value-p1.floor,'LineWidth',1.2); hold(ax,'off')
        xlim(ax,[0 p0.dMax]);
        legend(ax,{'Before','After'},'Location','northeast');
        % THE GLOW IS A PERCENTAGE AND THE VESSELS ARE NOT, and that is not an
        % inconsistency.  grad is an excess over the far field, so it is floor-free
        % and a ratio of two of them means what it says.  The in-vessel LEVEL of a
        % mean picture is not: it carries the camera offset and the field-wide
        % pedestal that no spatial filter can see, so dividing one by the other would
        % read as a vessel loss that is mostly pedestal.  The two levels are quoted
        % instead, in counts, and the retained-signal FRACTION is measured where it
        % is defined - on the dye alone, in the research document and in the test.
        title(ax,sprintf('%.0f%% of the glow left  ·  vessels %.0f → %.0f counts', ...
            100*p1.grad/p0.grad, p0.inVessel, p1.inVessel));
        xlabel(ax,'Distance from the nearest vessel, µm');
        ylabel(ax,'Intensity above the far field');
    else
        % Too little field between the vessels to measure a profile on.  Said out
        % loud on the page rather than drawn as an empty axes, because a blank
        % fourth tile reads as "nothing to report".
        text(ax,0.5,0.5,{'Not enough tissue between the vessels', ...
            'to measure the glow on this field'},'HorizontalAlignment','center');
        axis(ax,'off');
        title(ax,'The glow was not measured');
    end
catch ME
    delete(fh); rethrow(ME);        % reportSave never runs, so the figure goes here
end
reportSave(rep,fh,'background');    % reportSave titles the page and deletes fh
end

% =====================================================================
function showImage(ax,img)
%showImage  One panel, on ITS OWN percentile limits (see writeBackgroundPage).
imagesc(ax,img);
lo=prctile(img(:),1); hi=prctile(img(:),99);
if hi>lo, clim(ax,[lo hi]); end
axis(ax,'image'); colorbar(ax);
end

% =====================================================================
function s=hms(sec)
%hms  Seconds as minutes and seconds, for a line an operator reads while waiting.
s=sprintf('%d min %02d s', floor(sec/60), round(mod(sec,60)));
end

% =====================================================================
function n=shortName(f)
%shortName  Name and extension - what the user sees on disk.
[~,b,e]=fileparts(char(f));
n=[b e];
end
%------------- END OF CODE --------------
