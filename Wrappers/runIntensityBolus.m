%runIntensityBolus - Split a .cxd bolus recording into a BOLUS cube and an angiogram
%
%   THE THIRD ENTRY STEP OF THE FLUORESCENCE BRANCH, and the one a tracer injection is
%   for.  It divides one recording into two spans: a BOLUS span kept at the full frame
%   rate, which is what runCTTH measures the transit times on, and an ANGIOGRAM span
%   averaged into a single picture, which is what the segmentation is drawn on.  The
%   product is '<stem>_b_I_d/_r/_s.mat' - stage 'b', product 'I', branch 'bolus'.
%
%   THE RECORDING IS OPENED ONCE, EVER.  bfGetReader is 16 min 43 s on the 60 GB
%   reference recording, so a second open is a second wait; every path below goes
%   through the one reader this step opens per file.  What the version this replaces did
%   instead was worse than a second open: with the spans left empty it loaded THE WHOLE
%   STACK into memory to draw the picker - 129 GB as single for that recording - so the
%   interactive mode could not be used on the recordings it exists for.
%
%   SO THE PICKER IS FED BY A STREAM, AND THE CUBE IS READ AFTERWARDS.  With s.fBolus
%   empty the step makes one pass holding ONE plane at a time and keeping only three
%   numbers per frame (mean, min, max), draws the picker off those, and then reads the
%   chosen bolus span - through the same open reader - into the only large array this
%   step ever allocates.  With s.fBolus set there is no first pass at all: the bolus
%   span is read and the angiogram span accumulated in one traversal, which is what the
%   old step already did and is kept.
%
%   THE BOLUS SPAN IS NEVER DECIMATED.  A transit time is measured off the shape of the
%   leading edge and the whole physiological signal on the reference recording is
%   0.3-0.5 s wide; at 100 Hz that is thirty samples, and there is nothing to spend.
%
%   KEEP 25-30 SECONDS OF BOLUS SPAN.  Measured (research-ctth.md sections 5.2-5.3): the
%   mean transit time survives a 13.6 s span to within -6.6 to +5.6 %, but the WIDTH of
%   the transit-time distribution does not - on 13.6 s a true spread of 0.50 s comes back
%   as -0.62 s, and it takes about 20 s of noiseless record for that to be exact.  The
%   real tail is still rising at 0.45 % of the step per second, so allow margin.  This
%   step therefore SAYS SO - in the picker, and once as a warning when a shorter span is
%   chosen - and changes nobody's existing span silently.
%
%   AND THE PRE-BOLUS PERIOD IS PART OF THE SPAN.  Every marker runCTTH measures is
%   referenced to the level before the injection, and getBolusMetrics refuses by name
%   (getBolusMetrics:baselineNotFlat) a span that starts after the bolus has begun -
%   because a span that starts mid-rise derives a "baseline" out of the upslope, scores
%   1.00 on every confidence factor, and is wrong in every marker.  So the span must
%   begin at least 1.5 s BEFORE the rise, and the picker asks for it in those words.
%
%   IT RECORDS THE PIXEL SIZE, and it is the only step that can.  Every downstream step
%   that works in micrometres - background removal today, topology and wall motion later
%   - needs to know how many micrometres a pixel is, and nothing about a product says so
%   otherwise: runBackgroundRemoval ships with NO default and refuses by name rather than
%   guess (2.5 um/px is measured for one rig and a default would be a silent wrong answer
%   at another magnification).  This step sees the recording, so this is where the
%   operator's answer is recorded: it goes into results.pixelSize and into the settings
%   sidecar, and in the workbench it is a SHARED key, so the number typed here is the
%   number background removal is offered.  It is empty when nobody answered - one field,
%   two current states, and empty is what the refusal downstream is about.
%
%   THE REPORT PAGE SAYS WHETHER THE FIRST PASS SEPARATED.  It draws the angiogram, the
%   bolus trace, and across the trace the pre-injection level, the plateau level and the
%   UNRETURNED FRACTION - how much of the rise was still there when the recording ended.
%   On the reference recording that is 0.752: three quarters of the bolus never came back
%   out, which is what a non-clearing intravascular tracer does and is the single glance
%   that tells a user their first pass did not separate.  It is why there is no
%   first-pass area and no rCBV in this branch (author's decision, 07-Aug-2026), and why
%   the plateau step became the volume-like marker instead.
%
% Syntax:
%    runIntensityBolus(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%               .fBolus     [f1 f2] frame span kept at the full frame rate.  Leave
%                           empty to mark it on the streamed trace.
%               .fAngio     [f1 f2] frame span averaged into the angiogram.  Leave
%                           empty to mark it, or (with fBolus preset) to average every
%                           frame after the bolus span.
%               .pixelSize  micrometres per pixel.  Empty is allowed and is recorded as
%                           empty; the steps that need it refuse rather than guess.
%    fNames - cell array of .cxd file paths.
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.  This step can block on the picker, so cancel is only checked
%    between files.
%
% Outputs:
%    (none) - writes, through getResultsPath, beside the recording or in the project's
%             results folder:
%               <stem>_b_I_d.mat  SOURCE   source.data [Y x X x T], source.time [T x 1]
%               <stem>_b_I_r.mat  RESULTS  time [T x 1], timeStamp, imgI, pixelSize
%               <stem>_b_I_s.mat  SETTINGS settings.runIntensityBolus
%               <stem>_rep_bolus.jpg       angiogram + the bolus trace and its plateau
%
% Example:
%    s.libraryFolder = libraryFolder;
%    s.fBolus    = [];            % marked on the plot
%    s.fAngio    = [];            % everything after the bolus
%    s.pixelSize = 2.5;
%    fNames = getFileNamesList(rootFolder,'*BB*.cxd');   % BB names a bolus recording
%    runIntensityBolus(s, fNames(:));
%
% Dependencies: the Bio-Formats MATLAB library (bfGetReader, bfGetPlane).
% See also: runCTTH, runIntensity, runIntensityInternalCycle, runBackgroundRemoval,
%           getResultsPath, wbFileModel, wbStepRegistry
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL
% s.fBolus=[];     %frames kept at the full frame rate. Empty marks them on the plot.
%                  %Keep 25-30 s of it, starting at least 1.5 s before the injection
% s.fAngio=[];     %frames averaged into the picture. Empty takes everything after the bolus
% s.pixelSize=2.5; %micrometres per pixel on this microscope

%------------- BEGIN CODE --------------
function runIntensityBolus(s,fNames)

REPORTEVERYSEC = 30;    % how often a read that is running long says where it is
WANTBOLUSSEC   = 25;    % the bolus span the width marker needs - see the header
WANTPREBOLUSSEC= 1.5;   % quiet baseline getBolusMetrics needs before the rise

if ~all( cellfun(@(f) isempty(f) || contains(f,'.cxd'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".cxd".');
end

% The spans default to empty, i.e. "ask me".  pixelSize defaults to empty too and that
% is NOT the same kind of default: it is the recorded absence of an answer, and the
% steps that need micrometres refuse on it by name rather than assume a magnification.
if ~isfield(s,'fBolus'),    s.fBolus=[];    end
if ~isfield(s,'fAngio'),    s.fAngio=[];    end
if ~isfield(s,'pixelSize'), s.pixelSize=[]; end
checkPixelSize(s.pixelSize);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks are
% resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  This step can block
% on the span picker, so cancel is only checked between files.
rep=reportOpen(s,'Bolus',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end          % cooperative cancel between files
    if ~isempty(fNames{fidx})
        %THE SPANS ARE RESOLVED PER FILE FROM THE UNRESOLVED SETTING, never from the
        %previous file's answer.  s itself is not written back into, so marking a span
        %on file 1 cannot silently become file 2's span - which is what a shared s
        %mutated in the loop would have done.
        sf=s;
        sf.fName=char(fNames{fidx});
        reportFile(rep,fidx,sf.fName);
        clearvars results source settings data

        reader=bfGetReader(sf.fName);
        % The guard is what closes the recording when a refusal below throws.  It is an
        % ANONYMOUS handle capturing the reader by value: an onCleanup holding a nested
        % function handle would pin this workspace and never fire.
        guard=onCleanup(@() closeReader(reader));

        omeMeta=reader.getMetadataStore();
        sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
        sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
        sizeT=double(omeMeta.getPlaneCount(0));
        startT=omeMeta.getImageAcquisitionDate(0);
        timeStamp=round(posixtime(datetime(string(startT), ...
            'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC')))*1000; %ms
        dataType=matlabClass(char(omeMeta.getPixelsType(0).toString()));
        sampling=double(omeMeta.getPixelsTimeIncrement(0).value()); %s
        timeAll=((0:1:(sizeT-1))*sampling)';         % [T x 1] COLUMN, library-wide

        % ---- the spans -------------------------------------------------------
        if isempty(sf.fBolus)
            %PASS 1: the picker's trace, streamed.  One plane resident, three numbers
            %kept per frame - so the interactive mode costs the same on a 60 GB
            %recording as on a 2 GB one.
            [mn,mi,mx]=streamTrace(reader,sizeT,rep,REPORTEVERYSEC);
            [sf.fBolus,sf.fAngio]=pickSpans(timeAll,mn,mi,mx,sizeT,sampling, ...
                WANTBOLUSSEC,WANTPREBOLUSSEC,sf.fName);
        else
            sf.fBolus=checkSpan(sf.fBolus,sizeT,'fBolus',sf.fName);
        end
        if ~isempty(sf.fAngio)
            sf.fAngio=checkSpan(sf.fAngio,sizeT,'fAngio',sf.fName);
        end
        nBolus=sf.fBolus(2)-sf.fBolus(1)+1;
        warnIfSpanShort(nBolus*sampling,WANTBOLUSSEC,sf.fName);

        %REFUSED BEFORE THE READ, and told which span would have fitted.  An allocation
        %this size fails at the end of a twenty-minute traversal otherwise.
        checkCubeFits(sizeY,sizeX,nBolus,dataType,sampling,sf.fName);

        % ---- PASS 2: the bolus cube, and the angiogram accumulated beside it ---
        %ONE traversal of the frames that are wanted, in ascending plane order, through
        %the reader that is already open.  The angiogram is accumulated in DOUBLE - it
        %is a mean over thousands of frames and an integer accumulator would saturate.
        [aFrom,aTo]=angioSpan(sf.fAngio,sf.fBolus,sizeT,sf.fName);
        data=zeros(sizeY,sizeX,nBolus,dataType);
        imgI=zeros(sizeY,sizeX);
        nAngio=max(aTo-aFrom+1,0);
        want=unique([sf.fBolus(1):sf.fBolus(2), aFrom:aTo]);
        tRead=tic; tNext=REPORTEVERYSEC;
        for k=1:1:numel(want)
            t=want(k);
            plane=bfGetPlane(reader,t);
            if t>=sf.fBolus(1) && t<=sf.fBolus(2)
                data(:,:,t-sf.fBolus(1)+1)=plane;
            end
            if nAngio>0 && t>=aFrom && t<=aTo
                imgI=imgI+double(plane)./nAngio;
            end
            if toc(tRead)>=tNext
                el=toc(tRead);
                rep.emit(sprintf('  read %d of %d frames, %.0f%%, about %s left', ...
                    k,numel(want),100*k/numel(want),hms(el*(numel(want)-k)/k)));
                tNext=el+REPORTEVERYSEC;
            end
        end
        reader.close();

        %crop and rezero the time vector, and shift the absolute stamp with it
        time=timeAll(sf.fBolus(1):sf.fBolus(2));
        timeStamp=timeStamp+time(1)*1000;            % time is in s, timeStamp in ms
        time=time-time(1);

        results.time=time;
        results.timeStamp=timeStamp;
        results.imgI=imgI;
        %THE ONE FIELD THAT CLOSES A QUESTION NO LATER STEP CAN ANSWER FOR ITSELF.  See
        %the header: it is empty when the operator did not say, and empty is what the
        %micrometre steps refuse on.
        results.pixelSize=sf.pixelSize;
        settings.runIntensityBolus=reportSettings(sf);

        drawBolusPage(rep,results,data,time);

        %save the settings and results
        reportWriting(rep);
        source.data=data;
        source.time=time;
        %This is the ENTRY step, so this is where a raw path becomes a product path, and
        %the only place that has to know a project may keep its results apart from its
        %recordings.  sf.fName is NOT reassigned - it stays the raw recording.  With no
        %results folder set the name comes back verbatim.
        outBase=getResultsPath(sf.fName,sf);
        save(strrep(outBase,'.cxd','_b_I_d.mat'),'source','-v7.3','-nocompression');
        save(strrep(outBase,'.cxd','_b_I_r.mat'),'results','-v7.3','-nocompression');
        save(strrep(outBase,'.cxd','_b_I_s.mat'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
        %Both copies of the cube go before the next file allocates its own: a plain
        %re-assignment builds the new one while the old is still resident.
        clear data source
    end
end
reportClose(rep);

end

% =====================================================================
function [mn,mi,mx]=streamTrace(reader,sizeT,rep,everySec)
%streamTrace  The picker's three curves, ONE PLANE AT A TIME.
%   This is the whole reason the interactive mode is usable on a large recording: the
%   picker needs a trace over every frame, not the frames themselves, and a trace is
%   three numbers each.  Nothing here is kept but 3*sizeT doubles.
mn=zeros(sizeT,1); mi=zeros(sizeT,1); mx=zeros(sizeT,1);
tRead=tic; tNext=everySec;
for t=1:1:sizeT
    fd=double(bfGetPlane(reader,t));
    mn(t)=mean(fd(:)); mi(t)=min(fd(:)); mx(t)=max(fd(:));
    if toc(tRead)>=tNext
        el=toc(tRead);
        rep.emit(sprintf('  scanning %d of %d frames, %.0f%%, about %s left', ...
            t,sizeT,100*t/sizeT,hms(el*(sizeT-t)/t)));
        tNext=el+everySec;
    end
end
end

% =====================================================================
function [fBolus,fAngio]=pickSpans(timeAll,mn,mi,mx,sizeT,sampling,wantSec,preSec,fName)
%pickSpans  Mark the two spans on the streamed trace.
%
%   THE TITLES ARE THE STEP'S ONLY CHANCE TO SAY WHAT A GOOD SPAN IS, so they say it in
%   seconds and in the two ways it goes wrong: too short and the width of the transit
%   distribution cannot be measured, started too late and every marker is referenced to
%   a piece of the upslope.  Both are measured findings, not preferences.
h=figure;
h.WindowState='Maximize';
semilogy(timeAll,mn);
hold on
semilogy(timeAll,mi);
semilogy(timeAll,mx);
hold off
legend({'Mean','Min','Max'})
ylabel('Intensity')
xlabel('Time, s')
xlim([timeAll(1),timeAll(end)]);
ylim(logSpanOf([mi;mn;mx]));

[~,stem,ext]=fileparts(char(fName));
title(sprintf(['%s%s   Mark the BOLUS: drag a box over it, adjust, then double-click. ' ...
    'Keep about %g s of it, starting at least %g s BEFORE the rise'], ...
    stem,ext,wantSec,preSec));
rB=drawrectangle(gca,'Color','r','Label','Bolus','FaceAlpha',0.15);
wait(rB);
fBolus=boxToSpan(rB,sampling,sizeT,[1 sizeT]);

title(sprintf(['%s%s   Mark the ANGIOGRAM: the frames averaged into the picture the ' ...
    'vessels are found on'],stem,ext));
rA=drawrectangle(gca,'Color','g','Label','Angio','FaceAlpha',0.15);
wait(rA);
fAngio=boxToSpan(rA,sampling,sizeT,[1 sizeT]);
close(h);
end

% =====================================================================
function span=boxToSpan(r,sampling,sizeT,dflt)
%boxToSpan  The horizontal extent of a drawn box, as a frame span.
%   THE PICKER'S AXIS IS IN SECONDS and the span is in FRAMES, so the conversion happens
%   here and in one place.  The version this replaces plotted against the frame index,
%   which meant the operator judged a 25 s recommendation off an axis that did not show
%   seconds.  A box that was never drawn, or was cancelled, keeps everything.
if ~isvalid(r) || isempty(r.Position)
    span=dflt; return
end
t1=r.Position(1); t2=r.Position(1)+r.Position(3);
span=sort(round([t1 t2]./sampling)+1);
span=min(max(span,1),sizeT);
if span(2)<=span(1), span=dflt; end
end

% =====================================================================
function [aFrom,aTo]=angioSpan(fAngio,fBolus,sizeT,fName)
%angioSpan  The frames the picture is averaged over.
%   Empty means EVERYTHING AFTER THE BOLUS, which is what a bolus protocol usually
%   wants: by then the dye has filled and the vessels are at their most visible.  A
%   recording that stops with the bolus has nothing after it and says so.
if ~isempty(fAngio)
    aFrom=fAngio(1); aTo=fAngio(2); return
end
aFrom=fBolus(2)+1; aTo=sizeT;
if aTo<aFrom
    warning('runIntensityBolus:noAngiogramFrames', ...
        ['%s has no frames after the bolus span, so the picture the vessels are ' ...
         'found on would be blank.  Set s.fAngio to a span inside the recording.'], ...
        shortName(fName));
    aFrom=1; aTo=0;
end
end

% =====================================================================
function span=checkSpan(v,sizeT,name,fName)
%checkSpan  A pre-set span, checked against the recording BEFORE any plane is read.
span=double(v(:))';
if numel(span)~=2 || any(~isfinite(span)) || span(1)>span(2)
    error('runIntensityBolus:badSpan', ...
        's.%s must be [firstFrame lastFrame] with the first no later than the last; it is %s.', ...
        name,mat2str(v));
end
span=round(span);
if span(1)>sizeT
    error('runIntensityBolus:spanOutsideRecording', ...
        's.%s=[%d %d] starts after the %d frames of %s.', ...
        name,span(1),span(2),sizeT,shortName(fName));
end
if span(2)>sizeT
    warning('runIntensityBolus:spanClamped', ...
        's.%s(2)=%d exceeds the %d frames of %s; clamped.', ...
        name,span(2),sizeT,shortName(fName));
    span(2)=sizeT;
end
span(1)=max(span(1),1);
end

% =====================================================================
function warnIfSpanShort(spanSec,wantSec,fName)
%warnIfSpanShort  Say once that the width marker will be below its floor.
%
%   IT IS A WARNING AND IT CHANGES NOTHING.  Extending somebody's span by itself would
%   make this step's output depend on a number nobody typed, and shortening the analysis
%   is not this step's decision either.  What a short span costs is specific and
%   measured, so it is said in those terms: the transit-time SPREAD is what goes, and the
%   mean transit time survives.
if spanSec>=wantSec, return, end
warning('runIntensityBolus:shortSpan', ...
    ['The bolus span of %s is %.1f s.  The mean transit time survives that, but the ' ...
     'WIDTH of the transit-time distribution needs about %g s and will come back ' ...
     'below the floor this recording can resolve.'], ...
    shortName(fName),spanSec,wantSec);
end

% =====================================================================
function checkPixelSize(v)
%checkPixelSize  Empty is allowed and is an answer; a nonsense number is not.
if isempty(v), return, end
if ~isscalar(v) || ~isfinite(v) || v<=0
    error('runIntensityBolus:pixelSize', ...
        's.pixelSize must be a positive number of micrometres per pixel, or empty; it is %s.', ...
        mat2str(v));
end
end

% =====================================================================
function checkCubeFits(sizeY,sizeX,nBolus,cls,sampling,fName)
%checkCubeFits  REFUSE BEFORE THE READ, and say what would have fitted.
%   Half the free physical memory, which is the bound the other entry steps apply, and
%   for the same reason: the cube is not the only thing resident while it is written.
%   The number offered is a span in SECONDS, because that is what the operator draws.
frameBytes=sizeY*sizeX*bytesPerSample(cls);
cubeBytes =frameBytes*nBolus;
[~,mem]=memory;
avail=mem.PhysicalMemory.Available;
if cubeBytes<=avail/2, return; end
tFit=max(1,floor((avail/2)/frameBytes));
error('runIntensityBolus:cubeTooLarge', ...
    ['The %d frames of the bolus span of %s would be %.1f GB as %s and half the free ' ...
     'memory is %.1f GB.  A span of %d frames (%.1f s) would fit.'], ...
    nBolus,shortName(fName),cubeBytes/2^30,cls,avail/2/2^30,tFit,tFit*sampling);
end

% =====================================================================
function drawBolusPage(rep,results,data,time)
%drawBolusPage  '_rep_bolus.jpg' - the picture, the trace, and WHETHER THE FIRST PASS
%   SEPARATED.
%
%   THE UNRETURNED FRACTION IS THE POINT OF THIS PAGE.  (Fn - Bl) / (Pk - Bl): how much
%   of the rise was still there when the recording ended.  A tracer that clears returns
%   towards its baseline and the fraction is small; a non-clearing intravascular tracer
%   settles at a systemic level and it is large - 0.752 on the reference recording, which
%   is why this branch has no first-pass area and reads the PLATEAU as the volume-like
%   quantity instead.  A user who sees 0.75 here knows their first pass did not separate
%   before any transit time is computed.
%
%   ITS WINDOWS ARE THIS PAGE'S OWN AND ARE DELIBERATELY NOT runCTTH's.  This is an entry
%   step: it has no arterial input to derive a pre-bolus window from and must not grow a
%   settings group to acquire one.  The first tenth of the span and the last second are
%   enough for a page whose job is to show an order of magnitude.
%
%   A PAGE IS A BY-PRODUCT AND MUST NEVER KILL A RUN.  The product is already complete
%   when this is attempted.
fh=[];
try
    mn=squeeze(mean(data,[1,2]));
    mi=squeeze(min(data,[],[1,2]));
    mx=squeeze(max(data,[],[1,2]));
    nBl=max(3,round(0.1*numel(mn)));
    Bl=mean(double(mn(1:nBl)));
    Fn=mean(double(mn(time>=time(end)-1)));
    Pk=max(double(mn));
    unret=(Fn-Bl)/max(Pk-Bl,eps);

    fh=reportFigure(rep,'bolus');
    tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl);
    imagesc(ax,results.imgI)
    clim(ax,[prctile(results.imgI(:),1),prctile(results.imgI(:),99)])
    colorbar(ax)
    axis(ax,'image')
    title(ax,'The picture the vessels are found on')

    ax=nexttile(tl);
    semilogy(ax,time,mn);
    hold(ax,'on')
    semilogy(ax,time,mi);
    semilogy(ax,time,mx);
    %the two levels the fraction is read between - a number with no line beside it is a
    %number nobody can check
    yline(ax,Fn,'-','Label','plateau','Color',[0.15 0.45 0.75],'LineWidth',1.5, ...
        'LabelHorizontalAlignment','left');
    yline(ax,Bl,'--','Label','before the injection','Color',[0.45 0.45 0.45], ...
        'LabelHorizontalAlignment','left');
    hold(ax,'off')
    legend(ax,{'Mean','Min','Max'})
    ylabel(ax,'Intensity')
    xlabel(ax,'Time, s')
    xlim(ax,spanOf(time,[0 1]));
    ylim(ax,logSpanOf([mi(:);mn(:);mx(:);Bl;Fn]));
    title(ax,sprintf('%.0f s of bolus - %.3f of the rise had not come back by the end', ...
        time(end)-time(1),unret))

    reportSave(rep,fh,'bolus');         % reportSave titles the page and deletes fh
    fh=[];
catch
end
if ~isempty(fh) && isgraphics(fh), delete(fh); end
end

% =====================================================================
function cls=matlabClass(omeType)
%matlabClass  The OME pixel type as a class zeros() will accept.  OME says 'float'
%   where MATLAB says 'single', and nothing else in the vocabulary differs.
cls=char(omeType);
if strcmp(cls,'float'), cls='single'; end
end

% =====================================================================
function n=bytesPerSample(cls)
%bytesPerSample  Bytes one sample of a stored frame costs.
switch cls
    case {'uint8','int8'},            n=1;
    case {'uint16','int16'},          n=2;
    case {'single','uint32','int32'}, n=4;
    otherwise,                        n=8;      % double
end
end

% =====================================================================
function lim=spanOf(v,fallback)
%spanOf  [first last] of a vector, widened when the two are the same value.
lim=[v(1) v(end)];
if ~(lim(2)>lim(1)), lim=fallback; end
end

% =====================================================================
function lim=logSpanOf(v)
%logSpanOf  A y range a log axis can take: strictly positive and strictly increasing.
%   Handed EVERY trace on the axis rather than only the minima, because a dark corner
%   pixel makes every frame's minimum zero and the floor then has to come from the
%   traces that are not zero.
v=double(v(:));
a=min(v(v>0));
if isempty(a), a=1; end
b=max(v);
if ~(b>a), b=a*10; end
lim=[a b];
end

% =====================================================================
function s=hms(sec)
%hms  Seconds as minutes and seconds, for a line an operator reads while waiting.
s=sprintf('%d min %02d s',floor(sec/60),round(mod(sec,60)));
end

% =====================================================================
function n=shortName(f)
%shortName  Name and extension - what the user sees on disk.
[~,b,e]=fileparts(char(f));
n=[b e];
end

% =====================================================================
function closeReader(reader)
try
    reader.close();
catch
end
end
%------------- END OF CODE --------------
