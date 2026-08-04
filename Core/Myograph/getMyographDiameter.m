%getMyographDiameter  Outer, wall-centre and luminal diameter of a vessel in myograph video
%
%   intervals = getMyographDiameter(s,fName) reads the myograph recording
%   fName (any VideoReader-readable file, e.g. AVI) and, for every row Y of
%   every processed frame, locates the two dark vessel walls and returns the
%   left/right wall x-positions and the wall-to-wall diameter as
%   (frame x Y x measure) arrays.  The vessel is near-vertical: two dark walls sit
%   in a bright field and the lumen/centre and the outer field may be saturated
%   while the walls stay dark.
%
%   THREE MEASURES, ALWAYS.  Each wall is a dark BAND, not a line, so there is no
%   single correct wall position - there are three, and all three are returned from
%   the same per-row relative-contrast rule in one pass over the band:
%     • 'outer'  the band edge nearest the image border  -> outer-wall to outer-wall
%                diameter (the widest of the three),
%     • 'mid'    the darkness-weighted band centre       -> wall-centre to wall-centre
%                diameter (the default analysed measure),
%     • 'inner'  the band edge nearest the lumen         -> the lumen (the narrowest).
%   The third dimension of idxL/idxR/diameter is ordered exactly like .measures,
%   i.e. {'outer','mid','inner'}.  The three share one per-row validity test (a
%   genuine wall dip, contrast > s.wallContrast): a row with no wall has none of
%   the three.
%
%   A centreline is located from the average column-darkness profile (the midpoint
%   of the two wall peaks - the bright lumen when dilated, the middle of the dark
%   region when constricted); the left wall is then searched between the left edge
%   and the centreline and the right wall between the centreline and the right edge,
%   using the SAME relative-contrast rule so a faint/thick wall on one side is
%   measured identically to a strong/thin wall on the other.  The raw per-(frame,Y)
%   estimates are then made physically plausible, INDEPENDENTLY PER MEASURE: the
%   diameter is forced non-negative, estimates that jump unnaturally in time
%   (slow-change prior) or along Y (straight-vessel prior) are rejected, and
%   rejected/failed points are filled by interpolation between the nearest valid
%   points (nearest-repeat at the ends).  The mask marks measured points with 1
%   and interpolated points with 0.
%
%   TIME INTERVALS
%     If s.intervals is non-empty the diameter analysis is restricted to those
%     time windows and returned per interval; s.intervalNames labels them.  With
%     no intervals a single interval named 'default' spans the whole recording, or
%     the s.timeCrop window when one is given.
%
%     A CROP OR A SET OF INTERVALS MAKES THE RECORDING SHORTER, everywhere.  Frames
%     outside the analysed window are never read, never stored and never returned:
%     each interval carries its own time vector and its own diameter arrays, exactly
%     as long as that interval has frames, and the internal per-frame arrays are
%     sized by the analysed span rather than by the video.  Measuring 200 s of a
%     2-hour recording costs 200 s of memory.  Times stay ABSOLUTE seconds from the
%     start of the recording (they are not re-based to zero per interval), so an
%     interval can always be located back in the original footage.
%
%   INPUTS
%     s        parameter struct with fields (defaults filled if missing/empty):
%                • edgeMode      name of the DEFAULT ANALYSED measure, 'mid' (default)
%                                | 'outer' | 'inner'.  It is NO LONGER a selector: all
%                                three measures are always returned, and this field only
%                                tells the consumers (propagation, vasomotion, export,
%                                the GUI) which one to analyse and plot by default.
%                                'center' is accepted as a synonym of 'mid'.  'min'
%                                also means 'mid', but switches the CENTRE rule from the
%                                darkness-weighted centroid to the darkest point (with
%                                s.subpixel refinement) - it is an alternative rule for
%                                the centre measure only, and never affects outer/inner.
%                • wallContrast  min local (bright-dark) contrast for a row to count as
%                                a wall (relative, so faint walls are still measured)
%                • brightPrctile per-row percentile used as the bright level
%                • smoothSigma   Gaussian pre-smoothing, px
%                • dustRadius    vertical half-extent (px) of dust blobs removed by
%                                a vertical-line close; keep < wall height (0 = off)
%                • dustContrast  min bottom-hat response to treat a pixel as dust
%                • rowRange      [lo hi] rows to measure (others are filled)
%                • minWallGap    minimum valid diameter, px (enforces non-negativity)
%                • tSpan         temporal outlier window, frames (0 = off)
%                • ySpan         along-Y outlier window, rows (0 = off)
%                • outlierK      outlier threshold factor (scaled MADs)
%                • subpixel      parabolic sub-pixel refinement of the darkest-point
%                                centre rule, i.e. only when edgeMode='min' (true/false)
%                • smoothSpan    optional final along-Y smoothing span (0 = off)
%                • wallProm      a half's darkest peak must be >= this fraction of the
%                                other half's to count as a wall; below it the vessel is
%                                taken to have dilated beyond the FOV (frame marked invalid)
%                • tSmoothHz     temporal low-pass cutoff (Hz) enforcing the slow-change
%                                prior - the wall/diameter cannot jump faster than this
%                                (0 = off); window = frameRate / tSmoothHz frames
%                • intervals     [start end; ...] time windows, s ([] = whole recording)
%                • intervalNames cellstr of interval names ({} -> 'interval1',...)
%     fName    char/string path to the video file.
%
%   OUTPUT
%     intervals  struct array, one element per interval, with fields
%                • name      interval name
%                • measures  {'outer','mid','inner'} - the order of the 3rd dimension
%                • time      [nFrames x 1] frame times, s
%                • frames    [first last] THE INTERVAL'S FRAME RANGE IN THE ORIGINAL
%                            RECORDING, 1-based, empty when the interval holds no
%                            frames.  This function is the only one that can say it:
%                            it is the one that built the recording's frame axis, and
%                            once the intervals are all that is stored, nothing else
%                            can work out which frame of the file a row came from.
%                            It is what finds the window in the video again.
%                • idxL      [nFrames x nY x 3] left  wall position (x), px
%                • idxR      [nFrames x nY x 3] right wall position (x), px
%                • diameter  [nFrames x nY x 3] idxR - idxL, px (>= 0)
%                • mask      [nFrames x nY] 1 = measured, 0 = interpolated.  ONE mask
%                            for the three: it is the wall-centre ('mid') measure's,
%                            because the three come from the same per-row wall test
%                            and differ only in second-order outlier rejection.
%                • valid     [nFrames x 1] false = a wall was off-FOV (dilated out of
%                            view) so the diameter is only a lower bound (invalid).
%                            A property of the FRAME, so shared by the three measures.
%
%   DEPENDS ON
%     MATLAB base VideoReader, Image Processing Toolbox (imgaussfilt/imclose/
%     strel/im2double/rgb2gray) and Statistics Toolbox (prctile).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Wall detection
% s.edgeMode='mid';       % default analysed measure: 'mid' | 'outer' | 'inner'
%                         % ('min' = 'mid' measured as the darkest point instead)
% s.rowRange=[1 Inf];     % rows to measure (others are interpolated)
% %ADJUSTED IF NECESSARY - Robustness / smoothing
% s.smoothSigma=1;        % Gaussian pre-smoothing, px
% s.dustRadius=8;         % bottom-hat dust removal disk radius (0 = off)
% s.tSpan=25;             % temporal outlier window, frames
% s.ySpan=41;             % along-Y outlier window, rows
% s.outlierK=4;           % outlier threshold factor (scaled MADs)
% s.subpixel=true;        % parabolic sub-pixel refinement for 'min'
% %ADJUSTED IF NECESSARY - Time intervals
% s.intervals=[];         % [start end; ...] in s ([] = whole recording)
% s.intervalNames={};     % {'baseline','drug1',...} ({} -> 'interval1',...)

function intervals = getMyographDiameter(s,fName)

MEASURES={'outer','mid','inner'};   %the 3rd dimension of idxL/idxR/diameter, in this order
nM=numel(MEASURES);
iMid=2;                             %the wall-centre measure - the one whose mask is reported

% ---- default parameters (filled if missing or empty) ----
def.edgeMode='mid';       def.brightPrctile=90;
def.smoothSigma=1.2;      def.dustRadius=8;      def.dustContrast=0.06;
def.rowRange=[1 Inf];     def.minWallGap=3;
def.tSpan=25;             def.ySpan=31;          def.outlierK=3;
def.subpixel=true;        def.smoothSpan=15;     def.intervals=[];
def.intervalNames={};     def.timeCrop=[];       def.progressFcn=[];
def.cancelFcn=[];         def.searchWin=[];
def.wallContrast=0.05;    def.wallProm=0.25;
def.tSmoothHz=1;
fnd=fieldnames(def);
for i=1:numel(fnd)
    if ~isfield(s,fnd{i}) || isempty(s.(fnd{i}))
        s.(fnd{i})=def.(fnd{i});
    end
end

% ---- open video and build the frame time base ----
v=VideoReader(char(fName));
H=v.Height; frameRate=v.FrameRate;
try
    nFrames=v.NumFrames;
catch
    nFrames=[];
end
if isempty(nFrames) || ~isfinite(nFrames)
    nFrames=floor(v.Duration*frameRate);
end
time=(0:nFrames-1)'/frameRate;

% ---- interval windows and names ----
if isempty(s.intervals)
    if ~isempty(s.timeCrop) && numel(s.timeCrop)==2      % default interval = the crop
        ivT=[max(time(1),s.timeCrop(1)) min(time(end),s.timeCrop(2))];
    else
        ivT=[time(1) time(end)];
    end
    ivNames={'default'};
else
    ivT=s.intervals;
    ivNames=s.intervalNames;
    if isempty(ivNames)
        ivNames=arrayfun(@(k)sprintf('interval%d',k),(1:size(ivT,1))','UniformOutput',false);
    end
end
nIv=size(ivT,1);
if numel(ivNames)~=nIv
    error('numel(s.intervalNames) must equal the number of interval rows.');
end
inAny=false(nFrames,1);
for iv=1:nIv
    inAny=inAny | (time>=ivT(iv,1) & time<=ivT(iv,2));
end
%optional time crop: read/process nothing outside [t0 t1]
if ~isempty(s.timeCrop) && numel(s.timeCrop)==2
    inAny=inAny & (time>=s.timeCrop(1) & time<=s.timeCrop(2));
end

% ---- per-frame raw wall detection (only frames inside an interval/crop) ----
r0=max(1,round(s.rowRange(1)));
r1=min(H,round(s.rowRange(2)));
seDust=strel('line',2*max(round(s.dustRadius),1)+1,90); %vertical: removes short dark blobs (dust) but keeps tall vessel walls
%A crop (or a set of intervals) makes the recording SHORTER: only frames fStart..fEnd
%are ever read, so only those are stored.  The raw arrays are the biggest thing this
%function allocates - three positions per wall per row, single, 12 B per frame-row per
%wall - and sizing them by the whole video would cost the untouched footage as well.
%Row g of the raw arrays is absolute frame fStart+g-1 (rawRows does the mapping back).
fStart=find(inAny,1,'first'); fEnd=find(inAny,1,'last');
if isempty(fStart), fStart=1; fEnd=0; end
nProc=max(fEnd-fStart+1,0);                                % frames actually read
idxLraw=nan(nProc,H,nM,'single');
idxRraw=nan(nProc,H,nM,'single');
%the off-FOV flags stay on the ABSOLUTE frame index: 1 byte per frame is nothing, and
%it keeps the per-interval read a plain offLraw(fr).
offLraw=false(nFrames,1);                                  % left  wall dilated beyond the FOV this frame
offRraw=false(nFrames,1);                                  % right wall dilated beyond the FOV this frame
if nProc>0
    if fStart>1, v.CurrentTime=time(fStart); end          % seek past leading skipped footage
    f=fStart-1;
    while hasFrame(v) && f<fEnd
        f=f+1;
        frame=readFrame(v);
        if ~inAny(f), continue; end
        if ndims(frame)==3, I=im2double(rgb2gray(frame)); else, I=im2double(frame); end
        [L,R,offL,offR]=detectWalls(I,s,seDust,r0,r1);   % [H x nM] each
        idxLraw(f-fStart+1,:,:)=reshape(L,[1 H nM]);
        idxRraw(f-fStart+1,:,:)=reshape(R,[1 H nM]);
        offLraw(f)=offL; offRraw(f)=offR;
        if mod(f,500)==0
            myoProgress(s,f,fEnd);
            if ~isempty(s.cancelFcn) && s.cancelFcn()      % user asked to stop
                error('getMyographDiameter:cancelled','Diameter measurement stopped by user.');
            end
        end
    end
end

% ---- per-interval post-processing (INDEPENDENTLY PER MEASURE) and output ----
intervals=struct('name',{},'measures',{},'time',{},'frames',{},'idxL',{},'idxR',{}, ...
    'diameter',{},'mask',{},'valid',{});
[Lc,Rc,Dc,Mc]=deal(cell(1,nM));
for iv=1:nIv
    fr=find(time>=ivT(iv,1) & time<=ivT(iv,2));         % this interval's OWN frames
    [rawL,rawR]=rawRows(idxLraw,idxRraw,fr,fStart,nProc,H,nM);
    for m=1:nM
        %the slow-change-in-time prior, the along-Y prior, the minWallGap clamp and the
        %interpolation fill are all per measure - an outer edge and a luminal edge are
        %different physical series and must not be filtered as one.
        [Lc{m},Rc{m},Dc{m},Mc{m},valid]=postProcessWalls(rawL(:,:,m),rawR(:,:,m), ...
            offLraw(fr),offRraw(fr),s,frameRate);       % valid is per FRAME: identical for the three
    end
    intervals(iv).name=char(ivNames{iv});
    intervals(iv).measures=MEASURES;
    intervals(iv).time=time(fr);
    %WHERE THIS WINDOW SITS IN THE ORIGINAL RECORDING.  fr came from the file's own
    %frame axis, so it is contiguous and 1-based, and it is the only chance to record
    %it: after this the windows are all there is, and a stored row knows its absolute
    %TIME but not which frame of the .avi it was read from.
    if isempty(fr), intervals(iv).frames=[]; else, intervals(iv).frames=[fr(1) fr(end)]; end
    intervals(iv).idxL=cat(3,Lc{:});
    intervals(iv).idxR=cat(3,Rc{:});
    intervals(iv).diameter=cat(3,Dc{:});
    intervals(iv).mask=Mc{iMid};                        % the wall-centre measure's (see header)
    intervals(iv).valid=valid;                          % per-frame: false = a wall was off-FOV (diameter invalid)
end
end

% =====================================================================
function myoProgress(s,f,total)
%MYOPROGRESS  route frame progress to s.progressFcn (GUI) or the console.
%   Named myo* so it cannot SHADOW a Core/Reporting function: addpath(genpath(...))
%   puts the whole library on one path, and a local called report* would win inside
%   this file.  This hook is the Myograph subsystem's own, not the library seam.
if isfield(s,'progressFcn') && ~isempty(s.progressFcn)
    try, s.progressFcn(f,total); catch, end
else
    fprintf('  frame %d / %d\n',f,total);
end
end

% =====================================================================
function [A,B]=rawRows(idxLraw,idxRraw,fr,fStart,nProc,H,nM)
%RAWROWS  the raw wall positions of the ABSOLUTE frames fr, as [numel(fr) x H x nM]
%   The raw arrays hold only the frames that were read (absolute fStart..fStart+nProc-1).
%   An interval may legitimately name frames outside that window - it can reach past a
%   time crop, and the union of intervals can have gaps - and those frames were never
%   measured, so they come back NaN here and are filled by postProcessWalls exactly as
%   an unmeasurable frame inside the window is.
A=nan(numel(fr),H,nM,'single'); B=nan(numel(fr),H,nM,'single');
in = fr>=fStart & fr<=fStart+nProc-1;
if any(in)
    g=fr(in)-fStart+1;
    A(in,:,:)=idxLraw(g,:,:);
    B(in,:,:)=idxRraw(g,:,:);
end
end

% =====================================================================
function [L,R,offL,offR]=detectWalls(I,s,seDust,r0,r1)
%DETECTWALLS  raw left/right wall x-positions per row for one frame ([H x 3], NaN where absent)
%   The two vessel walls ALWAYS lie on opposite sides of the field-of-view centre
%   (the vessel is mounted centred), so the search is split at the FOV centre: the
%   left wall is the strongest column-darkness peak in the left half, the right wall
%   the strongest peak in the right half.  The darkness profile is CONTINUOUS (an
%   average over rows, not thresholded), so a FAINT wall still forms a peak and is
%   found, and a dust speck (dark in only a few rows) averages away.  If a half has
%   no real wall the vessel has dilated beyond the FOV on that side: the wall is
%   reported at the window edge and offL/offR flags that frame invalid.  Each wall is
%   then refined per row by the SAME relative-contrast rule within an adaptive band,
%   so a thick/faint wall is measured identically to a thin/strong one.  The refinement
%   returns THREE positions per row - outer edge, band centre, luminal edge - so L and R
%   are [H x 3], columns ordered {'outer','mid','inner'}.  offL/offR stay per FRAME:
%   being off-FOV is a property of the frame, not of the measure.
[H,W]=size(I);
L=nan(H,3); R=nan(H,3); offL=false; offR=false;
if s.smoothSigma>0, Is=imgaussfilt(I,s.smoothSigma); else, Is=I; end

%dust removal: bottom-hat finds small dark blobs -> inpaint to local background
if s.dustRadius>0
    clo=imclose(Is,seDust);
    dust=(clo-Is)>s.dustContrast;
    if any(dust(:)), Is(dust)=clo(dust); end
end

%illumination-/saturation-robust normalisation (bright outside AND, when dilated, inside)
rowBright=prctile(Is,s.brightPrctile,2);
Inorm=Is./max(rowBright,eps);

%--- one wall per HALF, split at the FOV centre (walls never share a side) ---
rr=false(H,1); rr(r0:r1)=true;
P=1-mean(Inorm(rr,:),1);                             % continuous darkness (no threshold)
P=P-min(P); if numel(P)>=5, P=smoothdata(P,'movmean',3); end
mg=max(2,round(0.02*W)); P([1:mg, W-mg+1:end])=0;    % ignore the frame margin
imgC=round(W/2);
if imgC<=mg+1 || (W-mg)<=imgC+1, return; end
[pL,xL]=max(P(mg:imgC));       xL=xL+mg-1;            % strongest wall left of the FOV centre
[pR,xR]=max(P(imgC+1:W-mg));   xR=xR+imgC;            % strongest wall right of the FOV centre
pk=max(pL,pR);
if pk<=0, offL=true; offR=true; return; end          % no vessel at all this frame
leftOn = pL>=s.wallProm*pk;                           % a real wall in this half? (else dilated beyond FOV)
rightOn= pR>=s.wallProm*pk;
offL=~leftOn; offR=~rightOn;

%--- LEFT wall: adaptive band around xL (clamped left of centre, padded toward lumen) ---
if leftOn
    [loL,hiL]=bandExtent(P,xL,0.4); padL=max(3,round(0.5*(hiL-loL)));
    loL=max(loL,mg); hiL=min(hiL+padL,imgC-1);
    if ~isempty(s.searchWin), loL=max(loL,xL-s.searchWin); hiL=min(hiL,xL+s.searchWin); end
    if hiL>loL, [iL,okL]=wallInWindow(Inorm,loL:hiL,s,'L'); L(okL,:)=iL(okL,:); end
else
    L(r0:r1,:)=1;                                     % dilated beyond FOV -> left window edge (invalid)
end

%--- RIGHT wall: adaptive band around xR (clamped right of centre, padded toward lumen) ---
if rightOn
    [loR,hiR]=bandExtent(P,xR,0.4); padR=max(3,round(0.5*(hiR-loR)));
    loR=max(loR-padR,imgC+1); hiR=min(hiR,W-mg);
    if ~isempty(s.searchWin), loR=max(loR,xR-s.searchWin); hiR=min(hiR,xR+s.searchWin); end
    if hiR>loR, [iR,okR]=wallInWindow(Inorm,loR:hiR,s,'R'); R(okR,:)=iR(okR,:); end
else
    R(r0:r1,:)=W;                                     % dilated beyond FOV -> right window edge (invalid)
end
L([1:r0-1,r1+1:H],:)=NaN; R([1:r0-1,r1+1:H],:)=NaN;
end

% =====================================================================
function [lo,hi]=bandExtent(prof,x0,frac)
%BANDEXTENT  columns of the dark band around x0 (where prof stays above frac*prof(x0))
th=frac*prof(x0); n=numel(prof);
lo=x0; while lo>1 && prof(lo-1)>th, lo=lo-1; end
hi=x0; while hi<n && prof(hi+1)>th, hi=hi+1; end
end

% =====================================================================
function [pos,ok]=wallInWindow(Inorm,cols,s,side)
%WALLINWINDOW  per-row OUTER/CENTRE/LUMINAL wall x-positions within a column band
%   The wall is located relative to the LOCAL bright level (lumen side) inside the
%   band, NOT an absolute threshold - so a faint wall (a shallow dip) is measured
%   exactly like a strong one, which is what makes the two walls symmetric.  A row
%   is valid (ok) only if its band shows a real dip: local (bright-dark) > wallContrast.
%
%   All three positions come out of ONE pass over seg, from that single rule: the
%   dark run (seg below the per-row half-drop level) has a first and a last column,
%   and which of the two is the OUTER edge depends only on the side; the CENTRE is
%   the darkness-weighted centroid of the same band.  No loop, and no per-measure
%   re-thresholding - hence the three cannot disagree about whether a wall is there.
%
%   pos is [H x 3], columns ordered {'outer','mid','inner'}; ok is [H x 1], shared.
seg=Inorm(:,cols);
segBright=max(seg,[],2); segDark=min(seg,[],2);
contrast=segBright-segDark;
ok=contrast>s.wallContrast;                          % a genuine wall dip is present in this row
thr=segBright-0.5*contrast;                          % per-row half-drop level (relative)

%--- the two band edges: first and last column of the dark run ---
dk=seg<thr;
[~,jF]=max(dk,[],2);          posFirst=cols(1)+jF-1;   posFirst=posFirst(:);
[~,jL]=max(fliplr(dk),[],2);  posLast =cols(end)-jL+1; posLast =posLast(:);
if side=='L'                                         % the border is on the low-column side
    outer=posFirst; inner=posLast;
else                                                 % the border is on the high-column side
    outer=posLast;  inner=posFirst;
end

%--- the band centre: darkness-weighted centroid, or the darkest point when asked ---
if strcmpi(s.edgeMode,'min')                         % alternative centre rule (sharp walls)
    [~,j]=min(seg,[],2); mid=(cols(1)+j-1); mid=mid(:);
    if s.subpixel, mid=subpix(Inorm,mid,ok); end
else                                                 % default: relative to the local bright level
    wdk=max(segBright-seg,0); wsum=sum(wdk,2);
    mid=(wdk*cols(:))./max(wsum,eps); mid=mid(:);
end

pos=[outer mid inner];
end

% =====================================================================
function idx=subpix(Inorm,idx,ok)
%SUBPIX  parabolic sub-pixel refinement of a per-row intensity minimum at column idx
[H,W]=size(Inorm);
rows=(1:H)';
interior=ok & idx>1 & idx<W;
r=rows(interior); c=idx(interior);
ym=Inorm(sub2ind([H W],r,c-1));
y0=Inorm(sub2ind([H W],r,c));
yp=Inorm(sub2ind([H W],r,c+1));
den=ym-2*y0+yp;
off=zeros(size(c));
g=abs(den)>eps;
off(g)=0.5*(ym(g)-yp(g))./den(g);
off=max(min(off,0.5),-0.5);
idx(interior)=c+off;
end

% =====================================================================
function [L,R,D,mask,valid]=postProcessWalls(L,R,offL,offR,s,frameRate)
%POSTPROCESSWALLS  enforce physical plausibility, filter in time/Y, fill, build masks
%   Runs on ONE measure at a time ([nFrames x nY] L and R); the caller loops the three.
%   valid (per frame) is false where a wall was off-FOV (dilated out of view); those
%   frames keep their edge-clamped position (a lower-bound diameter) and are excluded
%   from the temporal statistics/filtering so they cannot bias the valid frames.
[nF,nY]=size(L);
offL=logical(offL(:)); offR=logical(offR(:));
off=offL|offR;                                        % per-frame: a wall was beyond the FOV
L0=L; R0=R;                                           % keep the raw edge clamp of off-FOV frames

%enforce non-negative diameter / minimum wall separation
bad=~(isfinite(L)&isfinite(R)) | (R-L)<s.minWallGap;
L(bad)=NaN; R(bad)=NaN;

%hold off-FOV frames OUT of the time statistics/filtering (their clamp is restored later)
L(off,:)=NaN; R(off,:)=NaN;

%reject unnatural temporal jumps (slow-change prior) then along-Y outliers
if s.tSpan>=3 && nF>=5
    L(isoutlier(L,'movmedian',s.tSpan,1,'ThresholdFactor',s.outlierK))=NaN;
    R(isoutlier(R,'movmedian',s.tSpan,1,'ThresholdFactor',s.outlierK))=NaN;
end
if s.ySpan>=3 && nY>=5
    L(isoutlier(L,'movmedian',s.ySpan,2,'ThresholdFactor',s.outlierK))=NaN;
    R(isoutlier(R,'movmedian',s.ySpan,2,'ThresholdFactor',s.outlierK))=NaN;
end

%a point is valid only if both walls survived (and the frame was in-FOV)
measured=isfinite(L)&isfinite(R);
L(~measured)=NaN; R(~measured)=NaN;
mask=single(measured);

%fill failed points by interpolation (nearest-repeat at the ends)
L=fillWalls(L); R=fillWalls(R);

%enforce the smooth-wall prior along Y (walls change slowly row-to-row): a robust
%median first to kill residual spikes, then Savitzky-Golay for a <=~1 px/row curve
if s.smoothSpan>=3
    sp=min(s.smoothSpan,nY);
    L=smoothdata(L,2,'movmedian',sp); R=smoothdata(R,2,'movmedian',sp);
    L=smoothdata(L,2,'sgolay',sp);    R=smoothdata(R,2,'sgolay',sp);
end

%enforce the slow-change-in-TIME prior: the wall/diameter cannot jump faster than
%~tSmoothHz (no physiology drives sub-second back-and-forth), so low-pass along time
%with a zero-phase moving average (a robust median first removes residual 1-5 frame
%spikes).  win = frames per tSmoothHz cycle; skips off-FOV frames which are filled next.
if isfield(s,'tSmoothHz') && s.tSmoothHz>0 && frameRate>0 && nF>=5
    win=max(3,round(frameRate/s.tSmoothHz));
    if any(off), L(off,:)=NaN; R(off,:)=NaN; L=fillWalls(L); R=fillWalls(R); end
    L=smoothdata(L,1,'movmedian',min(win,nF)); R=smoothdata(R,1,'movmedian',min(win,nF));
    L=smoothdata(L,1,'movmean',  min(win,nF)); R=smoothdata(R,1,'movmean',  min(win,nF));
end

%restore the off-FOV edge clamp (these frames are invalid, shown as a lower bound)
if any(offL), L(offL,:)=L0(offL,:); end
if any(offR), R(offR,:)=R0(offR,:); end
valid=~off;                                           % per-frame validity

%guarantee a non-negative diameter
R=max(R,L);
D=R-L;
end

% =====================================================================
function X=fillWalls(X)
%FILLWALLS  fill NaN by linear interpolation along Y then frames (nearest at ends)
if all(isnan(X(:))), X(:)=0; return; end
X=fillmissing(X,'linear',2,'EndValues','nearest');
X=fillmissing(X,'linear',1,'EndValues','nearest');
if any(isnan(X(:)))
    X=fillmissing(X,'nearest',2);
    X=fillmissing(X,'nearest',1);
end
end
