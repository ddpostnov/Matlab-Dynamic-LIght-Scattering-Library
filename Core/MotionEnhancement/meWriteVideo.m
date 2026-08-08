%meWriteVideo - Write one product as an AVI, with everything it needs to be read on it.
%
%   A MOTION-MAGNIFIED MOVIE IS PERSUASIVE WHETHER OR NOT IT IS TRUE, and it travels
%   without its caption.  Given noise and a large enough alpha the pipeline will
%   produce convincing wall motion from a recording containing none, so a frame that
%   does not carry alpha, the band it was filtered at and the calibre it was
%   calibrated for is a misinformation hazard the moment it reaches a slide.  Every
%   parameter is therefore burnt into the picture rather than written beside it.
%
%   THE ORIGINAL TRAVELS WITH IT, ON THE SAME INTENSITY SCALE.  The panels are the
%   same pixels of the same frames, cropped identically by meEnhance, mapped through
%   one pair of limits.  A viewer comparing two pictures scaled differently is
%   comparing the scaling.
%
%   THE DIFFERENCE PANEL IS WHERE THE MAGNIFICATION ACTUALLY IS.  Magnified minus
%   original, multiplied up and drawn about mid-grey, so nothing is a flat grey
%   field.  It is the panel that separates a product from its null fastest: on a
%   product the difference sits on the vessel walls and moves with the beat, on a
%   null it is everywhere and it is noise.
%
%   A NULL SAYS SO ON EVERY FRAME.  A control that can be mistaken for a product at a
%   glance is worse than no control, so meta.isNull paints a red banner and names
%   what was destroyed.
%
%   A FRAME THE ANIMAL MOVED THROUGH SAYS SO ON EVERY FRAME.  The vasomotion run
%   spans the whole recording and 53 per cent of its decimated frames touch a
%   movement bout; they are MARKED, not dropped, because a gap a viewer cannot see is
%   worse than one they can.
%
%   THE CINE'S PER-BIN SEM IS ON THE FRAME, NOT IN A FOOTNOTE.  It is the only thing
%   in the picture that can contradict the movie: the strip draws the averaged cycle
%   over the wall against its own standard error, with a cursor at the bin being
%   shown.  On the reference recording the excursion is 9.5 times the SEM, and that
%   ratio is what makes the strip worth the space.
%
%   ONE CODEC FOR A PRODUCT AND ITS NULL, DELIBERATELY.  The plan asked for Motion
%   JPEG for the deliverables and an uncompressed control.  That is the wrong way
%   round: the requirement the whole comparison rests on is that the null gets the
%   IDENTICAL treatment, and two codecs are not identical.  s.videoCodec therefore
%   applies to both, and meRun measures what the codec costs by re-reading a written
%   frame rather than assuming it costs nothing.
%
%   VideoWriter APPENDS THE PROFILE'S EXTENSION, so an outPath of 'x.avi.partial'
%   produces 'x.avi.partial.avi'.  The path is given whole, with its .avi, and is
%   checked against what VideoWriter actually opened.
%
% Syntax:
%    v = meWriteVideo(s, stacks, meta)
%
% Inputs:
%    s      - settings.  Fields read:
%             .videoFps     frames per second written.
%             .videoLoops   times the stack is repeated - for a cine, which is one
%                           cardiac period and is meant to be looped.
%             .videoCodec   a VideoWriter profile name.
%             .videoQuality 1 to 100, for the profiles that have one.
%             .showDifference  add the third panel.
%             .diffGain     what the difference is multiplied by before it is drawn.
%             .videoPct     percentiles of the ORIGINAL mapped to black and white.
%             .videoFontSize
%             .scaleBarUm   length of the scale bar, micrometres.
%             .reportLevel  'quiet' silences the progress line.
%    stacks - .original   [rows columns frames], the unmagnified input
%             .magnified  [rows columns frames], the same frames magnified
%    meta   - .outPath    where to write, including .avi
%             .stem       the recording's name, burnt in
%             .product    what this run is, in words, burnt in
%             .alpha      the amplification, scalar or one per level
%             .lambdaFull the wavelengths those levels carry, recording pixels
%             .passband   [low high] of the temporal filter
%             .passbandUnit  'Hz' or 'cycles/cycle'
%             .differential  true when bulk motion was removed
%             .masked     true when amplification was confined to the wall mask
%             .calibre    one line saying which vessels the run is valid for
%             .pixelSize  micrometres per pixel OF THIS COPY
%             .axis       'time' or 'phase'
%             .axisValues [frames 1] seconds, or phase in [0,1)
%             .rateHz     the rhythm's frequency, for the cine's readout
%             .bad        [frames 1] logical, frames inside a movement bout
%             .isNull     true for a control
%             .nullText   what the null destroyed, one line
%             .controlName the matched control's file name, named on the product's
%                         own banner so the pair cannot be separated
%             .strip      optional .value .sem [nFrames 1] and .label - the cine's
%                         per-bin excursion and its standard error
%             .notes      optional {n 1} extra lines under the parameters
%
% Outputs:
%    v - .path .frames .size [height width] .bytes .seconds .codec .slowdown
%
% Example:
%    [mag, info] = meEnhance(s, pass);
%    meta = struct('outPath','rec_ME_puls_cine.avi','stem','rec','product','cardiac cine', ...
%        'alpha',s.alpha,'passband',info.passband,'passbandUnit','cycles/cycle', ...
%        'axis','phase','axisValues',(0:24)'/25,'pixelSize',info.geom.pixelSize);
%    v = meWriteVideo(s, struct('original',info.original,'magnified',mag), meta);
%
% Dependencies: Computer Vision Toolbox (insertText); VideoWriter.
% See also: meEnhance, meRun, meCine, meShift, meSettings
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function v = meWriteVideo(s, stacks, meta)

tAll = tic;
O = single(stacks.original);
M = single(stacks.magnified);
if ~isequal(size(O), size(M))
    error('meWriteVideo:size', ...
        ['The original is %s and the magnified stack is %s. meEnhance crops both by ' ...
         's.edgeCrop and returns them together for exactly this reason.'], ...
        mat2str(size(O)), mat2str(size(M)));
end
[ph, pw, nF] = size(O);
meta = fillMeta(meta, nF);

% ---- one pair of limits, from the original ----------------------------------
lim = prctile(O(1:max(1,floor(numel(O)/2e6)):end), s.videoPct);
if diff(lim) <= 0, lim = [min(O(:)) max(O(:))]; end

nPanel = 2 + double(s.showDifference);

% ---- the canvas -------------------------------------------------------------
g   = geometry(ph, pw, nPanel, s, meta);
tpl = template(g, s, meta);
barI= barTile(g);
stp = [];
if ~isempty(meta.strip)
    stp = stripBase(meta.strip, g.strip(4), g.strip(3), s.videoFontSize);
end

% ---- write ------------------------------------------------------------------
w = VideoWriter(meta.outPath, s.videoCodec);
w.FrameRate = s.videoFps;
if isprop(w,'Quality'), w.Quality = s.videoQuality; end
open(w);
cleanup = onCleanup(@() close(w));

nOut = nF*max(1,round(s.videoLoops));
for i = 1:nOut
    k = mod(i-1, nF) + 1;
    F = tpl;
    F = paste(F, g.panel(1,:), gray2rgb(O(:,:,k), lim));
    F = paste(F, g.panel(2,:), gray2rgb(M(:,:,k), lim));
    if nPanel == 3
        dif = (M(:,:,k)-O(:,:,k)).*s.diffGain + mean(lim);
        F   = paste(F, g.panel(3,:), gray2rgb(dif, lim));
    end
    if ~isempty(barI)
        F = paste(F, g.barBox, barI);        % AFTER the panels, or it is painted over
    end
    if ~isempty(stp)
        F = paste(F, g.strip, stripCursor(stp, k, numel(meta.strip.value)));
    end
    F = insertText(F, [g.pad g.readY], readout(meta, k, i, nF, s), ...
        'FontSize', s.videoFontSize, 'BoxOpacity', 0, 'TextColor', 'white');
    if meta.bad(k)
        F = insertText(F, [g.W-g.pad g.readY], 'ANIMAL MOVING', 'AnchorPoint','RightTop', ...
            'FontSize', s.videoFontSize, 'BoxOpacity', 0, 'TextColor', [255 90 70]);
        F = frameBorder(F, g, [255 90 70]);
    end
    writeVideo(w, F);
end
clear cleanup                                    % close before dir() reads the size

onDisk = dir(meta.outPath);
if isempty(onDisk)
    error('meWriteVideo:notWritten', ...
        ['VideoWriter did not produce %s. It appends the profile''s own extension, ' ...
         'so a path that already ends in the wrong one gets a second.'], meta.outPath);
end

v = struct('path', meta.outPath, 'frames', nOut, 'size', [g.H g.W], ...
    'bytes', onDisk.bytes, 'seconds', toc(tAll), 'codec', s.videoCodec, ...
    'fps', s.videoFps, 'slowdown', slowdown(meta, s));
if ~strcmpi(s.reportLevel,'quiet')
    fprintf('   wrote %s  %d frames  %dx%d  %.0f MB  %.1f s\n', ...
        nameOf(meta.outPath), nOut, g.W, g.H, onDisk.bytes/2^20, v.seconds);
end
end

% =====================================================================
function meta = fillMeta(meta, nF)
%fillMeta  The optional half of meta, so nothing below has to test for it.
def = struct('stem','','product','','alpha',[],'lambdaFull',[],'passband',[], ...
    'passbandUnit','Hz','differential',true,'masked',true,'calibre','', ...
    'pixelSize',NaN,'axis','time','axisValues',(0:nF-1)','rateHz',NaN, ...
    'bad',false(nF,1),'isNull',false,'nullText','','controlName','', ...
    'strip',[],'notes',{{}});
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(meta,f{i}) || (isempty(meta.(f{i})) && ~ischar(def.(f{i})))
        meta.(f{i}) = def.(f{i});
    end
end
meta.bad = logical(meta.bad(:));
if numel(meta.bad) ~= nF, meta.bad = false(nF,1); end
meta.axisValues = meta.axisValues(:);
end

% =====================================================================
function g = geometry(ph, pw, nPanel, s, meta)
%geometry  Where everything sits.  Even width and height, because some profiles
%   refuse an odd one and the failure is at the last frame.
%
%   THE CANVAS IS AS WIDE AS THE TEXT NEEDS, NOT ONLY AS WIDE AS THE PICTURES.  The
%   burnt-in parameters are the point of the frame, and a canvas sized off the panels
%   alone truncates them silently on a narrow product - which loses exactly the alpha
%   the plan insists must never travel separately.  The widths are MEASURED by
%   rendering each line once, not estimated from a character count.
g.pad = 10;  g.gap = 8;  g.line = round(s.videoFontSize*1.65);
g.banner = g.line + 6;                  % ALWAYS reserved - see bannerText
g.lines  = headerLines(meta);
g.hdrH   = numel(g.lines)*g.line;
g.labelH = g.line + 4;
g.stripH = 0;
if ~isempty(meta.strip), g.stripH = 62; end
g.footH  = g.line + 10;

wPanels = nPanel*pw + (nPanel-1)*g.gap;
wText   = 0;
probe   = [g.lines(:); {readoutTemplate(meta)}; {bannerText(meta)}];
for i = 1:numel(probe), wText = max(wText, textWidth(probe{i}, s.videoFontSize)); end
if any(meta.bad), wText = wText + textWidth('    ANIMAL MOVING', s.videoFontSize); end

g.W = 2*g.pad + max(wPanels, wText);
yTop      = g.pad + g.banner + g.hdrH + 4;
g.labelY  = yTop;
g.panelY  = yTop + g.labelH;
g.stripY  = g.panelY + ph + 6;
g.readY   = g.stripY + g.stripH + 6;
g.H       = g.readY + g.footH;
g.W = g.W + mod(g.W,2);
g.H = g.H + mod(g.H,2);

x0 = round((g.W - wPanels)/2);                       % centred when the text is wider
g.panel = zeros(nPanel,4);
for i = 1:nPanel
    g.panel(i,:) = [x0 + (i-1)*(pw+g.gap), g.panelY, pw, ph];
end
g.strip  = [x0 g.stripY wPanels g.stripH];
g.panelW = pw;  g.panelH = ph;  g.nPanel = nPanel;

% The scale bar sits over the bottom-left of the first picture, in a box of its own
% so it stays legible on a dark corner and on a bright one.  The box is as wide as
% the LONGER of the bar and its label - on the binned copy 200 um is 40 px and the
% label is ninety, and a box sized off the bar alone eats the last character.
g.fontSize = s.videoFontSize;
g.barText  = sprintf('%g um', s.scaleBarUm);
g.barLen = 0;  g.barBox = [];
if isfinite(meta.pixelSize) && meta.pixelSize > 0
    L  = round(s.scaleBarUm/meta.pixelSize);
    bw = max(L, textWidth(g.barText, s.videoFontSize)) + 14;
    bh = g.line + 12;
    if L > 4 && bw < pw-8 && bh < ph-8
        g.barLen = L;
        g.barBox = [g.panel(1,1)+6, g.panelY+ph-bh-6, bw, bh];
    end
end
end

% =====================================================================
function w = textWidth(txt, fs)
%textWidth  How wide insertText actually draws a line, measured rather than guessed.
if isempty(txt), w = 0; return; end
I = insertText(zeros(round(2.2*fs), 6000, 3, 'uint8'), [1 1], txt, ...
    'FontSize', fs, 'BoxOpacity', 0, 'TextColor', 'white');
c = find(any(any(I,3),1), 1, 'last');
if isempty(c), w = 0; else, w = double(c) + 4; end
end

% =====================================================================
function t = bannerText(meta)
%bannerText  What this video is, in the one place a viewer cannot miss.
if meta.isNull
    t = ['NULL CONTROL - ' meta.nullText];
elseif ~isempty(meta.controlName)
    t = ['PRODUCT - watch it beside its matched control, ' meta.controlName];
else
    t = 'PRODUCT';
end
end

% =====================================================================
function L = headerLines(meta)
%headerLines  The static text, in the order it is drawn.
L = {sprintf('%s   |   %s', meta.stem, meta.product), parameterLine(meta)};
if ~isempty(meta.calibre), L{end+1} = meta.calibre; end
L = [L, meta.notes(:)'];
L = L(:);
end

% =====================================================================
function B = barTile(g)
%barTile  The scale bar and its length, as one tile pasted over the picture.
B = [];
if isempty(g.barBox), return; end
w = g.barBox(3);  h = g.barBox(4);
B = zeros(h, w, 3, 'uint8') + uint8(20);
B = insertText(B, [7 2], g.barText, 'FontSize', g.fontSize, ...
    'BoxOpacity', 0, 'TextColor', 'white');
B(h-8:h-5, 7:7+g.barLen-1, :) = 255;
end

% =====================================================================
function tpl = template(g, s, meta)
%template  Everything that does not change from frame to frame, rendered once.
%   insertText is not free, and a 1 302-frame run would pay for the header 1 302
%   times over for text that is the same every time.
tpl = zeros(g.H, g.W, 3, 'uint8') + uint8(24);

% THE BANNER IS ALWAYS THERE, on a product as well as on a control, so that the two
% videos come out the SAME SIZE.  They are meant to be watched side by side, and two
% pictures of different heights are two pictures a viewer has to reconcile before
% they can compare them.
y   = g.pad;
col = [40 50 62];
if meta.isNull, col = [150 30 25]; end
for c = 1:3
    tpl(y:y+g.banner-1, g.pad:g.W-g.pad, c) = col(c);
end
tpl = insertText(tpl, [g.pad+6 y+3], bannerText(meta), ...
    'FontSize', s.videoFontSize, 'BoxOpacity', 0, 'TextColor', 'white');
y = y + g.banner;

for i = 1:numel(g.lines)
    col = 'white';
    if i > 1, col = [200 205 215]; end
    tpl = insertText(tpl, [g.pad y+(i-1)*g.line], g.lines{i}, ...
        'FontSize', s.videoFontSize, 'BoxOpacity', 0, 'TextColor', col);
end

caps = {'ORIGINAL','MAGNIFIED'};
if g.nPanel == 3, caps{3} = sprintf('DIFFERENCE x%g', s.diffGain); end
for i = 1:g.nPanel
    tpl = insertText(tpl, [g.panel(i,1) g.labelY], caps{i}, ...
        'FontSize', s.videoFontSize, 'BoxOpacity', 0, 'TextColor', [140 200 255]);
end
end

% =====================================================================
function t = parameterLine(meta)
%parameterLine  alpha, the scales it acts on, the band, and the two mode switches.
a = meta.alpha;
if isscalar(a), aTxt = sprintf('%.4g', a); else, aTxt = strtrim(sprintf('%.4g ', a)); end
t = sprintf('alpha %s', aTxt);
if ~isempty(meta.lambdaFull)
    t = [t sprintf('  on wavelengths %s px', strtrim(sprintf('%g ', meta.lambdaFull)))];
end
if ~isempty(meta.passband)
    t = [t sprintf('  |  band %.3g-%.3g %s', meta.passband(1), meta.passband(2), meta.passbandUnit)];
end
if meta.differential, t = [t '  |  differential']; else, t = [t '  |  absolute']; end
if meta.masked,       t = [t '  |  wall mask']; else, t = [t '  |  no mask']; end
end

% =====================================================================
function t = readout(meta, k, i, nF, s)
%readout  Where in the recording, or where in the beat, this frame is.
%   A NULL HAS NO CLOCK, and saying otherwise on the frame would be a caption about a
%   recording the stack no longer is: its frames are permuted, so 'axis' is 'index'
%   and the readout counts frames.
if strcmpi(meta.axis,'index')
    t = sprintf('frame %d of %d   |   %g fps written   |   no time axis: the frames are permuted', ...
        k, nF, s.videoFps);
elseif strcmpi(meta.axis,'phase')
    loop = floor((i-1)/nF) + 1;
    t = sprintf('cardiac phase %3.0f%%   |   beat %.1f ms shown in %.2f s (%.0fx slowed)   |   loop %d of %d', ...
        100*meta.axisValues(k), 1000/meta.rateHz, nF/s.videoFps, slowdown(meta,s), ...
        loop, max(1,round(s.videoLoops)));
else
    f = slowdown(meta, s);
    if f >= 1
        t = sprintf('t = %7.2f s   |   %g fps written, %.3gx slowed', ...
            meta.axisValues(k), s.videoFps, f);
    else
        t = sprintf('t = %7.2f s   |   %g fps written, %.3gx time-compressed', ...
            meta.axisValues(k), s.videoFps, 1/f);
    end
end
end

% =====================================================================
function t = readoutTemplate(meta)
%readoutTemplate  The widest the readout can get, for sizing the canvas.
if strcmpi(meta.axis,'index')
    t = 'frame 8888 of 8888   |   8888 fps written   |   no time axis: the frames are permuted';
elseif strcmpi(meta.axis,'phase')
    t = 'cardiac phase 100%   |   beat 888.8 ms shown in 88.88 s (8888x slowed)   |   loop 88 of 88';
else
    t = 't = 8888.88 s   |   8888 fps written, 8888x compressed';
end
end

% =====================================================================
function f = slowdown(meta, s)
%slowdown  How much slower than the recording the video plays.  Under one is
%   time compression, which is what a vasomotion run is for.
if strcmpi(meta.axis,'index')
    f = NaN;                                     % a permuted stack has no rate
elseif strcmpi(meta.axis,'phase')
    f = meta.rateHz*numel(meta.axisValues)/s.videoFps;
else
    dt = median(diff(meta.axisValues));
    if ~isfinite(dt) || dt <= 0, f = 1; return; end
    f = 1/(dt*s.videoFps);
end
end

% =====================================================================
function F = paste(F, box, I)
%paste  An RGB tile into the canvas at [x y w h].
F(box(2):box(2)+box(4)-1, box(1):box(1)+box(3)-1, :) = I;
end

% =====================================================================
function I = gray2rgb(X, lim)
%gray2rgb  One plane through one pair of limits, as RGB.
X = (double(X)-lim(1))./(lim(2)-lim(1));
X = uint8(255*min(max(X,0),1));
I = repmat(X, 1, 1, 3);
end

% =====================================================================
function S = stripBase(strip, H, W, fs)
%stripBase  The averaged cycle over the wall and its standard error, drawn once.
%   Straight into the pixels rather than through a figure: a figure per frame is
%   seconds, and this has to reach every frame of a 1 302-frame product.  Only the
%   cursor changes from frame to frame, and that is one column.
%
%   IT IS THE ONLY THING ON THE FRAME THAT CAN CONTRADICT THE MOVIE.  A magnified
%   movie is persuasive whether or not it is true; the standard error says which part
%   of the cycle is determined by the beats that went into it and which part is a
%   handful of them disagreeing.
S = zeros(H, W, 3, 'uint8') + uint8(38);
v = strip.value(:);  e = strip.sem(:);
n = numel(v);
if n < 2, return; end
top = 16;  bot = H-6;
% A CONTROL'S STRIP IS DRAWN ON THE PRODUCT'S SCALE.  Auto-scaling each would blow a
% flat control up until it filled the strip and looked like a cycle, which is the
% opposite of what the comparison is for.
if isfield(strip,'hi') && ~isempty(strip.hi)
    hi = strip.hi;
else
    hi = max([v+e; -(v-e)]);
end
if ~isfinite(hi) || hi <= 0, hi = 1; end
mid = round((top+bot)/2);
sc  = (bot-top)/2/(1.05*hi);
x   = round(linspace(1, W, n));

S(mid, :, :) = 80;                                          % the zero line
for j = 1:n-1                                               % the SEM ribbon
    xa = x(j):x(j+1);
    lo = round(mid - sc*interp1([x(j) x(j+1)],[v(j)-e(j) v(j+1)-e(j+1)],xa));
    up = round(mid - sc*interp1([x(j) x(j+1)],[v(j)+e(j) v(j+1)+e(j+1)],xa));
    for c = 1:numel(xa)
        r = max(min(up(c),lo(c)),top):min(max(up(c),lo(c)),bot);
        S(r, xa(c), 1) =  60;  S(r, xa(c), 2) = 100;  S(r, xa(c), 3) = 140;
    end
end
for j = 1:n-1                                               % the cycle itself
    xa = x(j):x(j+1);
    yv = round(mid - sc*interp1([x(j) x(j+1)],[v(j) v(j+1)],xa));
    for c = 1:numel(xa)
        r = max(min(yv(c),bot),top);
        S(r, xa(c), :) = 255;
    end
end
S = insertText(S, [4 0], strip.label, 'FontSize', max(9,round(0.78*fs)), ...
    'BoxOpacity', 0, 'TextColor', [170 190 210]);
end

% =====================================================================
function S = stripCursor(S, k, n)
%stripCursor  Where in the cycle the picture above is.
if n < 2, return; end
x  = round(linspace(1, size(S,2), n));
xc = x(min(max(k,1),n));
S(16:end-6, xc, 1) = 120;  S(16:end-6, xc, 2) = 240;  S(16:end-6, xc, 3) = 255;
end

% =====================================================================
function F = frameBorder(F, g, col)
%frameBorder  A border round the pictures, for a frame the animal moved through.
x0 = g.panel(1,1)-3;  x1 = g.panel(end,1)+g.panelW+2;
y0 = g.panelY-3;      y1 = g.panelY+g.panelH+2;
x0 = max(x0,1); y0 = max(y0,1); x1 = min(x1,size(F,2)); y1 = min(y1,size(F,1));
for c = 1:3
    F(y0:y0+2, x0:x1, c) = col(c);   F(y1-2:y1, x0:x1, c) = col(c);
    F(y0:y1, x0:x0+2, c) = col(c);   F(y0:y1, x1-2:x1, c) = col(c);
end
end

% =====================================================================
function n = nameOf(p)
[~,f,e] = fileparts(p);  n = [f e];
end
%------------- END OF CODE --------------
