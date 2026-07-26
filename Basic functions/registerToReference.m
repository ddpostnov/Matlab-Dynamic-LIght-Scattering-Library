%registerToReference - modality-agnostic cross-recording image registration engine
%
%   [tforms, diagOut] = registerToReference(images)
%   [tforms, diagOut] = registerToReference(images, s)
%
% DESCRIPTION
%   Registers a set of images to a common reference (images{1}) and returns the
%   estimated geometric transforms plus per-image diagnostics.  It is a pure
%   function - images in, transforms out - and knows nothing about any data
%   model, file naming or file I/O; a caller applies the returned transforms to
%   its own pristine data and writes its own reports.
%
%   For every follower it forms three candidate transforms - identity, an
%   intensity-based imregtform and a correlation-based imregcorr - and, in
%   'silent' mode, keeps whichever gives the smallest overlay difference (Delta)
%   to the reference (or a forced method).  In 'gui' mode it shows the classic
%   4-panel chooser (original / intensity / correlation / Start manual) and the
%   user decides, with the manual landmark tool as a fallback.
%
%   Estimation-time enhancement is OWNED here but never leaks out: with
%   'enhance' true each image is passed through enhanceForRegistration (using its
%   mask) purely to estimate and compare; the transforms are meant to be applied
%   by the caller to the ORIGINAL, un-enhanced data.  diagOut(k).warped holds the
%   warped ENHANCED proxy for the caller's report overlays only - never analysis
%   data.
%
% INPUT
%   images   cell array of 2-D images.  images{1} is the reference (fixed) and
%            must be non-empty.  Any images{k} may be [] -> tforms{k}=[], skipped.
%
% OPTIONAL s (struct; any missing field takes its default)
%   .mode          'silent' (default) | 'gui'
%   .method        'auto' (default) | 'intensity' | 'correlation'  (silent only).
%                  'auto' picks the smallest-Delta candidate.
%   .allowManual   logical, default true (gui only) - show the "Start manual" button.
%   .tFormType     imregtform type, default 'affine'.
%   .optimizer     imregtform optimizer, default from imregconfig('monomodal').
%   .metric        imregtform metric,    default from imregconfig('monomodal').
%   .rotationLimit deg; an intensity/correlation candidate rotating more than this
%                  is rejected (cannot win in silent mode, flagged in gui mode).
%                  [] or absent = no limit.
%   .masks         cell parallel to images: the region used by the internal
%                  enhancement (percentile clip).  Default all-true.  The overlay
%                  Delta is always computed over the whole output frame.
%   .enhance       logical, default true.  true = run enhanceForRegistration on
%                  each image internally.  false = images are already
%                  registration-ready proxies (the LSCI wrapper pre-enhances so it
%                  can reuse the proxies for its own report overlays).
%
% OUTPUT
%   tforms   cell parallel to images of MODERN 2-D transform objects
%            (affinetform2d / simtform2d / projtform2d).  tforms{1} = identity
%            affinetform2d(eye(3)).  Empty for skipped ([]) inputs.
%   diagOut  struct array parallel to images (the contract's "diag"):
%              .method   'original' | 'intensity' | 'correlation' | 'manual'
%              .reason   short human-readable selection reason
%              .delta    [dOrig dIntensity dCorrelation] whole-frame overlay diffs
%              .rejected logical(1,3); which candidates were over the rotation limit
%              .warped   the chosen warped ENHANCED image (diagnostics only).
%                        diagOut(1).warped is the enhanced reference, so a caller
%                        can build reference-vs-follower overlays.
%
% DEPENDS ON
%   enhanceForRegistration, manualByPointRegistration, and MATLAB's Image
%   Processing Toolbox (imregtform/imregconfig, imregcorr, imwarp, imref2d,
%   imgaussfilt, transformPointsForward, imshowpair).
%
% See also: enhanceForRegistration, manualByPointRegistration, registerRetinaLSCI
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 23-July-2026

function [tforms, diagOut] = registerToReference(images, s)

if nargin<2 || isempty(s), s = struct(); end
if ~iscell(images) || isempty(images)
    error('registerToReference:images','images must be a non-empty cell array.');
end
if isempty(images{1})
    error('registerToReference:noReference','images{1} (the reference) must be non-empty.');
end

% ---- resolve options -------------------------------------------------------
mode        = lower(getOpt(s,'mode','silent'));
method      = lower(getOpt(s,'method','auto'));
allowManual = logical(getOpt(s,'allowManual',true));
tFormType   = getOpt(s,'tFormType','affine');
enhance     = logical(getOpt(s,'enhance',true));
masks       = getOpt(s,'masks',{});
if isfield(s,'optimizer') && ~isempty(s.optimizer) && isfield(s,'metric') && ~isempty(s.metric)
    optimizer = s.optimizer; metric = s.metric;
else
    [optimizer,metric] = imregconfig('monomodal');
end
rotationLimit = [];
if isfield(s,'rotationLimit') && isnumeric(s.rotationLimit) && isscalar(s.rotationLimit)
    rotationLimit = abs(s.rotationLimit);
end
if ~any(strcmp(mode,{'silent','gui'}))
    error('registerToReference:mode','s.mode must be ''silent'' or ''gui''.');
end
if ~any(strcmp(method,{'auto','intensity','correlation'}))
    error('registerToReference:method','s.method must be ''auto'', ''intensity'' or ''correlation''.');
end

% ---- containers ------------------------------------------------------------
n       = numel(images);
tforms  = cell(size(images));
diagOut = repmat(struct('method','','reason','','delta',[NaN NaN NaN], ...
    'rejected',false(1,3),'warped',[]), size(images));

% ---- reference -------------------------------------------------------------
ref  = prepImage(images{1}, enhance, localMask(masks,1,size(images{1})));
Rout = imref2d(size(ref));
identity = affinetform2d(eye(3));
tforms{1}          = identity;
diagOut(1).method  = 'original';
diagOut(1).reason  = 'reference';
diagOut(1).delta   = [0 NaN NaN];
diagOut(1).warped  = ref;

% ---- followers -------------------------------------------------------------
for k = 2:n
    if isempty(images{k}), continue; end
    mov = prepImage(images{k}, enhance, localMask(masks,k,size(images{k})));

    % three candidates + their warped overlays and the whole-frame Delta
    fSize = floor(min(size(mov))/20)*2+1;
    [tform1,ok1] = tryIntensity(mov,ref,fSize,tFormType,optimizer,metric);
    [tform2,ok2] = tryCorrelation(mov,ref);
    w0 = imwarp(mov,identity,'OutputView',Rout,'FillValues',0);
    w1 = imwarp(mov,tform1,  'OutputView',Rout,'FillValues',0);
    w2 = imwarp(mov,tform2,  'OutputView',Rout,'FillValues',0);
    delta    = [sumAbs(ref,w0), sumAbs(ref,w1), sumAbs(ref,w2)];
    rejected = [false, ~ok1 || overRotated(tform1,rotationLimit), ...
                       ~ok2 || overRotated(tform2,rotationLimit)];
    cand      = {identity, tform1, tform2};
    warpedAll = {w0, w1, w2};
    nameShort = {'original','intensity','correlation'};

    if strcmp(mode,'gui')
        [best,reason] = guiSelect(ref,mov,w1,w2,rejected,allowManual);
        if best==4      % manual
            [mtf,mwarped] = manualByPointRegistration(ref,mov,'sideBySide');
            tforms{k}         = toModernTform(mtf);
            diagOut(k).method = 'manual';
            diagOut(k).reason = reason;
            diagOut(k).delta  = delta;
            diagOut(k).rejected = rejected;
            diagOut(k).warped = mwarped;
            continue
        end
    else
        % ---- silent auto-pick / forced method ------------------------------
        sel = delta; sel(rejected) = Inf;
        switch method
            case 'intensity',   best = 2; reason = 'forced';
            case 'correlation', best = 3; reason = 'forced';
            otherwise,          [~,best] = min(sel); reason = 'smallest \Delta';
        end
        if rejected(best)
            best = 1; reason = 'forced method over rotation limit -> original';
        end
    end

    tforms{k}           = toModernTform(cand{best});
    diagOut(k).method   = nameShort{best};
    diagOut(k).reason   = reason;
    diagOut(k).delta    = delta;
    diagOut(k).rejected = rejected;
    diagOut(k).warped   = warpedAll{best};
end
end

% ===========================================================================
function out = prepImage(img, enhance, mask)
%prepImage  Estimation-time proxy: optionally vessel-emphasise, always double.
if enhance
    out = double(enhanceForRegistration(img,'mask',mask));
else
    out = double(img);
end
end

% ===========================================================================
function m = localMask(masks,k,sz)
%localMask  Mask for image k (all-true when none supplied).
if iscell(masks) && numel(masks)>=k && ~isempty(masks{k})
    m = logical(masks{k});
else
    m = true(sz);
end
end

% ===========================================================================
function d = sumAbs(a,b)
%sumAbs  Whole-frame overlay difference (matches the LSCI wrapper's Delta).
d = sum(abs(double(a(:))-double(b(:))));
end

% ===========================================================================
function [tform,ok] = tryIntensity(mov,ref,fSize,tFormType,optimizer,metric)
%tryIntensity  Intensity-based imregtform on Gaussian-smoothed copies.  A failure
%   (divergence/error) returns identity with ok=false so it cannot win.
ok = true;
ws = warning('off','all');
try
    tform = imregtform(imgaussfilt(mov,fSize/5),imgaussfilt(ref,fSize/5), ...
        tFormType,optimizer,metric);
catch
    ok = false; tform = affinetform2d(eye(3));
end
warning(ws);
end

% ===========================================================================
function [tform,ok] = tryCorrelation(mov,ref)
%tryCorrelation  Phase-correlation imregcorr (similarity).  A failure returns
%   identity with ok=false so it cannot win.
ok = true;
ws = warning('off','all');
try
    tform = imregcorr(mov,ref,'similarity');
    if isempty(tform), ok = false; tform = affinetform2d(eye(3)); end
catch
    ok = false; tform = affinetform2d(eye(3));
end
warning(ws);
end

% ===========================================================================
function tf = overRotated(tform,rotationLimit)
%overRotated  True if a 2-D transform's |rotation| exceeds rotationLimit (deg).
%   Convention-agnostic: measures the rotation of the unit x-axis via
%   transformPointsForward, so it works for any geometric-transform object.
if isempty(rotationLimit)
    tf = false; return
end
[ox,oy] = transformPointsForward(tform,0,0);
[ux,uy] = transformPointsForward(tform,1,0);
tf = abs(atan2d(uy-oy,ux-ox)) > rotationLimit;
end

% ===========================================================================
function t = toModernTform(t)
%toModernTform  Normalise to the modern premultiply family.  imregtform/imregcorr/
%   estgeotform2d already return modern objects; only manualByPointRegistration's
%   identity fallback is a legacy affine2d, converted here (one-off .T read).
if isa(t,'affinetform2d') || isa(t,'simtform2d') || isa(t,'rigidtform2d') || ...
        isa(t,'transltform2d') || isa(t,'projtform2d')
    return
elseif isa(t,'affine2d')
    t = affinetform2d(t.T');
elseif isa(t,'projective2d')
    t = projtform2d(t.T');
else
    error('registerToReference:tformClass','Unsupported transform class %s.',class(t));
end
end

% ===========================================================================
function [best,reason] = guiSelect(ref,mov,w1,w2,rejected,allowManual)
%guiSelect  Classic 4-panel chooser.  Returns the picked candidate index
%   (1 original, 2 intensity, 3 correlation, 4 manual) and reason 'user'.
reason = 'user';
warn1 = ''; if rejected(2), warn1 = ' [ROT>LIM]'; end
warn2 = ''; if rejected(3), warn2 = ' [ROT>LIM]'; end

f = figure(1);
f.WindowState = 'Maximized';
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile
imshowpair(rot90(ref,size(ref,2)>size(ref,1)),rot90(mov,size(mov,2)>size(mov,1)))
axis image
if isequal(size(ref),size(mov))
    title(['Original, \Delta=',num2str(round(sumAbs(ref,mov)))])
else
    title('Original, size mismatch')
end
nexttile
imshowpair(rot90(ref,size(ref,2)>size(ref,1)),rot90(w1,size(w1,2)>size(w1,1)))
axis image
title(['Intensity registration, \Delta=',num2str(round(sumAbs(ref,w1))),warn1]);
nexttile
imshowpair(rot90(ref,size(ref,2)>size(ref,1)),rot90(w2,size(w2,2)>size(w2,1)))
axis image
title(['Correlation registration, \Delta=',num2str(round(sumAbs(ref,w2))),warn2]);

ax = nexttile; axis(ax,'off');
pnl = uipanel(f,'Units','normalized','Position',ax.Position,'BorderType','none');
vals = {'Use original','Use intensity registration','Use correlation registration'};
if allowManual, vals{end+1} = 'Start manual'; end
for kk = 1:numel(vals)
    label = vals{kk};
    uicontrol(pnl,'Style','pushbutton','String',label, ...
        'Units','normalized','Position',[0.1 1-0.2*kk 0.8 0.15], ...
        'Callback',@(src,evt)buttonDone(f,label));
end
uiwait(f);
choice = getappdata(f,'choice');
close(f);
best = find(strcmp(vals,choice));
if isempty(best), best = 1; end     % window closed without a choice -> original
end

% ===========================================================================
function buttonDone(f,sel)
%buttonDone  Store the chooser selection and release uiwait.
setappdata(f,'choice',sel);
uiresume(f);
end

% ===========================================================================
function v = getOpt(s,name,default)
%getOpt  Struct field with a default (empty field -> default).
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
