%getVesselTypeGuess - the artery/vein guess of setVesselTypes, without the GUI
%
%   [sld,ctx] = getVesselTypeGuess(results,fName)
%   guess     = getVesselTypeGuess('vote',sld,ctx)
%
% DESCRIPTION
%   The computational half of setVesselTypes, extracted for the same reason
%   getVascularTree was extracted from setVascularTree: the wrapper opens a paint GUI and
%   blocks on it, so anything left inside it can only be checked by a person looking at a
%   screen.  What is here is the part that decides WHAT the automatic guess is driven by
%   and what that guess comes out as; the wrapper keeps the widgets, the painting and the
%   file I/O.
%
%   THE FAILURE THIS EXISTS TO REMOVE IS A SILENT ONE.  Every driving slider is built
%   from a per-segment column looked up BY NAME, and an absent column simply produced no
%   slider - so a product that had never been through the step which writes those columns
%   opened a labelling GUI with nothing driving it, every segment came out "Uncertain",
%   and nothing anywhere said why.  That was already true on the speckle branch for a
%   contrast product that skipped pulsatility; a fluorescence product reproduced it for a
%   whole modality, because its columns are pv* or bs* and never ps*.  So this function
%   returns, beside the sliders it COULD build, the columns it looked for and did not
%   find and the name of the step that writes them - and the wrapper says so out loud.
%
%   WHICH COLUMNS DRIVE IT IS NOT DECIDED HERE.  getVascularCues answers that from the
%   product's name, once, and defaultFlowParams reads the same answer for the hierarchy,
%   so the labelling GUI and the tree editor cannot come to disagree about what orders
%   flow on a product.
%
%   THE SLIDERS, IN THE ORDER THEY ARE OFFERED
%     1  pulsatility: gradient   the pulsatility index at this segment's endpoint minus
%                                the next one's - the arterial/venous step ALONG a vessel
%     2  pulsatility: magnitude  the index itself; arteries are the pulsatile ones
%     3  BFI / diameter          flow per unit calibre.  SPECKLE ONLY - there is no BFI
%                                twin on the fluorescence branch [D10], so this is absent
%                                there because the quantity is, not because it went missing
%     4  std(BFI) / BFI          how variable that flow is.  Speckle only, same reason
%     5+ earlier <arrival>       one per arrival/transit column the product carries, each
%                                NEGATED so that the high end of every slider is the
%                                arterial end: earlier arrival, higher pulsatility, higher
%                                flow.  One direction convention across the whole panel
%
%   WHICH ONES START ON: the strongest cue this product has.  With a pulsatility index
%   that is sliders 1 and 2 and the rest start off, which is what the speckle branch has
%   always done; with only arrivals - a bolus - the arrival sliders start on instead,
%   because a panel where every slider is off is a guess that guesses nothing.
%
%   FIVE SLIDERS WERE RETIRED HERE, AND THEY HAD NOT BEEN CONNECTED TO ANYTHING FOR SOME
%   TIME.  The block used to end with t0B / t5B / tUB / t90 / tMB - bolus landmarks off
%   the step that runCTTH replaced.  Nothing in this library writes any of the five (the
%   rebuilt transit step DROPS both generations of its own columns from a product it is
%   handed), so they were five sliders that could never appear.  Their replacement is not
%   a rename: a bolus product is driven by bsDelay and bsMtt, which are a tracer arrival
%   and a mean transit time rather than a re-spelling of the landmark set.
%
% INPUTS
%   results  RESULTS struct of a segmented product.  Needs sMetrics with columns idx and
%            category; every driving column is optional and its absence is reported.
%   fName    the RESULTS member's path - what getVascularCues reads.
%
% OUTPUTS
%   sld      1xN struct array, ONE PER SLIDER THAT CAN BE BUILT (a slider whose source
%            column is absent is not returned at all, so the caller builds a widget per
%            element and never an inert one):
%              .name     the label the panel shows
%              .seg      [nRow x 1] per-segment value, NaN where it does not apply
%              .rng      [min max] full range
%              .lo,.hi   default thumb positions
%              .enabled  does it start switched on
%              .guesses  does it feed the guess (true for all of them today)
%              .source   cellstr of the sMetrics columns it reads
%   ctx      what the vote below needs, and what the wrapper says out loud:
%              .idxs .m5 .N     the segment index vector, the valid category-5 rows and
%                               the number of rows - the vote's own context
%              .driving         how many returned sliders are driven by a flow column
%              .missing         cellstr of flow columns this product should have carried
%                               and does not
%              .writtenBy       registry LABEL of the step that writes them ('' if none)
%              .what            one phrase naming what orders flow here ('' if nothing)
%              .cues            the getVascularCues answer, so a caller needs one call
%
%   guess    [nRow x 1] the automatic guess at the sliders' CURRENT thumbs: a segment
%            above a slider's high thumb scores +1 (artery), below its low thumb -1
%            (vein), between them 0; the scores of every enabled guess slider are summed
%            and spread over the segment's two endpoints.  ONE implementation, called
%            both by the wrapper's live recompute and by a test.
%
% See also: getVascularCues, defaultFlowParams, setVesselTypes, getVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function varargout = getVesselTypeGuess(varargin)

% the vote is reachable as its own entry point for the same reason validate and filter
% are on wbStepRegistry: it is the RULE the live GUI runs on every thumb drag, and a rule
% nothing can call is a rule nothing can check
if nargin>=1 && (ischar(varargin{1}) || isstring(varargin{1})) && strcmp(varargin{1},'vote')
    varargout{1} = voteGuess(varargin{2},varargin{3});
    return
end

[sld,ctx] = buildSliders(varargin{:});
varargout = {sld,ctx};
end

% =====================================================================
function [sld,ctx] = buildSliders(results,fName)

c = getVascularCues(fName);
M = results.sMetrics;
vn = M.Properties.VariableNames;

idxs = M{:,'idx'};
ctgs = M{:,'category'};
N    = numel(idxs);
% Valid category-5 rows the auto-guess may write to (need a 2-endpoint span
% idxs(i):idxs(i)+1 that stays inside the guess vector).
m5 = ~isnan(idxs) & idxs>=1 & (idxs+1)<=N & (ctgs==5);
ii = find(m5);                  % shared by every slider

sld     = emptySlider();
missing = {};

% ---- 1 and 2: the pulsatility index, as a gradient and as a magnitude -------------
% A bolus product has neither: a tracer says WHEN a vessel filled and nothing about how
% hard it is pulsing, so these two are absent there because the quantity is.
if ~isempty(c.pulsatility)
    if any(strcmp(c.pulsatility,vn))
        piCol = double(M{:,c.pulsatility});
        dv = nan(N,1);  mv = nan(N,1);
        dv(ii) = piCol(idxs(ii)) - piCol(idxs(ii)+1);   % this endpoint - next
        mv(ii) = piCol(idxs(ii));                        % the index here

        % slider 1: symmetric range spanning the endpoint differences (the category-5 vs
        % category-3 difference along a segment), thumbs at 0
        R1 = max(abs(dv),[],'omitnan');
        if ~isfinite(R1) || R1==0, R1 = 1; end
        sld(end+1) = oneSlider('pulsatility: gradient',dv,[-R1 R1],0,0,true,{c.pulsatility});

        % slider 2: range over the observed magnitudes, thumbs at the medians of the
        % category-5 index below / above a split point.  THE SPLIT POINT IS THE DIMMEST
        % DECILE'S median index where the product has a brightness column, and the
        % range's thirds where it has none - a fluorescence product has no per-segment
        % brightness at all [D10], so this is a stated default rather than a value that
        % silently became NaN.
        [r2,lo2,hi2] = magnitudeThumbs(mv,piCol,ctgs,M,vn,c);
        sld(end+1) = oneSlider('pulsatility: magnitude',mv,r2,lo2,hi2,true,{c.pulsatility});
    else
        missing{end+1} = c.pulsatility;
    end
end

% ---- 3: flow per unit calibre.  Both columns or neither --------------------------
if ~isempty(c.brightness) && numel(c.caliber)>=2 && ...
        any(strcmp(c.brightness,vn)) && any(strcmp('diameter',vn))
    b = double(M{:,c.brightness});  diam = double(M{:,'diameter'});
    r3 = nan(N,1);
    r3(ii) = b(idxs(ii)) ./ diam(idxs(ii));
    r3(~isfinite(r3)) = NaN;                 % guard divide-by-zero / missing
    [rr,lo,hi] = autoRangeThumbs(r3,m5);
    sld(end+1) = oneSlider([c.brightness ' / diameter'],r3,rr,lo,hi,false, ...
                           {c.brightness,'diameter'});
end

% ---- 4: how variable that flow is; thumbs at the parenchymal median --------------
if ~isempty(c.brightness) && ~isempty(c.brightnessStd) && ...
        any(strcmp(c.brightness,vn)) && any(strcmp(c.brightnessStd,vn))
    b = double(M{:,c.brightness});  sb = double(M{:,c.brightnessStd});
    r4 = nan(N,1);
    r4(ii) = sb(idxs(ii)) ./ b(idxs(ii));
    r4(~isfinite(r4)) = NaN;
    rP = sb(ctgs==1) ./ b(ctgs==1);          % ratio over parenchyma
    [rr,lo,hi] = thumbsAtValue(r4,m5,median(rP(isfinite(rP))));
    sld(end+1) = oneSlider([c.brightnessStd ' / ' c.brightness],r4,rr,lo,hi,false, ...
                           {c.brightness,c.brightnessStd});
end

% ---- 5+: one per arrival column, NEGATED so the high end is the arterial end ------
% THE STRONGEST AVAILABLE CUE IS THE ONE THAT STARTS ON, and that rule is what keeps a
% bolus product from opening with a guess that guesses nothing.  On a product with a
% pulsatility index the two sliders above are the primary discriminator and the arrival
% ones start off, exactly as the speckle branch has always behaved; on a bolus there IS
% no pulsatility column, so leaving the arrivals off as well would have opened the panel
% with every slider greyed and every vessel undecided - the same silent failure this
% function exists to remove, wearing a different hat.
arrivalOn = isempty(c.pulsatility) || ~any(strcmp(c.pulsatility,vn));
for k = 1:numel(c.arrival)
    col = c.arrival{k};
    if ~any(strcmp(col,vn)), missing{end+1} = col; continue; end %#ok<AGROW>
    xv = -double(M{:,col});                  % earlier arrival -> higher score
    sv = nan(N,1);  sv(ii) = xv(idxs(ii));
    sv(~isfinite(sv)) = NaN;
    mPar = xv(ctgs==1);                      % thumbs at the parenchymal median
    [rr,lo,hi] = thumbsAtValue(sv,m5,median(mPar(isfinite(mPar))));
    sld(end+1) = oneSlider(['earlier ' c.arrivalLabels{k}],sv,rr,lo,hi,arrivalOn,{col}); %#ok<AGROW>
end

% ---- what the wrapper has to be able to say --------------------------------------
driving = 0;
flowCols = [c.arrival(:); {c.pulsatility}];
flowCols = flowCols(~cellfun(@isempty,flowCols));
for k = 1:numel(sld)
    if any(ismember(sld(k).source,flowCols)), driving = driving+1; end
end

ctx = struct('idxs',idxs,'m5',m5,'N',N,'driving',driving,'missing',{missing}, ...
             'writtenBy',c.writtenBy,'what',c.what,'cues',c);
end

%% =========================  HELPERS (private)  ====================== %%
function g = voteGuess(sld,ctx)
%voteGuess  The automatic guess at the sliders' current thumbs.
%   value above a slider's HIGH thumb -> +1 (artery)
%   value below a slider's LOW  thumb -> -1 (vein)
%   value between the two             ->  0 (unknown)
%   Contributions from every enabled guess slider are summed and each segment's vote is
%   spread over its two endpoints idxs(i):idxs(i)+1.  A disabled slider contributes
%   nothing, so with none enabled every segment scores 0 and reads "Uncertain" - which is
%   correct, and is exactly why the wrapper has to say when that is all it could do.
votes = zeros(ctx.N,1);
for si = 1:numel(sld)
    if sld(si).guesses && sld(si).enabled && ~isempty(sld(si).seg)
        vv = sld(si).seg(:);
        votes = votes + double(vv > sld(si).hi) - double(vv < sld(si).lo);
    end
end
sel = ctx.m5;
g = accumarray(ctx.idxs(sel),   votes(sel), [ctx.N 1]) + ...
    accumarray(ctx.idxs(sel)+1, votes(sel), [ctx.N 1]);
end

% =====================================================================
function s = emptySlider()
s = struct('name',{},'seg',{},'rng',{},'lo',{},'hi',{},'enabled',{},'guesses',{},'source',{});
end

function s = oneSlider(name,seg,rng,lo,hi,enabled,source)
s = struct('name',name,'seg',seg,'rng',rng,'lo',lo,'hi',hi, ...
           'enabled',logical(enabled),'guesses',true,'source',{source});
end

% =====================================================================
function [rng,lo,hi] = magnitudeThumbs(mv,piCol,ctgs,M,vn,c)
% Range over the observed pulsatility magnitudes, thumbs at the medians of the
% category-5 index below / above a split point.  See the caller for what the split is.
thrsh = NaN;
if ~isempty(c.brightness) && any(strcmp(c.brightness,vn))
    b = double(M{:,c.brightness});
    thrsh = median(piCol(b<prctile(b(ctgs==5),10) & ctgs==5),'omitnan');
end
lo = median(piCol(ctgs==5 & piCol<thrsh),'omitnan');
hi = median(piCol(ctgs==5 & piCol>thrsh),'omitnan');
mvLo = min(mv,[],'omitnan');  mvHi = max(mv,[],'omitnan');
if ~isfinite(mvLo) || ~isfinite(mvHi) || mvLo==mvHi, mvLo=0; mvHi=1; end
if ~isfinite(lo), lo = mvLo + (mvHi-mvLo)/3;   end
if ~isfinite(hi), hi = mvLo + 2*(mvHi-mvLo)/3; end
pad = 0.02*(mvHi-mvLo);
rng = [min(mvLo,lo)-pad, max(mvHi,hi)+pad];
lo  = min(max(lo,rng(1)),rng(2));
hi  = min(max(hi,rng(1)),rng(2));
end

% =====================================================================
function [rng, lo, hi] = autoRangeThumbs(vals, vmask)
% Pick a range and default thumbs for a guess slider from the spread of its per-segment
% values (over the valid rows in vmask): range = observed min/max (padded), thumbs =
% medians of the lower / upper halves (roughly the quartiles).  Falls back to [0 1] split
% in thirds without enough data.
rr = vals(vmask);
rr = rr(isfinite(rr));
if isempty(rr)
    rng = [0 1]; lo = 1/3; hi = 2/3; return
end
rmin = min(rr);  rmax = max(rr);
if rmin==rmax, rmin = rmin-0.5; rmax = rmax+0.5; end
rmed = median(rr);
lo = median(rr(rr<rmed));            % ~ lower quartile
hi = median(rr(rr>rmed));            % ~ upper quartile
if ~isfinite(lo), lo = rmin + (rmax-rmin)/3;   end
if ~isfinite(hi), hi = rmin + 2*(rmax-rmin)/3; end
pd  = 0.02*(rmax-rmin);
rng = [rmin-pd, rmax+pd];
lo  = min(max(lo,rng(1)),rng(2));
hi  = min(max(hi,rng(1)),rng(2));
if hi<lo, hi=lo; end
end

% =====================================================================
function [rng, lo, hi] = thumbsAtValue(vals, vmask, defVal)
% Range spans the per-segment values over vmask (extended to include defVal), and BOTH
% thumbs start at defVal -- e.g. the parenchymal median, giving a single threshold there
% by default.  Falls back to the range midpoint / [0 1] when defVal or the data is
% unavailable.
rr = vals(vmask);
rr = rr(isfinite(rr));
haveData = ~isempty(rr);
haveDef  = isscalar(defVal) && isfinite(defVal);
if ~haveData && ~haveDef
    rng = [0 1]; lo = 0.5; hi = 0.5; return
end
if haveData, rmin = min(rr); rmax = max(rr); else, rmin = defVal; rmax = defVal; end
if haveDef,  rmin = min(rmin,defVal); rmax = max(rmax,defVal); end
if rmin==rmax, rmin = rmin-0.5; rmax = rmax+0.5; end
pd  = 0.02*(rmax-rmin);
rng = [rmin-pd, rmax+pd];
if haveDef, vc = defVal; else, vc = (rmin+rmax)/2; end
vc  = min(max(vc,rng(1)),rng(2));
lo  = vc;  hi = vc;
end
%------------- END OF CODE --------------
