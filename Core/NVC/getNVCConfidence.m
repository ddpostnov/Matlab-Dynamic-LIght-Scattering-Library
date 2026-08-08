%getNVCConfidence - The two confidence numbers of every (unit, epoch) pair
%
%   getNVCConfidence turns the marker block runNVC has just measured into two numbers
%   per (segment, epoch) - or per (pixel, epoch) - and the boolean the collapse and the
%   metrics tables are built on:
%
%       Conf     geometric mean of every live factor      [nUnit x nEp]
%       ConfMin  the WEAKEST factor                       [nUnit x nEp]
%       Trust    Conf >= s.nvcConfThreshold  AND  ConfMin >= s.nvcConfMinThreshold
%
%   Both are stored because neither answers the other's question.  With a zero factor
%   the geometric mean is zero and carries no information about the rest; ConfMin alone
%   cannot tell four mediocre factors from one bad one.  The author's own example is the
%   reason the second threshold exists: an epoch whose T10 is negative is not a response
%   for that segment however good everything else looks, so fT10 scores 0, ConfMin is 0
%   and the pair is cancelled.
%
%   IT IS A WHOLE-ARRAY OPERATION AND IT IS CALLED ONCE, AFTER THE MARKER LOOP.  Never
%   inside a parfor over segments and never per segment: the optional deviation group of
%   §6.2 measures each epoch against the OTHER epochs of the same unit, so one segment
%   in a worker does not carry enough of the problem to compute one.  That is the whole
%   reason this is a second core rather than part of getNVCMetrics.
%
%   THE DEFAULT SET JUDGES AN EPOCH ON ITS OWN.  Seven rules, one of them a family:
%
%     fSnrBl   |St| / (2*BlStd)                     min(1, stat)  BIGGER IS BETTER
%              "is there a response at all?"  Full marks once the stimulus-window
%              response exceeds two baseline SDs - the /2 IS the scale, there is no
%              coefficient.  The median segment of the reference recording scores
%              3.385/(2*0.598) = 2.83 and caps at 1; a non-responder at 0.2 baseline SDs
%              scores 0.1 and is cancelled by the 0.2 minimum threshold.  THIS is the
%              "only a part of the segments are responders" filter, applied per
%              (segment, epoch) with no ensemble statistics at all
%     fSnrFn   |St - (Fn - Bl)| / (2*FnStd)         min(1, stat)  BIGGER IS BETTER
%              "a response, or a step to a new level?"  A trace that rises and stays up
%              has St ~ Fn - Bl, so the numerator collapses and the pair is cancelled.
%              Nothing else in the set catches that.  Median segment: 1.92, capped at 1
%     fReturn  |Fn - Bl| / ((BlStd + FnStd)/2)      max(0, 1 - stat/nvcReturnScale)
%              "did the epoch come back?"  The median segment of the reference recording
%              sits at 1.02/0.608 = 1.68, so the scale is 5 rather than 3 - at 3 the
%              MEDIAN segment would score 0.44
%     fPeak    the time of the epoch's GLOBAL maximum, graded on the T50 ramp
%     fLoss    lost time as a FRACTION of the epoch span
%     fFirst   0 for epoch 1 when s.rejectFirstEpoch, else 1 - the only binary rule left
%     fT<p>    one per requested timing level, on the ramps below
%
%   THE FIRST THREE ARE RATIOS, NOT PRODUCTS.  They were first written as products, and
%   as products they carry BFI-squared and are non-monotonic in the wrong direction: a
%   large drift with a quiet baseline would score near zero - "good" - when it is the
%   most alarming case there is.  As ratios they are dimensionless and say what they were
%   meant to say.  fSnrBl and fSnrFn also run the OPPOSITE way from every other factor:
%   bigger is better, so they are min(1, stat) and NOT the max(0, 1 - stat/scale) shape.
%
%   THE TIMING RAMPS, in the author's words.  All times are on the stimulus clock: 0 is
%   the mark, D the stimulus end, Fs the finale window's start, E the epoch's last
%   sample.  All four anchors come from the layout, which is where the geometry lives.
%     below 50 % (T10)  0 before the mark.  Graded from the mark (1) to D (0)
%     at     50 % (T50)  0 before the mark.  1 anywhere inside the stimulus window, then
%                        graded from D (1) to Fs (0)
%     above  50 % (T90)  0 before the mark.  1 up to D, then graded from D (1) to E (0)
%   A NaN time - the trace never rose - gives 0, which cancels the epoch for that
%   segment.  That is correct and it has a consequence worth stating: a segment that
%   never responds has NO trusted epoch at all, its metrics row is the mean over the
%   recording's trusted set and its Conf is near zero.  Most of a field is not a
%   responder, and getNVCEpochTrust is built on exactly that fact.
%
%   THE T50 RAMP IS ONE FUNCTION CALLED WITH TWO TIMES - the 50 % level's own and the
%   peak's - rather than two pieces of code that have to agree.
%
%   THE DEVIATION GROUP IS OPTIONAL AND OFF BY DEFAULT (s.nvcDevRules).  Author's rule of
%   06-Aug-2026: "Move deviation rules as a separate group that is enabled optionally.
%   If enabled it acts as the rest of the rules."  So when it is on these five join the
%   set exactly - both the geometric mean and the minimum - and when it is off they are
%   not computed at all.  Each is a deviation along the EPOCH axis, per unit, on a robust
%   scale:
%
%       dev = |x(k) - median_k(x)| / (1.4826 * median_k|x - median_k(x)|)
%       f   = max(0, 1 - dev/s.nvcDevScale)
%
%     fBl     Bl        fFn   Fn        fEp   Ep, the TRUE whole-epoch mean
%     fAuc    AucRel    fNoise BlStd, ONE-SIDED: only a noisier-than-usual baseline
%                              counts, so dev uses max(0, x(k) - median)
%
%   The amplitude-outlier rule is on AucRel and not on Pk - the author's choice: AucRel
%   carries amplitude and persistence together, is the most reliable of the three
%   relatives, and is not read off a single sample of a noisy trace.
%
%   ONE SCALE, AND IT IS ARITHMETIC RATHER THAN TASTE.  s.nvcDevScale defaults to 3
%   because dev sits at ~0.8 for a perfectly ordinary epoch - that is simply the expected
%   absolute deviation of a normal sample.  The five retired reject*Coef settings worked
%   at 1 because they were a hard rejection at the 1-SD tail; used as the scale of a RAMP,
%   1 would score the median epoch 0.2 on all five factors and fail every epoch against a
%   0.6 threshold.  At 3 the median epoch scores 0.73, 1.5 SD out scores 0.5, 3 SD out
%   scores 0.  The robust scale is 1.4826*mad, not std, so the outlier under test cannot
%   inflate its own denominator; a zero spread means every epoch agrees, which is dev = 0
%   and f = 1, not a division to guard.  The deviations are taken over ALL epochs, never
%   over the trusted ones - a trust decision that fed its own inputs would not converge.
%
%   TWO SCALES, NOT ONE, AND THE REASON IS MEASURED.  s.nvcDevScale (3) and
%   s.nvcReturnScale (5) cannot be the same number: an ordinary epoch sits at 0.8 on a
%   deviation statistic and at 1.68 on the return statistic, so one scale that puts the
%   first at 0.73 puts the second at 0.44.
%
%   K.lossFrac = 0.02 IS THE ONE HARDCODED CONSTANT, and it is a FRACTION of the epoch
%   span rather than an absolute half-second, because a number nobody can tune has to
%   hold on a 5 s epoch as well as a 30 s one.  0.02 of 30 s is 0.6 s, the old default to
%   within a frame.
%
%   fLoss and fFirst are properties of the EPOCH, so they are the same for every unit.
%   That is what makes "unless the entire image quality was really bad" fall out for
%   free: a lost half-second zeroes every segment's fLoss, so no area responds, and
%   getNVCEpochTrust rejects the epoch recording-wide without a rule of its own.
%
%   A NON-FINITE FACTOR IS A ZERO, NEVER A MISSING ONE.  An invalid trace is NaN in every
%   marker, so every factor of that pair is 0, Conf is 0 and the pair is untrusted.  A
%   segment with no trustworthy response has no trusted epoch, which is the correct answer
%   and not a gap in the table.
%
%   CONF IS THE GEOMETRIC MEAN OVER HOWEVER MANY FACTORS ARE LIVE - nine with the default
%   three timing levels and the deviation group off, fourteen with the group on.  So it is
%   comparable between recordings processed the same way and NOT between one with the
%   group on and one with it off.  The saved settings are what records which; that is one
%   more reason the settings file is the record.
%
%   THE INDIVIDUAL FACTORS ARE RETURNED AND NOT MEANT TO BE SAVED.  Eight to fourteen more
%   [nSeg x nEp] arrays would answer a question ConfMin plus the confidence page already
%   answer.  The pages and the epoch editor draw them from here; the producer does not
%   write them into the tree.
%
%   AND EVERY FACTOR CARRIES ITS OWN SENTENCE.  C.factorWords is the same length and the
%   same order as C.factorNames and says, in words a biologist reads, what a LOW score on
%   that factor means: 'fT10' is 'the response did not start with the stimulus'.  It lives
%   here rather than in the two windows that show it (the epoch editor and guiResponse)
%   for the reason this module keeps repeating about the count: THE FACTOR SET IS NOT A
%   CONSTANT.  It is nine or fourteen, the timing family is derived from s.nvcAreaPcts,
%   and a phrasebook kept anywhere but beside the loop that builds the factors is a second
%   list that goes stale the first time a factor is added.  Written here it cannot: the
%   words are produced from the same field names, in the same pass.
%
% Syntax:
%    C = getNVCConfidence(M, names, timeLoss, layout, s)
%
% Inputs:
%    M        - [nUnit x nMetric x nEp] the marker block the producer built from
%               layout.blockNames - segments or pixels, single or double.
%    names    - {1 x nMetric} the metric names of M, i.e. layout.blockNames.  It must
%               carry the COMPLETE marker set: the producer measures with a full-set
%               layout and gates only what reaches the tree.
%    timeLoss - [1 x nEp] seconds of dropped frames in each epoch; [] reads as none.
%    layout   - the struct getNVCMetrics SETUP returned: .aStimEnd .aFinale .aLast
%               .span and .pcts are what is read.
%    s        - parameter struct.  Fields:
%                 nvcConfThreshold     geometric-mean threshold  (default 0.6)
%                 nvcConfMinThreshold  weakest-factor threshold  (default 0.2)
%                 nvcReturnScale       fReturn's scale, in noise units (default 5)
%                 nvcDevRules          add the five deviation factors (default false)
%                 nvcDevScale          their shared scale, in robust SDs (default 3)
%                 rejectFirstEpoch     zero epoch 1's fFirst (default true)
%
% Outputs:
%    C - struct with
%          .conf        [nUnit x nEp] single, the geometric mean of every live factor
%          .confMin     [nUnit x nEp] single, the weakest factor
%          .trust       [nUnit x nEp] logical, both thresholds cleared
%          .factors     struct of [nUnit x nEp] single arrays, one per factor
%          .factorNames {1 x nF} their order, which is the order of the minimum
%          .factorWords {1 x nF} what a LOW score on each of them means, in words, in
%                       the same order - the phrasebook the two windows show
%
% Example:
%    C = getNVCConfidence(M,layout.blockNames,timeLoss,layout,s);
%    fprintf('%.0f %% of pairs trusted\n',100*mean(C.trust(:)));
%
% See also: getNVCMetrics, getNVCEpochTrust, runNVC
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 06-August-2026

%------------- BEGIN CODE --------------
function C = getNVCConfidence(M,names,timeLoss,layout,s)

% THE ONE HARDCODED CONSTANT.  Every other number in this file is a setting; this one is
% not, by the author's ruling of 06-Aug-2026 ("remove coef completely, user will have no
% control over those"), and it is a fraction so that it holds at any epoch length.
K.lossFrac=0.02;

if nargin<5
    error('getNVCConfidence:nargin', ...
        'Use C=getNVCConfidence(M,names,timeLoss,layout,s).');
end
if ~isnumeric(M) || ndims(M)>3
    error('getNVCConfidence:badBlock', ...
        'M must be a numeric [nUnit x nMetric x nEp] marker block.');
end
if ischar(names)||isstring(names), names=cellstr(names); end
if ~iscellstr(names) || numel(names)~=size(M,2)
    error('getNVCConfidence:badNames', ...
        ['names must be the %d metric names of M - pass layout.blockNames, which is ' ...
         'what the producer built the block from.'],size(M,2));
end
nUnit=size(M,1);
nEp  =size(M,3);

% ---- settings, resolved in ONE place so the caller can record what ran -----
cThr   =optScalar(s,'nvcConfThreshold'   ,0.6  ,'the geometric-mean confidence threshold');
cMinThr=optScalar(s,'nvcConfMinThreshold',0.2  ,'the weakest-factor threshold');
retScl =optScalar(s,'nvcReturnScale'     ,5    ,'fReturn''s scale, in noise units');
devOn  =optFlag  (s,'nvcDevRules'        ,false,'whether the deviation group is on');
devScl =optScalar(s,'nvcDevScale'        ,3    ,'the deviation scale, in robust SDs');
first  =optFlag  (s,'rejectFirstEpoch'   ,true ,'whether epoch 1 is rejected');
if ~(retScl>0), error('getNVCConfidence:nvcReturnScale', ...
        's.nvcReturnScale must be positive, not %g.',retScl); end
if ~(devScl>0), error('getNVCConfidence:nvcDevScale', ...
        's.nvcDevScale must be positive, not %g.',devScl); end

% ---- the marker planes -----------------------------------------------------
% Every factor is written on [nUnit x nEp] arrays, so the whole file is one expression
% per rule and there is no loop over segments anywhere in it.
Bl   =plane(M,names,'Bl'    ,nUnit,nEp);
BlStd=plane(M,names,'BlStd' ,nUnit,nEp);
Fn   =plane(M,names,'Fn'    ,nUnit,nEp);
FnStd=plane(M,names,'FnStd' ,nUnit,nEp);
St   =plane(M,names,'St'    ,nUnit,nEp);
tPeak=plane(M,names,'tPeak' ,nUnit,nEp);

% ---- the seven always-on rules --------------------------------------------
f=struct();
f.fSnrBl =min(1,abs(St)./(2*BlStd));                       % bigger is better
f.fSnrFn =min(1,abs(St-(Fn-Bl))./(2*FnStd));               % bigger is better
f.fReturn=max(0,1-(abs(Fn-Bl)./((BlStd+FnStd)/2))/retScl);
f.fPeak  =timeRamp(tPeak,layout.aStimEnd,layout.aFinale);  % THE T50 RAMP, second call

% fLoss is a fraction of the epoch SPAN, never an absolute number of seconds.  Turning it
% back into seconds because a fraction looks surprising would make the rule mean two
% different things on a 30 s epoch and a 5 s one.
tl=double(timeLoss(:))';
if isempty(tl), tl=zeros(1,nEp); end
if numel(tl)~=nEp
    error('getNVCConfidence:timeLoss', ...
        'timeLoss must have one entry per epoch (%d), not %d.',nEp,numel(tl));
end
f.fLoss=repmat(max(0,1-(tl/layout.span)/K.lossFrac),nUnit,1);

f.fFirst=ones(nUnit,nEp);
if first && nEp>=1, f.fFirst(:,1)=0; end

% ---- one timing factor per requested level, on its own ramp ---------------
for k=1:numel(layout.pcts)
    p=layout.pcts(k);
    T=plane(M,names,sprintf('T%d',p),nUnit,nEp);
    if p<50
        r=timeRamp(T,layout.aMark   ,layout.aStimEnd);     % 1 at the mark, 0 at D
    elseif p==50
        r=timeRamp(T,layout.aStimEnd,layout.aFinale);      % 1 to D, 0 at the finale
    else
        r=timeRamp(T,layout.aStimEnd,layout.aLast);        % 1 to D, 0 at the last sample
    end
    f.(sprintf('fT%d',p))=r;
end

% ---- the optional deviation group -----------------------------------------
% When it is off these are not computed at all, which is what "off by default" has to
% mean for a factor set whose size decides what Conf means.
if devOn
    f.fBl   =devFactor(plane(M,names,'Bl'    ,nUnit,nEp),devScl,false);
    f.fFn   =devFactor(plane(M,names,'Fn'    ,nUnit,nEp),devScl,false);
    f.fEp   =devFactor(plane(M,names,'Ep'    ,nUnit,nEp),devScl,false);
    f.fAuc  =devFactor(plane(M,names,'AucRel',nUnit,nEp),devScl,false);
    f.fNoise=devFactor(plane(M,names,'BlStd' ,nUnit,nEp),devScl,true);   % one-sided
end

% ---- the two numbers -------------------------------------------------------
% A non-finite factor is a zero: an invalid trace is NaN in every marker, and a pair with
% no trustworthy response must fail rather than go missing.
fn=fieldnames(f)';
F =zeros(nUnit,nEp,numel(fn));
for k=1:numel(fn)
    v=f.(fn{k});
    v(~isfinite(v))=0;
    v=min(1,max(0,v));
    f.(fn{k})=single(v);
    F(:,:,k)=v;
end
% log(0) is -Inf and exp(-Inf) is 0, so one zero factor zeroes Conf exactly as a product
% would - and the mean of the logs is what keeps fourteen factors from underflowing.
conf   =exp(mean(log(F),3));
confMin=min(F,[],3);

C.conf       =single(conf);
C.confMin    =single(confMin);
C.trust      =conf>=cThr & confMin>=cMinThr;
C.factors    =f;
C.factorNames=fn;
C.factorWords=cellfun(@factorWord,fn,'UniformOutput',false);
end

% =====================================================================
function w=factorWord(name)
%factorWord  What a LOW score on one factor means, in a sentence.
%
%   THE ONE PLACE THIS LIBRARY SAYS WHAT A CHECK IS ABOUT.  The epoch editor puts it in
%   its table and guiResponse puts it in its numbers panel, and neither of them owns a
%   list of its own: a factor added above without a sentence here reads as its own name,
%   which is visible the first time anybody looks, rather than as somebody else's
%   sentence, which is not.
%
%   The voice is the library's: no factor ids on screen, no bracketed explanations, and
%   the sentence describes the OBSERVATION rather than the arithmetic - "the rise did not
%   follow the stimulus", not "T10 sits late on its ramp".
%
%   EACH ONE IS A CLAUSE WITH ITS OWN SUBJECT, AND SHORT ENOUGH TO BE A TABLE CELL.  It
%   has to read on its own in the epoch editor's "weakest check" column and inside a
%   sentence in the same window's summary line - "its weakest check, at 0.26, is that the
%   peak fell outside the window" - so it carries a subject and stays under about forty
%   characters.  A second, shorter list for the column would be the drift this whole
%   arrangement exists to prevent.
%
%   THE TIMING FAMILY IS PARSED, NOT LISTED, because its members are named from
%   s.nvcAreaPcts and a protocol asking for [5 25 75 95] has four of them.  Which side of
%   50 % a level sits on is what decides its ramp above, so it is what decides its
%   sentence here.
switch name
    case 'fSnrBl'
        w='the response was no bigger than the noise';
    case 'fSnrFn'
        w='flow settled at a new level';
    case 'fReturn'
        w='flow had not come back to rest';
    case 'fPeak'
        w='the peak fell outside the window';
    case 'fLoss'
        w='the camera dropped frames';
    case 'fFirst'
        w='the first repetition is never used';
    case 'fBl'
        w='the resting level was unlike the others';
    case 'fFn'
        w='the end level was unlike the others';
    case 'fEp'
        w='the overall level was unlike the others';
    case 'fAuc'
        w='the response size was unlike the others';
    case 'fNoise'
        w='the baseline was noisier than usual';
    otherwise
        %AN EMPTY MATCH IS AN EMPTY CELL, and `if isnan([])` is false rather than true -
        %so the miss is tested on the token itself and never on a number derived from it.
        tok=regexp(name,'^fT(\d+)$','tokens','once');
        if isempty(tok)
            w=name;                 % a factor added above without a sentence here
            return
        end
        p=str2double(tok{1});
        if p<50
            w='the rise did not follow the stimulus';
        elseif p==50
            w='half the response arrived late';
        else
            w='it was still building at the end';
        end
end
end

% =====================================================================
function f=timeRamp(x,xOne,xZero)
%timeRamp  The one graded rule every timing factor and fPeak are written on.
%   1 at xOne and everywhere before it, 0 at xZero and everywhere after, linear between,
%   and 0 for ANY time before the stimulus mark or missing altogether.  The pre-mark
%   clause is not the clamp: at a negative time the linear piece would read above 1, and
%   a response that had already happened is not a good response, it is drift.
if xZero<=xOne
    f=double(x<=xOne);                    % a degenerate geometry is a step, not a NaN
else
    f=min(1,max(0,(xZero-x)/(xZero-xOne)));
end
f(~(x>=0))=0;                             % before the mark, or NaN
end

% =====================================================================
function f=devFactor(x,scale,oneSided)
%devFactor  |x - median| in robust SDs of the across-EPOCH spread, mapped to 1..0.
%   The scale is 1.4826*mad rather than std so the outlier under test cannot inflate its
%   own denominator, and it is computed inline rather than through mad() so the core
%   needs no toolbox.  A ZERO SPREAD IS AGREEMENT, NOT A DIVISION TO GUARD: every epoch
%   identical gives dev = 0 and f = 1.
med=median(x,2,'omitnan');
rs =1.4826*median(abs(x-med),2,'omitnan');
d  =x-med;
if oneSided, d=max(0,d); else, d=abs(d); end
dev=zeros(size(x));
ok =rs>0;
dev(ok,:)=d(ok,:)./rs(ok);
f=max(0,1-dev/scale);
f(isnan(x))=NaN;                          % a NaN marker is a zero factor, not a 1
end

% =====================================================================
function P=plane(M,names,want,nUnit,nEp)
%plane  One metric of the block as a [nUnit x nEp] double, named in the error if absent.
%   A missing marker means the producer measured with a GATED layout; the confidence
%   needs the complete set, so say that rather than fail on an empty index.
k=find(strcmp(names,want),1);
if isempty(k)
    error('getNVCConfidence:missingMetric', ...
        ['The marker block has no ''%s''. getNVCConfidence needs the COMPLETE marker ' ...
         'set - measure with a full-set layout and gate only what reaches the tree.'], ...
        want);
end
P=reshape(double(M(:,k,:)),nUnit,nEp);
end

% =====================================================================
function v=optScalar(s,name,dflt,what)
%optScalar  An optional finite numeric setting with a documented default.
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=double(s.(name)); end
if ~(isscalar(v)&&isfinite(v))
    error('getNVCConfidence:badSetting', ...
        's.%s must be a finite scalar (%s).',name,what);
end
end

% =====================================================================
function v=optFlag(s,name,dflt,what)
%optFlag  An optional logical setting with a documented default.
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=s.(name); end
if ~(isscalar(v)&&(islogical(v)||(isnumeric(v)&&isfinite(v))))
    error('getNVCConfidence:badSetting', ...
        's.%s must be true or false (%s).',name,what);
end
v=logical(v);
end
