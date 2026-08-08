%getBolusConfidence - The two confidence numbers of every bolus unit, and their factors
%
%   getNVCConfidence's twin for a bolus, and deliberately the same shape:
%
%       Conf     geometric mean of every live factor      [nUnit x 1]
%       ConfMin  the WEAKEST factor                       [nUnit x 1]
%       Trust    Conf >= s.ctthConfThreshold  AND  ConfMin >= s.ctthConfMinThreshold
%
%   Both are stored because neither answers the other's question.  With a zero factor the
%   geometric mean is zero and says nothing about the rest; ConfMin alone cannot tell four
%   mediocre factors from one fatal one.
%
%   IT IS ONE-DIMENSIONAL WHERE THE NVC ONE IS TWO.  A recording has ONE bolus, so there
%   is no epoch axis, and with no epoch axis there is no deviation group: every factor
%   here judges a unit against the recording's own geometry and its own noise, never
%   against the other units.  That is not a simplification, it is the same principle -
%   "only a part of the segments are responders" - reaching its natural form when there
%   is nothing to average over.
%
%   THE SEVEN FACTORS
%     fSnr     Step / (2*BlStd)                     min(1,stat)   BIGGER IS BETTER
%              "did a bolus arrive at all?"  The /2 IS the scale: full marks once the
%              plateau step clears two baseline SDs.  This is the mask - a pixel outside
%              the vasculature scores near zero and is cancelled by the minimum threshold,
%              with no threshold on the image and no morphology
%     fSettle  the tail SLOPE, as a share of Step per record length  ramp on
%              ctthSettleFrac.  "the curve had not settled".  A slow unit still filling at
%              the record end has a plateau step that is not its plateau, and every marker
%              divided by it is wrong by the same factor.
%              A SLOPE AND NOT A SPREAD: the first version read the SD inside the plateau
%              window, which is drift on a segment average and PHOTON NOISE on a single
%              pixel - measured, it ran 0.639-0.934 per segment and 0.000-0.496 per pixel
%              and bound 99.9 % of pixels, i.e. it had quietly become a second, badly
%              scaled copy of fSnr
%     fRise    riseClean - the share of the leading edge whose slope keeps up with the
%              rise's own mean slope.  "the rise was not clean" - motion, a partial
%              injection, a shoulder.  It counts SLOW samples rather than reversals,
%              because a shoulder does not reverse a rise, it pauses one
%     fArrive  T10 - t(1), against the pre-bolus window's own width      ramp
%              "the recording started after the bolus".  THE ONE FACTOR THAT BLAMES THE
%              ENTRY STEP RATHER THAN THE PREPARATION, and it is worth a factor of its own
%              because the fix is a different frame span, not a different mouse
%     fPeakIn  t(end) - TPeak, against two plateau windows               ramp
%              "the bolus had not peaked yet"
%     fOrder   share of adjacent level times that increase
%              "the timing came out inconsistent"
%     fInput   max(0, -Delay), against ctthInputScale                    ramp
%              "it filled before the artery did".  A unit ahead of its own input is either
%              part of the input or a mistake, and either way its transit markers are not
%              transit markers
%
%   Cth IS NOT GATED BY ANY OF THEM, AND THAT IS THE POINT.  Its failure is systematic -
%   see getBolusCthFloor - so it is gated once for the whole recording by a calibration
%   through the recording's own input, not per unit by a statistic that cannot see it.  A
%   factor that appeared to check Cth would be worse than none.
%
%   A NON-FINITE FACTOR IS A ZERO, NEVER A MISSING ONE, exactly as in getNVCConfidence: an
%   invalid trace is NaN in every marker, so every factor of that unit is 0 and the unit is
%   untrusted rather than absent from the table.
%
%   AND EVERY FACTOR CARRIES ITS OWN SENTENCE.  C.factorWords is the same length and order
%   as C.factorNames and says what a LOW score means, in words a biologist reads.  It lives
%   here, beside the loop that builds the factors, for the reason getNVCConfidence gives:
%   the factor set is not a constant (the timing family is named from s.ctthLevelPcts), and
%   a phrasebook kept anywhere else is a second list that goes stale the first time a
%   factor is added.
%
%   A PER-UNIT STATISTIC MUST BE CHECKED AT BOTH UNIT SIZES.  A pixel is not a small
%   segment - it has a different noise budget - and fSettle is the worked example: it was
%   calibrated on segments and inverted its meaning per pixel without changing its median.
%   Anything added here owes the same two-unit comparison.
%
% Syntax:
%    C = getBolusConfidence(M, names, L, s)
%
% Inputs:
%    M     - [nUnit x nMetric] the marker block built from L.blockNames.
%    names - {1 x nMetric} its metric names, i.e. L.blockNames.  The COMPLETE set.
%    L     - the layout getBolusMetrics SETUP returned.
%    s     - parameter struct.  Fields:
%              ctthConfThreshold     geometric-mean threshold  (default 0.6)
%              ctthConfMinThreshold  weakest-factor threshold  (default 0.2)
%              ctthSettleFrac        fSettle's scale, as a share of Step per record
%                                    length (default 1)
%              ctthInputScale        fInput's scale, seconds (default 0.5)
%
% Outputs:
%    C - struct with .conf .confMin (single [nUnit x 1]), .trust (logical),
%        .factors (struct of single [nUnit x 1]), .factorNames, .factorWords.
%
% Example:
%    C = getBolusConfidence(M,L.blockNames,L,s);
%    fprintf('%.0f %% of segments trusted\n',100*mean(C.trust));
%
% See also: getBolusMetrics, getBolusCthFloor, getBolusInput, getNVCConfidence, runCTTH
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function C = getBolusConfidence(M,names,L,s)

if nargin<4
    error('getBolusConfidence:nargin','Use C=getBolusConfidence(M,names,L,s).');
end
if ischar(names)||isstring(names), names=cellstr(names); end
if ~iscellstr(names) || numel(names)~=size(M,2)
    error('getBolusConfidence:badNames', ...
        'names must be the %d metric names of M - pass L.blockNames.',size(M,2));
end
n=size(M,1);
cThr   =opt(s,'ctthConfThreshold'   ,0.6 );
cMinThr=opt(s,'ctthConfMinThreshold',0.2 );
% ctthSettleFrac = 1 is an ANCHOR, not a taste: a unit scores 0 when it would still move
% a WHOLE plateau step over another record length, i.e. when it has barely begun to
% settle.  MEASURED, the reference recording's median unit sits at 0.092 - it is still
% filling at 9 % of its own step per record length - and scores 0.91.
setFrac=opt(s,'ctthSettleFrac'      ,1.0 );
inScl  =opt(s,'ctthInputScale'      ,0.5 );

col=@(w) reshape(double(M(:,idx(names,w))),n,1);
Step =col('Step');  BlStd=col('BlStd');
TPeak=col('TPeak'); Delay=col('Delay');
T10  =col(sprintf('T%d',L.pcts(1)));

f=struct();
f.fSnr   =min(1,Step./(2*BlStd));                              % bigger is better

% fSettle reads the tail SLOPE, not the tail spread.  See getBolusMetrics: the spread
% version was measured and it is a noise statistic per pixel, not a drift one.
f.fSettle=max(0,1-col('tailSlope')/setFrac);

% riseClean and orderOk are two of the FIVE LOWERCASE fields getBolusMetrics carries for
% exactly this - shape checks that are properties of the trace and could never be tree
% columns.
f.fRise  =col('riseClean');
f.fArrive=min(1,max(0,(T10-L.time(1))/max(L.guardSec,eps)));
f.fPeakIn=min(1,max(0,(L.time(end)-TPeak)/max(2*L.plateauSec,eps)));
f.fOrder =col('orderOk');
f.fInput =max(0,1-max(0,-Delay)/inScl);

fn=fieldnames(f)';
F=zeros(n,numel(fn));
for k=1:numel(fn)
    v=f.(fn{k}); v(~isfinite(v))=0; v=min(1,max(0,v));
    f.(fn{k})=single(v); F(:,k)=v;
end
conf   =exp(mean(log(F),2));      % log(0) = -Inf, exp(-Inf) = 0: one zero zeroes Conf
confMin=min(F,[],2);

C.conf=single(conf); C.confMin=single(confMin);
C.trust=conf>=cThr & confMin>=cMinThr;
C.factors=f; C.factorNames=fn;
C.factorWords=cellfun(@word,fn,'UniformOutput',false);
end

% =====================================================================
function w=word(name)
%word  What a LOW score on one factor means, in a sentence a biologist reads.
%   The voice is the library's: no factor ids on screen, no bracketed explanations, and
%   the sentence describes the OBSERVATION rather than the arithmetic.  Each is a clause
%   with its own subject and under about forty characters, so one list serves a table cell
%   and a sentence alike.  A factor added above without a sentence here reads as its own
%   name, which is visible the first time anybody looks.
switch name
    case 'fSnr',    w='no bolus arrived above the noise';
    case 'fSettle', w='the curve had not settled';
    case 'fRise',   w='the rise was not clean';
    case 'fArrive', w='the recording started after the bolus';
    case 'fPeakIn', w='the bolus had not peaked yet';
    case 'fOrder',  w='the timing came out inconsistent';
    case 'fInput',  w='it filled before the artery did';
    otherwise,      w=name;
end
end

% =====================================================================
function k=idx(names,want)
k=find(strcmp(names,want),1);
if isempty(k)
    error('getBolusConfidence:missingMetric', ...
        ['The marker block has no ''%s''. getBolusConfidence needs the COMPLETE ' ...
         'marker set - measure with a full-set layout and gate only what reaches the ' ...
         'tree.'],want);
end
end

% =====================================================================
function v=opt(s,name,dflt)
v=dflt;
if isfield(s,name)&&~isempty(s.(name)), v=double(s.(name)); end
end
