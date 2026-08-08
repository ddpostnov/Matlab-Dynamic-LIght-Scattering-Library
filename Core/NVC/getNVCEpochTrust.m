%getNVCEpochTrust - Which epochs the RECORDING trusts, from responding area
%
%   getNVCEpochTrust answers the one question the per-(segment, epoch) confidence
%   cannot: which epochs of this recording are good enough to average.  It is the
%   author's rule of 06-Aug-2026 and it has no precedent in the old step, which judged
%   an epoch by an average over its segments:
%
%   > "only a part of the segments are responders.  So it might be that only 10 % of the
%   > segments get accepted for the epoch but that would still be correct.  The user
%   > should use % of total area as responding marker, possibly it should be % of
%   > interconnected responding area - so individual segments spatially distributed that
%   > got accepted just by noisy chance are also excluded.  But the point is it cannot be
%   > simply judged on average, unless the entire image quality was really bad."
%
%   AN EPOCH IS JUDGED BY AREA, NEVER BY A COUNT AND NEVER BY AN AVERAGE.
%
%       areaFrac(k) = area of the segments trusted in epoch k / area of all segments
%       epochTrust  = areaFrac >= s.nvcEpochAreaFrac
%
%   Areas are pixel counts read off results.sMap, which is the only spatial map the
%   product carries - so this is a segmented-flow rule, and the caller passes the trust
%   of the segmented flow signal (sData, else dvsData, else dvsDiameter) and says on the
%   page which one it used.
%
%   WITH s.nvcEpochConnected THE NUMERATOR IS THE LARGEST CONNECTED COMPONENT of the
%   trusted mask painted back through sMap, 8-connected.  That is the author's "%
%   interconnected responding area": a scatter of individually-trusted segments spread
%   across the field is what noise looks like when it clears a threshold by chance, and
%   it must not add up to a response.  A single coherent patch covering a tenth of the
%   field is a response; the same area in forty specks is not.
%
%   THE DENOMINATOR IS THE AREA OF EVERY SEGMENT THAT HAS A ROW, not the pixel count of
%   the whole image.  Background is not a segment that failed to respond; it is not a
%   segment.  Counting it would make the threshold mean something different on every
%   recording, because it would be reading the segmentation's coverage rather than the
%   response's extent.
%
%   THE DEFAULT 0.10 IS A GUESS AND IS MEANT TO BE MEASURED.  It is written down as such
%   in the plan: the producer session reports the measured areaFrac of every epoch of the
%   reference recording, connected and unconnected, and does NOT tune the default to
%   reproduce what the old rules accepted.
%
%   IT GATES THE COLLAPSE AND NOTHING ELSE.  The per-segment metrics table averages over
%   this recording-wide set too - the author's ruling, so that every row of the table
%   describes the same epochs - but the per-segment Trust is what feeds this rule.  So a
%   segment's own confidence chooses which epochs the RECORDING trusts, and the
%   recording's choice is what its table row averages over.
%
%   DEPENDS ON
%     MATLAB Image Processing Toolbox (bwlabel), and only on the connected path.  sMap
%     comes from runSegmentation, which needs the same toolbox, so this adds nothing.
%
% Syntax:
%    [epochTrust, areaFrac] = getNVCEpochTrust(trust, sMap, s)
%
% Inputs:
%    trust - [nSeg x nEp] logical, the per-(segment, epoch) trust getNVCConfidence
%            returned for the segmented flow signal.
%    sMap  - [Y x X] segment label map; 0 is background, labels are 1..nSeg and index
%            the rows of `trust`.
%    s     - parameter struct.  Fields:
%              nvcEpochAreaFrac  the responding-area share an epoch needs (default 0.10)
%              nvcEpochConnected count only the largest connected trusted region
%                                (default true)
%
% Outputs:
%    epochTrust - [nEp x 1] logical, the recording-wide trusted set.
%    areaFrac   - [nEp x 1] double, the share the decision was made on.  Double rather
%                 than single because it is a handful of numbers and the producer casts
%                 it into the tree exactly as it casts epochStart.
%
% Example:
%    [epochTrust,areaFrac] = getNVCEpochTrust(C.trust,results.sMap,s);
%    fprintf('%d of %d epochs trusted, best %.0f %% of the field\n', ...
%        sum(epochTrust),numel(epochTrust),100*max(areaFrac));
%
% See also: getNVCConfidence, getNVCMetrics, runNVC, bwlabel
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 06-August-2026

%------------- BEGIN CODE --------------
function [epochTrust,areaFrac] = getNVCEpochTrust(trust,sMap,s)

if nargin<3
    error('getNVCEpochTrust:nargin', ...
        'Use [epochTrust,areaFrac]=getNVCEpochTrust(trust,sMap,s).');
end
if ~(islogical(trust)||isnumeric(trust)) || ndims(trust)>2 %#ok<ISMAT>
    error('getNVCEpochTrust:badTrust', ...
        'trust must be an [nSeg x nEp] logical array.');
end
trust=logical(trust);
nSeg=size(trust,1);
nEp =size(trust,2);

if ~isnumeric(sMap) || ndims(sMap)>2 %#ok<ISMAT>
    error('getNVCEpochTrust:badMap', ...
        'sMap must be a [Y x X] numeric label map (0 is background).');
end
lab=double(sMap);
if any(lab(:)<0) || any(lab(:)~=round(lab(:)))
    error('getNVCEpochTrust:badMap', ...
        'sMap must hold whole, non-negative segment labels.');
end
maxLab=max([0;lab(:)]);
if maxLab>nSeg
    error('getNVCEpochTrust:segmentCount', ...
        ['sMap holds %d segments but the trust array has %d rows; they must be the ' ...
         'same segmentation.'],maxLab,nSeg);
end

areaThr  =optScalar(s,'nvcEpochAreaFrac',0.10,'the responding-area share an epoch needs');
connected=optFlag  (s,'nvcEpochConnected',true,'whether only the largest region counts');
if ~(areaThr>=0 && areaThr<=1)
    error('getNVCEpochTrust:nvcEpochAreaFrac', ...
        's.nvcEpochAreaFrac must lie between 0 and 1, not %g.',areaThr);
end

% ---- the denominator: every segment that has a row -------------------------
area=zeros(nSeg,1);
if maxLab>0
    counted=accumarray(lab(lab>0),1,[maxLab 1]);
    area(1:maxLab)=counted;
end
totalArea=sum(area);

areaFrac=zeros(nEp,1);
if totalArea==0
    epochTrust=false(nEp,1);
    return
end

for k=1:nEp
    sel=trust(:,k);
    if ~any(sel), continue, end
    if connected
        % Paint the trusted segments back through sMap with a label lookup - O(Y*X)
        % rather than an ismember over every pixel - and keep only the biggest piece.
        keep=[false;sel];                       % index 1 is background (label 0)
        mask=reshape(keep(lab(:)+1),size(lab));
        cc=bwlabel(mask,8);
        if any(cc(:))
            areaFrac(k)=max(accumarray(cc(cc>0),1))/totalArea;
        end
    else
        areaFrac(k)=sum(area(sel))/totalArea;
    end
end
epochTrust=areaFrac>=areaThr;
end

% =====================================================================
function v=optScalar(s,name,dflt,what)
%optScalar  An optional finite numeric setting with a documented default.
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=double(s.(name)); end
if ~(isscalar(v)&&isfinite(v))
    error('getNVCEpochTrust:badSetting', ...
        's.%s must be a finite scalar (%s).',name,what);
end
end

% =====================================================================
function v=optFlag(s,name,dflt,what)
%optFlag  An optional logical setting with a documented default.
v=dflt;
if isfield(s,name) && ~isempty(s.(name)), v=s.(name); end
if ~(isscalar(v)&&(islogical(v)||(isnumeric(v)&&isfinite(v))))
    error('getNVCEpochTrust:badSetting', ...
        's.%s must be true or false (%s).',name,what);
end
v=logical(v);
end
