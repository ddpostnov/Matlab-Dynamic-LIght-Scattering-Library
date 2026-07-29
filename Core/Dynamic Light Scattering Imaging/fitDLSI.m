%fitDLSI - Per-pixel fit of g2 to Dynamic Light Scattering Imaging models.
%
% Syntax:
%    model = fitDLSI(g2, lags, iniTau)
%    model = fitDLSI(g2, lags, iniTau, 'Name', Value, ...)
%
% Description:
%    Fits the normalized intensity autocorrelation g2 (from getNormalizedG2)
%    pixel by pixel to one of the DLSI models, returning maps of the fitted
%    parameters and goodness of fit. For each pixel the usable lag range is
%    trimmed to the monotonically decaying part of g2, the model is fitted by
%    nonlinear least squares (Curve Fitting Toolbox), and, for the multi-model
%    types, the best-fitting motion regime is selected per pixel.
%
%    Pixels are fitted in the order given by 'order' (longest decorrelation
%    time first by default). When 'isAdaptive' is true the beta and p bounds
%    are tightened on the fly from the statistics of already-fitted, well-fit
%    pixels, which stabilises the fit of the faster-decorrelating pixels.
%
%    Model types (n is the motion exponent in exp(-(x/tau)^n)):
%       'MDSN' - mixed dynamics + static + noise (ordered n=1 and 2, or
%                unordered n=0.5 and 1 combined; best of the two mixtures).
%       'DSN'  - single dynamics + static + noise (best of n = 2, 1, 0.5).
%       'DN'   - single dynamics + noise (best of n = 2, 1, 0.5).
%       'D'    - single dynamics only (best of n = 2, 1, 0.5).
%
%    Coherence-loss form ('cLoss') sets the coefficient of the 2p(1-p) static/
%    dynamic cross-term in g2 (Liu et al., Biomed. Opt. Express 15(2):579,
%    2024, Eq. 15/16): 'sqrtbeta' uses sqrt(beta) (static field not fluctuating
%    in time), 'beta' uses full beta (coherence loss on both fields). Passing a
%    cell {'sqrtbeta','beta'} fits each form fully independently - its own
%    adaptive bounds and fit history, so the mutually exclusive forms never
%    co-adapt - and keeps the higher-R^2 form per pixel. Only MDSN and DSN
%    carry this cross-term; DN and D are unaffected.
%
% Inputs:
%    g2     - Normalized intensity autocorrelation, 3-D [Y, X, nLags].
%    lags   - Lag-time vector, length nLags (e.g. (0:nLags-1)/fps, seconds).
%    iniTau - Initial decorrelation-time map [Y, X] (from getTauC), also used
%             to order the pixels.
%
% Optional Name-Value Pair Arguments:
%    'type'       - (String) Model type, see above. Default: 'MDSN'.
%    'cLoss'      - (Char or cell of chars) Coherence-loss form of the cross-
%                   term: 'sqrtbeta' (default), 'beta', or both as a cell
%                   {'sqrtbeta','beta'}. Affects MDSN and DSN only.
%    'pointsMin'  - (Numeric) Minimum number of lag points used per pixel.
%                   Default: 5.
%    'order'      - (Numeric) Linear-index fitting order, numel(iniTau) long.
%                   Default: [] (iniTau sorted descending).
%    'betaMin'    - (Numeric) Lower beta bound, scalar or [Y, X]. Default: []
%                   (derived as prctile(g2(:,:,1),0.01)-1).
%    'betaMax'    - (Numeric) Upper beta bound, scalar or [Y, X]. Default: []
%                   (1 everywhere).
%    'bs'         - (Numeric) Beta scaling (speckle sampling) factor.
%                   Default: 1.
%    'mask'       - (Logical/Numeric) Pixels to fit, [Y, X]. Default: [] (all).
%    'isAdaptive' - (Logical) Adapt beta/p bounds during fitting.
%                   Default: true.
%    'progressFcn'- (Function handle) Progress sink, progressFcn(frac,label).
%                   Default: [] - the per-batch percentage is printed to the
%                   command window exactly as it always has been, so every
%                   existing caller is unaffected.  Supplied, one bar spans all
%                   the coherence-loss passes instead.
%
% Outputs:
%    model - Lightweight fit struct (does NOT store g2) with fields:
%       .settings - fit configuration used: type, cLoss, pointsMin, lags,
%                   lagStart, order, isAdaptive, fitSet.
%       .range    - per-pixel lag range, artifact/low-point flags and the
%                   (winning form's) fitting bounds.
%       .fit      - fitted maps [Y, X]: beta, tau, p, d (regime), c (noise),
%                   goodness of fit r2 and sse, and cLoss (coherence-loss form
%                   picked per pixel: 1 = sqrtbeta, 2 = beta, NaN for DN/D).
%    The heavy g2 is intentionally not duplicated in the model: it belongs in
%    the source (*_g_d.mat) data struct. Store one or more returned models in
%    a results struct (e.g. results.mMDSN) for the *_g_r.mat results file.
%
% Example:
%    tauC  = getTauC(g2, 'fps', fps);
%    lags  = (0:size(g2,3)-1)/fps;
%    model = fitDLSI(g2, lags, tauC, 'type', 'MDSN');
%    imagesc(model.fit.tau); axis image; colorbar
%
% Dependencies: Curve Fitting Toolbox (fit/fittype/fitoptions); Parallel
%    Computing Toolbox (parfor, optional).
% See also: getNormalizedG2, getTauC
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 23-July-2026

%------------- BEGIN CODE --------------
function model = fitDLSI(g2, lags, iniTau, varargin)

p = inputParser;
p.KeepUnmatched = false;
addRequired(p, 'g2',     @(x) isnumeric(x) && ndims(x)==3);
addRequired(p, 'lags',   @(x) isnumeric(x) && isvector(x));
addRequired(p, 'iniTau', @(x) isnumeric(x) && ismatrix(x));
addParameter(p, 'type', 'MDSN', @(x) any(validatestring(x, {'MDSN','DSN','DN','D'})));
addParameter(p, 'cLoss', 'sqrtbeta', @(x) ischar(x) || iscell(x));
addParameter(p, 'pointsMin', 5, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p, 'order', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'betaMin', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'betaMax', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'bs', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'mask', [], @(x) isempty(x) || isnumeric(x) || islogical(x));
addParameter(p, 'isAdaptive', true, @(x) islogical(x) || ismember(x,[0 1]));
addParameter(p, 'progressFcn', [], @(x) isempty(x) || isa(x,'function_handle'));
parse(p, g2, lags, iniTau, varargin{:});

% Absent (the default) the per-batch percentage is printed exactly as it always
% has been, so the DLSI launchers and everything in Drafts/ are unaffected.
progressFcn = p.Results.progressFcn;

type       = validatestring(p.Results.type, {'MDSN','DSN','DN','D'});
pointsMin  = p.Results.pointsMin;
order      = p.Results.order;
bs         = p.Results.bs;
isAdaptive = logical(p.Results.isAdaptive);

% Coherence-loss form(s) for the p(1-p) cross-term (Liu et al. 2024, Eq. 15/16):
%   'sqrtbeta' -> sqrt(beta) cross-term (static field not fluctuating in time);
%   'beta'     -> full beta   cross-term (coherence loss on both fields).
% Only MDSN and DSN carry this cross-term; DN and D are unaffected. A cell of
% both forms fits each independently and keeps the higher-R^2 result per pixel.
cLoss = p.Results.cLoss;
if ischar(cLoss), cLoss = {cLoss}; end
if ~iscell(cLoss) || ~all(cellfun(@ischar, cLoss(:)))
    error('fitDLSI:cLoss','cLoss must be a char or a cell array of chars.');
end
cLoss = cellfun(@(c) validatestring(c,{'beta','sqrtbeta'}), cLoss(:)', 'UniformOutput', false);
cLoss = unique(cLoss,'stable');
cLossCoef = cell(size(cLoss));      % cross-term coefficient string per variant
cLossCode = zeros(size(cLoss));     % code stored in fit.cLoss: 1 = sqrtbeta, 2 = beta
for ci = 1:numel(cLoss)
    if strcmp(cLoss{ci},'sqrtbeta')
        cLossCoef{ci} = 'sqrt(beta)'; cLossCode(ci) = 1;
    else
        cLossCoef{ci} = 'beta';       cLossCode(ci) = 2;
    end
end

size2D = size(iniTau);
size1D = size(iniTau,1).*size(iniTau,2);

% Defaults that absorb the former example-side boilerplate ------------------
if isempty(order)
    [~,order] = sort(iniTau(:), 'descend');
end
mask = p.Results.mask;
if isempty(mask)
    mask = ones(size2D);
end
betaMax = p.Results.betaMax;
if isempty(betaMax)
    betaMax = ones(size2D);
elseif isscalar(betaMax)
    betaMax = betaMax*ones(size2D);
end
betaMin = p.Results.betaMin;
if isempty(betaMin)
    g20     = g2(:,:,1);
    betaMin = ones(size2D)*(prctile(g20(:),0.01)-1);   % beta >= g2_0-1
elseif isscalar(betaMin)
    betaMin = betaMin*ones(size2D);
end
% Beta scaling and a safety gap between the bounds (as in the original setup).
btMax = min(betaMax./bs, 1);
btMin = min(max(betaMin./bs, betaMin), btMax-2e-7);

% Lightweight settings record. The heavy g2 is NOT stored here: it lives once
% in the source (*_g_d.mat) data struct and must not be duplicated in the fit.
model.settings.type=type;
model.settings.pointsMin=pointsMin;
model.settings.lags=lags;
model.settings.lagStart=2;
model.settings.order=order;
model.settings.isAdaptive=isAdaptive;
model.settings.cLoss=cLoss;
model.settings.fitSet={'Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off'};
fitSet=model.settings.fitSet;
lagStart=model.settings.lagStart;

updStep=floor(size1D/20);

%points range
model.range.tailPoints=2;
model.range.idx=zeros(size(g2,1),size(g2,2),2);
model.range.g2Low=zeros(size(g2,1),size(g2,2));
model.range.lowPoints=zeros(size(g2,1),size(g2,2));
model.range.artifacts=squeeze(g2(:,:,1)>2) | squeeze(g2(:,:,2)<1);

model.range.idx(:,:,1)=lagStart;
for i=1:1:size(g2,1)
    for ii=1:1:size(g2,2)
        ts=squeeze(g2(i,ii,:));
        model.range.idx(i,ii,2)=length(ts);
        idx=find((ts(1:end-1)-ts(2:end))<=0,1,'first');
        if ~isempty(idx)
            model.range.idx(i,ii,2)=idx;
        end
        model.range.idx(i,ii,2)=min(model.range.idx(i,ii,2)+model.range.tailPoints,size(g2,3));
        model.range.pointsN=model.range.idx(i,ii,2)-model.range.idx(i,ii,1)-model.range.tailPoints;
        if (model.range.idx(i,ii,2)-model.range.idx(i,ii,1))-model.range.tailPoints<pointsMin-1
            model.range.idx(i,ii,2)=model.range.idx(i,ii,1)+pointsMin-1+model.range.tailPoints;
            model.range.lowPoints(i,ii)=1;
        end
    end
end
model.range.lowPointsTDC=round(iniTau/(lags(2)-lags(1)))<=pointsMin;

%Reshape to 1D and sort by fitting order
iniTau=iniTau(order);
mask=mask(order);

g2=reshape(g2,size1D,size(g2,3));
artifacts=reshape(model.range.artifacts,size1D,1);
lowPointsTDC=reshape(model.range.lowPointsTDC,size1D,1);
rangeIdxs=reshape(model.range.idx,size1D,2);
g2=g2(order,:);
btMin=btMin(order);
btMax=btMax(order);
lowPointsTDC=lowPointsTDC(order);
artifacts=artifacts(order);
rangeIdxs=rangeIdxs(order,:);
[~,idxRev]=sort(order);

%set initial boundaries
tauMin=0;
tauMax=lags(squeeze(rangeIdxs(:,2)));
pMin=zeros(size2D);
pMax=ones(size2D);

%temporary arrays for parfor: per-pixel trimmed g2 and lag vectors.
g2Sub=cell(size(iniTau));
lagsSub=cell(size(iniTau));
lagIdxAll=1:1:size(g2,2);                       % fixed lag-index set
for i=1:1:size(g2,1)
    if artifacts(i)==0 && mask(i)
        sel=lagIdxAll(lagIdxAll>rangeIdxs(i,1) & lagIdxAll<rangeIdxs(i,2));
        g2Sub{i}=g2(i,sel)';
        lagsSub{i}=lags(sel);
    end
end

% Fit each requested coherence-loss form INDEPENDENTLY (its own adaptive bounds
% and its own fit history) - the forms are mutually exclusive hypotheses and
% must never co-adapt - then keep the higher-R^2 form per pixel. DN and D have
% no static/dynamic cross-term, so cLoss does not apply and a single pass runs.
if any(strcmp(type,{'MDSN','DSN'}))
    nForms = numel(cLoss);
else
    nForms = 1;
end
validMask = artifacts==0 & mask~=0;
btMin0 = btMin; btMax0 = btMax; pMin0 = pMin; pMax0 = pMax;   % shared pre-adaptation start

% best-across-forms results (seeded by the first form, then improved per pixel)
fitBeta=zeros(size(iniTau)); fitTau=zeros(size(iniTau)); fitR2=zeros(size(iniTau));
fitSse=zeros(size(iniTau));  fitD=zeros(size(iniTau));   fitC=zeros(size(iniTau));
fitP=ones(size(iniTau));     fitCLoss=nan(size(iniTau));  % 1 = sqrtbeta, 2 = beta, NaN = N/A
btMinBest=btMin0; btMaxBest=btMax0; pMinBest=pMin0; pMaxBest=pMax0;
poolSize=64;

for fi=1:nForms
    cc=cLossCoef{fi}; clCode=cLossCode(fi);              % coherence-loss form for this pass
    btMin=btMin0; btMax=btMax0; pMin=pMin0; pMax=pMax0;  % reset the independent adaptive state
    betaUpdated=0;
    fBeta=zeros(size(iniTau)); fTau=zeros(size(iniTau)); fR2=zeros(size(iniTau));
    fSse=zeros(size(iniTau));  fD=zeros(size(iniTau));   fC=zeros(size(iniTau));
    fP=ones(size(iniTau));
    for ii=0:poolSize:size(g2,1)-1
        if isempty(progressFcn)
            if nForms>1
                disp([cc,': ',num2str(round(ii/size(g2,1)*100)),'%'])
            else
                disp([num2str(round(ii/size(g2,1)*100)),'%'])
            end
        else
            % One bar across all the coherence-loss passes, not one per pass.
            progressFcn((fi-1+ii/size(g2,1))/nForms, 'DLSI fit');
        end
        iEnd=min(ii+poolSize,size(g2,1));           % last batch covers the remainder
        if isAdaptive && ii>updStep
            tsR2=fR2((ii-updStep+1):ii);
            tsBeta=fBeta((ii-updStep+1):ii);
            validPixels= artifacts((ii-updStep+1):ii)==0 & mask((ii-updStep+1):ii)==1 & tsR2>median(tsR2) & lowPointsTDC((ii-updStep+1):ii)==0;
            tsBeta=tsBeta(validPixels);
            tsP=fP((ii-updStep+1):ii);
            tsP=tsP(validPixels);

            if sum(validPixels)>updStep/3
                updValMax=btMin(ii:end)+2e-7;
                updValMin=btMax(ii:end)-2e-7;
                tmpMax=median(tsBeta)+2*std(tsBeta);
                tmpMin=median(tsBeta)-6*std(tsBeta);
                updValMax(updValMax<tmpMax)=tmpMax;
                updValMin(updValMin>tmpMin)=tmpMin;

                tmpMax=btMax(ii:end);
                updValMax(updValMax>tmpMax)=tmpMax(updValMax>tmpMax);
                updValMax(updValMax>1)=1;
                btMax(ii:end)=updValMax;

                if betaUpdated==0
                    tmpMin=btMin(ii:end);
                    updValMin(updValMin<tmpMin)=tmpMin(updValMin<tmpMin);
                    updValMin(lowPointsTDC(ii:end)==1)=tmpMin(lowPointsTDC(ii:end)==1);
                    updValMin(updValMin<0)=0;
                    btMin(ii:end)=updValMin;
                    betaUpdated=1;
                end

                updVal=pMax(ii:end)-2e-7;
                tmp=median(tsP)-4*std(tsP);
                updVal(updVal>tmp)=tmp;
                tmp=pMin(ii:end);
                updVal(updVal<tmp)=tmp(updVal<tmp);
                pMin(ii:end)=updVal;
            end
        end
        if strcmp(type,'MDSN')
            parfor i=ii+1:iEnd
                if artifacts(i)==0 && mask(i)==1
                    ts=double(g2Sub{i});
                    xq=double(lagsSub{i});

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,0,pMin(i),tauMin],...
                        'Upper',[btMax(i),1,1,pMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]); %#ok<PFBNS>
                    fitType = fittype(['1+',num2str(bs),'*( beta*((p^2)*(d*d*exp(-2*((x/tau)))+(1-d)*(1-d)*exp(-2*((x/tau)^2)) + 2*d*(1-d)*exp(-((x/tau)))*exp(-((x/tau)^2))))+',cc,'*2*p*(1-p)*(d*exp(-((x/tau)))+(1-d)*exp(-((x/tau)^2))))+c'],'options',fitOpt);
                    [foSOSUMO,gofSOSUMO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,0,pMin(i),tauMin],...
                        'Upper',[btMax(i),1,1,pMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*( beta*((p^2)*(d*d*exp(-2*((x/tau)))+(1-d)*(1-d)*exp(-2*((x/tau)^0.5)) +2*d*(1-d)*exp(-((x/tau)^0.5))*exp(-((x/tau)))))+',cc,'*2*p*(1-p)*((1-d)*exp(-((x/tau)^0.5))+d*exp(-((x/tau)))))+c'],'options',fitOpt);
                    [foSUMOMU,gofSUMOMU] = fit(xq',ts,fitType);

                    if gofSOSUMO.rsquare>gofSUMOMU.rsquare
                        fBeta(i)=foSOSUMO.beta; fTau(i)=foSOSUMO.tau; fD(i)=2-foSOSUMO.d; fC(i)=foSOSUMO.c; fP(i)=foSOSUMO.p; fR2(i)=gofSOSUMO.rsquare; fSse(i)=gofSOSUMO.sse;
                    else
                        fBeta(i)=foSUMOMU.beta; fTau(i)=foSUMOMU.tau; fD(i)=foSUMOMU.d; fC(i)=foSUMOMU.c; fP(i)=foSUMOMU.p; fR2(i)=gofSUMOMU.rsquare; fSse(i)=gofSUMOMU.sse;
                    end
                end
            end

        elseif strcmp(type,'DSN')
            parfor i=ii+1:iEnd
                if artifacts(i)==0  && mask(i)==1
                    ts=double(g2Sub{i});
                    xq=double(lagsSub{i});
                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,pMin(i),tauMin],...
                        'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]); %#ok<PFBNS>
                    fitType = fittype(['1+',num2str(bs),'*(beta*(p^2)*(exp(-2*((x/tau)^2)))+',cc,'*2*p*(1-p)*(exp(-((x/tau)^2))))+c'],'options',fitOpt);
                    [foSO,gofSO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,pMin(i),tauMin],...
                        'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*(beta*(p^2)*(exp(-2*((x/tau)^1)))+',cc,'*2*p*(1-p)*(exp(-((x/tau)^1))))+c'],'options',fitOpt);
                    [foSUMO,gofSUMO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,pMin(i),tauMin],...
                        'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*(beta*(p^2)*(exp(-2*((x/tau)^0.5)))+',cc,'*2*p*(1-p)*(exp(-((x/tau)^0.5))))+c'],'options',fitOpt);
                    [foMU,gofMU] = fit(xq',ts,fitType);

                    if gofSO.rsquare>max(gofSUMO.rsquare,gofMU.rsquare)
                        fBeta(i)=foSO.beta; fTau(i)=foSO.tau; fD(i)=2; fC(i)=foSO.c; fP(i)=foSO.p; fR2(i)=gofSO.rsquare; fSse(i)=gofSO.sse;
                    elseif gofSUMO.rsquare>max(gofSO.rsquare,gofMU.rsquare)
                        fBeta(i)=foSUMO.beta; fTau(i)=foSUMO.tau; fD(i)=1; fC(i)=foSUMO.c; fP(i)=foSUMO.p; fR2(i)=gofSUMO.rsquare; fSse(i)=gofSUMO.sse;
                    else
                        fBeta(i)=foMU.beta; fTau(i)=foMU.tau; fD(i)=0; fC(i)=foMU.c; fP(i)=foMU.p; fR2(i)=gofMU.rsquare; fSse(i)=gofMU.sse;
                    end
                end
            end
        elseif strcmp(type,'DN')
            parfor i=ii+1:iEnd
                if artifacts(i)==0  && mask(i)==1
                    ts=double(g2Sub{i});
                    xq=double(lagsSub{i});

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,tauMin],...
                        'Upper',[btMax(i),1,tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,iniTau(i)]); %#ok<PFBNS>
                    fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)^2))))+c'],'options',fitOpt);
                    [foSO,gofSO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,tauMin],...
                        'Upper',[btMax(i),1,tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)))))+c'],'options',fitOpt);
                    [foSUMO,gofSUMO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),-1,tauMin],...
                        'Upper',[btMax(i),1,tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)^0.5))))+c'],'options',fitOpt);
                    [foMU,gofMU] = fit(xq',ts,fitType);

                    if gofSO.rsquare>max(gofSUMO.rsquare,gofMU.rsquare)
                        fBeta(i)=foSO.beta; fTau(i)=foSO.tau; fC(i)=foSO.c; fP(i)=1; fR2(i)=gofSO.rsquare; fSse(i)=gofSO.sse; fD(i)=2;
                    elseif gofSUMO.rsquare>max(gofSO.rsquare,gofMU.rsquare)
                        fBeta(i)=foSUMO.beta; fTau(i)=foSUMO.tau; fC(i)=foSUMO.c; fP(i)=1; fR2(i)=gofSUMO.rsquare; fSse(i)=gofSUMO.sse; fD(i)=1;
                    else
                        fBeta(i)=foMU.beta; fTau(i)=foMU.tau; fC(i)=foMU.c; fP(i)=1; fR2(i)=gofMU.rsquare; fSse(i)=gofMU.sse; fD(i)=0;
                    end
                end
            end
        elseif strcmp(type,'D')
            parfor i=ii+1:iEnd
                if artifacts(i)==0  && mask(i)==1
                    ts=double(g2Sub{i});
                    xq=double(lagsSub{i});

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),tauMin],...
                        'Upper',[btMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]); %#ok<PFBNS>
                    fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)^2))))'],'options',fitOpt);
                    [foSO,gofSO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),tauMin],...
                        'Upper',[btMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)))))'],'options',fitOpt);
                    [foSUMO,gofSUMO] = fit(xq',ts,fitType);

                    fitOpt = fitoptions(fitSet{:},...
                        'Lower',[btMin(i),tauMin],...
                        'Upper',[btMax(i),tauMax(i)],...
                        'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]);
                    fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)^0.5))))'],'options',fitOpt);
                    [foMU,gofMU] = fit(xq',ts,fitType);

                    if gofSO.rsquare>max(gofSUMO.rsquare,gofMU.rsquare)
                        fBeta(i)=foSO.beta; fTau(i)=foSO.tau; fC(i)=0; fP(i)=1; fR2(i)=gofSO.rsquare; fSse(i)=gofSO.sse; fD(i)=2;
                    elseif gofSUMO.rsquare>max(gofSO.rsquare,gofMU.rsquare)
                        fBeta(i)=foSUMO.beta; fTau(i)=foSUMO.tau; fC(i)=0; fP(i)=1; fR2(i)=gofSUMO.rsquare; fSse(i)=gofSUMO.sse; fD(i)=1;
                    else
                        fBeta(i)=foMU.beta; fTau(i)=foMU.tau; fC(i)=0; fP(i)=1; fR2(i)=gofMU.rsquare; fSse(i)=gofMU.sse; fD(i)=0;
                    end
                end
            end
        end
    end
    % Keep the higher-R^2 form per pixel; the first form seeds the best arrays.
    if fi==1
        fitBeta=fBeta; fitTau=fTau; fitD=fD; fitC=fC; fitP=fP; fitR2=fR2; fitSse=fSse;
        btMinBest=btMin; btMaxBest=btMax; pMinBest=pMin; pMaxBest=pMax;
        if any(strcmp(type,{'MDSN','DSN'}))
            fitCLoss(validMask)=clCode;
        end
    else
        adopt = fR2>fitR2;
        fitBeta(adopt)=fBeta(adopt); fitTau(adopt)=fTau(adopt); fitD(adopt)=fD(adopt);
        fitC(adopt)=fC(adopt); fitP(adopt)=fP(adopt); fitR2(adopt)=fR2(adopt); fitSse(adopt)=fSse(adopt);
        fitCLoss(adopt)=clCode;
        btMinBest(adopt)=btMin(adopt); btMaxBest(adopt)=btMax(adopt);
        pMinBest(adopt)=pMin(adopt); pMaxBest(adopt)=pMax(adopt);
    end
end
%get data from temporary arrays and reshape them to 2d
model.range.tauMax=reshape(tauMax(idxRev),size2D);
model.range.tauMin=tauMin;
model.range.btMax=reshape(btMaxBest(idxRev),size2D);
model.range.btMin=reshape(btMinBest(idxRev),size2D);
model.range.pMax=reshape(pMaxBest(idxRev),size2D);
model.range.pMin=reshape(pMinBest(idxRev),size2D);
model.fit.c=reshape(fitC(idxRev),size2D);
model.fit.beta=reshape(fitBeta(idxRev),size2D);
model.fit.p=reshape(fitP(idxRev),size2D);
model.fit.d=reshape(fitD(idxRev),size2D);
model.fit.tau=reshape(fitTau(idxRev),size2D);
model.fit.r2=reshape(fitR2(idxRev),size2D);
model.fit.sse=reshape(fitSse(idxRev),size2D);
model.fit.cLoss=reshape(fitCLoss(idxRev),size2D);   % 1 = sqrtbeta, 2 = beta (NaN if N/A)
end
