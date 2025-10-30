% fitNormalizedG2 - fits g2 to DLSI models allowing mixed dynamics
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

%------------- BEGIN CODE --------------

function model=fitDLSI_basic(g2,lags,pointsMin,order,type,iniTau,btMin,btMax,bs,mask,isAdaptive)

size2D=size(iniTau);
size1D=size(iniTau,1).*size(iniTau,2);

model.input.g2=g2;
model.input.iniTau=iniTau;
model.input.type=type;
model.input.pointsMin=pointsMin;
model.input.lags=lags;
model.input.lagStart=2;
model.input.order=order;
model.input.isAdaptive=isAdaptive;
model.input.fitSet={'Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off'}; %'Robust','Bisquare',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
fitSet=model.input.fitSet;
lagStart=model.input.lagStart;

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

%Reshape to 1D and sort
%iniTau=medfilt2(iniTau,[7,7]).*mask;
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
%pMin=(g2(:,1)-1)./btMax;
%pMax=(g2(:,1)-1)./btMin;
pMin=zeros(size2D);
pMax=ones(size2D);



%temporary arrays for parfor
g2Sub=cell(size(iniTau));
lagsSub=cell(size(iniTau));
%idxs=unique(round(logspace(0,log10(size(g2,2)),100)));
idxs=1:1:size(g2,2);
for i=1:1:size(g2,1)
    if artifacts(i)==0 && mask(i)
        idxs=idxs(idxs>rangeIdxs(i,1)& idxs<rangeIdxs(i,2));
        g2Sub{i}=squeeze(g2(i,idxs))';
        lagsSub{i}=lags(idxs);
    end
end

fitBeta=zeros(size(iniTau));
fitTau=zeros(size(iniTau));
fitR2=zeros(size(iniTau));
fitSse=zeros(size(iniTau));
fitD=zeros(size(iniTau));
fitC=zeros(size(iniTau));
fitP=ones(size(iniTau));
figure
poolSize=64;
betaUpdated=0;
for ii=0:poolSize:size(g2,1)-poolSize
    disp([num2str(ii/size(g2,1)*100),'%'])
    tic
    ii
    if isAdaptive && ii>updStep
        tsR2=fitR2((ii-updStep+1):ii);
        tsBeta=fitBeta((ii-updStep+1):ii);
        validPixels= artifacts((ii-updStep+1):ii)==0 & mask((ii-updStep+1):ii)==1 & tsR2>median(tsR2) & lowPointsTDC((ii-updStep+1):ii)==0;
        tsBeta=tsBeta(validPixels);
        tsP=fitP((ii-updStep+1):ii);
        tsP=tsP(validPixels);

        if sum(validPixels)>updStep/3
                        
                            hold on
                            updValMax=btMin(ii:end)+2e-7;
                            updValMin=btMax(ii:end)-2e-7;
                            tmpMax=median(tsBeta)+2*std(tsBeta);
                            tmpMin=median(tsBeta)-6*std(tsBeta);
                            errorbar(ii,median(tsBeta),std(tsBeta),'r')
                            scatter(ii,tmpMax,'>r')
                            scatter(ii,tmpMin,'<r')
                            updValMax(updValMax<tmpMax)=tmpMax;
                            updValMin(updValMin>tmpMin)=tmpMin;
            
                            tmpMax=btMax(ii:end);
                            updValMax(updValMax>tmpMax)=tmpMax(updValMax>tmpMax);
                            updValMax(updValMax>1)=1;
                            btMax(ii:end)=updValMax;
                            scatter(ii,max(updValMax(lowPointsTDC(ii:end)==0)),'or')

                            if betaUpdated==0
                                tmpMin=btMin(ii:end);
                                updValMin(updValMin<tmpMin)=tmpMin(updValMin<tmpMin);
                                updValMin(lowPointsTDC(ii:end)==1)=tmpMin(lowPointsTDC(ii:end)==1);
                                updValMin(updValMin<0)=0;
                                btMin(ii:end)=updValMin;
                                scatter(ii,min(updValMin(lowPointsTDC(ii:end)==0)),'or')
                                betaUpdated=1;
                            end
            
            
                        hold on
                        updVal=pMax(ii:end)-2e-7;
                        tmp=median(tsP)-4*std(tsP);
                        scatter(ii,tmp,'<b')
                        errorbar(ii,median(tsP),std(tsP),'b')
                        updVal(updVal>tmp)=tmp;
                        tmp=pMin(ii:end);
                        updVal(updVal<tmp)=tmp(updVal<tmp);
                        pMin(ii:end)=updVal;
                        scatter(ii,min(updVal),'ob')
                        hold off
                        pause(0.1)

        end
    end
    if strcmp(type,'MDSN')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0 && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});

                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),-1,0,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*( beta*((p^2)*(d*d*exp(-2*((x/tau)))+(1-d)*(1-d)*exp(-2*((x/tau)^2)) + 2*d*(1-d)*exp(-((x/tau)))*exp(-((x/tau)^2))))+sqrt(beta)*2*p*(1-p)*(d*exp(-((x/tau)))+(1-d)*exp(-((x/tau)^2))))+c'],'options',fitOpt);
                [foSOSUMO,gofSOSUMO] = fit(xq',ts,fitType);


                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),-1,0,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*( beta*((p^2)*(d*d*exp(-2*((x/tau)))+(1-d)*(1-d)*exp(-2*((x/tau)^0.5)) +2*d*(1-d)*exp(-((x/tau)^0.5))*exp(-((x/tau)))))+sqrt(beta)*2*p*(1-p)*((1-d)*exp(-((x/tau)^0.5))+d*exp(-((x/tau)))))+c'],'options',fitOpt);
                [foSUMOMU,gofSUMOMU] = fit(xq',ts,fitType);

                if gofSOSUMO.rsquare>gofSUMOMU.rsquare
                    fitBeta(i)=foSOSUMO.beta;
                    fitTau(i)=foSOSUMO.tau;
                    fitD(i)=2-foSOSUMO.d;
                    fitC(i)=foSOSUMO.c;
                    fitP(i)=foSOSUMO.p;
                    fitR2(i)=gofSOSUMO.rsquare;
                    fitSse(i)=gofSOSUMO.sse;

                else
                    fitBeta(i)=foSUMOMU.beta;
                    fitTau(i)=foSUMOMU.tau;
                    fitD(i)=foSUMOMU.d;
                    fitC(i)=foSUMOMU.c;
                    fitP(i)=foSUMOMU.p;
                    fitR2(i)=gofSUMOMU.rsquare;
                    fitSse(i)=gofSUMOMU.sse;

                end
            end
        end
        toc

    elseif strcmp(type,'DSN')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),-1,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*(beta*(p^2)*(exp(-2*((x/tau)^2)))+(beta^0.5)*2*p*(1-p)*(exp(-((x/tau)^2))))+c'],'options',fitOpt);
                [foSO,gofSO] = fit(xq',ts,fitType);

                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),-1,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*(beta*(p^2)*(exp(-2*((x/tau)^1)))+(beta^0.5)*2*p*(1-p)*(exp(-((x/tau)^1))))+c'],'options',fitOpt);
                [foSUMO,gofSUMO] = fit(xq',ts,fitType);

                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),-1,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*(beta*(p^2)*(exp(-2*((x/tau)^0.5)))+(beta^0.5)*2*p*(1-p)*(exp(-((x/tau)^0.5))))+c'],'options',fitOpt);
                [foMU,gofMU] = fit(xq',ts,fitType);

                if gofSO.rsquare>max(gofSUMO.rsquare,gofMU.rsquare)
                    fitBeta(i)=foSO.beta;
                    fitTau(i)=foSO.tau;
                    fitD(i)=2;
                    fitC(i)=foSO.c;
                    fitP(i)=foSO.p;
                    fitR2(i)=gofSO.rsquare;
                    fitSse(i)=gofSO.sse;
                elseif gofSUMO.rsquare>max(gofSO.rsquare,gofMU.rsquare)
                    fitBeta(i)=foSUMO.beta;
                    fitTau(i)=foSUMO.tau;
                    fitC(i)=foSUMO.c;
                    fitP(i)=foSUMO.p;
                    fitR2(i)=gofSUMO.rsquare;
                    fitSse(i)=gofSUMO.sse;
                    fitD(i)=1;
                else %MU
                    fitBeta(i)=foMU.beta;
                    fitTau(i)=foMU.tau;
                    fitC(i)=foMU.c;
                    fitP(i)=foMU.p;
                    fitR2(i)=gofMU.rsquare;
                    fitSse(i)=gofMU.sse;
                    fitD(i)=0;
                end
            end
        end
        toc
    elseif strcmp(type,'DN')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});

                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),-1,tauMin],...
                    'Upper',[btMax(i),1,tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,iniTau(i)]);
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
                    fitBeta(i)=foSO.beta;
                    fitTau(i)=foSO.tau;
                    fitC(i)=foSO.c;
                    fitP(i)=1;
                    fitR2(i)=gofSO.rsquare;
                    fitSse(i)=gofSO.sse;
                    fitD(i)=2;
                elseif gofSUMO.rsquare>max(gofSO.rsquare,gofMU.rsquare)
                    fitBeta(i)=foSUMO.beta;
                    fitTau(i)=foSUMO.tau;
                    fitC(i)=foSUMO.c;
                    fitP(i)=1;
                    fitR2(i)=gofSUMO.rsquare;
                    fitSse(i)=gofSUMO.sse;
                    fitD(i)=1;
                else %MU
                    fitBeta(i)=foMU.beta;
                    fitTau(i)=foMU.tau;
                    fitC(i)=foMU.c;
                    fitP(i)=1;
                    fitR2(i)=gofMU.rsquare;
                    fitSse(i)=gofMU.sse;
                    fitD(i)=0;
                end
            end
        end
        toc
    elseif strcmp(type,'D')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});

                fitOpt = fitoptions(fitSet{:},...
                    'Lower',[btMin(i),tauMin],...
                    'Upper',[btMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]);
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
                    fitBeta(i)=foSO.beta;
                    fitTau(i)=foSO.tau;
                    fitC(i)=0;
                    fitP(i)=1;
                    fitR2(i)=gofSO.rsquare;
                    fitSse(i)=gofSO.sse;
                    fitD(i)=2;
                elseif gofSUMO.rsquare>max(gofSO.rsquare,gofMU.rsquare)
                    fitBeta(i)=foSUMO.beta;
                    fitTau(i)=foSUMO.tau;
                    fitC(i)=0;
                    fitP(i)=1;
                    fitR2(i)=gofSUMO.rsquare;
                    fitSse(i)=gofSUMO.sse;
                    fitD(i)=1;
                else %MU
                    fitBeta(i)=foMU.beta;
                    fitTau(i)=foMU.tau;
                    fitC(i)=0;
                    fitP(i)=1;
                    fitR2(i)=gofMU.rsquare;
                    fitSse(i)=gofMU.sse;
                    fitD(i)=0;
                end
            end
        end
        toc

    end
end
%get data from temporary arrays and reshape them to 2d
model.range.tauMax=reshape(tauMax(idxRev),size2D);
model.range.tauMin=tauMin;
model.range.btMax=reshape(btMax(idxRev),size2D);
model.range.btMin=reshape(btMin(idxRev),size2D);
model.range.pMax=reshape(pMax(idxRev),size2D);
model.range.pMin=reshape(pMin(idxRev),size2D);
model.fit.c=reshape(fitC(idxRev),size2D);
model.fit.beta=reshape(fitBeta(idxRev),size2D);
model.fit.p=reshape(fitP(idxRev),size2D);
model.fit.d=reshape(fitD(idxRev),size2D);
model.fit.tau=reshape(fitTau(idxRev),size2D);
model.fit.r2=reshape(fitR2(idxRev),size2D);
model.fit.sse=reshape(fitSse(idxRev),size2D);
end