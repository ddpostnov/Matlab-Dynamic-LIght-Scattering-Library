% fitNormalizedG2 - fits g2 to DLSI models allowing mixed dynamics
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

%------------- BEGIN CODE --------------

function model=fitDLSI_old(g2,lags,pointsMin,order,type,iniTau,btMin,btMax,bs,mask)

%Model types are
% SOSUMO - single-ordered to single-unordered/multiple-ordered mixed model
% SUMOMU - single-unordered/multiple-ordered to multiple-unordered mixed model

size2D=size(iniTau);
size1D=size(iniTau,1).*size(iniTau,2);

model.input.g2=g2;
model.input.iniTau=iniTau;
model.input.type=type;
model.input.pointsMin=pointsMin;
model.input.lags=lags;
model.input.lagStart=1;
model.input.order=order;
lagStart=model.input.lagStart;

updStep=floor(size1D/100);

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
        
        model.range.idx(i,ii,2)=find(ts<=max(ts(end),double(1+ts(2)/100)),1,'first');
        model.range.idx(i,ii,2)=min(model.range.idx(i,ii,2)+model.range.tailPoints,size(g2,3));
        if (model.range.idx(i,ii,2)-model.range.idx(i,ii,1))-model.range.tailPoints<pointsMin-1
            model.range.idx(i,ii,2)=model.range.idx(i,ii,1)+pointsMin-1+model.range.tailPoints;
            model.range.lowPoints(i,ii)=1;
        end
        %if  model.range.idx(i,ii,2)-model.range.idx(i,ii,1)>model.input.pointsMin*2
        %    model.range.idx(i,ii,1)=model.range.idx(i,ii,1)+model.input.pointsMin;
        %end
    end
end

%Reshape to 1D and sort
%iniTau=medfilt2(iniTau,[7,7]).*mask;
iniTau=iniTau(order);
mask=mask(order);

g2=reshape(g2,size1D,size(g2,3));
artifacts=reshape(model.range.artifacts,size1D,1);
rangeIdxs=reshape(model.range.idx,size1D,2);
g2=g2(order,:);
btMin=btMin(order);
btMax=btMax(order);
artifacts=artifacts(order);
rangeIdxs=rangeIdxs(order,:);
[~,idxRev]=sort(order);

%set fixed boundaries
tauMin=0;
tauMax=lags(squeeze(rangeIdxs(:,2)));
%pMin=(g2(:,1)-1)./btMax;
%pMax=(g2(:,1)-1)./btMin;
pMin=zeros(size2D);
pMax=ones(size2D);



%temporary arrays for parfor
g2Sub=cell(size(iniTau));
lagsSub=cell(size(iniTau));
idxs=unique(round(logspace(0,log10(size(g2,2)),100)));
for i=1:1:size(g2,1)
    if artifacts(i)==0 && mask(i)
        idxs(idxs>rangeIdxs(i,1)& idxs<rangeIdxs(i,2));
        g2Sub{i}=squeeze(g2(i,idxs(idxs>rangeIdxs(i,1)& idxs<rangeIdxs(i,2))))';
        lagsSub{i}=lags(idxs(idxs>rangeIdxs(i,1)& idxs<rangeIdxs(i,2)));
    end
end

fitBeta=zeros(size(iniTau));
fitTau=zeros(size(iniTau));
fitR2=zeros(size(iniTau));
fitSse=zeros(size(iniTau));
fitM=zeros(size(iniTau));
fitD=zeros(size(iniTau));
fitC=zeros(size(iniTau));
fitP=ones(size(iniTau));

poolSize=64;
for ii=0:poolSize:size(g2,1)-poolSize
    disp([num2str(ii/size(g2,1)*100),'%'])
    tic
    if ii>updStep
        tsBeta=fitBeta((ii-updStep+1):ii);
        validPixels= artifacts((ii-updStep+1):ii)==0 & mask((ii-updStep+1):ii);
        tsBeta=tsBeta(validPixels);

        tsP=fitP((ii-updStep+1):ii);        
        tsP=tsP(validPixels);

        if sum(validPixels)>100
            updVal=btMin(ii:end)+2e-7;
            tmp=median(tsBeta)+2*std(tsBeta);
            updVal(updVal<tmp)=tmp;
            tmp=btMax(ii:end);
            updVal(updVal>tmp)=tmp(updVal>tmp);
            btMax(ii:end)=updVal;


            updVal=pMin(ii:end)+2e-7;
            tmp=median(tsP)-2*std(tsP);
            updVal(updVal>tmp)=tmp;
            tmp=pMin(ii:end);
            updVal(updVal<tmp)=tmp(updVal<tmp);
            pMin(ii:end)=updVal;
        end
    end
    if strcmp(type,'SOSUMO')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0 && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off',...%'Robust','Bisquare',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),-1,0,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*( beta*((p^2)*(m*m*exp(-2*((x/tau)))+(1-m)*(1-m)*exp(-2*((x/tau)^2)) + 2*m*(1-m)*exp(-((x/tau)))*exp(-((x/tau)^2))))+sqrt(beta)*2*p*(1-p)*(m*exp(-((x/tau)))+(1-m)*exp(-((x/tau)^2))))+c'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitM(i)=fitobject.m;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                
            end
        end
        toc
    elseif strcmp(type,'SUMOMU')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off',...%'Robust','Bisquare',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),-1,0,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*( beta*((p^2)*(d*d*exp(-2*((x/tau)^0.5))+(1-d)*(1-d)*exp(-2*((x/tau))) +2*d*(1-d)*exp(-((x/tau)^0.5))*exp(-((x/tau)))))+sqrt(beta)*2*p*(1-p)*(d*exp(-((x/tau)^0.5))+(1-d)*exp(-((x/tau)))))+c'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitD(i)=fitobject.d;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                
            end
        end
        toc

    elseif strcmp(type,'DCS n=2')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
               
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off',...%'Robust','Bisquare',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),-1,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*(beta*( (p^2)*(exp(-2*((x/tau)^2)))) + (beta^0.5)*2*p*(1-p)*(exp(-((x/tau)^2))))+c'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitM(i)=0;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                
            end
        end
        toc

    elseif strcmp(type,'DCS n=1')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off',...%'Robust','Bisquare',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),-1,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*(beta*( (p^2)*(exp(-2*((x/tau))))) + (beta^0.5)*2*p*(1-p)*(exp(-((x/tau)))))+c'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitM(i)=1;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                
            end
        end
        toc

    elseif strcmp(type,'DCS n=0.5')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off',...%'Robust','Bisquare',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),-1,pMin(i),tauMin],...
                    'Upper',[btMax(i),1,pMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),0,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*(beta*( (p^2)*(exp(-2*((x/tau)^0.5)))) + (beta^0.5)*2*p*(1-p)*(exp(-((x/tau)^0.5))))+c'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitD(i)=1;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;               
            end
        end
        toc
    elseif strcmp(type,'DCO n=2')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off','MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),tauMin],...
                    'Upper',[btMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)^2))))'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitM(i)=0;
                fitC(i)=0;
                fitP(i)=1;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                
            end
        end
        toc

    elseif strcmp(type,'DCO n=1')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off','MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),tauMin],...
                    'Upper',[btMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)))))'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitM(i)=1;
                fitC(i)=0;
                fitP(i)=1;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                
            end
        end
        toc

    elseif strcmp(type,'DCO n=0.5')
        parfor i=ii+1:ii+poolSize
            if artifacts(i)==0  && mask(i)==1
                ts=double(g2Sub{i});
                xq=double(lagsSub{i});
                
                fitOpt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off','MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[btMin(i),tauMin],...
                    'Upper',[btMax(i),tauMax(i)],...
                    'StartPoint',[btMax(i)-0.9*(btMax(i)-btMin(i)),iniTau(i)]);
                fitType = fittype(['1+',num2str(bs),'*beta*( (exp(-2*((x/tau)^0.5))))'],'options',fitOpt);
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitD(i)=1;
                fitC(i)=0;
                fitP(i)=1;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;               
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
model.fit.m=reshape(fitM(idxRev),size2D);
model.fit.c=reshape(fitC(idxRev),size2D);
model.fit.beta=reshape(fitBeta(idxRev),size2D);
model.fit.p=reshape(fitP(idxRev),size2D);
model.fit.d=reshape(fitD(idxRev),size2D);
model.fit.tau=reshape(fitTau(idxRev),size2D);
model.fit.r2=reshape(fitR2(idxRev),size2D);
model.fit.sse=reshape(fitSse(idxRev),size2D);
end