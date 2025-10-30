% fitNormalizedG2 - fits g2 to DLSI models allowing mixed dynamics
%
%
% Authors: DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 2018

%------------- BEGIN CODE --------------

function model=fitNormalizedG2(g2,fps,iniTau,pointsMin,type,betaMin,betaMax,mask)

%Model types are
% SOSUMO - single-ordered to single-unordered/multiple-ordered mixed model
% SUMOMU - single-unordered/multiple-ordered to multiple-unordered mixed model

size2D=size(iniTau);
size1D=size(iniTau,1).*size(iniTau,2);

model.input.g2=g2;
model.input.iniTau=iniTau;
model.input.type=type;
model.input.pointsMin=pointsMin;
model.input.fps=fps;


%Reshape to 1D and sort
iniTau=medfilt2(iniTau,[7,7]).*mask;
[~,idxSort]=sort(iniTau(:),'descend');
iniTau=iniTau(idxSort(:));
mask=mask(idxSort(:));
g2=reshape(g2,size1D,size(g2,3));
for i=1:1:size(g2,2)
    ts=squeeze(g2(:,i));
    ts=ts(idxSort);
    g2(:,i)=ts;
end
parAdaptStep=floor(sum(mask(:))/4);


betaMin=betaMin(idxSort(:));
betaMax=betaMax(idxSort(:));

[~,idxRev]=sort(idxSort(:));
%points range
model.range.g2Coef=0.1;
model.range.idx=zeros(size(g2,1),size(g2,2),2);
model.range.g2Low=zeros(size(g2,1),size(g2,2));
model.range.lowPoints=zeros(size(g2,1),size(g2,2));

rangeArtifacts=squeeze(g2(:,2)>2) | squeeze(g2(:,2)<1);

startLagIdx=2;
figure
for i=1:1:size(g2,1)
    if rangeArtifacts(i)==0 && mask(i)
    ts=squeeze(g2(i,:));
    model.range.g2Low(i)=mean(ts(end-25:end));
    model.range.idx(i,1)=startLagIdx;
    
    val=find(ts<=(model.range.g2Low(i)+(ts(1)-model.range.g2Low(i)).*model.range.g2Coef),1);
    if isempty(val)
    rangeArtifacts(i)=1;   
    mask(i)=0;
    else
    model.range.idx(i,2)=val;
    
    if (model.range.idx(i,2)-model.range.idx(i,1))<pointsMin-1
        model.range.idx(i,2)=model.range.idx(i,1)+pointsMin-1;
        model.range.lowPoints(i)=1;
    end
    end
    end
end

%temporary arrays for parfor
g2Sub=cell(size(iniTau));
for i=1:1:size(g2,1)
     if rangeArtifacts(i)==0 && mask(i)
    g2Sub{i}=squeeze(g2(i,model.range.idx(i,1):model.range.idx(i,2)))';
     end
end

tauMin=0+(1./fps)./5;
tauMax=(size(g2,2)-1)./fps;

pMin=zeros(size(iniTau));
pMax=ones(size(iniTau));

fitBeta=zeros(size(iniTau));
fitTau=zeros(size(iniTau));
fitR2=zeros(size(iniTau));
fitSse=zeros(size(iniTau));
fitM=zeros(size(iniTau));
fitD=zeros(size(iniTau));
fitC=zeros(size(iniTau));
fitP=ones(size(iniTau));

poolSize=64;
figure
if strcmp(type,'SOSUMO')
    for ii=0:poolSize:size(g2,1)-poolSize
        disp([num2str(ii/size(g2,1)*100),'%'])
        tic
        if ii>parAdaptStep
            validPixels= rangeArtifacts(ii-parAdaptStep+1:ii)==0 & mask(ii-parAdaptStep+1:ii);
            tsBeta=fitBeta(ii-parAdaptStep+1:ii);
            
            tsP=fitP(ii-parAdaptStep+1:ii);
            tsP=tsP(tsBeta>=betaMin(ii-1));
            validPixels=validPixels(tsBeta>=betaMin(ii-1));
            tsBeta=tsBeta(tsBeta>=betaMin(ii-1));
            
            tsBeta=tsBeta(validPixels);
            tsP=tsP(validPixels);
            
            %make it so max adapted only based on mean, but then std is
            %added when fitted
            
            if sum(validPixels)>100
                betaMax(ii:end)=max(min(mean(tsBeta)+2*std(tsBeta),max(tsBeta)),betaMin(ii:end)+2e-7);
                pMin(ii:end)=min(max(mean(tsP)-2*std(tsP),min(tsP)),pMax(ii-1)-2e-7);
            end
        end
        parfor i=ii+1:1:ii+poolSize
            if rangeArtifacts(i)==0 && mask(i)==1
                fitOpt = fitoptions('Method','NonlinearLeastSquares',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[betaMin(i),-1,0,pMin(i),tauMin],...
                    'Upper',[betaMax(i),1,1,pMax(i),tauMax],...
                    'StartPoint',[betaMax(i)-0.9*(betaMax(i)-betaMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype('1+beta*( (p^2)*(m*m*exp(-2*((x/tau)))+(1-m)*(1-m)*exp(-2*((x/tau)^2)) + 2*m*(1-m)*exp(-((x/tau)))*exp(-((x/tau)^2)) )) + (beta^0.5)*2*p*(1-p)*(m*exp(-((x/tau)))+(1-m)*exp(-((x/tau)^2)))+c','options',fitOpt);
                ts=double(g2Sub{i});
                xq=(startLagIdx:1:(length(ts)+1))-1;
                xq=xq./fps;
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitM(i)=fitobject.m;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                %plot(fitobject,xq,ts);
                %title(['b=',num2str(fitBeta(i)),' c=',num2str(fitC(i)),' p=',num2str(fitP(i)),' bMax=',num2str(betaMax(i)),' pMin=',num2str(pMin(i))])
            end
        end
        toc
    end
elseif strcmp(type,'SUMOMU')
    for ii=0:poolSize:size(g2,1)-poolSize
        disp([num2str(ii/size(g2,1)*100),'%'])
        tic
        if ii>parAdaptStep
             validPixels= rangeArtifacts(ii-parAdaptStep+1:ii)==0 & mask(ii-parAdaptStep+1:ii);
            tsBeta=fitBeta(ii-parAdaptStep+1:ii);
            
            tsP=fitP(ii-parAdaptStep+1:ii);
            tsP=tsP(tsBeta>=betaMin(ii-1));
            validPixels=validPixels(tsBeta>=betaMin(ii-1));
            tsBeta=tsBeta(tsBeta>=betaMin(ii-1));
            
            
            tsBeta=tsBeta(validPixels);
            tsP=tsP(validPixels);
            
            
            
            if sum(validPixels)>100
                betaMax(ii:end)=max(min(mean(tsBeta)+2*std(tsBeta),max(tsBeta)),betaMin(ii:end)+2e-7);
                pMin(ii:end)=min(max(mean(tsP)-2*std(tsP),min(tsP)),pMax(ii-1)-2e-7);
            end
            % mean(fitBeta(ii-poolSize+1:ii))
        end
        parfor i=ii+1:1:ii+poolSize
            if rangeArtifacts(i)==0  && mask(i)==1
                fitOpt = fitoptions('Method','NonlinearLeastSquares',...%'MaxFunEvals',3000,'DiffMinChange',10^(-12),'MaxIter',1500,'TolFun',10^(-12),'TolX',10^(-12),...
                    'Lower',[betaMin(i),-1,0,pMin(i),tauMin],...
                    'Upper',[betaMax(i),1,1,pMax(i),tauMax],...
                    'StartPoint',[betaMax(i)-0.9*(betaMax(i)-betaMin(i)),0,0.1,pMax(i)-(pMax(i)-pMin(i))/10,iniTau(i)]);
                fitType = fittype('1+beta*( (p^2)*(d*d*exp(-2*((x/tau)^0.5))+(1-d)*(1-d)*exp(-2*((x/tau))) +2*d*(1-d)*exp(-((x/tau)^0.5))*exp(-((x/tau))))) + (beta^0.5)*2*p*(1-p)*(d*exp(-((x/tau)^0.5))+(1-d)*exp(-((x/tau))))+c','options',fitOpt);
                ts=double(g2Sub{i});
                xq=(startLagIdx:1:(length(ts)+1))-1;
                xq=xq./fps;
                [fitobject,gof] = fit(xq',ts,fitType);
                fitBeta(i)=fitobject.beta;
                fitTau(i)=fitobject.tau;
                fitD(i)=fitobject.d;
                fitC(i)=fitobject.c;
                fitP(i)=fitobject.p;
                fitR2(i)=gof.rsquare;
                fitSse(i)=gof.sse;
                %plot(fitobject,xq,ts);
                %title(['b=',num2str(fitBeta(i)),' c=',num2str(fitC(i)),' p=',num2str(fitP(i)),' d=',num2str(fitD(i)),' m=',num2str(fitM(i)),' d> ',num2str(fitD(i)<10e-10),' m>',num2str(fitM(i)<10e-10)])
            end
        end
        toc
    end
end

%get data from temporary arrays and reshape them to 2d
model.range.tauMax=tauMax;
model.range.tauMin=tauMin;
model.range.betaMax=reshape(betaMax(idxRev),size2D);
model.range.betaMin=reshape(betaMin(idxRev),size2D);
model.range.pMax=reshape(pMax(idxRev),size2D);
model.range.pMin=reshape(pMin(idxRev),size2D);
model.range.artifacts=reshape(rangeArtifacts(idxRev),size2D);
model.fit.m=reshape(fitM(idxRev),size2D);
model.fit.c=reshape(fitC(idxRev),size2D);
model.fit.beta=reshape(fitBeta(idxRev),size2D);
model.fit.p=reshape(fitP(idxRev),size2D);
model.fit.d=reshape(fitD(idxRev),size2D);
model.fit.tau=reshape(fitTau(idxRev),size2D);
model.fit.r2=reshape(fitR2(idxRev),size2D);
model.fit.sse=reshape(fitSse(idxRev),size2D);
end