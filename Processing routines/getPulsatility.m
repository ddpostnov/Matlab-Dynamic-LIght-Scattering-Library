%getPulsatility  Derive pulsatility indices and optional harmonic fits
%
%   getPulsatility(s,fNames) loads every *_BFI_d.mat file in fNames (the
%   “cycle” datasets from getExternalCycle / getInternalCycle), calculates
%   classical pulsatility metrics for each spatial ROI and for every pixel,
%   and—if s.fitData is "regions" or "all"—fits the user-supplied harmonic
%   model (s.fitEquation) to the BFI waveform:
%
%        y(x) = Σ aₙ·sin(2π·n·x + bₙ)  +  c      ,  n = 1‥5
%
%   Fit coefficients and R² are appended to sMetrics, dvsMetrics, and
%   (when "all") stored as 3-D maps in RESULTS.fCoeffs.  Updated MAT-files
%   overwrite the originals; SOURCE is re-saved only when pixel-wise
%   fitting is requested.
%
%   INPUTS
%     s        parameter struct with fields  
%                • fitData      "none" | "regions" | "all"  
%                • fitEquation  fittype object (see example)  
%                • fitOptions   fitoptions handle  
%                • coefNames    cell/array of coefficient labels  
%     fNames   cell array of *_BFI_d.mat paths.
%
%   OUTPUT FILES (side-effects)
%       *_BFI_d.mat   SOURCE   – BFI cube (unchanged or refitted)  
%       *_BFI_r.mat   RESULTS  – extended ROI metrics + imgPI/imgRI layers  
%       *_BFI_s.mat   SETTINGS – field settings.pulsatilityCalculation added
%
%   EXAMPLE
%     p = defaultPulsatilityParams();        % fills fit* fields
%     files = dir(fullfile(dataRoot,'*_BFI_d.mat'));
%     getPulsatility(p, fullfile({files.folder}',{files.name}'));
%
%   DEPENDS ON
%     MATLAB Curve Fitting Toolbox for ‘fit’; core LSCI library utilities.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: Header generated with ChatGPT; minor inconsistencies may remain—
%   please verify before distribution.
%   ----------------------------------------------------------------------

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Waveform fitting
% s.fitData="all"; % "regions" or "all". Any other option e.g. "none" will not perform data fitting
% %ADJUSTED IF NECESSARY - Waveform fitting
% fitSettings={'Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off'};%,'Robust','Bisquare','MaxFunEvals',1200,'MaxIter',1000};
% fitLimits={'Lower',[0,0,0,0,0,-pi,-pi,-pi,-pi,-pi,0],...
%     'Upper',[Inf,Inf,Inf,Inf,Inf,pi,pi,pi,pi,pi,Inf],...
%     'StartPoint',[0.9,0.8,0.7,0.6,0.5,0.1,0.2,0.3,0.4,0.5,100]};
% s.fitOptions = fitoptions(fitSettings{:},fitLimits{:});
% s.fitEquation = fittype('a1*sin(2*pi*x+b1)+a2*sin(4*pi*x+b2)+a3*sin(6*pi*x+b3)+a4*sin(8*pi*x+b4)+a5*sin(10*pi*x+b5)+c','options',s.fitOptions);
% s.coefNames={'a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c','PR2'};

function getPulsatility(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_BFI_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_BFI_d.mat".');
end

for fidx=1:1:numel(fNames)

    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    load(s.fName)
    load(strrep(fNames{fidx},'_d.mat','_s.mat'));
    load(strrep(fNames{fidx},'_d.mat','_r.mat'));

    if strcmp(s.fitData,"all") || strcmp(s.fitData,"regions")
        eq=s.fitEquation;
        data=single(results.sData);
        fCoeffs=zeros(size(data,2),numel(s.coefNames),'single');        
        x=linspace(0,5,size(data,1).*5);
        xx=linspace(0,1,size(data,1));
        for ii=1:size(data,2)
            ts=double(repmat(squeeze(data(:,ii)),5,1));
            if all(~isnan(ts) & ts~=inf)
                [f,gof] = fit(x',ts,eq);
                data(:,ii)=feval(f,xx);
                fCoeffs(ii,:)=[coeffvalues(f),gof.rsquare];
            else
                data(:,ii)=NaN;
                fCoeffs(ii,:)=NaN;
            end
        end
        for i=6:1:10
            tmp=fCoeffs(:,i);
            tmp(tmp<0)=tmp(tmp<0)+2*pi;
            tmp=unwrap(tmp,[],3);
            if (sum(tmp<0)>0)
                tmp=tmp+2*pi;
            end
            fCoeffs(:,i)=tmp;
        end
        results.sData=data;
        for i=1:numel(s.coefNames)
            results.sMetrics.(s.coefNames{i})=fCoeffs(:,i);
        end

        if isfield(results,"dvsData")
            data=single(results.dvsData);
            fCoeffs=zeros(size(data,2),numel(s.coefNames),'single');   
            x=linspace(0,5,size(data,1).*5);
            xx=linspace(0,1,size(data,1));
            for ii=1:size(data,2)
                ts=double(repmat(squeeze(data(:,ii)),5,1));
                if all(~isnan(ts) & ts~=inf)
                    [f,gof] = fit(x',ts,eq);
                    data(:,ii)=feval(f,xx);
                    fCoeffs(ii,:)=[coeffvalues(f),gof.rsquare];
                    else
                data(:,ii)=NaN;
                fCoeffs(ii,:)=NaN;
                end
            end
            for i=6:1:10
                tmp=fCoeffs(:,i);
                tmp(tmp<0)=tmp(tmp<0)+2*pi;
                tmp=unwrap(tmp,[],3);
                if (sum(tmp<0)>0)
                    tmp=tmp+2*pi;
                end
                fCoeffs(:,i)=tmp;
            end
            results.dvsData=data;
            for i=1:numel(s.coefNames)
                results.dvsMetrics.(s.coefNames{i})=fCoeffs(:,i);
            end

            data=single(results.dvsDiameter);
            fCoeffs=zeros(size(data,2),numel(s.coefNames),'single');   
            x=linspace(0,5,size(data,1).*5);
            xx=linspace(0,1,size(data,1));
            for ii=1:size(data,2)
                ts=double(repmat(squeeze(data(:,ii)),5,1));
                if all(~isnan(ts) & ts~=inf)
                    [f,gof] = fit(x',ts,eq);
                    data(:,ii)=feval(f,xx);
                    fCoeffs(ii,:)=[coeffvalues(f),gof.rsquare];
                    else
                data(:,ii)=NaN;
                fCoeffs(ii,:)=NaN;
                end
            end
            for i=6:1:10
                tmp=fCoeffs(:,i);
                tmp(tmp<0)=tmp(tmp<0)+2*pi;
                tmp=unwrap(tmp,[],3);
                if (sum(tmp<0)>0)
                    tmp=tmp+2*pi;
                end
                fCoeffs(:,i)=tmp;
            end
            results.dvsDiameter=data;
            for i=1:numel(s.coefNames)
                results.dvsMetrics.([s.coefNames{i},'D'])=fCoeffs(:,i);
            end
        end
    end

    [maxVal,maxTime]=max(results.sData,[],1);
    [minVal,minTime]=min(results.sData,[],1);
    maxTime=results.time(maxTime);
    minTime=results.time(minTime);
    idxs=minTime>maxTime;
    minTime(idxs)=minTime(idxs)-results.time(end);
    results.sMetrics.('BFI')=mean(results.sData,1,'omitnan')';
    results.sMetrics.('std(BFI)')=std(results.sData,0,1,'omitnan')';
    results.sMetrics.('totalTime')=results.time(end).*ones(height(results.sMetrics),1);
    results.sMetrics.('PI(BFI)')=((maxVal-minVal)./mean(results.sData,1,'omitnan'))';
    results.sMetrics.('RI(BFI)')=((maxVal-minVal)./maxVal)';
    results.sMetrics.('minBFI')=minVal';
    results.sMetrics.('maxBFI')=maxVal';
    results.sMetrics.('timeMaxBFI')=maxTime';
    results.sMetrics.('timeMinBFI')=minTime';    
    results.sMetrics.('integralBFI')=trapz(results.time,results.sData,1)';
    results.sMetrics.('integralRBFI')=trapz(results.time,results.sData./minVal,1)';
    results.sMetrics.('integralDBFI')=trapz(results.time,results.sData-minVal,1)';


    if isfield(results,"dvsData")
        [maxVal,maxTime]=max(results.dvsData,[],1);
        [minVal,minTime]=min(results.dvsData,[],1);
        maxTime=results.time(maxTime);
        minTime=results.time(minTime);
        idxs=minTime>maxTime;
        minTime(idxs)=minTime(idxs)-results.time(end);
        results.dvsMetrics.('BFI')=mean(results.dvsData,1,'omitnan')';
        results.dvsMetrics.('std(BFI)')=std(results.dvsData,0,1,'omitnan')';
        results.dvsMetrics.('totalTime')=results.time(end).*ones(height(results.dvsMetrics),1);
        results.dvsMetrics.('PI(BFI)')=((maxVal-minVal)./mean(results.dvsData,1,'omitnan'))';
        results.dvsMetrics.('RI(BFI)')=((maxVal-minVal)./maxVal)';
        results.dvsMetrics.('minBFI')=minVal';
        results.dvsMetrics.('maxBFI')=maxVal';
        results.dvsMetrics.('timeMaxBFI')=maxTime';
        results.dvsMetrics.('timeMinBFI')=minTime';        
        results.dvsMetrics.('integralBFI')=trapz(results.time,results.dvsData,1)';
        results.dvsMetrics.('integralRBFI')=trapz(results.time,results.dvsData./minVal,1)';
        results.dvsMetrics.('integralDBFI')=trapz(results.time,results.dvsData-minVal,1)';

        [maxVal,maxTime]=max(results.dvsDiameter,[],1);
        [minVal,minTime]=min(results.dvsDiameter,[],1);
        maxTime=results.time(maxTime);
        minTime=results.time(minTime);
        idxs=minTime>maxTime;
        minTime(idxs)=minTime(idxs)-results.time(end);
        results.dvsMetrics.('D')=mean(results.dvsDiameter,1,'omitnan')';
        results.dvsMetrics.('PI(D)')=((maxVal-minVal)./mean(results.dvsDiameter,1,'omitnan'))';
        results.dvsMetrics.('RI(D)')=((maxVal-minVal)./maxVal)';
        results.dvsMetrics.('minD')=minVal';
        results.dvsMetrics.('maxD')=maxVal';
        results.dvsMetrics.('timeMaxD')=maxTime';
        results.dvsMetrics.('timeMinD')=minTime';
        results.dvsMetrics.('integralD')=trapz(results.time,results.dvsDiameter,1)';
        results.dvsMetrics.('integralRD')=trapz(results.time,results.dvsDiameter./minVal,1)';
        results.dvsMetrics.('integralDD')=trapz(results.time,results.dvsDiameter-minVal,1)';
    end



    if strcmp(s.fitData,"all")
        data=single(source.data);
        fCoeffs=zeros(size(data,1),size(data,2),numel(s.coefNames),'single');   
        x=linspace(0,5,size(data,3).*5);
        xx=linspace(0,1,size(data,3));
        mask=results.cMask>0;
        for i=1:1:size(data,1)
            subdata=squeeze(data(i,:,:));
            subcoeffs=zeros(size(data,2),numel(s.coefNames));
            parfor ii=1:size(data,2)
                if mask(i,ii)==1
                    ts=double(repmat(squeeze(subdata(ii,:)),1,5));
                    if all(~isnan(ts) & ts~=inf)
                    [f,gof] = fit(x',ts',eq);
                    subdata(ii,:)=feval(f,xx);
                    subcoeffs(ii,:)=[coeffvalues(f),gof.rsquare];
                    else
                subdata(ii,:)=NaN;
                subcoeffs(ii,:)=NaN;
                    end
                end
            end
            if any(i==round(linspace(1,size(data,1),11)))
                fprintf('\rFitting the source data %3d%%',round(100*i./(size(data,1))));
                drawnow limitrate
            end
            data(i,:,:)=subdata;
            fCoeffs(i,:,:)=subcoeffs;
        end
        fprintf('\n');


        for i=6:1:10
            img=fCoeffs(:,:,i);
            img(img<0)=img(img<0)+2*pi;
            img=unwrap(img,[],3);
            if (sum(img<0)>0)
                img=img+2*pi;
            end
            fCoeffs(:,:,i)=img;
        end
        source.data=data;
        results.fCoeffs=fCoeffs;
    end

    [maxVal,maxTime]=max(source.data,[],3);
    [minVal,minTime]=min(source.data,[],3);
    maxTime=source.time(maxTime);
    minTime=source.time(minTime);
    idxs=minTime>maxTime;
    minTime(idxs(:))=minTime(idxs(:))-results.time(end);
    results.imgPI=(maxVal-minVal)./mean(source.data,3,'omitnan');
    results.imgRI=(maxVal-minVal)./maxVal;
    results.extendedMetrics.imgMinBFI=minVal;
    results.extendedMetrics.imgMaxBFI=maxVal;
    results.extendedMetrics.imgTimeMaxBFI=maxTime;
    results.extendedMetrics.imgTimeMinBFI=minTime;
    results.extendedMetrics.imgIntegralBFI=trapz(source.time,source.data,3);
    results.extendedMetrics.imgIntegralRBFI=trapz(source.time,source.data./minVal,3);
    results.extendedMetrics.imgIntegralDBFI=trapz(source.time,source.data-minVal,3);

    settings.pulsatilityCalculation=s;
    disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    if strcmp(s.fitData,"all")
        save(fNames{fidx},'source','-v7.3');
    end
    save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
    save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');
end
end
