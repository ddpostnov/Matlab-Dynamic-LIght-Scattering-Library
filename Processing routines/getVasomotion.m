function getVasomotion(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_BFI_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_BFI_d.mat".');
end

for fidx=1:1:numel(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    if s.analysePerPixel
        load(s.fName,'source')
    end
    load(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings');
    load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

    %estimate sampling frequency
    fs=squeeze(1./mean(results.time(2:end)-results.time(1:end-1)));

    if s.keepSpectrum && ~isempty(s.tgtFS) && s.tgtFS>fs
        error('Target sampling frequency cannot be higher than source sampling frequency')
    end

    %define cwt filter bank
    fb = cwtfilterbank('SignalLength', numel(results.time), ...
        'SamplingFrequency', fs, ...
        'Wavelet', 'amor', ...
        'FrequencyLimits',s.wFR,...
        'VoicesPerOctave', s.wVPO);


    [~,fwt,coi]=wt(fb,results.time);
    f=cat(1,fwt,s.vFR',s.cFR');
    f=sort(unique(f),'descend');
    idxsVFR=f>=s.vFR(1) & f<=s.vFR(2);
    idxsCFR=f>=s.cFR(1) & f<=s.cFR(2);
    timeVSM=results.time(coi<0.05);
    avgN=[];


    fn = fieldnames(results);
    for k=1:numel(fn)
        if strcmp(fn{k}, 'sData') || strcmp(fn{k}, 'dvsData') || strcmp(fn{k}, 'dvsDiameter')
            subDataIn=(results.(fn{k})-mean(results.(fn{k}),1))./mean(results.(fn{k}),1);
            if numel(subDataIn)>0
                sz=size(subDataIn);
                if s.reconstructData
                    rData=zeros(numel(timeVSM),sz(2));
                end
                spctPCTS=zeros(sz(2),numel(f),numel(s.pcts)-1,2);
                spctMean=zeros(sz(2),numel(f));
                spctSTD=zeros(sz(2),numel(f));
                spctSKW=zeros(sz(2),numel(f));
                spctKRT=zeros(sz(2),numel(f));
                spctFLR=zeros(sz(2),numel(f),2);
                spctSLC=zeros(sz(2),numel(f),2);
                statFLR=zeros(sz(2),4);
                statSLC=zeros(sz(2),4);

                adVFR=zeros(sz(2),2);
                adFVFR=zeros(sz(2),2);
                adSVFR=zeros(sz(2),2);
                adCFR=zeros(sz(2),2);
                adFCFR=zeros(sz(2),2);
                adSCFR=zeros(sz(2),2);

                if s.keepClustering
                    vTS=zeros(numel(timeVSM),sz(2));
                    vFlare=zeros(numel(timeVSM),sz(2));
                end

                if s.keepSpectrum
                    if ~isempty(s.tgtFS)
                        avgN=floor(fs./s.tgtFS/2)*2+1;
                        tWT=movmean(timeVSM,avgN,'Endpoints','discard');
                        tWT=tWT(1:avgN:end);
                    else
                        tWT=timeVSM;
                    end
                    sWT=zeros(sz(2),numel(f),numel(tWT),'single');
                end

                keepSpectrum=s.keepSpectrum;
                tgtFS=s.tgtFS;
                pcts=s.pcts;
                vFR=s.vFR;
                cFR=s.cFR;
                maxN = s.otsuMaxN;
                otsuElbow =s.otsuElbow;
                reconstructData=s.reconstructData;
                keepClustering=s.keepClustering;

                parfor i=1:sz(2)
                    if sum(isnan(subDataIn(:,i)),'all')==0
                        wtts=wt(fb,squeeze(subDataIn(:,i)));

                        %Synchrosqueezing is worth considering to use in the future
                        %wtts = wsst(squeeze(subDataIn(:,i)), fs, 'amor', 'VoicesPerOctave', s.wVPO);
                        %wtts=abs(wtts(idxsWFR,coi<0.05));

                        if reconstructData
                            tmp=icwt(wtts,'amor',fwt,vFR);
                            rData(:,i)=tmp(coi<0.05);
                        end

                        wtts=abs(wtts(:,coi<0.05));
                        wtts=interp1(fwt,wtts,f);
                        spctMean(i,:)=mean(wtts,2);

                        spctSTD(i,:)=std(wtts,0,2);
                        spctSKW(i,:)=skewness(wtts,0,2);
                        spctKRT(i,:)=kurtosis(wtts,0,2);

                        if keepSpectrum
                            if ~isempty(tgtFS)
                                tmp=movmean(wtts,avgN,2,'Endpoints','discard');
                                sWT(i,:,:)=tmp(:,1:avgN:end);
                            else
                                sWT(i,:,:)=wtts;
                            end
                        end

                        % get spectrums for time percentiles.
                        ts=-trapz(f(idxsVFR),wtts(idxsVFR,:),1);
                        ts=ts(:)./(vFR(2)-vFR(1));
                        adVFR(i,:)=[mean(ts),std(ts)];

                        tsc=-trapz(f(idxsCFR),wtts(idxsCFR,:),1);
                        tsc=tsc(:)./(cFR(2)-cFR(1));
                        adCFR(i,:)=[mean(tsc),std(tsc)];


                        vpcts=prctile(ts,pcts);
                        tmp=zeros(numel(f),numel(pcts)-1,2);
                        for ii=1:1:numel(vpcts)-1
                            tmp(:,ii,1)=mean(wtts(:,ts>=vpcts(ii) & ts<(vpcts(ii+1)+eps)),2);
                            tmp(:,ii,2)=std(wtts(:,ts>=vpcts(ii) & ts<(vpcts(ii+1)+eps)),0,2);
                        end
                        spctPCTS(i,:,:,:)=tmp;


                        % %worth considering flattening the signal in the future
                        % intervalGroups = discretize((1:numel(ts))', unique([1;find(islocalmin(ts));numel(ts)]));
                        % intervalMeans = splitapply(@mean, ts(~isnan(intervalGroups)), intervalGroups(~isnan(intervalGroups)));
                        % tmp=ts;
                        % tmp(~isnan(intervalGroups)) = intervalMeans(intervalGroups(~isnan(intervalGroups)));


                        %Evaluate optimal number for otsu clasters (elbow method)
                        optN = zeros(maxN, 1);
                        for n = 1:maxN, [~, optN(n)] = multithresh(ts, n); end
                        optN = diff(optN);
                        optN = find(optN < otsuElbow, 1, 'first');
                        if isempty(optN)
                            flare=zeros(size(ts));
                        else
                            minT=ceil((fs)./min(vFR)); %minimum flare duration
                            flare = multithresh(ts,optN);
                            flare=ts>flare(1);
                            tmp=padarray(flare,[minT,0],"replicate");
                            tmp=bwareaopen(tmp,minT)>0;
                            flare=tmp(minT+1:end-minT);
                        end


                        if any(flare)
                            spctFLR(i,:,:)=cat(2,mean(wtts(:,flare),2),std(wtts(:,flare),0,2));
                            statFLR(i,:)=[sum(flare)/fs,sum(flare)/fs,0,NaN];
                            statSLC(i,:)=[0,0,0,NaN];

                            adFVFR(i,:)=[mean(ts(flare)),std(ts(flare))];
                            adFCFR(i,:)=[mean(tsc(flare)),std(tsc(flare))];
                        end
                        if any(~flare)
                            spctSLC(i,:,:)=cat(2,mean(wtts(:,~flare),2),std(wtts(:,~flare),0,2));
                            statFLR(i,:)=[0,0,0,NaN];
                            statSLC(i,:)=[sum(~flare)/fs,sum(~flare)/fs,0,NaN];

                            adSVFR(i,:)=[mean(ts(~flare)),std(ts(flare))];
                            adSCFR(i,:)=[mean(tsc(~flare)),std(tsc(flare))];
                        end
                        if any(flare) && any(~flare)
                            toggleIdx = find(diff(flare(:)'));
                            internalLengths = diff(toggleIdx);
                            internalStates = flare(toggleIdx(1:end-1) + 1);
                            durFLR = internalLengths(internalStates == 1);
                            durSLC = internalLengths(internalStates == 0);
                            statFLR(i,:)=[sum(flare)/fs,mean(durFLR),std(durFLR),median(durFLR)]/fs;
                            statSLC(i,:)=[sum(~flare)/fs,mean(durSLC),std(durSLC),median(durSLC)]/fs;
                        end


                        if keepClustering
                            vTS(:,i)=ts;
                            vFlare(:,i)=flare;
                        end
                    end
                end

                toc

                results.vsm.(fn{k}).adFVFR=adFVFR;
                results.vsm.(fn{k}).adSVFR=adSVFR;
                results.vsm.(fn{k}).adFCFR=adFCFR;
                results.vsm.(fn{k}).adSCFR=adSCFR;


                results.vsm.(fn{k}).statSLC=statSLC;
                results.vsm.(fn{k}).statFLR=statFLR;
                results.vsm.(fn{k}).spctFLR=spctFLR;
                results.vsm.(fn{k}).spctSLC=spctSLC;
                results.vsm.(fn{k}).spctSKW=spctSKW;
                results.vsm.(fn{k}).spctKRT=spctKRT;
                results.vsm.(fn{k}).spctMean=spctMean;
                results.vsm.(fn{k}).spctSTD=spctSTD;
                results.vsm.(fn{k}).f=f;
                results.vsm.(fn{k}).spctPCTS=spctPCTS;
                results.vsm.(fn{k}).coordPCTS=(s.pcts(1:end-1)+s.pcts(2:end))./2;
                if s.keepClustering
                    results.vsm.(fn{k}).vTS=vTS;
                    results.vsm.(fn{k}).vFlare=vFlare;
                end
                if reconstructData
                    results.vsm.(fn{k}).rData=rData.*mean(results.(fn{k}),1) + mean(results.(fn{k}),1);
                end
                if s.keepSpectrum
                    results.vsm.(fn{k}).tWT=tWT;
                    results.vsm.(fn{k}).sWT=sWT;
                end

                switch fn{k}
                    case 'sData'
                        results.sMetrics.('adVFR')=adVFR(:,1);
                        results.sMetrics.('std(adVFR)')=adVFR(:,2);
                        results.sMetrics.('adCFR')=adCFR(:,1);
                        results.sMetrics.('std(adCFR)')=adCFR(:,2);
                    case 'dvsData'
                         results.dvsMetrics.('adVFR')=adVFR(:,1);
                        results.dvsMetrics.('std(adVFR)')=adVFR(:,2);
                        results.dvsMetrics.('adCFR')=adCFR(:,1);
                        results.dvsMetrics.('std(adCFR)')=adCFR(:,2);
                    case 'dvsDiameter'
                         results.dvsMetrics.('adDVFR')=adVFR(:,1);
                        results.dvsMetrics.('std(adDVFR)')=adVFR(:,2);
                        results.dvsMetrics.('adDCFR')=adCFR(:,1);
                        results.dvsMetrics.('std(adDCFR)')=adCFR(:,2);
                end
            end
        end
    end


    % % For future - delay estimate based on correlation of reconstructed data
    % figure
    % for k=1:1:size(rData,3)
    % delay=zeros(size(rData,2),1);
    % cc=zeros(size(rData,2),1);
    % for i=1:1:size(rData,2)
    % [xc, lags] = xcorr(rData(:,i,k),rData(:,14,k), 400, 'normalized');
    % [cc(i), max_idx] = max(xc);
    % delay(i) = lags(max_idx) / fs;
    % end
    % subplot(1,2,1)
    % imagesc(delay(results.sMap))
    % axis image
    % subplot(1,2,2)
    % imagesc(cc(results.sMap))
    % axis image
    % title(num2str(rFR(k)))
    % pause(0.1)
    % end

    if s.analysePerPixel
        sz=size(source.data);
        results.vsm.ppxAD=zeros(sz(1),sz(2),2,sum(coi<0.05),'single');
        tic
        for i=1:sz(1)
            subDataIn=squeeze((source.data(i,:,:)-mean(source.data(i,:,:),3))./mean(source.data(i,:,:),3));
            subDataOut=squeeze(results.vsm.rwt(i,:,:,:));
            parfor j=1:sz(2)
                wtts=abs(wt(fb,squeeze(subDataIn(j,:))));
                wtts=abs(wtts(:,coi<0.05));
                wtts=interp1(fwt,wtts,f);
                ts=-trapz(f(idxsVFR),wtts(idxsVFR,:),1);
                ts=ts(:)./(vFR(2)-vFR(1));
                tsc=-trapz(f(idxsCFR),wtts(idxsCFR,:),1);
                tsc=tsc(:)./(cFR(2)-cFR(1));
                subDataOut(j,:,:)=cat(1,ts',tsc');
            end
            results.vsm.ppxAD(i,:,:,:)=subDataOut;
            fprintf('Batch %d/%d processed. Elapsed: %.2fs\n', i, sz(1), toc);
        end
    end

    settings.getVasomotions=s;
    disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    save(strrep(fNames{fidx},'_d.mat','_r.mat'),'results','-v7.3');
    save(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings','-v7.3');

end