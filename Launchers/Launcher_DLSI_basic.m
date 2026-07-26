% Launcher_DLSI_basic - DLSI pipeline: raw .mraw recordings to a per-pixel fit.
%
%   Batch launcher for Dynamic Light Scattering Imaging.  STEP 1 turns one or
%   more Photron .mraw high-speed recordings into the normalized intensity
%   autocorrelation g2; STEP 2 fits every pixel of each g2 stack to a DLSI model
%   with fitDLSI; STEP 3 shows the fitted parameter maps.  Run STEP 0 once per
%   MATLAB session, then run the step cells (%%) in order.
%
%   Everything is written using the library's three-file data model, one triplet
%   per recording (the heavy g2 is stored once, never duplicated into the fit):
%       *_g_d.mat  SOURCE   - source.g2 (heavy), source.data (binned), source.lags
%       *_g_r.mat  RESULTS  - contrast maps, results.iniTau, results.mask and the
%                             model fit(s) (results.mMDSN, results.mDSN, ...)
%       *_g_s.mat  SETTINGS - settings.getNormalizedG2, settings.fitDLSI
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.  Header generation and
% script formatting were done with Claude Code.

%% STEP 0 - RUN EVERY TIME YOU RESTART MATLAB - THEN PROCEED TO THE STEP YOU STOPPED AT (BY DEFAULT STEP 1)
%LIBRARY PATH - add YOUR path manually here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
rootFolder = 'C:\Dropbox\Work\Data'; %root folder for the raw .mraw files lookup


%% STEP 1 Process raw .mraw recordings into the normalized autocorrelation g2
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - AUTOCORRELATION
s.lagMax=151; %number of lags (g2 has lagMax+1 points spanning lags 0..lagMax)

%ADJUSTED IF NECESSARY - PERFORMANCE
s.procType='gpu'; %'gpu' if a high-end GPU is available, 'cpu' otherwise (auto-falls back to cpu)
s.bandsN=5; %horizontal bands read at a time (raise it for recordings larger than RAM)

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*.mraw'); %the companion .cihx must sit next to each .mraw
%fNames{1}='C:\Dropbox\Work\Data\HSLS_recording.mraw'; %or list the files manually

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    fName=fNames{i,1};
    fprintf('STEP 1 - g2 for %d/%d: %s\n',i,size(fNames,1),fName);

    %read metadata (first frame only) and derive the band tiling / time binning
    [frame1,fps,~,sizeT]=readMRAW(fName,0,1);
    sizeY=size(frame1,1); sizeX=size(frame1,2);
    bandSize=ceil(sizeY/s.bandsN);
    nTimeBins=floor(sizeT/s.lagMax);
    framesUsed=nTimeBins*s.lagMax;

    %process the recording band by band (keeps memory bounded)
    source=struct();
    source.g2=zeros(sizeY,sizeX,s.lagMax+1,'single');
    source.data=zeros(sizeY,sizeX,nTimeBins,'single'); %temporally binned intensity
    for b=1:s.bandsN
        yRange=[(b-1)*bandSize+1,min(b*bandSize,sizeY)];
        fprintf('   band %d/%d (rows %d-%d)\n',b,s.bandsN,yRange(1),yRange(2));
        data=readMRAW(fName,0,[],[yRange;1,sizeX]);

        %temporally binned intensity stack (average every lagMax frames)
        tmp=data(:,:,1:framesUsed);
        tmp=reshape(tmp,size(tmp,1),size(tmp,2),s.lagMax,nTimeBins);
        source.data(yRange(1):yRange(2),:,:)=squeeze(mean(tmp,3));

        %normalized intensity autocorrelation for this band
        source.g2(yRange(1):yRange(2),:,:)=getNormalizedG2(data,'lagMax',s.lagMax,'procType',s.procType);
    end
    source.lags=(0:s.lagMax)/fps; %seconds

    %companion contrast maps (results).  getK has no GPU fallback, and the binned
    %stack is small, so the spatial contrast is computed on the CPU.
    results=struct();
    results.tLSCI=std(source.data,0,3)./mean(source.data,3); %temporal contrast
    sLSCI=getK(source.data,'spatial','kernelSize',5,'procType','cpu');
    results.sLSCI=mean(sLSCI,3); %spatial contrast
    results.betaMax=max(results.sLSCI(:)).^2; %coherence upper bound

    %save the three-file triplet (base name = recording without the .mraw extension)
    s.fps=fps; s.fileName=fName;
    settings=struct(); settings.getNormalizedG2=s;
    base=regexprep(fName,'\.mraw$','','ignorecase');
    save([base '_g_d.mat'],'source','-v7.3');   %heavy data
    save([base '_g_r.mat'],'results','-v7.3');  %maps
    save([base '_g_s.mat'],'settings','-v7.3'); %settings
end
disp('STEP 1 complete');

%% STEP 2 Fit every pixel of each g2 stack to a DLSI model
close all
clearvars -except fNames libraryFolder rootFolder
s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - MODEL
s.type='MDSN'; %'MDSN' (mixed dynamics+static+noise), 'DSN', 'DN' or 'D'
s.cLoss='sqrtbeta'; %coherence-loss form of the static/dynamic cross-term:
% 'sqrtbeta' (default), 'beta', or {'sqrtbeta','beta'} to fit both and keep the
% better one per pixel (the choice is stored in model.fit.cLoss).  MDSN/DSN only.

%ADJUSTED IF NECESSARY - FIT
s.pointsMin=5;      %minimum number of lag points used per pixel
s.isAdaptive=true;  %adapt the beta/p bounds on the fly from already-fitted pixels
s.spatialDS=1;      %>1 (e.g. 2 or 4) spatially downsamples for a quick test - fitting a full frame can take many minutes

%SET FILE NAMES HERE
fNames=getFileNamesList(rootFolder,'*_g_d.mat'); %the g2 sources written by STEP 1

%RUN THE PROCESSING ROUTINE
for i=1:1:size(fNames,1)
    dName=fNames{i,1};
    fprintf('STEP 2 - fit %d/%d: %s\n',i,size(fNames,1),dName);

    %load the g2 source and the existing results/settings (STEP 1) to augment them
    load(dName,'source');
    base=regexprep(dName,'_g_d\.mat$','');
    rName=[base '_g_r.mat']; sName=[base '_g_s.mat'];
    if isfile(rName), load(rName,'results');  else, results  = struct(); end
    if isfile(sName), load(sName,'settings'); else, settings = struct(); end
    fps=1/(source.lags(2)-source.lags(1));

    %(optional) spatial downsampling for a quick test
    g2=source.g2;
    if s.spatialDS>1
        k=ones(s.spatialDS,s.spatialDS,1);
        g2=convn(g2,k,'same')./convn(ones(size(g2)),k,'same');
        g2=g2(1:s.spatialDS:end,1:s.spatialDS:end,:);
    end
    lags=(0:size(g2,3)-1)/fps; %seconds

    %initial decorrelation-time map (seeds and orders the fit) and valid-pixel mask
    results.iniTau=getTauC(g2,'method','decay','fps',fps);
    results.mask=~((g2(:,:,1)>2)|(g2(:,:,2)<1));

    %fit the model - a lightweight struct stored as results.m<type> (e.g. mMDSN).
    %The heavy g2 is NOT copied into the results; it stays in *_g_d.mat.
    mName=['m',s.type];
    results.(mName)=fitDLSI(g2,lags,results.iniTau,'type',s.type, ...
        'pointsMin',s.pointsMin,'isAdaptive',s.isAdaptive,'cLoss',s.cLoss,'mask',results.mask);
    % To store more models side by side, uncomment (each is a lightweight struct):
    % results.mDSN=fitDLSI(g2,lags,results.iniTau,'type','DSN','pointsMin',s.pointsMin,'mask',results.mask);
    % results.mDN =fitDLSI(g2,lags,results.iniTau,'type','DN', 'pointsMin',s.pointsMin,'mask',results.mask);
    settings.fitDLSI=s;

    %save results (g_r) and settings (g_s); the g2 stays in g_d and is not rewritten
    save(rName,'results','-v7.3');
    save(sName,'settings','-v7.3');
end
disp('STEP 2 complete');

%% STEP 3 See the results - initial and fitted parameter maps for one recording
close all
clearvars -except fNames libraryFolder rootFolder

%SET FILE NAME HERE (defaults to the first g2 source found)
fNames=getFileNamesList(rootFolder,'*_g_d.mat');
dName=fNames{1,1};
base=regexprep(dName,'_g_d\.mat$','');
load([base '_g_r.mat'],'results');
load([base '_g_s.mat'],'settings');
mName=['m',settings.fitDLSI.type];   %the model fitted in STEP 2
model=results.(mName);

f=figure('Color','w','Name','DLSI - fitted maps'); f.WindowState='maximized';
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
imagesc(results.iniTau*1e6); axis image off; colormap(gca,'jet');
clim(prctile(results.iniTau(results.iniTau>0)*1e6,[1 99])); colorbar;
title('Initial \tau_c (\mus)');

nexttile;
imagesc(model.fit.tau*1e6); axis image off; colormap(gca,'jet');
clim(prctile(model.fit.tau(model.fit.tau>0)*1e6,[1 99])); colorbar;
title(sprintf('%s fitted \\tau (\\mus)',settings.fitDLSI.type));

nexttile;
imagesc(model.fit.beta); axis image off; colormap(gca,'jet');
clim(prctile(model.fit.beta(model.fit.beta>0),[1 99])); colorbar;
title('Fitted \beta');
drawnow
