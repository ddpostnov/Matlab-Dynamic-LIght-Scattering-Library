% Either used to compare steady-state perfusion or responses to intervention.
% Use Launcher_slowDynamics.m. Ensure to edit the file names and settings according to the project needs.


%% STEP 1 Process .rls files to get the temoiral contrast for segmentation and vasomotion analysis
%LIBRARY PATH - add YOUR path manualy here:
libraryFolder = 'C:\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
libraryFolder = 'C:\Users\AU707705\Dropbox\Work\GitHub\Matlab-Dynamic-LIght-Scattering-Library';
addpath(genpath(libraryFolder));
s.libraryFolder=libraryFolder;

s.fBolus=[301,1500];
s.fAngio=[2101,9800];

%SET FILE NAMES HERE
%OPTION 1 - AUTOMATIC LOOKUP USING REGULAR EXPRESSIONS
s.reCalculate=true; %rewrites the files if existing
rootFolder = 'C:\Dropbox\Work\Data'; %root folder for the files lookup
files      = dir(fullfile(rootFolder,'**','*BB.cxd')); %<--- use regexp to define the files of interest
fNames     = fullfile({files.folder}', {files.name}');
if ~s.reCalculate
    outNames  = regexprep(fNames, '\.cxd$', '_I_d.mat');
    hasOut    = cellfun(@(f) isfile(f), outNames);
    fNames    = fNames(~hasOut);
end

getBolus(s,fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 2 Define pixel categories 
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - CONTRAST CALCULATION
s.maxK=0.5; % Maximum valid contrast - helps with initial masking
s.minK=0.0001; % Minimum valid contrast
s.regionsN=1; %Numer of regions for manual selection. 0 if using entire window.
s.lSizeN=121; % Odd, approximately 2 times larger than the largest vessel
s.sSizeN=5; % Odd, approximately 2 times larger than small vessels diameter
s.sens=0.3; % Segmentation sensitivity - increase if missing vessels, decrease to minimize segmentation noise
s.deSens=1;

%ADJUSTED IF NECESSARY - SEGMENTATION ADJUSTEMNTS
s.lThinN=2; % Large vessels thinning (appears as internal edges)
s.imOpen=0; % Small vessels thinning (appears as internal edges)
s.iniSizeN=7; % Odd number equal or larger than the spatial contrast kernel

%DO NOT CHANGE - META DATA
s.categories={'background','parenchyma','unsegmented','outerEdge','innerEdge','lumen'}; %CATEGORIES

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*_b_I_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

getCategories(s,fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 4 Assign same categories to the internal cycle (pulsatility) data 
% (OPTIONAL: used in combination with categories extration based on temporal contrast analysis only)

close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
files      = dir(fullfile(rootFolder,'**','*c_K_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');
fNamesRef=regexprep(fNames, '\_c_K_d.mat$', '_t_K_d.mat');

assignCategories(fNames,fNamesRef); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 5 (OPTIONAL. Only use if 1 or more regions are defined in step 2) Split the regions. 
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
s.deleteOriginal=false; %true or false. USE TRUE IF YOU DO NOT PLAN TO RE-DEFINE REGIONS

files      = dir(fullfile(rootFolder,'**','*_K_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

%USES THE SAME FILE NAMES AS ABOVE as STEP 2
splitRegions(s,fNames); %LAUNCHES THE UTILITY ROUTINE

%% STEP 6 Perform segmentation
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED (OR VERIFIED) PER PROTOCOL - BASIC PARAMETERS
s.attmemptDS=false; %attempt to perform automated dynamic segmentation or not
s.sMinL=15; % Minimum length for segments
s.simR=0.3;
s.difR=0.4;
s.correctNodes=true;

s.prchNSize=90; % Parenchymal pixels neighbourhoud.

%ADJUSTED (OR VERIFIED) PER PROTOCOL - DYNAMIC SEGMENTATION
s.sMinP2R2=0.95; %Min accepted R2 of 3-degree polynom fit
s.sMaxLBI=(1/5)./s.sMinL; %Max local bending (0 to pi per pixel)
s.sMaxCLR=1.3; %Maximum accepted CLR of the segment 1 perfectly straight, 1.5 - slow bend, 2 - coil
s.sMaxDK=0.2; %Max accepted std/mean for the initial diameter estimation
s.sMaxKK=0.3; %Max accepted std/mean for the initial contrast estimation
s.iniNSize=7; % Odd number equal or larger than the spatial contrast kernel
s.sMaxP2D=3; %Max accepted deviation of the fit from center estimate

%ADJUSTED IF NECESSARY - QUALITY CHECK AND INTERPOALTION
s.minOverlapMask=0.6; %minimum overlap between the initial center line and segmentation mask present in each frame
s.minOverlapSelf=0.2; %minimum size of segmented area compared to the initial ROI
s.pInterpF=10; % leave as is

%SET FILE NAMES HERE - IF REGIONS SPLIT WAS PERFORMED, OTHERWISE KEEP THE
%SAME NAMES
files      = dir(fullfile(rootFolder,'**','*_b_I_d.mat')); 
fNames     = fullfile({files.folder}', {files.name}');

getSegmentation(s, fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 7 (OPTIONAL. Use if multiple recordings of the same field of view have to be compared to each other) Register LSCI files to the first file in the list
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;

%ADJUSTED IF NECESSARY - REGISTRATION SETTINGS
[s.optimizer,s.metric] = imregconfig("monomodal");
s.optimizer.MaximumIterations=500;
s.tFormType='affine';
s.matchSegmentation=true;
s.prchNSize=30; % Parenchymal pixels neighbourhoud.

files      = dir(fullfile(rootFolder,'**','*t_K_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_K_d.mat$', '_c_K_d.mat');
fNames=cat(1,fNamesPLS,fNamesVSM);
registerLSCItoLSCI(s,fNames); %LAUNCHES THE UTILITY ROUTINE

%% STEP 8 Convert contrast to blood flow index
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
s.deleteOriginal=true; %true or false
%ADJUSTED IF NECESSARY - CONVERSION METHOD
s.method="basic"; %only "basic" is avaliable

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*_K_d.mat'));
fNames     = fullfile({files.folder}', {files.name}');

getBFI(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 9 Perform vasomotion analysis
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
s.vFR=[0.05,0.25];
s.cFR=[0.4,0.6];
s.wFR=[0.01,1];
s.wVPO=10;
s.tgtFS=1; %Hz
s.pcts=0:10:100;

s.otsuMaxN=5;
s.otsuElbow= 0.05;

s.analysePerPixel  = false;
s.keepSpectrum=false;
s.keepClustering=true;
s.reconstructData=true;

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*t_BFI_d.mat'));
fNames     = fullfile({files.folder}', {files.name}');

getVasomotion(s,fNames);  %LAUNCHES THE PROCESSING ROUTINE

%% STEP 10 Perform pulsatility analysis (strictly after conversion to BFI)
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
%ADJUSTED (OR VERIFIED) PER PROTOCOL - Waveform fitting
s.fitData="regions"; %"regions" - for all segmented data, "all" - for segmented data AND pixel by pixel fits, or "none"
%ADJUSTED IF NECESSARY - Waveform fitting
fitSettings={'Method','NonlinearLeastSquares','Algorithm','Trust-Region','Display','off'};%,'Robust','Bisquare','MaxFunEvals',1200,'MaxIter',1000};
fitLimits={'Lower',[0,0,0,0,0,-pi,-pi,-pi,-pi,-pi,0],...
    'Upper',[Inf,Inf,Inf,Inf,Inf,pi,pi,pi,pi,pi,Inf],...
    'StartPoint',[0.9,0.8,0.7,0.6,0.5,0.1,0.2,0.3,0.4,0.5,100]};
s.fitOptions = fitoptions(fitSettings{:},fitLimits{:});
s.fitEquation = fittype('a1*sin(2*pi*x+b1)+a2*sin(4*pi*x+b2)+a3*sin(6*pi*x+b3)+a4*sin(8*pi*x+b4)+a5*sin(10*pi*x+b5)+c','options',s.fitOptions);
s.coefNames={'a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c','PR2'};

%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*c_BFI_d.mat')); %<---ALWAYS REFER TO "c_BFI_d.mat" files, but you may use regexp to define specific "c_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

getPulsatility(s, fNames); %LAUNCHES THE PROCESSING ROUTINE

%% STEP 11 Assign vessel types and regions of interest
close all
clearvars -except fNames libraryFolder rootFolder

s.libraryFolder=libraryFolder;
s.useReference=true; 
s.libraryFolder=libraryFolder;

files      = dir(fullfile(rootFolder,'**','*t_BFI_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_BFI_d.mat$', '_c_BFI_d.mat');
fNames=cat(1,fNamesPLS,fNamesVSM);
s.refFName=fNames{1};
setTypes(s,fNames);

%% STEP 11 (OPTIONAL) Export key results to an excel table
%SET FILE NAMES HERE
files      = dir(fullfile(rootFolder,'**','*_BFI_d.mat')); %<---ALWAYS REFER TO "_BFI_d.mat" files, but you may use regexp to define specific "_BFI_d.mat" files of interest
fNames     = fullfile({files.folder}', {files.name}');

exportToExcel(fNames); %LAUNCHES THE UTILITY ROUTINE

%%

figure
plot(squeeze(max(source.data,[],[1,2])))
hold on
plot(squeeze(mean(source.data,[1,2])))
hold off

%%

dataBL=source.data(:,:,1:575);
dataRP=source.data(:,:,576:end-100);
dataFN=source.data(:,:,end-99:end);

figure
cval=[mean(dataRP(:,:,1),'all'),prctile(max(dataRP,[],3),99,'all')];
for i=1:1:size(dataRP,3)
imagesc(dataRP(:,:,i))
clim(cval)
axis image
colormap gray
title(num2str(i/100))
pause(0.05)
end


dataRPN=single(dataRP)-mean(single(dataBL),3);%
dataRPN=mat2gray(dataRPN);
fSize=151;
parfor i=1:1:size(dataRPN,3)
i
dataRPN(:,:,i)=dataRPN(:,:,i)-imopen(medfilt2(dataRPN(:,:,i),[15,15],'symmetric'),strel('disk',fSize));
end

figure
cval=[mean(dataRPN(:,:,1),'all'),prctile(max(dataRPN,[],3),99,'all')];
for i=1:1:size(dataRPN,3)
imagesc(dataRPN(:,:,i))
clim(cval)
axis image
colormap gray
title(num2str(i/100))
pause(0.5)
end
%%
img=img-imopen(medfilt2(img,[15,15],'symmetric'),strel('disk',fSize));
data=applyDirectionalFilter(single(dataRP), img);
[fitP, fitQ, fitG,bMin,bMax] = fitGammaVariatePerPixel(data, time);
%%
figure
imagesc(mean(data,3))
mask=roipoly;


%%
figure
img=squeeze(MTT);
imagesc(img)
clim(prctile(img(mask(:)),[0.1,99.9]))
colorbar
axis image
colormap jet
title('MTT')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);

%%
figure
img=squeeze(fitP(:,:,3));
imagesc(img)
clim([0,2.5])
colorbar
axis image
colormap jet
title('Arrival time')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);

figure
img=squeeze(fitP(:,:,5));
imagesc(img)
clim(prctile(img(mask(:)),[0.1,99.9]))
colorbar
axis image
colormap jet
title('Peak time')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);

figure
img=squeeze(fitP(:,:,4));
imagesc(img)
clim([1,4])
colorbar
axis image
colormap jet
title('Alpha')
set(gcf,'Color','w');
xticklabels([]);
yticklabels([]);


%%
files      = dir(fullfile(rootFolder,'**','*t_BFI_d.mat')); %<---ALWAYS REFER TO "_K_d.mat" files, but you may use regexp to define specific "_K_d.mat" files of interest
fNamesVSM  = fullfile({files.folder}', {files.name}');
fNamesPLS=regexprep(fNamesVSM, '\_t_BFI_d.mat$', '_c_BFI_d.mat');

clearvars vsm
for i=1:1:numel(fNamesVSM)
load(strrep(fNamesVSM{i},'_d.mat','_r.mat'),'results');
vsm(i)=results;
clearvars results
end

clearvars pls
for i=1:1:numel(fNamesPLS)
load(strrep(fNamesPLS{i},'_d.mat','_r.mat'),'results');
pls(i)=results;
clearvars results
end

f=vsm(1).vsm.sData.f;
coordPCTS=vsm(1).vsm.sData.coordPCTS;
types={'Arteries','Parenchyma','Veins'};
conds={'AWK','ISO','K/X'};

spctPCTS=zeros(3,3,numel(f),numel(coordPCTS));

adVFR=zeros(numel(vsm),3);
adCFR=zeros(numel(vsm),3);


for i=1:1:numel(vsm)
    adVFR(i,1)= mean(vsm(i).sMetrics.adVFR(strcmp(vsm(i).sMetrics.type,"Artery")));
    adVFR(i,2)=mean(vsm(i).sMetrics.adVFR(vsm(i).sMetrics.category==1));
    adVFR(i,3)= mean(vsm(i).sMetrics.adVFR(strcmp(vsm(i).sMetrics.type,"Vein")));

    adCFR(i,1)= mean(vsm(i).sMetrics.adCFR(strcmp(vsm(i).sMetrics.type,"Artery")));
    adCFR(i,2)=mean(vsm(i).sMetrics.adCFR(vsm(i).sMetrics.category==1));
    adCFR(i,3)= mean(vsm(i).sMetrics.adCFR(strcmp(vsm(i).sMetrics.type,"Vein")));

    spctPCTS(i,1,:,:)= mean(vsm(i).vsm.sData.spctPCTS(strcmp(vsm(i).sMetrics.type,"Artery"),:,:,1),1);
    spctPCTS(i,2,:,:)= mean(vsm(i).vsm.sData.spctPCTS(vsm(i).sMetrics.category==1,:,:,1),1);
    spctPCTS(i,3,:,:)= mean(vsm(i).vsm.sData.spctPCTS(strcmp(vsm(i).sMetrics.type,"Vein"),:,:,1),1);
end

figure
for i=1:1:3
for j=1:1:3
subplot(3,3,(i-1)*3+j)
plot(f,squeeze(spctPCTS(i,j,:,:)))
title([conds{i},', ',types{j}])
ylim([0,0.2])
xlabel('Frequency, Hz')
ylabel('WT Amplitude')
end
end
set(gcf,'Color','w');

figure
for i=1:1:3
subplot(1,3,i)
plot(squeeze(adVFR(i,:)))
hold on
plot(squeeze(adCFR(i,:)))
hold off
xlim([0,4])
ylabel('Amplitude density')
legend({'VSM: 0.05-0.25Hz','CTR: 0.4-0.6Hz'})
xticks([1,2,3]);
xticklabels(types)
title(conds{i})
ylim([0.01,0.07])
end
set(gcf,'Color','w');
