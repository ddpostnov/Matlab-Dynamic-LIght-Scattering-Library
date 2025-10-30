close all
clear

%% Set the file names in a strctured way to allow multi-file comparison

%If you were smart with your filenaming:
rootFolder = 'C:\Dropbox\Work\Data\LH2';
for ai=1:1 %loop over animals
for ri=1:2 %loop over rois
files      = dir(fullfile(rootFolder,'**',['*_cycle_sLSCIMM_roi',num2str(ri),'*_pi.mat']));
files     = fullfile({files.folder}', {files.name}');
fNames(1:numel(files),ai,ri)=files;
end
end

%If not - do it manually.
% fNames{1,1,1}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";
% fNames{1,1,2}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";
% fNames{1,1,3}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";
% fNames{1,2,1}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";
% fNames{1,2,2}="O:\HE_VSM\ChristianStaehr\data\Cardiac_arrest_CBF_mouse\LSCI\ID251\20241217_ID251_2hfollowup.rls";

%% Load the first file and inspect what is inside
load(fNames{1,1,1})
r=results;

%You can plot BFI and PI images as well as the masks obtained in previous
%steps
figure
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
t1=nexttile;
imagesc(r.imgBFI)
clim(prctile(r.imgBFI(r.cMask(:)>0),[5,99]));
title('BFI')
colorbar
axis image
t2=nexttile;
imagesc(r.imgPI)
clim(prctile(r.imgPI(r.cMask(:)>0),[5,99]));
title('PI')
axis image
t3=nexttile;
imagesc(r.cMask)
title('Categorical mask')
axis image
t4=nexttile;
imagesc(r.mapPI)
title('Vessel types')
axis image

cmap=ones(4,3);
cmap(2,:)=[1,0,0];
cmap(3,:)=[0,0,1];
cmap(4,:)=[0,1,0];
colormap(t1,"parula");
colormap(t2,"parula");
colormap(t3,"parula");
colormap(t4,cmap);

%Using cMask you can also plot average BFI and PI distributions depending
%on vessel type
figure
tiledlayout(2,3,'TileSpacing','tight','Padding','tight')
nexttile
histogram(r.imgPI(r.cMask(:)==1)) 
title('Parenchymal pulsatility')
nexttile
histogram(r.imgBFI(r.cMask(:)==1)) 
title('Parenchymal BFI')
nexttile
histogram(r.imgPI(r.cMask(:)>3 & r.mapPI(:)==1)) 
title('Arterial pulsatility')
nexttile
histogram(r.imgBFI(r.cMask(:)>3 & r.mapPI(:)==1)) 
title('Arterial BFI')
nexttile
histogram(r.imgPI(r.cMask(:)>3 & r.mapPI(:)==2)) 
title('Venous pulsatility')
nexttile
histogram(r.imgBFI(r.cMask(:)>3 & r.mapPI(:)==2)) 
title('Venous BFI')


%Or you can use vsMetrics and vsData to plot features of vascular segments
%identified without using dynamic diameter segmetnation.
vClass=r.vsMetrics{:,"Type"}; %get the vessel types

figure
tiledlayout(2,3,'TileSpacing','tight','Padding','tight')
for i=[2,3,4,9,12]
    nexttile
    vName=r.vsMetrics.Properties.VariableNames(i);
    vData=r.vsMetrics{:,i};
    
    boxplot(vData(vClass==1 | vClass==2),vClass(vClass==1 | vClass==2)) %Only choose the vessels identified as arteries (1) or veins (2)
    title(vName)
end

figure
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')
nexttile
plot(r.time,r.vsData(:,vClass==1)) % all identified arteries
title('Arteries BFI over cycle')
xlabel('Time, s')
ylabel('BFI')
grid on
axis tight

nexttile
plot(r.time,r.vsData(:,vClass==2)) % all identified veins
title('Veins BFI over cycle')
xlabel('Time, s')
ylabel('BFI')
grid on
axis tight


%You can do the same for dynamic vessels segmentation if any
vClass=r.dvsMetrics{:,"Type"}; %get the vessel types
figure
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
nexttile
plot(r.time,r.dvsData(:,vClass==1)) % all identified arteries
title('Arteries BFI over cycle')
xlabel('Time, s')
ylabel('BFI')
grid on
axis tight

nexttile
plot(r.time,r.dvsData(:,vClass==2)) % all identified veins
title('Veins BFI over cycle')
xlabel('Time, s')
ylabel('BFI')
grid on
axis tight

nexttile
plot(r.time,r.dvsDiameter(:,vClass==1)) % all identified arteries
title('Arteries diameter over cycle')
xlabel('Time, s')
ylabel('Diameter, px')
grid on
axis tight

nexttile
plot(r.time,r.dvsDiameter(:,vClass==2)) % all identified veins
title('Veins diameter over cycle')
xlabel('Time, s')
ylabel('Diameter, px')
grid on
axis tight

%%Make a simple animation of flow dynamics
%First get color limits as percentiles over the entire data set but limited to
%the mask
data=reshape(r.data,[numel(r.cMask),size(r.data,3)]);
colorLim=prctile(data(r.cMask(:)>0),[5,99]);

%re-use the data array for normlised blood flow dynamics visualisation
data=(r.data-prctile(r.data,1,3))./(prctile(r.data,99,3)-prctile(r.data,1,3));
data(data<0)=0;
data(data>1)=1;

figure
t=tiledlayout(1,2)
for i=1:1:size(r.data,3)
    nexttile(1);
    imagesc(r.data(:,:,i))
    clim(colorLim)
    axis image
    title('BFI')
    nexttile(2);
    imagesc(data(:,:,i))
    clim([0,1])
    axis image
    title('Normalised BFI')   

    nexttile(1);

    title(t,[num2str(round(r.time(i)*1000)),' ms'])
    pause(0.1)
end

%%Interrogate a region in the dataset (we re-use the normalised data
%%calculated before)

colors=turbo(10);
f=figure;
imagesc(r.imgBFI)
clim(prctile(r.imgBFI(r.cMask(:)>0),[5,99]));
ROIs=getAndShowROIS(f,colors,[],[]);
data=reshape(r.data,[numel(r.cMask),size(r.data,3)]); %reshape the data so it would be easier to get ROI pixels from it
f=figure;
tiledlayout(1,2,"TileSpacing","tight","Padding","tight")
nexttile
imagesc(r.imgBFI)
axis image
clim(prctile(r.imgBFI(r.cMask(:)>0),[5,99]));
getAndShowROIS(f,colors,[],ROIs);
nexttile
hold on
for i=1:size(ROIs,3)
    mask=ROIs(:,:,i).*(r.cMask>0); % you could also set it to only vessels or only arteries or whatever following examples above
    ts=mean(data(mask(:)==1,:),1);
    plot(r.time,ts,'Color',colors(i,:))
end
title('ROI plots')
hold off
ylabel('BFI')
xlabel('Time, s')
grid on



%% Load variables of interest or the entire file
for fidx1=size(fNames,1)
    for fidx2=size(fNames,2)
        for fidx3=size(fNames,3)
            load(fNames{fidx1,fidx2,fidx3})
            
        end
    end
end
