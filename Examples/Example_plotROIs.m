libraryFolder = 'C:\Dropbox\Work\GitHub\DDPLab-private\Dynamic LIght Scattering Library v2.0';
addpath(genpath(libraryFolder));

fName='C:\Dropbox\Work\Data\Mia\LSCI_20250317_1WTTM02BN_t_BFI_d.mat';
load(fName,'source')
load(strrep(fName,'_d.mat','_r.mat'),'results');

img=mean(source.data,3);
if isfield(results,"cMask")
    mask=results.cMask>0;
else
    mask=~isnan(img) & ~isinf(img) & img>0;
end

h=figure;
imagesc(img)
clim(prctile(img(mask(:)),[5,99]));
axis image;
colors=lines(999);
ROIs=getAndShowROIS(h,colors,[],[]);
ts=zeros(numel(source.time),size(ROIs,3));
data=reshape(source.data,[],numel(source.time));
for i=1:1:size(ROIs,3)
roi=squeeze(ROIs(:,:,i))==1;
    ts(:,i)=squeeze(mean(data(roi(:),:),'omitnan'));
end
h=figure;
tiledlayout(2,4,"TileSpacing","tight","Padding","tight");
nexttile([2,2]);
imagesc(img);
clim(prctile(img(mask(:)),[5,99]));
axis image;
getAndShowROIS(h,colors,[],ROIs);
nexttile([1,2]);
plot(source.time, ts);
xlabel('Time, s');
ylabel('BFI');
axis tight;

nexttile([1,2]);
plot(source.time, ts./mean(ts,1));
xlabel('Time, s');
ylabel('Mean normalised BFI');
axis tight;

set(gcf,'Color','w');





