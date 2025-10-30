%Pipeline_registerLSCItoFluorescence - pipeline for LSCI registration to Fluorescences mask/angiogram

% Copyright © : Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 27-Feb-2025

%------------- BEGIN CODE --------------
%% Pipeline settingsuration
%Clear the workspace to avoid potential conflicts.
clear
close all
%Add all LSCI library folders to the matlab path
lsciLibraryFolder = 'C:\Dropbox\Work\GitHub\DDPLab-private';
addpath(genpath(lsciLibraryFolder));


%Set file-names - option 1 - automatic lookup within the root folder
% reCalc=1; %set to 1 if the existing processed file should be rewritten
% rootFolder='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer';
% files=dir([rootFolder,'\**\*_cycle_sLSCIMM.mat']); %configure the file name that contains LSCI data
% count=1;
% fileNames={};
% for fidx=1:1:size(files,1)
%     fileName=[files(fidx).folder,'\',files(fidx).name];
%     fileNames{count,1}=fileName;
%     fileNames{count,2}=strrep(strrep(strrep(fileName,'_cycle_sLSCIMM','_mask'),'BP','BB'),'LSCI','EPFL');
%     if ~exist(fileNames{count,2},"file")
%         error("Mask file not found, execution aborted")
%     end
%     count=count+1;
%end

%Path to rls files (change indexes correspondingly to the number of files)
fileNames{1,1}="O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\awake\r2pm\LSCI_25112024_1PTFM02BP_cycle_sLSCIMM.mat";
fileNames{2,1}="O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\awake\r3pm\LSCI_25112024_1PTFM03BP_cycle_sLSCIMM.mat";
fileNames{3,1}="O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\awake\r5pf\LSCI_25112024_1PTFF05BP_cycle_sLSCIMM.mat";
fileNames{4,1}="O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\iso\iso\r2pm\LSCI_12122024_1PTFM02BP_cycle_sLSCIMM.mat";
fileNames{5,1}="O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\iso\iso\r3pm\LSCI_12122024_1PTFM03BP_cycle_sLSCIMM.mat";
fileNames{6,1}="O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\iso\iso\r5pf\LSCI_12122024_1PTFF05BP_cycle_sLSCIMM.mat";

fileNames{1,2}='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\awake\r2pm\EPFL_25112024_1PTFM02BB_mask.mat';
fileNames{2,2}='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\awake\r3pm\EPFL_25112024_1PTFM03BB1_mask.mat';
fileNames{3,2}='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\awake\r5pf\EPFL_25112024_1PTFF05BB1_mask.mat';
fileNames{4,2}='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\iso\iso\r2pm\EPFL_12122024_1PTFM02BB1_mask.mat';
fileNames{5,2}='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\iso\iso\r3pm\EPFL_12122024_1PTFM03BB1_mask.mat';
fileNames{6,2}='O:\HE_BFI-DATA\Rasmus\pulsetransfer\pulsetransfer\iso\iso\r5pf\EPFL_12122024_1PTFF05BB1_mask.mat';

for fidx=1:1:size(fileNames,1)
    disp(['File ',num2str(fidx),' out of ',num2str(size(fileNames,1))])
    fileName=char(fileNames{fidx,2});
    load(fileName);
    maskRef=results.mask;
    refImg=results.imgMean;    
    fileName=char(fileNames{fidx,1});
    load(fileName)
    img=1./results.imgK;


    tform=manualByPointRegistration(refImg,img,'sideBySide');
    ov=affineOutputView(size(refImg),tform,"BoundsStyle","SameAsInput");
    fn = fieldnames(results);
    for k=1:numel(fn)
        results.(fn{k})=results.(fn{k});

        if size(results.(fn{k}),1)==size(img,1) & size(results.(fn{k}),2)==size(img,2)
            tmp2=zeros(size(refImg,1),size(refImg,2),size(results.(fn{k}),3),class(results.(fn{k})));
            for i=1:1:size(results.(fn{k}),3)
                tmp2(:,:,i)=imwarp(results.(fn{k})(:,:,i),tform,"OutputView",affineOutputView(size(refImg),tform,"BoundsStyle","SameAsInput"));
            end
            results.(fn{k})=tmp2;

        end
    end

    results.mask=logical(results.mask.*(maskRef>0));
    results.maskType=results.mask.*maskRef;
    settings.registration.tform=tform;
    settings.registration.refImg=refImg;
    settings.registration.ov=ov;

    save(fileName,'results','settings','-v7.3');
end