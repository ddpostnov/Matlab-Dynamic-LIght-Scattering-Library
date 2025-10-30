

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;

function getGuidedContrast(s,fNames,fNamesRaw)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".mat".');
end

for fidx=1:1:length(fNames)
    %set file name to load data
    s.fName=char(fNames{fidx});
    s.fNameRaw=char(fNamesRaw{fidx});
    clearvars results source settings
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');

    sMap=results.sMap;
    pixelIndices = find(sMap > 0);
    objectLabelsForPixels = double(sMap(pixelIndices));
    uniqueLabels = unique(objectLabelsForPixels);
    numObjects = numel(uniqueLabels);
    pixelsPerObject = accumarray(objectLabelsForPixels, 1, [max(uniqueLabels) 1]);
    pixelsPerObject = pixelsPerObject(uniqueLabels)'; 


    [fPointer,cfg]=getPointerRLS(s.fNameRaw);
    timeStamps=zeros(cfg.sizeT,1,'int64');
    results.sDataHD = zeros(cfg.sizeT, numObjects,'double');

    tic
    for i=1:1:cfg.sizeT
        if mod(i,100)==0; disp(['Processed frame ',num2str(i),' out of ',num2str(cfg.sizeT),'. Time elapsed ',num2str(toc)]); end
        timeStamps(i)=fread(fPointer,1,'*uint64');
        frame=fread(fPointer,[cfg.sizeY,cfg.sizeX],['*',cfg.dataType]);
        frame=double(frame(pixelIndices));
        sumX = accumarray(objectLabelsForPixels, frame, [], @sum, 0);
        sumX2 = accumarray(objectLabelsForPixels, frame.^2, [], @sum, 0);
        meanVal = sumX(uniqueLabels)' ./ pixelsPerObject;
        varVal = (sumX2(uniqueLabels)' ./ pixelsPerObject - meanVal.^2) .* (pixelsPerObject ./ (pixelsPerObject - 1));
        stdVal = sqrt(varVal);
        results.sDataHD(i, :) = single(stdVal ./ meanVal);
    end
    fclose(fPointer);
    cfg.fName=s.fName;
    cfg.fNameRaw=s.fNameRaw;
    cfg.libraryFolder=s.libraryFolder;
    results.timeHD=(timeStamps-timeStamps(1)).*1000; %conversion to seconds.

    % Save the settings and results
    settings.guidedContrast=cfg;
     disp(['Saving the results. Elapsed time ',num2str(round(toc)),'s']);
    save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
    save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
    disp('Saving complete');
end
end