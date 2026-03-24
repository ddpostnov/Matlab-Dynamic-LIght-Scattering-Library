%getGuidedIntensity

function getGuidedIntensity(fNamesRef,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.cxd'), fNames(:)) & cellfun(@(s) isempty(s) || contains(s,'_r.mat'), fNamesRef(:)) )
    error('One or more *non-empty* entries do not contain ".cxd" or "_r.mat".');
end

for fidx=1:1:length(fNames)
    %set file name to load data
    s.fName=char(fNames{fidx});
    s.fNameRef=char(fNamesRef{fidx});
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    clearvars results source settings
    load(s.fNameRef,'results');

    if ~isfield(results, 'sMap') || ~isfield(results, 'cMask')
        error('Reference file is missing sMap or cMask variables')
    end

    reader = bfGetReader(s.fName);
    omeMeta = reader.getMetadataStore();
    sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
    sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
    sizeT=double(omeMeta.getPlaneCount(0)); 
    startT=omeMeta.getImageAcquisitionDate(0);
    timeStamp = round(posixtime(datetime(string(startT), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', 'TimeZone', 'UTC')))*1000;
    dataType=char(omeMeta.getPixelsType(0).toString());
    sampling=double(omeMeta.getPixelsTimeIncrement(0).value()); %s
    time=(0:1:(sizeT-1))*sampling;
    time=time(:);

    frame=zeros(sizeX,sizeY,dF,dataType);

    %read data and close the file
    for t=1:1:sizeT
        for j=1:1:dF
            frame(:,:,j) = bfGetPlane(reader, (t-1).*dF+j);
        end
        data(:,:,t)=mean(frame,3);
    end
    reader.close();

    results.time=time;
    results.timeStamp=timeStamp;
    results.imgI=mean(data,3);
    settings.getIntensity=s;

    h=figure;
    h.WindowState='Maximize';
    subplot(1,2,1)
    imagesc(results.imgI)
    clim([prctile(results.imgI(:),1),prctile(results.imgI(:),99)])
    colorbar
    axis image
    subplot(1,2,2)
    semilogy(time,squeeze(mean(data,[1,2])));
        hold on
        semilogy(time,squeeze(min(data,[],[1,2])));
        semilogy(time,squeeze(max(data,[],[1,2])));
        hold off
        legend({'Mean','Min','Max'})
        ylabel('Intensity')
        xlabel('Time, s')
        xlim([time(1),time(end)]);
        ylim([min(data(:)),max(data(:))])   
        fNameshort=split(s.fName,'\');
    fNameshort=fNameshort(end);
    sgtitle(strrep(fNameshort,'_',' '));
    drawnow
    print(h,strrep(s.fName,'.cxd','_I.jpg'), '-djpeg', '-r300');

    % Save the settings and results
    disp('Saving the results');
    if s.saveSource
        source.data=data;
        source.time=time;
        save(strrep(s.fName,'.cxd','_I_d.mat'),'source','-v7.3');
    end
    save(strrep(s.fName,'.cxd','_I_r.mat'),'results','-v7.3');
    save(strrep(s.fName,'.cxd','_I_s.mat'),'settings','-v7.3');
end
end