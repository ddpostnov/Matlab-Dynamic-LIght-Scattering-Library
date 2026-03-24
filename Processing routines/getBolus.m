%getBolus
% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% s.fBolus=[1,2000];
% s.fAngio=[2001,30000];

function getBolus(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.cxd'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".cxd".');
end

for fidx=1:1:length(fNames)
    %set file name to load data
    s.fName=char(fNames{fidx});
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    clearvars results source settings

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
    if isempty(s.fBolus)
        data=zeros(sizeX,sizeY,sizeT,dataType);
    else
        data=zeros(sizeX,sizeY,s.fBolus(2)-s.fBolus(1)+1,dataType);
    end
    
    imgI=zeros(sizeX,sizeY,'double');

    if isempty(s.fBolus)
        for t=1:1:sizeT
            data(:,:,t) = bfGetPlane(reader, t);
        end
        h=figure;
        h.WindowState='Maximize';
        semilogy(squeeze(mean(data,[1,2])));
        hold on
        semilogy(squeeze(min(data,[],[1,2])));
        semilogy(squeeze(max(data,[],[1,2])));
        hold off
        legend({'Mean','Min','Max'})
        ylabel('Intensity')
        xlabel('Frames')
        xlim([1,size(data,3)]);
        ylim([min(data(:)),max(data(:))])
        s.fBolus=input("Enter the desired bolus frames span as [f1,f2] vector:");
        s.fAngio=input("Enter the desired angiogram frames span as [f1,f2] vector:");
        if isempty(s.fBolus)
            s.fBolus=[1,size(data,3)];
        end
        if isempty(s.fAngio)
            s.fAngio=[1,size(data,3)];
        end
        imgI=squeeze(mean(data(:,:,s.fAngio(1):s.fAngio(2)),3));
        data=data(:,:,s.fBolus(1):s.fBolus(2));

    else
        for t=s.fBolus(1):s.fBolus(2)
            data(:,:,t-s.fBolus(1)+1) = bfGetPlane(reader, t);
        end

        if isempty(s.fAngio)
            for t=(s.fBolus(2)+1):sizeT
                imgI=imgI+ double(bfGetPlane(reader, t))./(sizeT-s.fBolus(2)+1);
            end
        else
            for t=s.fAngio(1):s.fAngio(2)
                imgI=imgI+ double(bfGetPlane(reader, t))./(s.fAngio(2)-s.fAngio(1)+1);
            end
        end
    end
    reader.close();
    time=time(s.fBolus(1):s.fBolus(2));
    timeStamp=timeStamp+time(1);
    time=time-time(1);

    results.time=time;
    results.timeStamp=timeStamp;
    results.imgI=imgI;
    settings.getBolus=s;

    h=figure;
    h.WindowState='Maximize';
    subplot(1,2,1)
    imagesc(results.imgI)
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
    print(h,strrep(s.fName,'.cxd','_b_I.jpg'), '-djpeg', '-r300');

    % Save the settings and results
    disp('Saving the results');
    source.data=data;
    source.time=time;
    save(strrep(s.fName,'.cxd','_b_I_d.mat'),'source','-v7.3');
    save(strrep(s.fName,'.cxd','_b_I_r.mat'),'results','-v7.3');
    save(strrep(s.fName,'.cxd','_b_I_s.mat'),'settings','-v7.3');
end
end