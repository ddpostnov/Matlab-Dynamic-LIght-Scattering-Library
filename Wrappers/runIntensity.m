%runIntensity - Read .cxd recordings, compute the mean-intensity image and save.
%
%   For each .cxd file, reads the frames (Bio-Formats), optionally averages every
%   decimFactor frames, and stores the intensity stack, its time-mean image and a
%   diagnostic JPG.  Results are written next to each file as _I_d/_r/_s.mat.
%
% Syntax:
%    runIntensity(s, fNames)
%
% Inputs:
%    s      - parameter struct:
%               .decimFactor  frames averaged per output frame (1 = none).
%               .dataTypeOut  output class, e.g. 'single' ([] = source type).
%               .saveSource   logical; also save the full intensity stack.
%    fNames - cell array of .cxd file paths.
%
% Outputs:
%    (none) - writes <file>_I_d.mat (source), _I_r.mat (results), _I_s.mat
%             (settings) and <file>_I.jpg next to each input.
%
% Dependencies: the Bio-Formats MATLAB library (bfGetReader, bfGetPlane).
% See also: runContrastFromRLS, readCXD, readDV
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% s.decimFactor=1;
% s.dataTypeOut='single';
% s.saveSource=true;

function runIntensity(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.cxd'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".cxd".');
end

for fidx=1:1:length(fNames)
    if ~isempty(fNames{fidx})
        %set file name to load data
        s.fName=char(fNames{fidx});
        disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        clearvars results source settings
        dF=s.decimFactor;

        reader = bfGetReader(s.fName);
        omeMeta = reader.getMetadataStore();
        sizeX=double(omeMeta.getPixelsSizeX(0).getValue());
        sizeY=double(omeMeta.getPixelsSizeY(0).getValue());
        sizeTRaw=double(omeMeta.getPlaneCount(0));
        sizeT=floor(sizeTRaw/dF);
        startT=omeMeta.getImageAcquisitionDate(0);
        timeStamp = round(posixtime(datetime(string(startT), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss', 'TimeZone', 'UTC')))*1000;
        dataType=char(omeMeta.getPixelsType(0).toString());
        sampling=double(omeMeta.getPixelsTimeIncrement(0).value()); %s
        time=(0:1:(sizeT-1))*sampling*dF;
        time=time(:);

        if isempty(s.dataTypeOut)
            data=zeros(sizeY,sizeX,sizeT,dataType);
        else
            data=zeros(sizeY,sizeX,sizeT,s.dataTypeOut);
        end

        frame=zeros(sizeY,sizeX,dF,dataType);

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
        settings.runIntensity=s;

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
end