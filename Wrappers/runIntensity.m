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
%    Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%    s.cancelFcn()->tf.
%
% Outputs:
%    (none) - writes <file>_I_d.mat (source), _I_r.mat (results), _I_s.mat
%             (settings) and <file>_rep_intensity.jpg next to each input.
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

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Intensity',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        %set file name to load data
        s.fName=char(fNames{fidx});
        reportFile(rep,fidx,s.fName);
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
        settings.runIntensity=reportSettings(s);

        fh=reportFigure(rep,'intensity');
        tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
        ax=nexttile(tl);
        imagesc(ax,results.imgI)
        clim(ax,[prctile(results.imgI(:),1),prctile(results.imgI(:),99)])
        colorbar(ax)
        axis(ax,'image')
        ax=nexttile(tl);
        semilogy(ax,time,squeeze(mean(data,[1,2])));
        hold(ax,'on')
        semilogy(ax,time,squeeze(min(data,[],[1,2])));
        semilogy(ax,time,squeeze(max(data,[],[1,2])));
        hold(ax,'off')
        legend(ax,{'Mean','Min','Max'})
        ylabel(ax,'Intensity')
        xlabel(ax,'Time, s')
        xlim(ax,[time(1),time(end)]);
        ylim(ax,[min(data(:)),max(data(:))])
        reportSave(rep,fh,'intensity');     % reportSave titles the page

        % Save the settings and results
        reportWriting(rep);
        if s.saveSource
            source.data=data;
            source.time=time;
            save(strrep(s.fName,'.cxd','_I_d.mat'),'source','-v7.3','-nocompression');
        end
        save(strrep(s.fName,'.cxd','_I_r.mat'),'results','-v7.3','-nocompression');
        save(strrep(s.fName,'.cxd','_I_s.mat'),'settings','-v7.3','-nocompression');
        reportSaved(rep);
    end
end
reportClose(rep);

end