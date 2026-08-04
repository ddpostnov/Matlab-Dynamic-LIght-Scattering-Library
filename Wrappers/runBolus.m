%runBolus  Split a *.cxd bolus recording into a bolus cube and an angiogram
%
%   runBolus(s,fNames) reads each Bio-Formats *.cxd acquisition, then divides
%   the recording into two sections: a "bolus" span that is kept at full
%   frame rate, and an "angiogram" span that is averaged into a single
%   mean-intensity image.  When the frame spans are left empty the function
%   plots the intensity trace and lets the user mark each span by dragging a
%   box on it; otherwise it uses the pre-set s.fBolus / s.fAngio values.
%   Three MAT-files and a JPEG
%   preview are written per recording.
%
%   INPUTS
%     s        parameter structure
%              ├─ fBolus  [f1,f2] frame span kept at full frame rate.
%              │           Leave empty ([]) to mark it with a box on the plot.
%              └─ fAngio  [f1,f2] frame span averaged into the angiogram.
%                          Leave empty ([]) to mark it with a box on the plot,
%                          or (when fBolus is preset) to average every frame
%                          after the bolus span.
%     fNames   cell array of char vectors / strings with full paths to
%              *.cxd files produced by the acquisition software.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUT SIDE-EFFECTS (per file)
%     <name>_b_I_d.mat   SOURCE  – bolus cube (source.data) + time (source.time, s)
%     <name>_b_I_r.mat   RESULTS – angiogram (imgI), time, timeStamp
%     <name>_b_I_s.mat   SETTINGS – copy of the parameter structure s
%     <name>_rep_bolus.jpg  report page: angiogram + intensity trace
%
%   EXAMPLE
%     s.libraryFolder = libraryFolder;
%     s.fBolus = [1,2000];
%     s.fAngio = [2001,30000];
%     D = dir(fullfile(dataRoot,'*.cxd'));
%     runBolus(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     Bio-Formats (bfGetReader, bfGetPlane), plus core functions in the
%     LSCI processing library.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% s.fBolus=[1,2000];
% s.fAngio=[2001,30000];

function runBolus(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.cxd'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".cxd".');
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.  This step can
% block on interactive span selection, so cancel is only checked between files.
rep=reportOpen(s,'Bolus',fNames);

for fidx=1:1:numel(fNames)
     if reportCancelled(rep), break; end        % cooperative cancel between files
     if ~isempty(fNames{fidx})
    % The blanket close-every-figure call that used to start each iteration is gone.
    % It existed only to mop up the report figure this wrapper leaked, and it closed
    % UIFIGURES too - it could take the Processing Workbench window down mid-run.  The
    % span picker below is closed explicitly and reportSave deletes the report page on
    % every path, so there is nothing left to sweep up.
    s.fName=char(fNames{fidx});
    reportFile(rep,fidx,s.fName);
    clearvars results source settings

    %read the file meta data
    reader    = bfGetReader(s.fName);
    omeMeta   = reader.getMetadataStore();
    sizeX     = double(omeMeta.getPixelsSizeX(0).getValue());
    sizeY     = double(omeMeta.getPixelsSizeY(0).getValue());
    sizeT     = double(omeMeta.getPlaneCount(0));
    startT    = omeMeta.getImageAcquisitionDate(0);
    timeStamp = round(posixtime(datetime(string(startT),'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone','UTC')))*1000; %ms
    dataType  = char(omeMeta.getPixelsType(0).toString());
    if strcmp(dataType,'float'); dataType='single'; end %OME reports 'float'; MATLAB needs 'single'
    sampling  = double(omeMeta.getPixelsTimeIncrement(0).value()); %s
    time      = ((0:1:(sizeT-1))*sampling)';

    imgI=zeros(sizeY,sizeX,'double');

    if isempty(s.fBolus)
        %load the whole recording, then ask for the spans
        data=zeros(sizeY,sizeX,sizeT,dataType);
        for t=1:1:sizeT
            data(:,:,t)=bfGetPlane(reader,t);
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
        ylim([double(min(data(:))),double(max(data(:)))])

        %mark the bolus span by dragging a box over it (only the horizontal
        %extent is used); adjust, then double-click the box to confirm
        title('Mark the BOLUS frames: drag a box over them, adjust, then double-click to confirm');
        rBolus=drawrectangle(gca,'Color','r','Label','Bolus','FaceAlpha',0.15);
        wait(rBolus);
        if isvalid(rBolus) && ~isempty(rBolus.Position)
            s.fBolus=round([rBolus.Position(1),rBolus.Position(1)+rBolus.Position(3)]);
        else
            s.fBolus=[1,size(data,3)];
        end

        %mark the angiogram span the same way
        title('Mark the ANGIOGRAM frames: drag a box over them, adjust, then double-click to confirm');
        rAngio=drawrectangle(gca,'Color','g','Label','Angio','FaceAlpha',0.15);
        wait(rAngio);
        if isvalid(rAngio) && ~isempty(rAngio.Position)
            s.fAngio=round([rAngio.Position(1),rAngio.Position(1)+rAngio.Position(3)]);
        else
            s.fAngio=[1,size(data,3)];
        end

        %clamp the selected spans to the available frame range
        s.fBolus=min(max(s.fBolus,1),size(data,3));
        s.fAngio=min(max(s.fAngio,1),size(data,3));
        close(h);

        imgI=squeeze(mean(data(:,:,s.fAngio(1):s.fAngio(2)),3));
        data=data(:,:,s.fBolus(1):s.fBolus(2));
    else
        %validate the pre-set bolus span against the available frame count
        if s.fBolus(1)>sizeT || s.fBolus(1)>s.fBolus(2)
            warning('s.fBolus=[%d,%d] is outside the %d available frames; skipping file.',s.fBolus(1),s.fBolus(2),sizeT);
            reader.close();
            continue;
        end
        if s.fBolus(2)>sizeT
            warning('s.fBolus(2)=%d exceeds the %d available frames; clamped.',s.fBolus(2),sizeT);
            s.fBolus(2)=sizeT;
        end
        s.fBolus(1)=max(s.fBolus(1),1);

        %load only the bolus span at full frame rate
        data=zeros(sizeY,sizeX,s.fBolus(2)-s.fBolus(1)+1,dataType);
        for t=s.fBolus(1):s.fBolus(2)
            data(:,:,t-s.fBolus(1)+1)=bfGetPlane(reader,t);
        end

        %average the angiogram span (or everything after the bolus)
        if isempty(s.fAngio)
            nAngio=sizeT-s.fBolus(2);
            if nAngio<1
                warning('No frames after the bolus span; angiogram left blank.');
            end
            for t=(s.fBolus(2)+1):sizeT
                imgI=imgI+double(bfGetPlane(reader,t))./nAngio;
            end
        else
            %validate the pre-set angiogram span against the available frame count
            if s.fAngio(1)>sizeT || s.fAngio(1)>s.fAngio(2)
                warning('s.fAngio=[%d,%d] is outside the %d available frames; skipping file.',s.fAngio(1),s.fAngio(2),sizeT);
                reader.close();
                continue;
            end
            if s.fAngio(2)>sizeT
                warning('s.fAngio(2)=%d exceeds the %d available frames; clamped.',s.fAngio(2),sizeT);
                s.fAngio(2)=sizeT;
            end
            s.fAngio(1)=max(s.fAngio(1),1);
            nAngio=s.fAngio(2)-s.fAngio(1)+1;
            for t=s.fAngio(1):s.fAngio(2)
                imgI=imgI+double(bfGetPlane(reader,t))./nAngio;
            end
        end
    end
    reader.close();

    %crop and rezero the time vector, shift the absolute time stamp
    time=time(s.fBolus(1):s.fBolus(2));
    timeStamp=timeStamp+time(1)*1000; %time is in s, timeStamp in ms
    time=time-time(1);

    results.time=time;
    results.timeStamp=timeStamp;
    results.imgI=imgI;
    settings.runBolus=reportSettings(s);

    fh=reportFigure(rep,'bolus');
    tl=tiledlayout(fh,1,2,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl);
    imagesc(ax,results.imgI)
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
    ylim(ax,[double(min(data(:))),double(max(data(:)))])
    reportSave(rep,fh,'bolus');             % reportSave titles the page

    %save the settings and results
    reportWriting(rep);
    source.data=data;
    source.time=time;
    %This is the ENTRY step, so this is where a raw path becomes a product path, and the
    %only place that has to know a project may keep its results apart from its
    %recordings.  s.fName is NOT reassigned - it stays the raw recording.  With no
    %results folder set the name comes back verbatim.
    outBase=getResultsPath(s.fName,s);
    save(strrep(outBase,'.cxd','_b_I_d.mat'),'source','-v7.3','-nocompression');
    save(strrep(outBase,'.cxd','_b_I_r.mat'),'results','-v7.3','-nocompression');
    save(strrep(outBase,'.cxd','_b_I_s.mat'),'settings','-v7.3','-nocompression');
    reportSaved(rep);
     end
end
reportClose(rep);

end