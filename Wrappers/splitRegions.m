%splitRegions  Crop LSCI result files into separate ROI-specific datasets
%
%   splitRegions(s,fNames) looks for a binary/label mask called
%   results.regionsMask in every *_r.mat file listed in fNames.  For each
%   non-zero region ID the function:
%       • crops every image-sized variable in RESULTS and SOURCE
%       • writes three new MAT-files whose names are prefixed with 'RoiN_'
%         (N = region index) and keep the original *_d / *_r / *_s suffixes
%       • copies the input parameter structure *s* to SETTINGS and stamps
%         the sub-field settings.splitRegions
%
%   Optionally (s.deleteOriginal == true) the original trio of files is
%   deleted after all ROIs are extracted.
%
%   INPUTS
%     s        parameter structure with at least
%                • deleteOriginal   logical flag (true/false)
%                • libraryFolder    path to toolbox root (not used here but
%                                   stored in SETTINGS for traceability)
%     fNames   cell array of full paths to *_r.mat files produced by the
%              LSCI pipelines.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUTS
%     None – function acts via side-effects:
%       RoiN_<name>_d.mat   SOURCE  (cropped)
%       RoiN_<name>_r.mat   RESULTS (cropped)
%       RoiN_<name>_s.mat   SETTINGS (updated)
%
%   EXAMPLE
%     p.deleteOriginal = false;
%     D = dir(fullfile(dataRoot,'*_K_r.mat'));
%     splitRegions(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     Only MATLAB built-ins and data structures produced by the LSCI
%     processing library.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026



% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
% s.deleteOriginal=true; %true or false

function splitRegions(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_K_r.mat'').']);
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Split regions',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end         % cooperative cancel between files
    if ~isempty(fNames{fidx})
        s.fName=fNames{fidx};
        reportFile(rep,fidx,s.fName);
        clearvars results source settings
        load(getProductPath(s.fName,'d'),'source')
        load(getProductPath(s.fName,'s'),'settings');
        load(s.fName,'results');

        resultsIni=results;
        sourceIni=source;

        if isfield(results,'regionsMask')
            regionsMask=resultsIni.regionsMask;
        elseif isfield(results,'cMask')
            regionsMask=resultsIni.cMask>0;
        else
            % No regions defined (setRegions skipped or no ROI drawn) and nothing
            % segmented yet: fall back to the whole window (all-true mask), so the file
            % is copied out as a single Roi1_ dataset instead of erroring.
            regionsMask=true(size(sourceIni.data,1),size(sourceIni.data,2));
        end



        % Every region is cropped and written out in the same pass, so the whole
        % loop below IS the writing of this recording's results.
        nRoi=double(max(regionsMask(:)));
        reportWriting(rep);
        for ridx=1:1:nRoi
            results=resultsIni;
            source=sourceIni;
            [y,x] = find(regionsMask==ridx);
            y=[min(y),max(y)];
            x=[min(x),max(x)];

            fn = fieldnames(results);
            for k=1:numel(fn)
                if size(results.(fn{k}),1)==size(regionsMask,1) & size(results.(fn{k}),2)==size(regionsMask,2)
                    tmp=zeros(y(2)-y(1)+1,x(2)-x(1)+1,size(results.(fn{k}),3),class(results.(fn{k})));
                    for i=1:1:size(results.(fn{k}),3)
                        tmp(:,:,i)=results.(fn{k})(y(1):y(2),x(1):x(2),i);
                    end
                    results.(fn{k})=tmp;
                end
            end

            fn = fieldnames(source);
            for k=1:numel(fn)
                if size(source.(fn{k}),1)==size(regionsMask,1) & size(source.(fn{k}),2)==size(regionsMask,2)
                    tmp=zeros(y(2)-y(1)+1,x(2)-x(1)+1,size(source.(fn{k}),3),class(source.(fn{k})));
                    for i=1:1:size(source.(fn{k}),3)
                        tmp(:,:,i)=source.(fn{k})(y(1):y(2),x(1):x(2),i);
                    end
                    source.(fn{k})=tmp;
                end
            end

            fn = fieldnames(settings);
            for k=1:numel(fn)
                if size(settings.(fn{k}),1)==size(regionsMask,1) & size(settings.(fn{k}),2)==size(regionsMask,2)
                    tmp=zeros(y(2)-y(1)+1,x(2)-x(1)+1,size(settings.(fn{k}),3),class(settings.(fn{k})));
                    for i=1:1:size(settings.(fn{k}),3)
                        tmp(:,:,i)=settings.(fn{k})(y(1):y(2),x(1):x(2),i);
                    end
                    settings.(fn{k})=tmp;
                end
            end

            settings.splitRegions=reportSettings(s);
            [path,name,extension]=fileparts(fNames{fidx});
            fName = fullfile(path,['Roi' num2str(ridx) '_' name extension]);
            save(getProductPath(fName,'d'),'source','-v7.3','-nocompression');
            save(fName,'results','-v7.3','-nocompression');
            save(getProductPath(fName,'s'),'settings','-v7.3','-nocompression');
        end
        reportSaved(rep);
        if s.deleteOriginal
            delete(getProductPath(s.fName,'d'));
            delete(getProductPath(s.fName,'s'));
            delete(s.fName);
        end
    end
end
reportClose(rep);

end