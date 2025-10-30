%splitRegions  Crop LSCI result files into separate ROI-specific datasets
%
%   splitRegions(s,fNames) looks for a binary/label mask called
%   results.regionsMask in every *_d.mat file listed in fNames.  For each
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
%     fNames   cell array of full paths to *_d.mat files produced by the
%              LSCI pipelines.
%
%   OUTPUTS
%     None – function acts via side-effects:
%       RoiN_<name>_d.mat   SOURCE  (cropped)
%       RoiN_<name>_r.mat   RESULTS (cropped)
%       RoiN_<name>_s.mat   SETTINGS (updated)
%
%   EXAMPLE
%     p.deleteOriginal = false;
%     D = dir(fullfile(dataRoot,'*_K_d.mat'));
%     splitRegions(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     Only MATLAB built-ins and data structures produced by the LSCI
%     processing library.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: Header generated with ChatGPT; minor inconsistencies may remain—
%   please verify before distribution.
%   ----------------------------------------------------------------------



% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - DELETE THE ORIGINAL FILES
% s.deleteOriginal=true; %true or false

function splitRegions(s,fNames)

if ~all( cellfun(@(s) isempty(s) || contains(s,'.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain ".mat".');
end

for fidx=1:1:numel(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    load(s.fName,'source')
    load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
    load(strrep(s.fName,'_d.mat','_r.mat'),'results');

    resultsIni=results;
    sourceIni=source;

    if ~isfield(results,'regionsMask')
        error('No regionsMask detected in the file.');
    end

    regionsMask=resultsIni.regionsMask;
    for ridx=1:1:max(regionsMask(:))
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

        disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames)),'. Region ',num2str(ridx),' out of ',num2str(max(regionsMask(:)))])

        settings.splitRegions=s;
        [path,name,extension]=fileparts(fNames{fidx});
        fName = fullfile(path,['Roi' num2str(ridx) '_' name extension]);
        save(fName,'source','-v7.3');
        save(strrep(fName,'_d.mat','_r.mat'),'results','-v7.3');
        save(strrep(fName,'_d.mat','_s.mat'),'settings','-v7.3');
    end
    if s.deleteOriginal
        delete(fNames{fidx});
        delete(strrep(fNames{fidx},'_d.mat','_s.mat'));
        delete(strrep(fNames{fidx},'_d.mat','_r.mat'));
    end
end
end