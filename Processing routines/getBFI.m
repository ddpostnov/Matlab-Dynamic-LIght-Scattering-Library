%getBFI  Convert contrast cubes (*. _K_d.mat) to blood-flow index (BFI)
%
%   getBFI(s,fNames) scans every *_K_d.mat file in fNames, replaces every
%   field whose name contains the substring ‘data’ with its blood-flow
%   index, updates the corresponding metrics, and writes three new MAT-files
%   per dataset:
%
%       *_BFI_d.mat   SOURCE   – all ‘data’ variables converted to BFI
%       *_BFI_r.mat   RESULTS  – BFI images, <std(BFI)>, metrics updated
%       *_BFI_s.mat   SETTINGS – original parameter struct plus the
%                                sub-field settings.calculateBFI
%
%   By default (s.deleteOriginal == true) the original *_K_d / *_K_r /
%   *_K_s triplet is deleted after successful conversion.
%
%   INPUTS
%     s        parameter structure  
%                • deleteOriginal   true / false  
%                • method           currently only "basic" (=1/K²)
%     fNames   cell array of full paths to *_K_d.mat files.
%
%   OUTPUTS
%     None – the routine operates by overwriting / writing files on disk.
%
%   EXAMPLE
%     p.deleteOriginal = false;
%     p.method         = "basic";
%     D = dir(fullfile(dataRoot,'*_K_d.mat'));
%     getBFI(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     No external code beyond the LSCI processing library.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 05-Aug-2025
%
%   Note: This header was generated with ChatGPT and may contain minor
%   inconsistencies—please verify before distribution.
%   ----------------------------------------------------------------------


% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
% s.deleteOriginal=true; %true or false
% %ADJUSTED IF NECESSARY - CONVERSION METHOD
% s.method="basic"; %only "basic" is avaliable

function getBFI(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

for fidx=1:1:numel(fNames)
    tic
    disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    s.fName=fNames{fidx};
    clearvars results source settings
    load(s.fName,'source')
    load(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings');
    load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

    fn = fieldnames(source);
    for k=1:numel(fn)
        if contains(fn{k}, 'data', 'IgnoreCase', true )
            source.(fn{k})=calculateBFI(source.(fn{k}),s.method);
            if strcmp(fn{k},'data')
                results.imgBFI=mean(source.data,3,'omitnan');
                results.extendedMetrics.imgStdBFI=std(source.data,0,3,'omitnan');

            end
            disp(['Variable source.',(fn{k}),' has been converted'])
        end
    end

    fn = fieldnames(results);
    for k=1:numel(fn)
        if contains(fn{k}, 'data', 'IgnoreCase', true )
            results.(fn{k})=calculateBFI(results.(fn{k}),s.method);
            if strcmp(fn{k},'sData')
                results.sMetrics.('BFI')=squeeze(mean(results.sData,1,'omitnan'))';
                results.sMetrics.('std(BFI)')=squeeze(std(results.sData,0,1,'omitnan'))';
            elseif strcmp(fn{k},'dvsData')
                results.dvsMetrics.('BFI')=mean(results.dvsData,1,'omitnan')';
                results.dvsMetrics.('std(BFI)')=std(results.dvsData,0,1,'omitnan')';
            end
            disp(['Variable results.',(fn{k}),' has been converted'])
        end
    end

    disp(['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))])
    settings.calculateBFI=s;
    save(strrep(fNames{fidx},'_K_d.mat','_BFI_d.mat'),'source','-v7.3');
    save(strrep(fNames{fidx},'_K_d.mat','_BFI_r.mat'),'results','-v7.3');
    save(strrep(fNames{fidx},'_K_d.mat','_BFI_s.mat'),'settings','-v7.3');

    if s.deleteOriginal
        delete(fNames{fidx});
        delete(strrep(fNames{fidx},'_d.mat','_s.mat'));
        delete(strrep(fNames{fidx},'_d.mat','_r.mat'));
    end
end

    function data=calculateBFI(data,method)
        switch method
            case "basic"
                data=(1./data.^2);
            otherwise
                error("Method not recognised")
        end
    end
end