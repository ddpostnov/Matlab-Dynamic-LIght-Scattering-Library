%runBFI  Convert contrast cubes (*_K_r.mat) to blood-flow index (BFI)
%
%   runBFI(s,fNames) scans every *_K_r.mat file in fNames, replaces every
%   field whose name contains the substring ‘data’ with its blood-flow
%   index, updates the corresponding metrics, and writes three new MAT-files
%   per dataset:
%
%       *_BFI_d.mat   SOURCE   – all ‘data’ variables converted to BFI
%       *_BFI_r.mat   RESULTS  – BFI images, <std(BFI)>, metrics updated
%       *_BFI_s.mat   SETTINGS – original parameter struct plus the
%                                sub-field settings.calculateBFI
%
%   THE SWEEP TAKES results.gsData WITH THE REST, and that is correct: this step
%   runs on a *_K_ product, whose guided traces are contrast.  The product name
%   is the whole answer to what a field holds - a *_BFI_* file is blood-flow
%   index throughout, a *_K_* file is contrast throughout.
%
%   By default (s.deleteOriginal == true) the original *_K_d / *_K_r /
%   *_K_s triplet is deleted after successful conversion.
%
%   INPUTS
%     s        parameter structure  
%                • deleteOriginal   true / false  
%                • method           currently only "basic" (=1/K²)
%     fNames   cell array of full paths to *_K_r.mat files - the RESULTS member
%              of each contrast product.  The SOURCE cube and the SETTINGS are
%              named from it (getProductPath).
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUTS
%     None – the routine operates by overwriting / writing files on disk.
%
%   EXAMPLE
%     p.deleteOriginal = false;
%     p.method         = "basic";
%     D = dir(fullfile(dataRoot,'*_K_r.mat'));
%     runBFI(p, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     No external code beyond the LSCI processing library.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026


% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED IF NECESSARY - DELETE ORIGINAL FILES
% s.deleteOriginal=true; %true or false
% %ADJUSTED IF NECESSARY - CONVERSION METHOD
% s.method="basic"; %only "basic" is avaliable

function runBFI(s,fNames)
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_K_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_K_r.mat'').']);
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'BFI',fNames);

for fidx=1:1:numel(fNames)
     if reportCancelled(rep), break; end        % cooperative cancel between files
     if ~isempty(fNames{fidx})
    s.fName=fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results source settings
    load(getProductPath(s.fName,'d'),'source')
    load(getProductPath(s.fName,'s'),'settings');
    load(s.fName,'results');

    fn = fieldnames(source);
    for k=1:numel(fn)
        if contains(fn{k}, 'data', 'IgnoreCase', true )
            source.(fn{k})=getBFI(source.(fn{k}),s.method);
            if strcmp(fn{k},'data')
                results.imgBFI=mean(source.data,3,'omitnan');
                results.extendedMetrics.imgStdBFI=std(source.data,0,3,'omitnan');

            end
        end
    end

    fn = fieldnames(results);
    for k=1:numel(fn)
        if contains(fn{k}, 'data', 'IgnoreCase', true )
            results.(fn{k})=getBFI(results.(fn{k}),s.method);
            if strcmp(fn{k},'sData')
                results.sMetrics.('BFI')=squeeze(mean(results.sData,1,'omitnan'))';
                results.sMetrics.('std(BFI)')=squeeze(std(results.sData,0,1,'omitnan'))';
            elseif strcmp(fn{k},'dvsData')
                results.dvsMetrics.('BFI')=mean(results.dvsData,1,'omitnan')';
                results.dvsMetrics.('std(BFI)')=std(results.dvsData,0,1,'omitnan')';
            end
        end
    end

    reportWriting(rep);
    settings.calculateBFI=reportSettings(s);
    % the BFI triplet SUBSTITUTES the product token of the input's own name, so
    % the branch flag it carries ('_t', '_c', ...) is preserved by construction
    bfiName=strrep(s.fName,'_K_r.mat','_BFI_r.mat');
    save(getProductPath(bfiName,'d'),'source','-v7.3','-nocompression');
    save(bfiName,'results','-v7.3','-nocompression');
    save(getProductPath(bfiName,'s'),'settings','-v7.3','-nocompression');

    if s.deleteOriginal
        delete(getProductPath(s.fName,'d'));
        delete(getProductPath(s.fName,'s'));
        delete(s.fName);
    end
    reportSaved(rep);
     end
end
reportClose(rep);

end