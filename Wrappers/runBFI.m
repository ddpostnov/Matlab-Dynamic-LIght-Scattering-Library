%runBFI  Convert contrast cubes (*. _K_d.mat) to blood-flow index (BFI)
%
%   runBFI(s,fNames) scans every *_K_d.mat file in fNames, replaces every
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
%     Optional workbench hooks in s (no-op when absent): s.progressFcn(frac,label),
%     s.stageFcn(stage,detail), s.cancelFcn()->tf.
%
%   OUTPUTS
%     None – the routine operates by overwriting / writing files on disk.
%
%   EXAMPLE
%     p.deleteOriginal = false;
%     p.method         = "basic";
%     D = dir(fullfile(dataRoot,'*_K_d.mat'));
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
if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

% Optional workbench hooks resolved to no-ops when absent (see header); s is never
% mutated and the hooks are stripped from the settings before saving.
[progressFcn,stageFcn,cancelFcn]=resolveHooks(s);

for fidx=1:1:numel(fNames)
     if cancelFcn(), break; end                 % cooperative cancel between files
     if ~isempty(fNames{fidx})
    tic
    msg=['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))];
    disp(msg); stageFcn('runBFI',msg);
    s.fName=fNames{fidx};
    clearvars results source settings
    load(s.fName,'source')
    load(strrep(fNames{fidx},'_d.mat','_s.mat'),'settings');
    load(strrep(fNames{fidx},'_d.mat','_r.mat'),'results');

    fn = fieldnames(source);
    for k=1:numel(fn)
        if contains(fn{k}, 'data', 'IgnoreCase', true )
            source.(fn{k})=getBFI(source.(fn{k}),s.method);
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
            results.(fn{k})=getBFI(results.(fn{k}),s.method);
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

    msgSave=['Saving file ',num2str(fidx),' out of ',num2str(numel(fNames))];
    disp(msgSave); stageFcn('runBFI',msgSave);
    settings.calculateBFI=stripHooks(s);
    save(strrep(fNames{fidx},'_K_d.mat','_BFI_d.mat'),'source','-v7.3');
    save(strrep(fNames{fidx},'_K_d.mat','_BFI_r.mat'),'results','-v7.3');
    save(strrep(fNames{fidx},'_K_d.mat','_BFI_s.mat'),'settings','-v7.3');

    if s.deleteOriginal
        delete(fNames{fidx});
        delete(strrep(fNames{fidx},'_d.mat','_s.mat'));
        delete(strrep(fNames{fidx},'_d.mat','_r.mat'));
    end
    progressFcn(fidx/numel(fNames),msg);        % coarse per-file progress
     end
end

end

% =====================================================================
function [progressFcn,stageFcn,cancelFcn]=resolveHooks(s)
%resolveHooks  Optional workbench callbacks, defaulted to no-ops when absent (progress/
%   stage take any args and do nothing; cancel returns false).  See the header.
progressFcn=@(varargin)[]; stageFcn=@(varargin)[]; cancelFcn=@()false;
if isfield(s,'progressFcn')&&~isempty(s.progressFcn), progressFcn=s.progressFcn; end
if isfield(s,'stageFcn')  &&~isempty(s.stageFcn),   stageFcn  =s.stageFcn;   end
if isfield(s,'cancelFcn') &&~isempty(s.cancelFcn),  cancelFcn =s.cancelFcn;  end
end

% =====================================================================
function s=stripHooks(s)
%stripHooks  Drop the transport callbacks before s is written to a settings file
%   (no-op when absent, so a hook-free call saves a byte-identical settings struct).
for h={'progressFcn','stageFcn','cancelFcn'}
    if isfield(s,h{1}), s=rmfield(s,h{1}); end
end
end