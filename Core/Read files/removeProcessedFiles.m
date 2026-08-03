%removeProcessedFiles - Drop files whose analysis step has already been done.
%
%   Inspects every entry of the cell array fNames (1-D or 2-D) and empties the
%   files that have already been through the requested processing step, so they
%   are not analysed twice.  For each file the token extStr in the name is
%   replaced by extReplStr to locate the matching SETTINGS *.mat (holding a
%   "settings" struct); if settings.(currentStep) exists the step is considered
%   done and that entry is set to [] (the array keeps its size and shape).
%
% Syntax:
%    fNames = removeProcessedFiles(fNames, extStr, extReplStr, currentStep)
%    fNames = removeProcessedFiles(fNames, extStr, extReplStr, currentStep, keepFirst)
%
% Inputs:
%    fNames      - cell array of file names (may be 2-D); empty entries are left
%                  in place.
%    extStr      - token/extension to replace, e.g. '.rls' or '_r.mat'.  It must
%                  really occur in the names being filtered: a raw '.rls' list
%                  filtered on '_r.mat' matches nothing and drops nobody.
%    extReplStr  - replacement giving the SETTINGS file, e.g. '_c_K_s.mat'.
%    currentStep - settings field name to test for, e.g. 'registration'.
%    keepFirst   - (optional, default false) when true the first column
%                  fNames(:,1) is always kept (the registration reference).
%
% Outputs:
%    fNames - the input list, same size/shape, with every already-processed
%             entry set to [] (test downstream with isempty).
%
% Example:
%    D = dir(fullfile(dataRoot,'*.rls'));
%    fNames = fullfile({D.folder}',{D.name}');
%    fNames = removeProcessedFiles(fNames,'.rls','_c_K_s.mat','registration');
%
% Dependencies: base MATLAB; SETTINGS *.mat files saved by the LSCI workflow.
% See also: getFileNamesList, getProductPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------

function fNames = removeProcessedFiles(fNames,extStr,extReplStr,currentStep,keepFirst)
if ~iscell(fNames)
    error('fNames must be a cell array of file names.');
end
if nargin<5 || isempty(keepFirst)
    keepFirst=false;
end

done=false(size(fNames));
for idx=1:numel(fNames)
    if isempty(fNames{idx})
        continue                                   % no file here -> nothing to check
    end
    settingsFile=strrep(fNames{idx},extStr,extReplStr);
    if ~endsWith(settingsFile,'.mat','IgnoreCase',true) || ~isfile(settingsFile)
        continue                                   % no SETTINGS *.mat yet -> keep for analysis
    end
    S=load(settingsFile,'settings');
    done(idx)= isfield(S,'settings') && isfield(S.settings,currentStep);
end

if keepFirst && ~isempty(done)
    if isvector(done)
        done(1)=false;                             % plain vector -> keep only the first entry
    else
        done(:,1)=false;                           % 2-D array -> keep the whole first column
    end
end

fNames(done)={[]};                                 % empty the already-processed entries (shape preserved)
disp([num2str(nnz(done)),' of ',num2str(numel(done)),' file(s) already processed for "',currentStep,'" and set to empty.'])
end