%removeProcessedFiles  Drop files whose analysis step has already been done
%
%   fNames = checkIfProcessed(fNames,extStr,extReplStr,currentStep,keepFirst)
%   inspects every entry of the cell array fNames (1-D or 2-D) and empties the
%   files that have already been through the requested processing step, so they
%   are not analysed a second time.  For each file the token extStr in the name
%   is replaced by extReplStr to obtain the matching SETTINGS file—a *.mat
%   holding a "settings" structure—which is then loaded.  If
%   settings.(currentStep) exists the step is considered done and that entry is
%   set to [] (the array keeps its original size and shape).
%
%   INPUTS
%     fNames       cell array of file names, may be 2-D.  Empty entries are
%                  ignored and left in place.
%     extStr       token/extension to replace, e.g. '.rls' or '_d.mat'
%     extReplStr   replacement giving the SETTINGS file, e.g. '_c_K_s.mat'
%     currentStep  settings field name to test for, e.g. 'registration'
%     keepFirst    logical (optional, default false).  When true the first
%                  column fNames(:,1) is always kept, even if already
%                  processed, as it is typically the registration reference.
%
%   OUTPUT
%     fNames       the input list with the same size and shape, but every
%                  already-processed entry set to [] so it can be treated as
%                  non-existent with an isempty check downstream.
%
%   EXAMPLE
%     D = dir(fullfile(dataRoot,'*.rls'));
%     fNames = fullfile({D.folder}',{D.name}');
%     fNames = checkIfProcessed(fNames,'.rls','_c_K_s.mat','registration');
%     % fNames now lists only the recordings still awaiting registration
%
%   DEPENDS ON
%     MATLAB base only (strrep, endsWith, isfile, load, isfield).  Files must
%     have been saved by the LSCI workflow so the SETTINGS *.mat contains a
%     "settings" structure.
%
%   ----------------------------------------------------------------------
%   Copyright © 2025 Dmitry D Postnov, Aarhus University
%   e-mail: dpostnov@cfin.au.dk
%   Last revision: 02-Jul-2026
%   ----------------------------------------------------------------------

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