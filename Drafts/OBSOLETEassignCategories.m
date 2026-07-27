%assignCategories - Copy a reference category mask (cMask) onto other datasets.
%
%   For each target/reference pair of *_K_d.mat datasets, copies the reference
%   file's category mask (results.cMask) and the runCategories edge size into the
%   target's RESULTS and SETTINGS.  Use it to reuse one categorisation across
%   co-registered recordings instead of re-running runCategories.
%
% Syntax:
%    assignCategories(fNames, fNamesRef)
%
% Inputs:
%    fNames    - cell array of target *_K_d.mat paths.
%    fNamesRef - cell array of reference *_K_d.mat paths (same size as fNames);
%                each provides the cMask copied onto the matching target.
%
% Outputs:
%    (none) - updates each target's *_K_r.mat (results.cMask) and *_K_s.mat
%             (settings.runCategories) in place.
%
% See also: runCategories, setVesselTypes, runRegistration
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function assignCategories(fNames,fNamesRef)

if ~all( cellfun(@(s) isempty(s) || contains(s,'_K_d.mat'), fNames(:)) )
    error('One or more *non-empty* entries do not contain "_K_d.mat".');
end

if numel(fNames)~=numel(fNamesRef)
    error('Number of reference and target files should match');
end

for fidx=1:1:numel(fNames)
    if ~isempty(fNames{fidx}) && ~isempty(fNamesRef{fidx})
        tic
        disp(['Processing file ',num2str(fidx),' out of ',num2str(numel(fNames))])
        s.fName=fNames{fidx};
        s.fNameCRef=fNamesRef{fidx};
        clearvars results source settings
        load(strrep(s.fNameCRef,'_d.mat','_r.mat'),'results');
        load(strrep(s.fNameCRef,'_d.mat','_s.mat'),'settings');
        s.edgeSize=settings.runCategories.edgeSize;
        cMask=results.cMask;
        clearvars results

        load(strrep(s.fName,'_d.mat','_s.mat'),'settings');
        load(strrep(s.fName,'_d.mat','_r.mat'),'results');

        results.cMask=cMask;
        settings.runCategories=s;
        %Save the data
        disp(['Saving the results. Elapsed time ',num2str(round(toc)),'s']);
        save(strrep(s.fName,'_d.mat','_s.mat'),'settings','-v7.3');
        save(strrep(s.fName,'_d.mat','_r.mat'),'results','-v7.3');
        disp('Saving complete');
    end
end
end