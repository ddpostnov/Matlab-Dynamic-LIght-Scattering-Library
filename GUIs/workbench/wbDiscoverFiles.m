%wbDiscoverFiles - Three loaders producing one grouped, reference-first file set.
%
%   The workbench discovers files three ways, all returning the SAME grouped
%   structure (rows = animal/group, column 1 = the reference/template when a
%   reference rule is given):
%     'structured' - thin wrapper over getFileNamesList (reused verbatim).
%     'folder'     - recurse a folder for a glob, group by an animal regexp (or
%                    one group when no regexp is given).
%     'manual'     - group an explicit path list by an animal regexp (or one
%                    group).  No disk scan.
%
% Syntax:
%    disc = wbDiscoverFiles('structured', root, glob, animalRegexp, refRegexp)
%    disc = wbDiscoverFiles('folder',     root, glob, animalRegexp)
%    disc = wbDiscoverFiles('manual',     pathList, animalRegexp)
%
% Inputs:
%    root         - folder searched recursively (structured/folder).
%    glob         - dir() glob, e.g. '*_c_K_d.mat'.
%    animalRegexp - regexp extracting the group id (e.g. '[A-Z]+\d+'); empty ->
%                   a single group.
%    refRegexp    - regexp pinning the reference file into column 1 (structured).
%    pathList     - cell array of full paths (manual).
%
% Outputs:
%    disc - struct with fields:
%       fNames        2-D cell (group rows, reference-first) - the execution grid.
%       models        2-D cell of wbFileModel aligned to fNames ([] where empty).
%       flat          1xM struct array of the non-empty models.
%       groups        struct array: .name and .rowIndex.
%       referenceMode logical (column 1 is a reference/template).
%
% See also: getFileNamesList, wbFileModel, removeProcessedFiles
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function disc = wbDiscoverFiles(mode, varargin)

referenceMode = false;
switch mode
    case 'structured'
        [root, glob, animalRegexp, refRegexp] = parseArgs(varargin, 4);
        fNames = getFileNamesList(root, glob, animalRegexp, refRegexp);
        referenceMode = ~isempty(refRegexp);

    case 'folder'
        [root, glob, animalRegexp] = parseArgs(varargin, 3);
        if isempty(animalRegexp)
            flat = getFileNamesList(root, glob);        % N-by-1 flat list
            fNames = reshape(flat, 1, []);              % one group (a single row)
        else
            fNames = getFileNamesList(root, glob, animalRegexp);
        end

    case 'manual'
        [pathList, animalRegexp] = parseArgs(varargin, 2);
        pathList = pathList(~cellfun(@isempty, pathList));
        pathList = pathList(:);
        if isempty(animalRegexp)
            fNames = reshape(pathList, 1, []);          % one group
        else
            fNames = groupManual(pathList, animalRegexp);
        end

    otherwise
        error('wbDiscoverFiles:badMode','Unknown mode "%s".',mode);
end

% ---- build the parallel model grid ------------------------------------------
[nr, nc] = size(fNames);
models = cell(nr, nc);
flat   = struct('path',{},'folder',{},'name',{},'ext',{},'modality',{}, ...
                'roi',{},'roiPrefix',{},'stem',{},'identity',{},'flags',{}, ...
                'temporal',{},'modifier',{},'product',{},'role',{}, ...
                'isRaw',{},'isReference',{},'group',{});
groups = struct('name',{},'rowIndex',{});
for r = 1:nr
    gname = groupNameForRow(fNames(r,:));
    groups(r) = struct('name',gname,'rowIndex',r);
    for c = 1:nc
        if isempty(fNames{r,c}), continue; end
        m = wbFileModel(fNames{r,c});
        m.group = gname;
        m.isReference = referenceMode && c==1;
        models{r,c} = m;
        flat(end+1) = m; %#ok<AGROW>
    end
end

disc = struct('fNames',{fNames},'models',{models},'flat',flat, ...
              'groups',groups,'referenceMode',referenceMode);
end

% =====================================================================
function varargout = parseArgs(args, n)
%parseArgs  Positional args padded with '' so optional trailing ones may be omitted.
out = repmat({''},1,n);
for i = 1:min(n,numel(args)), out{i} = args{i}; end
varargout = out;
end

% =====================================================================
function fNames = groupManual(pathList, animalRegexp)
%groupManual  Bucket an explicit path list into group rows by an animal regexp.
names = cell(size(pathList));
for i = 1:numel(pathList)
    [~,nm,ex] = fileparts(pathList{i});
    names{i} = [nm ex];
end
ids = regexp(names, animalRegexp, 'match', 'once');
ids(cellfun(@isempty, ids)) = {'(ungrouped)'};
uniq = unique(ids, 'stable');
rows = cell(numel(uniq),1);
maxc = 0;
for g = 1:numel(uniq)
    rows{g} = pathList(strcmp(ids, uniq{g}));
    maxc = max(maxc, numel(rows{g}));
end
fNames = cell(numel(uniq), maxc);
for g = 1:numel(uniq)
    fNames(g,1:numel(rows{g})) = reshape(rows{g},1,[]);
end
end

% =====================================================================
function gname = groupNameForRow(rowPaths)
%groupNameForRow  Best-effort group label from the first non-empty file's stem.
gname = '';
for c = 1:numel(rowPaths)
    if ~isempty(rowPaths{c})
        m = wbFileModel(rowPaths{c});
        gname = m.stem;
        return
    end
end
end
