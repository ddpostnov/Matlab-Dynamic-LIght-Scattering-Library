%wbDiscoverFiles - Four loaders producing one per-animal, reference-first file set.
%
%   The workbench discovers files four ways, all returning the SAME grouped
%   structure (rows = ANIMAL, column 1 = the reference/template when a reference
%   rule is given):
%     'structured' - thin wrapper over getFileNamesList (reused verbatim).
%     'folder'     - recurse a folder for a glob, bucket by an animal regexp (or
%                    one animal row when no regexp is given).
%     'manual'     - bucket an explicit path list by an animal regexp (or one
%                    animal row).  No disk scan.
%     'curated'    - re-grid an already-labelled path list: bucket by the CURATED
%                    animal label (hand overrides included) and pin each animal's
%                    reference recording into column 1.  No disk scan, no regexp -
%                    this is the loop the Files tab runs after every edit.
%
%   The grid rows are the ANIMAL axis (getFileNamesList's animalIdentifier), the
%   scope of registration / vessel typing / the reference recording.  It is NOT
%   the experimental group - the two are independent per-file labels; see
%   wbTypeModel, which derives all four label axes and which this function calls
%   to stamp .animal / .type / .index / .expGroup onto every model.
%
% Syntax:
%    disc = wbDiscoverFiles('structured', root, glob, animalRegexp, refRegexp, typeRegexp, groupRegexp)
%    disc = wbDiscoverFiles('folder',     root, glob, animalRegexp, typeRegexp, groupRegexp)
%    disc = wbDiscoverFiles('manual',     pathList, animalRegexp, typeRegexp, groupRegexp)
%    disc = wbDiscoverFiles('curated',    pathList, labels, animalRefMap)
%
% Inputs:
%    root         - folder searched recursively (structured/folder).
%    glob         - dir() glob, e.g. '*_c_K_r.mat'.
%    animalRegexp - regexp extracting the animal id (e.g. '[A-Z]+\d+'); empty ->
%                   a single animal row.
%    refRegexp    - regexp pinning the reference file into column 1 (structured).
%    typeRegexp   - (optional) regexp labelling the recording TYPE; no match ->
%                   wbTypeModel's '(untyped)' bucket.  Never validated against a
%                   list - the vocabulary is whatever the names contain.
%    groupRegexp  - (optional) regexp labelling the EXPERIMENTAL group; no match
%                   -> '(ungrouped)'.
%    pathList     - cell array of full paths (manual/curated).
%    labels       - (curated) a wbTypeModel labels struct aligned to pathList.
%    animalRefMap - (curated, optional) containers.Map animal -> the reference
%                   RECORDING IDENTITY; that recording's files lead its row.
%
% Outputs:
%    disc - struct with fields:
%       fNames        2-D cell (animal rows, reference-first) - the execution grid.
%       models        2-D cell of wbFileModel aligned to fNames ([] where empty).
%       flat          1xM struct array of the non-empty models.
%       animals       struct array: .name and .rowIndex.
%       referenceMode logical (column 1 is a reference/template).
%       patterns      struct of the regexps used (animal/type/index/expGroup/ref).
%
% Notes:
%    Adding typeRegexp/groupRegexp changes NOTHING about the grid: fNames, the
%    animals rows and every non-label model field are what the same call without
%    them produces.  The regexps only fill in label fields that would otherwise
%    hold their no-match bucket.
%
% See also: getFileNamesList, wbFileModel, wbTypeModel, removeProcessedFiles
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function disc = wbDiscoverFiles(mode, varargin)

referenceMode = false;
patterns      = wbTypeModel('emptyPatterns');
patterns.ref  = '';
labels        = [];
refMap        = containers.Map('KeyType','char','ValueType','char');

switch mode
    case 'structured'
        [root, glob, animalRegexp, refRegexp, typeRegexp, groupRegexp] = parseArgs(varargin, 6);
        fNames = getFileNamesList(root, glob, animalRegexp, refRegexp);
        referenceMode = ~isempty(refRegexp);
        patterns = setPatterns(patterns, animalRegexp, typeRegexp, groupRegexp, refRegexp);

    case 'folder'
        [root, glob, animalRegexp, typeRegexp, groupRegexp] = parseArgs(varargin, 5);
        if isempty(animalRegexp)
            flat = getFileNamesList(root, glob);        % N-by-1 flat list
            fNames = reshape(flat, 1, []);              % one animal row
        else
            fNames = getFileNamesList(root, glob, animalRegexp);
        end
        patterns = setPatterns(patterns, animalRegexp, typeRegexp, groupRegexp, '');

    case 'manual'
        [pathList, animalRegexp, typeRegexp, groupRegexp] = parseArgs(varargin, 4);
        pathList = cleanPaths(pathList);
        if isempty(animalRegexp)
            fNames = reshape(pathList, 1, []);          % one animal row
        else
            fNames = bucketByAnimal(pathList, animalRegexp);
        end
        patterns = setPatterns(patterns, animalRegexp, typeRegexp, groupRegexp, '');

    case 'curated'
        [pathList, labels, refMapIn] = parseArgs(varargin, 3);
        pathList = cleanPaths(pathList);
        if isa(refMapIn,'containers.Map'), refMap = refMapIn; end
        [fNames, labels, referenceMode] = gridFromLabels(pathList, labels, refMap);

    otherwise
        error('wbDiscoverFiles:badMode','Unknown mode "%s".',mode);
end

% ---- derive the per-file labels (curated mode arrives with them) -------------
if isempty(labels)
    labels = wbTypeModel('derive', reshape(fNames(~cellfun(@isempty,fNames)),1,[]), patterns);
end
%   A caller that lets the user correct the parsed MODALITY (the Files tab does)
%   hands it in as labels.modality, aligned like the regexp axes; without it the
%   parser's guess stands.
hasModality = isfield(labels,'modality');
labelOf = containers.Map('KeyType','char','ValueType','any');
for i = 1:labels.n
    L = struct('animal',labels.animal{i},'type',labels.type{i}, ...
               'index',labels.index{i},'expGroup',labels.expGroup{i},'modality','');
    if hasModality, L.modality = labels.modality{i}; end
    labelOf(labels.path{i}) = L;
end

% ---- build the parallel model grid ------------------------------------------
[nr, nc] = size(fNames);
models = cell(nr, nc);
flat   = struct('path',{},'folder',{},'name',{},'ext',{},'modality',{}, ...
                'roi',{},'roiPrefix',{},'stem',{},'identity',{},'flags',{}, ...
                'stage',{},'branch',{},'product',{},'role',{}, ...
                'isRaw',{},'isReference',{},'animal',{},'type',{},'index',{},'expGroup',{});
animals = struct('name',{},'rowIndex',{});
labelled = ~isempty(patterns.animal) || strcmp(mode,'curated');   % is .animal meaningful?
for r = 1:nr
    aname = animalNameForRow(fNames(r,:), labelOf, labelled);
    animals(r) = struct('name',aname,'rowIndex',r);
    for c = 1:nc
        if isempty(fNames{r,c}), continue; end
        m = wbFileModel(fNames{r,c});
        L = labelOf(fNames{r,c});
        m.animal      = L.animal;
        m.type        = L.type;
        m.index       = L.index;
        m.expGroup    = L.expGroup;
        if ~isempty(L.modality), m.modality = L.modality; end
        m.isReference = referenceMode && c==1;
        models{r,c} = m;
        flat(end+1) = m; %#ok<AGROW>
    end
end

disc = struct('fNames',{fNames},'models',{models},'flat',flat, ...
              'animals',animals,'referenceMode',referenceMode,'patterns',patterns);
end

% =====================================================================
function varargout = parseArgs(args, n)
%parseArgs  Positional args padded with '' so optional trailing ones may be omitted.
out = repmat({''},1,n);
for i = 1:min(n,numel(args)), out{i} = args{i}; end
varargout = out;
end

% =====================================================================
function p = setPatterns(p, animalRegexp, typeRegexp, groupRegexp, refRegexp)
%setPatterns  Record the regexps the caller asked for (index has no loader arg).
p.animal   = charOf(animalRegexp);
p.type     = charOf(typeRegexp);
p.expGroup = charOf(groupRegexp);
p.ref      = charOf(refRegexp);
end
function s = charOf(v), if isempty(v), s = ''; else, s = char(v); end, end

% =====================================================================
function pathList = cleanPaths(pathList)
%cleanPaths  Drop empties and return a column cellstr.
if isempty(pathList), pathList = {}; return; end
if ~iscell(pathList), pathList = {pathList}; end
pathList = pathList(~cellfun(@isempty, pathList));
pathList = reshape(pathList, [], 1);
end

% =====================================================================
function fNames = bucketByAnimal(pathList, animalRegexp)
%bucketByAnimal  Bucket an explicit path list into animal rows by an animal regexp.
names = cell(size(pathList));
for i = 1:numel(pathList)
    [~,nm,ex] = fileparts(pathList{i});
    names{i} = [nm ex];
end
ids = regexp(names, animalRegexp, 'match', 'once');
ids(cellfun(@isempty, ids)) = {wbTypeModel('default','animal')};
fNames = gridFromIds(pathList, ids);
end
% =====================================================================
function id = identityOf(pth)
%identityOf  The recording identity of one path (branch flags stripped).
m = wbFileModel(pth); id = m.identity;
end

% =====================================================================
function [fNames, labels, pinned] = gridFromLabels(pathList, labels, refMap)
%gridFromLabels  Curated grid: rows = the curated animal label, reference first.
%   Every file of the animal's reference RECORDING (identity match) leads the row,
%   so column 1 is that recording whichever branch product happens to be listed.
if isempty(labels) || ~isstruct(labels)
    labels = wbTypeModel('derive', pathList, wbTypeModel('emptyPatterns'));
end
pinned = false;
if isempty(pathList), fNames = {}; return; end
byPath = containers.Map('KeyType','char','ValueType','char');
for i = 1:labels.n, byPath(labels.path{i}) = labels.animal{i}; end
ids = cell(1,numel(pathList));
for i = 1:numel(pathList)
    if isKey(byPath, pathList{i}), ids{i} = byPath(pathList{i});
    else,                          ids{i} = wbTypeModel('default','animal'); end
end
[fNames, pinned] = gridFromIds(pathList, ids, refMap);
end

% =====================================================================
function [fNames, pinned] = gridFromIds(pathList, ids, refMap)
%gridFromIds  One row per unique id (first-appearance order), reference first.
if nargin<3, refMap = containers.Map('KeyType','char','ValueType','char'); end
pinned = false;
uniq = unique(ids, 'stable');
rows = cell(numel(uniq),1);
maxc = 0;
for g = 1:numel(uniq)
    member = pathList(strcmp(ids, uniq{g}));
    if isKey(refMap, uniq{g})
        refId  = refMap(uniq{g});
        isRef  = cellfun(@(p) strcmp(identityOf(p), refId), member);
        if any(isRef)
            member = [member(isRef); member(~isRef)];
            pinned = true;
        end
    end
    rows{g} = member;
    maxc = max(maxc, numel(member));
end
fNames = cell(numel(uniq), maxc);
for g = 1:numel(uniq)
    fNames(g,1:numel(rows{g})) = reshape(rows{g},1,[]);
end
end

% =====================================================================
function aname = animalNameForRow(rowPaths, labelOf, labelled)
%animalNameForRow  The row's animal id: the derived/curated ANIMAL label of its
%   first file when there is an animal rule at all, else - no rule, so the label
%   would be a meaningless bucket - the first file's stem, as before.
aname = '';
for c = 1:numel(rowPaths)
    if isempty(rowPaths{c}), continue; end
    if labelled && isKey(labelOf, rowPaths{c})
        L = labelOf(rowPaths{c});
        if ~isempty(L.animal), aname = L.animal; return; end
    end
    m = wbFileModel(rowPaths{c});
    aname = m.stem;
    return
end
end
