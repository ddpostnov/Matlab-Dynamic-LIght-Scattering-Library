%getFileNamesList - Collect (and optionally group/sort) file paths under a folder.
%
%   Searches rootFolder recursively for files matching a glob and, optionally,
%   groups them by an "animal" id (one row per animal), ordering each group's
%   columns either by a set of regexp patterns or by pinning a reference file
%   first.  Full paths are always returned.
%
% Syntax:
%    fNames = getFileNamesList(rootFolder, fileTypeIdentifier)
%    fNames = getFileNamesList(rootFolder, fileTypeIdentifier, animalIdentifier)
%    fNames = getFileNamesList(rootFolder, fileTypeIdentifier, animalIdentifier, sortingPatterns)
%
% Inputs:
%    rootFolder         - folder searched recursively.
%    fileTypeIdentifier - dir() glob, e.g. '*c_K_d.mat' ('*' is allowed).
%    animalIdentifier   - (optional) regexp extracting the animal id, e.g.
%                         '[A-Z]+\d+'.  Empty/omitted -> no grouping (flat list).
%    sortingPatterns    - (optional) either a cell array of regexps (one output
%                         column per pattern) or a single "reference" regexp
%                         whose matching file is forced into the first column.
%                         Passed to regexp - use regex syntax, not globs.
%
% Outputs:
%    fNames - full paths.  N-by-1 flat list when no animalIdentifier is given;
%             otherwise one row per animal, columns set by the sorting rule.
%
% Examples:
%    fNames = getFileNamesList(root, '*c_K_d.mat');
%    fNames = getFileNamesList(root, '*c_K_d.mat', '[A-Z]+\d+', '1BP\.mat');
%    fNames = getFileNamesList(root, '*c_K_d.mat', '[A-Z]+\d+', {'BV\.mat$','BP\.mat$'});
%
% See also: dir, regexp, removeProcessedFiles
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-July-2026

%------------- BEGIN CODE --------------
function fNames = getFileNamesList(rootFolder, fileTypeIdentifier, animalIdentifier, sortingPatterns)

    % ---- defaults ------------------------------------------------------
    if nargin < 1 || isempty(rootFolder)
        error('getFileNamesList:noRoot', 'rootFolder is required.');
    end
    if nargin < 2 || isempty(fileTypeIdentifier)
        error('getFileNamesList:noFileType', 'fileTypeIdentifier is required.');
    end
    if nargin < 3,                              animalIdentifier = '';  end
    if nargin < 4 || isempty(sortingPatterns),  sortingPatterns  = {};  end

    % ---- validate / normalise inputs ----------------------------------
    rootFolder = char(rootFolder);
    if ~isfolder(rootFolder)
        error('getFileNamesList:badRoot', '"%s" is not a valid folder.', rootFolder);
    end
    fileTypeIdentifier = char(fileTypeIdentifier);

    if isstring(animalIdentifier)
        animalIdentifier = char(animalIdentifier);
    elseif ~ischar(animalIdentifier) && ~isempty(animalIdentifier)
        error('getFileNamesList:badAnimalId', 'animalIdentifier must be char or string.');
    end

    % The 4th argument is EITHER:
    %   - a cell array of regexps   -> one output column per pattern (sorting), or
    %   - a single char/string      -> a "reference" regexp; the file matching it
    %                                  is forced to be the first column for each
    %                                  animal, the rest follow in default order.
    referencePattern = '';
    if ischar(sortingPatterns) || (isstring(sortingPatterns) && isscalar(sortingPatterns))
        referencePattern = char(sortingPatterns);   % reference mode (may be '')
        sortingPatterns  = {};
    elseif isstring(sortingPatterns)                 % multi-element string -> patterns
        sortingPatterns = cellstr(sortingPatterns);
    elseif iscell(sortingPatterns)
        sortingPatterns = cellfun(@char, sortingPatterns, 'UniformOutput', false);
    else
        error('getFileNamesList:badSorting', ...
            ['sortingPatterns must be a cell array of regexps (one per column) ', ...
             'or a single reference regexp string.']);
    end
    sortingPatterns  = sortingPatterns(:).';
    referenceMode    = ~isempty(referencePattern);

    % ---- find files (files only, deterministic order) ------------------
    d = dir(fullfile(rootFolder, '**', fileTypeIdentifier));
    d = d(~[d.isdir]);
    if isempty(d)
        warning('getFileNamesList:noFiles', ...
            'No files matching "%s" found under "%s".', fileTypeIdentifier, rootFolder);
        fNames = {};
        return;
    end
    fullPaths = fullfile({d.folder}, {d.name});
    [fullPaths, order] = sort(fullPaths);
    files = {d.name};
    files = files(order);

    % ---- case 1: no animal grouping -> flat list -----------------------
    if isempty(animalIdentifier)
        fNames = fullPaths(:);
        return;
    end

    % ---- group by animal ----------------------------------------------
    try
        expIdentifiers = regexp(files, animalIdentifier, 'match', 'once');
    catch ME
        error('getFileNamesList:badAnimalRegex', ...
            'animalIdentifier "%s" is not a valid regular expression:\n%s', ...
            animalIdentifier, ME.message);
    end
    uniqueGroups = unique(expIdentifiers(~cellfun(@isempty, expIdentifiers)));
    groupCount   = numel(uniqueGroups);

    if groupCount == 0
        warning('getFileNamesList:noGroups', ...
            'No matches for animalIdentifier "%s"; returning a flat list.', animalIdentifier);
        fNames = fullPaths(:);
        return;
    end

    % ---- reference mode: pin the matching file first, rest in default order
    if referenceMode
        grouped  = cell(groupCount, 1);
        maxFiles = 0;
        for rowIdx = 1:groupCount
            inGroup    = strcmp(expIdentifiers, uniqueGroups{rowIdx});
            groupPaths = fullPaths(inGroup);
            groupFiles = files(inGroup);
            try
                isRef = ~cellfun(@isempty, regexp(groupFiles, referencePattern, 'once'));
            catch ME
                error('getFileNamesList:badRefRegex', ...
                    'reference expression "%s" is not a valid regular expression:\n%s', ...
                    referencePattern, ME.message);
            end
            % reference match(es) first, then everything else in default order
            grouped{rowIdx} = [groupPaths(isRef), groupPaths(~isRef)];
            maxFiles = max(maxFiles, numel(groupPaths));
        end
        fNames = cell(groupCount, maxFiles);
        for rowIdx = 1:groupCount
            n = numel(grouped{rowIdx});
            fNames(rowIdx, 1:n) = grouped{rowIdx};
        end
        return;
    end

    % ---- case 2: grouping, no sorting patterns -------------------------
    if isempty(sortingPatterns)
        grouped  = cell(groupCount, 1);
        maxFiles = 0;
        for rowIdx = 1:groupCount
            inGroup = strcmp(expIdentifiers, uniqueGroups{rowIdx});  % exact match
            grouped{rowIdx} = fullPaths(inGroup);
            maxFiles = max(maxFiles, nnz(inGroup));
        end
        fNames = cell(groupCount, maxFiles);
        for rowIdx = 1:groupCount
            n = numel(grouped{rowIdx});
            fNames(rowIdx, 1:n) = grouped{rowIdx};
        end
        return;
    end

    % ---- case 3: grouping + sorting patterns ---------------------------
    nCols  = numel(sortingPatterns);
    fNames = cell(groupCount, nCols);
    for rowIdx = 1:groupCount
        inGroup    = strcmp(expIdentifiers, uniqueGroups{rowIdx});
        groupPaths = fullPaths(inGroup);
        groupFiles = files(inGroup);
        for colIdx = 1:nCols
            try
                matchIndex = ~cellfun(@isempty, ...
                    regexp(groupFiles, sortingPatterns{colIdx}, 'once'));
            catch ME
                error('getFileNamesList:badSortRegex', ...
                    'sortingPatterns{%d} = "%s" is not a valid regular expression:\n%s', ...
                    colIdx, sortingPatterns{colIdx}, ME.message);
            end
            if any(matchIndex)
                if nnz(matchIndex) > 1
                    warning('getFileNamesList:multipleMatches', ...
                        'Animal "%s", column %d: %d files matched "%s"; using the first.', ...
                        uniqueGroups{rowIdx}, colIdx, nnz(matchIndex), sortingPatterns{colIdx});
                    firstHit = find(matchIndex, 1, 'first');
                    matchIndex(:) = false;
                    matchIndex(firstHit) = true;
                end
                fNames(rowIdx, colIdx) = groupPaths(matchIndex);
            end
        end
    end
end
