function fNames = setFileNamesList(rootFolder, fileTypeIdentifier, animalIdentifier, sortingPatterns)
%SETFILENAMESLIST  Collect (and optionally group/sort) file paths under a folder.
%
%   fNames = setFileNamesList(rootFolder, fileTypeIdentifier)
%   fNames = setFileNamesList(rootFolder, fileTypeIdentifier, animalIdentifier)
%   fNames = setFileNamesList(rootFolder, fileTypeIdentifier, animalIdentifier, sortingPatterns)
%
%   Inputs:
%     rootFolder         - char/string. Folder searched recursively.
%     fileTypeIdentifier - (required) dir() wildcard, e.g. '*c_K_d.mat'.
%                          NOTE: this is a glob, so '*' is allowed here.
%     animalIdentifier   - (optional) regexp extracting the animal id, e.g.
%                          '[A-Z]+\d+'. If empty/omitted, no grouping is done.
%     sortingPatterns    - (optional) defines the second dimension. Either:
%                          * a cell array of regexps -> one output column per
%                            pattern (e.g. {'BV\.mat$','BP\.mat$'}), or
%                          * a single regexp string -> "reference" file; the
%                            file matching it is placed first for each animal,
%                            the remaining files follow in default name order.
%                          NOTE: these go to regexp, so use regex syntax
%                          (e.g. 'BV\.mat$' / '1BP\.mat'), not globs.
%
%   Output (full paths in every case):
%     - no animalIdentifier                       : N-by-1 flat list of files.
%     - animalIdentifier, no sortingPatterns      : one row per animal, columns
%                                                   are that animal's files in
%                                                   name order.
%     - animalIdentifier + reference string       : as above, but the file
%                                                   matching the reference is
%                                                   forced into the first column.
%     - animalIdentifier + sortingPatterns cell   : one row per animal, one
%                                                   column per pattern.
%
%   Examples:
%     fNames = setFileNamesList(rootFolder, '*c_K_d.mat');
%     fNames = setFileNamesList(rootFolder, '*c_K_d.mat', '[A-Z]+\d+', '1BP\.mat');
%     fNames = setFileNamesList(rootFolder, '*c_K_d.mat', '[A-Z]+\d+', ...
%                               {'BV\.mat$', 'BP\.mat$'});

    % ---- defaults ------------------------------------------------------
    if nargin < 1 || isempty(rootFolder)
        error('setFileNamesList:noRoot', 'rootFolder is required.');
    end
    if nargin < 2 || isempty(fileTypeIdentifier)
        error('setFileNamesList:noFileType', 'fileTypeIdentifier is required.');
    end
    if nargin < 3,                              animalIdentifier = '';  end
    if nargin < 4 || isempty(sortingPatterns),  sortingPatterns  = {};  end

    % ---- validate / normalise inputs ----------------------------------
    rootFolder = char(rootFolder);
    if ~isfolder(rootFolder)
        error('setFileNamesList:badRoot', '"%s" is not a valid folder.', rootFolder);
    end
    fileTypeIdentifier = char(fileTypeIdentifier);

    if isstring(animalIdentifier)
        animalIdentifier = char(animalIdentifier);
    elseif ~ischar(animalIdentifier) && ~isempty(animalIdentifier)
        error('setFileNamesList:badAnimalId', 'animalIdentifier must be char or string.');
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
        error('setFileNamesList:badSorting', ...
            ['sortingPatterns must be a cell array of regexps (one per column) ', ...
             'or a single reference regexp string.']);
    end
    sortingPatterns  = sortingPatterns(:).';
    referenceMode    = ~isempty(referencePattern);

    % ---- find files (files only, deterministic order) ------------------
    d = dir(fullfile(rootFolder, '**', fileTypeIdentifier));
    d = d(~[d.isdir]);
    if isempty(d)
        warning('setFileNamesList:noFiles', ...
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
        error('setFileNamesList:badAnimalRegex', ...
            'animalIdentifier "%s" is not a valid regular expression:\n%s', ...
            animalIdentifier, ME.message);
    end
    uniqueGroups = unique(expIdentifiers(~cellfun(@isempty, expIdentifiers)));
    groupCount   = numel(uniqueGroups);

    if groupCount == 0
        warning('setFileNamesList:noGroups', ...
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
                error('setFileNamesList:badRefRegex', ...
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
                error('setFileNamesList:badSortRegex', ...
                    'sortingPatterns{%d} = "%s" is not a valid regular expression:\n%s', ...
                    colIdx, sortingPatterns{colIdx}, ME.message);
            end
            if any(matchIndex)
                if nnz(matchIndex) > 1
                    warning('setFileNamesList:multipleMatches', ...
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
