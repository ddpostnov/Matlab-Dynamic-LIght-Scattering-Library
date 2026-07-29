%wbTypeModel - Per-file label layer: animal / type / index / experimental group.
%
%   The workbench labels every discovered file on four INDEPENDENT axes, each
%   derived by its own regexp over the file NAME and each overridable by hand:
%
%     animal   - the subject.  Same brain, same FOV, co-registerable.  Owns ONE
%                reference recording, valid for all of its files.  The scope of
%                registration and vessel typing.
%     type     - the recording's experimental role; owns the processing
%                configuration.  Global across animals.
%     index    - the recording index inside an animal (a stratifier).
%     expGroup - the EXPERIMENTAL group, a comparison label.  Used by Export and
%                Explore only; processing ignores it.
%
%   The axes DO NOT NEST: a group may span many animals and one animal may span
%   several groups, so no axis is ever derived from another.  One mechanism
%   serves all four because they are all just per-file labels.
%
%   THE VOCABULARIES ARE OPEN.  There is no built-in list of animals, types or
%   groups anywhere in the workbench.  Every token is whatever the user's regexp
%   matched or typed; the set is discovered at scan time and changes when a label
%   is edited.  Nothing here (or downstream) may enumerate, validate against, or
%   branch on a token's value - 'BV'/'BN'/'BP', 'ctrl'/'stim'/'washout' and
%   '1'/'2'/'3' all behave identically.  Only the AXES are fixed; their values
%   are not.  A file whose name does not match an axis's regexp is not an error -
%   it lands in that axis's no-match bucket ('(untyped)', '(ungrouped)', ...),
%   which is an ordinary value like any other.
%
%   Pure data: no graphics, no file I/O, no disk access - unit-testable.
%
% Syntax:
%    labels = wbTypeModel('derive', models, patterns)
%    labels = wbTypeModel('applyOverrides', labels, overrides)   % or (.., models, ovr)
%    v      = wbTypeModel('values', labels, axis)
%    ax     = wbTypeModel('axes')
%    tok    = wbTypeModel('default', axis)
%    ovr    = wbTypeModel('emptyOverrides')
%    p      = wbTypeModel('emptyPatterns')
%
% Inputs:
%    models    - struct array of wbFileModel structs, OR a cellstr of paths.
%    patterns  - struct with any of the fields animal/type/index/expGroup, each a
%                regexp matched against the file NAME ('match','once' - the same
%                rule getFileNamesList uses for its animalIdentifier).  A missing
%                or empty pattern sends every file to that axis's no-match bucket.
%    overrides - struct with the same axis fields, each a containers.Map keyed by
%                the file PATH -> the hand-assigned value.  Overrides win over the
%                regexp and, being keyed by path, survive a re-scan (the mechanism
%                guiExplore>applyGroupOverrides has proven).
%    axis      - 'animal' | 'type' | 'index' | 'expGroup'.
%
% Outputs:
%    labels - struct with one 1xN cellstr per axis plus .path (1xN, the alignment
%             key) and .n.  N = numel(models).
%    v      - ordered unique values of that axis (order of first appearance).
%
% See also: wbFileModel, wbDiscoverFiles, guiWorkbench, getFileNamesList
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbTypeModel(action, varargin)

switch action
    case 'axes',           out = axisList();
    case 'default',        out = defaultToken(varargin{1});
    case 'emptyOverrides', out = emptyOverrides();
    case 'emptyPatterns',  out = emptyPatterns();
    case 'derive',         out = derive(varargin{:});
    case 'applyOverrides', out = applyOverrides(varargin{:});
    case 'values',         out = axisValues(varargin{:});
    otherwise
        error('wbTypeModel:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function ax = axisList()
%axisList  The FIXED set of label axes (their values are not fixed - see above).
ax = {'animal','type','index','expGroup'};
end

% =====================================================================
function tok = defaultToken(axis)
%defaultToken  The no-match bucket of one axis - an ordinary value, not an error.
switch axis
    case 'animal',   tok = '(unassigned)';
    case 'type',     tok = '(untyped)';
    case 'index',    tok = '1';            % "every file is the same index"
    case 'expGroup', tok = '(ungrouped)';
    otherwise,       tok = '(none)';
end
end

% =====================================================================
function p = emptyPatterns()
%emptyPatterns  A pattern struct with every axis switched off.
ax = axisList();
p = struct();
for i = 1:numel(ax), p.(ax{i}) = ''; end
end

% =====================================================================
function o = emptyOverrides()
%emptyOverrides  One empty path->value Map per axis.
ax = axisList();
o = struct();
for i = 1:numel(ax)
    o.(ax{i}) = containers.Map('KeyType','char','ValueType','char');
end
end

% =====================================================================
function labels = derive(models, patterns)
%derive  Regexp-match every axis against each file NAME.
if nargin<2 || isempty(patterns), patterns = emptyPatterns(); end
[names, paths] = namesAndPaths(models);
n  = numel(names);
ax = axisList();

labels = struct('n',n,'path',{paths});
for i = 1:numel(ax)
    a   = ax{i};
    pat = '';
    if isstruct(patterns) && isfield(patterns,a), pat = char(patterns.(a)); end
    dflt = defaultToken(a);
    v = repmat({dflt},1,n);
    if ~isempty(pat) && n > 0
        try
            tok = regexp(names, pat, 'match', 'once');
        catch ME
            error('wbTypeModel:badRegexp', ...
                'The %s pattern "%s" is not a valid regular expression:\n%s', a, pat, ME.message);
        end
        hit = ~cellfun(@isempty, tok);
        v(hit) = tok(hit);                       % no match -> stays in the bucket
    end
    labels.(a) = v;
end
end

% =====================================================================
function labels = applyOverrides(labels, varargin)
%applyOverrides  Hand assignments beat the regexp, per axis, keyed by path.
%   The alignment key lives in labels.path, so the models are optional - both
%   applyOverrides(labels,overrides) and applyOverrides(labels,models,overrides)
%   are accepted (the models argument is ignored).
overrides = varargin{end};
if isempty(overrides) || ~isstruct(overrides), return; end
ax = axisList();
for i = 1:numel(ax)
    a = ax{i};
    if ~isfield(overrides,a), continue; end
    m = overrides.(a);
    if ~isa(m,'containers.Map') || m.Count==0, continue; end
    for k = 1:labels.n
        if isKey(m, labels.path{k}), labels.(a){k} = char(m(labels.path{k})); end
    end
end
end

% =====================================================================
function v = axisValues(labels, axis)
%axisValues  The ordered unique values of one axis (first-appearance order).
if ~isfield(labels,axis), v = {}; return; end
v = unique(labels.(axis),'stable');
v = reshape(v,1,[]);
end

% =====================================================================
function [names, paths] = namesAndPaths(models)
%namesAndPaths  Accept a wbFileModel array or a plain cellstr of paths.
if isempty(models), names = {}; paths = {}; return; end
if iscell(models)
    paths = reshape(cellfun(@char, models, 'UniformOutput', false), 1, []);
    names = cell(1,numel(paths));
    for i = 1:numel(paths)
        [~,nm,ex] = fileparts(paths{i});
        names{i} = [nm ex];
    end
else
    paths = reshape({models.path},1,[]);
    names = reshape({models.name},1,[]);
end
end
