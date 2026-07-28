%wbSession - Save / load a Processing-Workbench session sidecar (.mat).
%
%   Round-trips the workbench's in-memory session state to a plain .mat file so
%   a long processing job survives a MATLAB restart.  The session captures the
%   loaded file set (as the discovery grid, so it re-derives its disk baseline
%   on load), the group/order view, the settings model (bag + step/file
%   overrides), the preset reference, and the per-cell overlay (which cells the
%   user checked and which an in-session edit pushed stale).  Nothing on disk is
%   re-derived here - wbSession only (de)serialises; the caller re-runs
%   wbStateEngine after load to overlay saved state on the fresh disk picture.
%
%   The three containers.Map members (fileOverrides, checked, staleOverlay) are
%   flattened to parallel key/value cell arrays so the file is Map-free on disk
%   (the same trick wbSettingsModel uses for its preset).
%
% Syntax:
%    wbSession('save', path, session)      % write the sidecar
%    session = wbSession('load', path)     % read it back (Maps rebuilt)
%
% Inputs:
%    path    - full path of the .mat sidecar.
%    session - struct with fields:
%       fNames        2-D cell of char paths (the discovery grid).
%       referenceMode logical (column 1 is a reference/template).
%       groupNames    cellstr, one per group row.
%       modality      char, the active modality driving the columns.
%       rowOrder      double vector, the display order of the flat rows ([]=natural).
%       bag           struct, the shared settings bag (wbSettingsModel).
%       stepOverrides struct of structs, per-step overrides.
%       fileOverrides containers.Map 'identity||stepId' -> struct.
%       checked       containers.Map 'identity||stepId' -> logical true.
%       staleOverlay  containers.Map 'identity||stepId' -> logical true.
%       presetRef     char, path of the last-used preset ('' if none).
%
% Outputs:
%    session - the same struct shape, with the three Maps reconstructed.
%
% See also: wbSettingsModel, wbStateEngine, wbDiscoverFiles, guiWorkbench
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 28-July-2026

%------------- BEGIN CODE --------------
function out = wbSession(action, varargin)

switch action
    case 'save', saveSession(varargin{:}); out = [];
    case 'load', out = loadSession(varargin{:});
    otherwise,   error('wbSession:badAction','Unknown action "%s".',action);
end
end

% =====================================================================
function saveSession(pth, session)
%saveSession  Flatten the Maps and write a plain struct as 'wbSessionData'.
wbSessionData = struct();
wbSessionData.fNames        = session.fNames;
wbSessionData.referenceMode = session.referenceMode;
wbSessionData.groupNames    = session.groupNames;
wbSessionData.modality      = session.modality;
wbSessionData.rowOrder      = session.rowOrder;
wbSessionData.bag           = session.bag;
wbSessionData.stepOverrides = session.stepOverrides;
wbSessionData.presetRef     = session.presetRef;

[wbSessionData.foKeys,     wbSessionData.foVals]     = mapToCells(session.fileOverrides);
[wbSessionData.chkKeys,    ~]                        = mapToCells(session.checked);
[wbSessionData.staleKeys,  ~]                        = mapToCells(session.staleOverlay);

save(pth,'wbSessionData','-v7.3');
end

% =====================================================================
function session = loadSession(pth)
%loadSession  Read the sidecar and rebuild the containers.Map members.
S = load(pth,'wbSessionData');
p = S.wbSessionData;

session = struct();
session.fNames        = p.fNames;
session.referenceMode = p.referenceMode;
session.groupNames    = p.groupNames;
session.modality      = p.modality;
session.rowOrder      = p.rowOrder;
session.bag           = p.bag;
session.stepOverrides = p.stepOverrides;
session.presetRef     = p.presetRef;

session.fileOverrides = cellsToMap(p.foKeys, p.foVals);
session.checked       = cellsToMap(p.chkKeys, trueVals(p.chkKeys));
session.staleOverlay  = cellsToMap(p.staleKeys, trueVals(p.staleKeys));
end

% =====================================================================
function [k, v] = mapToCells(m)
%mapToCells  containers.Map -> parallel cell arrays (Map-free on disk).
if isa(m,'containers.Map') && m.Count > 0
    k = keys(m);
    v = values(m);
else
    k = {};
    v = {};
end
end

% =====================================================================
function m = cellsToMap(k, v)
%cellsToMap  Parallel cell arrays -> containers.Map (char keys, any values).
m = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(k)
    m(k{i}) = v{i};
end
end

% =====================================================================
function v = trueVals(k)
%trueVals  A logical-true value for each key (checked/stale are boolean sets).
v = num2cell(true(1,numel(k)));
end
