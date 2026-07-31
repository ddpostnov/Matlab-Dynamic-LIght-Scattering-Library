%reportOpen - Open the reporting context for ONE wrapper call.
%
%   Reporting used to have no owner: every wrapper carried its own private copy of
%   the hook-resolving helper, so the contract had 36 definitions and no single
%   place to change it.  reportOpen is that place.  It resolves the two optional
%   workbench callbacks that may ride in s (s.stageFcn, s.cancelFcn) to no-ops
%   when they are absent, starts the image ledger and the call timer, and hands
%   back a context that every later report* call is passed.
%
%   THREE LINES PER FILE, AND NOTHING ELSE.  A wrapper says what it started, that
%   it is writing, and what it finished:
%
%       Starting Contrast on Mouse1.rls
%       Writing results
%       Finished Contrast on Mouse1.rls, time elapsed: 12.4 s
%
%   There is no progress axis and no stage detail.  A long loop is silent between
%   its two lines - that is the contract, not an omission: the command window of a
%   sixty-file batch has to stay readable, and a percentage that rewrites itself
%   cannot be read back afterwards anyway.  ONLY WRAPPERS REPORT; a core function
%   neither prints nor receives a reporting argument.
%
%   MUTABILITY WITHOUT A CLASSDEF.  rep is a plain struct, but the parts of it that
%   change during a call - the image ledger, the current file, that file's timer -
%   live in a containers.Map held inside it.  A Map is a handle, so reportFile and
%   reportSave mutate it in place and no call site has to write
%   rep = reportSave(rep,...).
%
% Syntax:
%    rep = reportOpen(s, procLabel, fNames)
%
% Inputs:
%    s         - the wrapper's parameter struct.  Optional fields, all defaulted:
%                  .stageFcn(stage,detail)   host log hook
%                  .cancelFcn()->tf          host cancel hook
%                  .reportLevel      'normal' (default) or 'quiet' (nothing is
%                                    printed and nothing reaches the host log)
%    procLabel - human procedure name for this call ('Contrast', 'Vasomotion'),
%                NEVER a function name; it is what the operator reads, and it is
%                the same string in the Starting and the Finished line.
%    fNames    - the file list the call was given (cellstr, string array, or the
%                2-D per-animal cell whose ROWS are files).
%
% Outputs:
%    rep - the reporting context: .procLabel, .nFiles, .level, the two resolved
%          hooks .stageFcn / .cancelFcn, the .emit(text) sink and the mutable
%          .state Map.
%
% Example:
%    rep = reportOpen(s,'Contrast',fNames);
%    for fidx = 1:numel(fNames)
%        if reportCancelled(rep), break; end
%        reportFile(rep,fidx,fNames{fidx});
%        ...
%        reportWriting(rep);
%        save(dName,'source','-v7.3');
%        reportSaved(rep);
%    end
%    reportClose(rep);
%
% See also: reportFile, reportWriting, reportSaved, reportFigure, reportSave,
%           reportClose, reportSettings, reportCancelled
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 31-July-2026

%------------- BEGIN CODE --------------
function rep = reportOpen(s, procLabel, fNames)

if nargin<2, procLabel = ''; end
if nargin<3, fNames = {}; end

% ---- the hook seam, resolved ONCE -----------------------------------------
% Character for character what the 18 private resolveHooks copies did, so a
% hook-free call behaves and saves exactly as it did before this module existed.
stageFcn=@(varargin)[]; cancelFcn=@()false;
if isfield(s,'stageFcn')  &&~isempty(s.stageFcn),   stageFcn  =s.stageFcn;   end
if isfield(s,'cancelFcn') &&~isempty(s.cancelFcn),  cancelFcn =s.cancelFcn;  end

% ---- knobs ----------------------------------------------------------------
level = 'normal';
if isfield(s,'reportLevel') && ~isempty(s.reportLevel)
    v = lower(char(string(s.reportLevel)));
    if any(strcmp(v,{'normal','quiet'})), level = v; end
end

rep             = struct();
rep.procLabel   = char(string(procLabel));
rep.nFiles      = countFiles(fNames);
rep.level       = level;
rep.stageFcn    = stageFcn;
rep.cancelFcn   = cancelFcn;

% ---- the mutable half ------------------------------------------------------
st = containers.Map('KeyType','char','ValueType','any');
st('t0')     = tic;          % the call timer, read by reportClose
st('images') = cell(1,0);    % the ledger, in page order
st('fIdx')   = 0;
st('fName')  = '';
rep.state    = st;

% The routing policy is captured BY VALUE, so emit stays valid even though it is
% built before rep is finished.
pol      = struct('level',level,'tag',rep.procLabel,'stageFcn',stageFcn);
rep.emit = @(txt) emitLine(pol,txt);
end

% =====================================================================
function emitLine(pol,txt)
%emitLine  The ONE sink.  Every line this module emits identifies a recording, so
%   every line belongs in both places: the command window a launcher run is
%   watched in, and the host log the workbench paints beside its matrix.  There is
%   no class table any more because there is no longer more than one kind of line.
if strcmp(pol.level,'quiet'), return; end
txt = char(txt);
disp(txt);
pol.stageFcn(pol.tag,txt);
end

% =====================================================================
function n = countFiles(fNames)
%countFiles  How many files this call was given.  A per-animal wrapper is handed a
%   2-D cell whose ROWS are the files (column 1 is the template), so numel would
%   over-count it.
if isempty(fNames), n = 0; return; end
if ischar(fNames), n = 1; return; end
if ~isvector(fNames) && size(fNames,2)>1
    n = size(fNames,1);
else
    n = numel(fNames);
end
end
