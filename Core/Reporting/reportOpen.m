%reportOpen - Open the reporting context for ONE wrapper call.
%
%   Reporting used to have no owner: every wrapper carried its own private copy of
%   the hook-resolving helper, so the contract had 36 definitions and no single
%   place to change it.  reportOpen is that place.  It resolves the three optional
%   workbench callbacks that may ride in s (s.progressFcn, s.stageFcn, s.cancelFcn)
%   to no-ops when they are absent, reads the reporting knobs, starts the image
%   ledger and the call timer, and hands back a context that every later report*
%   call is passed.
%
%   THE SINK DECISION IS MADE HERE AND ONLY HERE.  Every emission carries a class
%   ('file', 'stage', 'progress', 'save', 'report', 'warn', 'error') and this
%   function decides once which classes reach the command window and which are
%   forwarded to the host through s.stageFcn.  No call site ever asks whether a
%   hook exists.
%
%   MUTABILITY WITHOUT A CLASSDEF.  rep is a plain struct, but the parts of it that
%   change during a call - the image ledger, the progress throttle, the current
%   file - live in a containers.Map held inside it.  A Map is a handle, so
%   reportFile and reportSave mutate it in place and no call site has to write
%   rep = reportSave(rep,...).
%
% Syntax:
%    rep = reportOpen(s, procLabel, fNames)
%
% Inputs:
%    s         - the wrapper's parameter struct.  Optional fields, all defaulted:
%                  .progressFcn(frac,label)  host progress hook
%                  .stageFcn(stage,detail)   host log hook
%                  .cancelFcn()->tf          host cancel hook
%                  .reportLevel      'normal' (default) or 'quiet' (warnings and
%                                    errors only)
%                  .reportPdf        assemble a PDF in reportClose (default true)
%                  .reportKeepImages keep the report images after the PDF is
%                                    written (default true)
%    procLabel - human procedure name for this call ('Speckle contrast'), NEVER a
%                function name; it is what the operator reads.
%    fNames    - the file list the call was given (cellstr, string array, or the
%                2-D per-animal cell whose ROWS are files).  Used for the [k/N]
%                banner only.
%
% Outputs:
%    rep - the reporting context: .procLabel, .nFiles, .level, .pdf, .keepImages,
%          the three resolved hooks .progressFcn / .stageFcn / .cancelFcn, the
%          class-routing .emit(class,text), and the mutable .state Map.
%
% Example:
%    rep = reportOpen(s,'Speckle contrast',fNames);
%    for fidx = 1:numel(fNames)
%        if reportCancelled(rep), break; end
%        reportFile(rep,fidx,fNames{fidx});
%        reportStage(rep,'Speckle contrast');
%    end
%    reportClose(rep);
%
% See also: reportFile, reportStage, reportProgress, reportSaved, reportFigure,
%           reportSave, reportClose, reportSettings, reportCancelled
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function rep = reportOpen(s, procLabel, fNames)

if nargin<2, procLabel = ''; end
if nargin<3, fNames = {}; end

% ---- the hook seam, resolved ONCE -----------------------------------------
% Character for character what the 18 private resolveHooks copies did, so a
% hook-free call behaves and saves exactly as it did before this module existed.
progressFcn=@(varargin)[]; stageFcn=@(varargin)[]; cancelFcn=@()false;
if isfield(s,'progressFcn')&&~isempty(s.progressFcn), progressFcn=s.progressFcn; end
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
rep.pdf         = optFlag(s,'reportPdf',true);
rep.keepImages  = optFlag(s,'reportKeepImages',true);
rep.progressFcn = progressFcn;
rep.stageFcn    = stageFcn;
rep.cancelFcn   = cancelFcn;
% True when a GUI is listening, so reportProgress knows whether a tick has a
% display to repaint.  A launcher run must not pay for drawnow at all.
rep.hosted      = isfield(s,'progressFcn') && ~isempty(s.progressFcn);

% ---- the mutable half ------------------------------------------------------
st = containers.Map('KeyType','char','ValueType','any');
st('t0')        = tic;          % the call timer, read by reportProgress/reportClose
st('images')    = cell(1,0);    % the ledger, in page order
st('fIdx')      = 0;
st('fName')     = '';
st('lastTickT') = -inf;         % last progress tick: seconds into the call ...
st('lastTickF') = -inf;         % ... and the fraction it reported
st('nTicks')    = 0;
st('openTick')  = false;        % a \r progress line is on screen, unterminated
rep.state       = st;

% The routing policy is captured BY VALUE, so emit stays valid even though it is
% built before rep is finished.  st rides in it too: a Map is a handle, so this is
% the same object rep.state points at, not a copy.
pol      = struct('level',level,'tag',rep.procLabel,'stageFcn',stageFcn,'state',st);
rep.emit = @(cls,txt) emitLine(pol,cls,txt);
end

% =====================================================================
function emitLine(pol,cls,txt)
%emitLine  The ONE place that decides which class of line reaches which sink.
%   'file' and 'save' identify a recording and belong in both the command window
%   and the host log; 'stage' and 'report' are command-window detail the host log
%   would drown in (it already shows the step and lists the images); 'warn' and
%   'error' go everywhere and survive 'quiet'.  'progress' never comes through
%   here - reportProgress renders it in place and pushes it to the host hook.
%
%   A progress tick leaves its line unterminated so the next tick can overwrite
%   it with \r.  Any ordinary line therefore has to close that line first, or it
%   would land on top of the percentage instead of under it.
if strcmp(pol.level,'quiet') && ~any(strcmp(cls,{'warn','error'})), return; end
st = pol.state;
if st('openTick'), fprintf('\n'); st('openTick') = false; end %#ok<NASGU> Map is a handle
txt = char(txt);
switch cls
    case {'file','save','warn','error'}
        disp(txt); pol.stageFcn(pol.tag,txt);
    otherwise
        disp(txt);
end
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

% =====================================================================
function tf = optFlag(s,name,dflt)
%optFlag  An optional scalar true/false field, absent or empty meaning the default.
tf = dflt;
if isfield(s,name) && ~isempty(s.(name)) && isscalar(s.(name))
    tf = logical(s.(name));
end
end
