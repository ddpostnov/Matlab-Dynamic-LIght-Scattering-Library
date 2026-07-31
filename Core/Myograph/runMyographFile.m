%runMyographFile  Headless end-to-end processing of one myograph video
%
%   out = runMyographFile(s,fName) runs the whole analysis for a single
%   recording and (by default) writes <video>_myograph.mat.  It is the headless
%   core the GUI is a thin layer over: it only orchestrates the existing
%   functions, so every step stays independently callable.
%
%   PIPELINE
%     1. Diameter.  Default is "detect once, then slice" (Path B): the diameter is
%        measured over the whole recording or the time crop with a SINGLE
%        post-processing pass (getMyographDiameter with no intervals), and any
%        predefined s.intervals then act as LABELS via cutMyographIntervals.  Set
%        s.detectPerInterval=true for the old per-interval detection (Path A).
%     2. Vasomotion.  getMyographVasomotion per interval (s.runVasomotion, default true).
%     3. Propagation. getMyographPropagation per interval (s.runPropagation, default true).
%     4. Save.  intervals + s + meta to a *_myograph.mat (-v7.3), s.saveResult.
%
%   Per-file failure is isolated: any error is caught and returned as a failed
%   status (so a batch can continue).  s.progressFcn(f,total) receives frame
%   progress; s.stageFcn(stage,detail) receives stage messages (pool startup,
%   per-interval progress) for the GUI.
%
%   INPUTS
%     s      parameter struct - all getMyographDiameter / getMyographVasomotion
%            / getMyographPropagation fields, plus (defaults filled):
%              • intervals/intervalNames  predefined windows (labels in Path B)
%              • timeCrop        [t0 t1] s ([] = whole recording)
%              • pixelSize       µm per px ([] or 0 = uncalibrated -> px)
%              • detectPerInterval  false (Path B) | true (Path A)
%              • runVasomotion / runPropagation  (default true)
%              • saveResult      write the .mat (default true)
%              • outPath         output file ([] -> <video>_myograph.mat)
%              • progressFcn / stageFcn  GUI callbacks ([])
%     fName  path to the video.
%
%   OUTPUT
%     out    struct: status ('done'|'failed'), message, path, nIntervals,
%            meta, timings (detect/vasomotion/propagation/save, s), intervals.
%
%   DEPENDS ON  getMyographDiameter, cutMyographIntervals,
%   getMyographVasomotion, getMyographPropagation (all in Core/Myograph/).
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

function out = runMyographFile(s,fName)

def.intervals=[]; def.intervalNames={}; def.timeCrop=[]; def.pixelSize=[];
def.detectPerInterval=false; def.runVasomotion=true; def.runPropagation=true;
def.saveResult=true; def.outPath=[]; def.progressFcn=[]; def.stageFcn=[];
fnd=fieldnames(def);
for i=1:numel(fnd)
    if ~isfield(s,fnd{i}) || (isempty(s.(fnd{i})) && ~islogical(def.(fnd{i}))), s.(fnd{i})=def.(fnd{i}); end
end

out=struct('status','failed','message','','path','','nIntervals',0, ...
    'meta',struct(),'timings',struct(),'intervals',[]);
tim=struct('detect',0,'vasomotion',0,'propagation',0,'save',0);

try
    % ---- 1. diameter ----
    myoStage(s,'detect',sprintf('measuring diameter (%s)',ternary(s.detectPerInterval,'per interval','whole/crop')));
    t0=tic;
    if s.detectPerInterval && ~isempty(s.intervals)
        intervals=getMyographDiameter(s,fName);                  % Path A
    else
        sB=s; sB.intervals=[]; sB.intervalNames={};
        big=getMyographDiameter(sB,fName);                       % detect once (whole/crop)
        if ~isempty(s.intervals)
            intervals=cutMyographIntervals(big(1),s.intervals,s.intervalNames); % labels
        else
            intervals=big;
        end
    end
    tim.detect=toc(t0);
    nIv=numel(intervals);

    % ---- 2+3. vasomotion + propagation (shared with the GUI's Vasomotion stage) ----
    if s.runVasomotion || s.runPropagation
        myoStage(s,'vasomotion',sprintf('%d interval(s) - starting parallel pool if needed',nIv));
        [intervals,at]=analyzeMyographIntervals(s,intervals,s.stageFcn);
        tim.vasomotion=at.vasomotion; tim.propagation=at.propagation;
    end

    % ---- meta + save ----
    meta=buildMeta(s,fName,intervals);
    outPath=s.outPath;
    if isempty(outPath)
        [pp,nn]=fileparts(char(fName)); outPath=fullfile(pp,[nn,'_myograph.mat']);
    end
    if s.saveResult
        myoStage(s,'save',outPath);
        s=stripFcnHandles(s);                            % drop GUI callbacks so no figure is serialised into the .mat
        t0=tic; save(outPath,'intervals','s','meta','-v7.3'); tim.save=toc(t0);
    end

    out.status='done'; out.message='ok'; out.path=outPath;
    out.nIntervals=nIv; out.meta=meta; out.timings=tim; out.intervals=intervals;
catch ME
    if contains(lower(ME.identifier),'cancel') || contains(lower(ME.message),'stopped by user')
        out.status='cancelled'; out.message='Stopped by user.';
    else
        out.status='failed'; out.message=ME.message;
    end
    out.timings=tim; myoStage(s,ternary(strcmp(out.status,'cancelled'),'cancel','error'),out.message);
end
end

% =====================================================================
function meta=buildMeta(s,fName,intervals)
meta.formatVersion=2;
meta.fName=char(fName);
meta.frameRate=NaN;
if ~isempty(intervals) && numel(intervals(1).time)>1
    meta.frameRate=1/median(diff(double(intervals(1).time)),'omitnan');
end
meta.pixelSize=s.pixelSize;                          % []/0 -> uncalibrated (px)
meta.timeCrop=s.timeCrop;
if isfield(s,'rowRange'), meta.rowRange=s.rowRange; else, meta.rowRange=[1 Inf]; end
meta.detectPerInterval=s.detectPerInterval;
meta.codeVersion=codeVersion();
meta.createdTimestamp=datestr(now,'yyyy-mm-ddTHH:MM:SS'); %#ok<TNOW1,DATST>
end

% =====================================================================
function v=codeVersion()
%CODEVERSION  short git hash of the Myograph folder if available, else 'unknown'
v='unknown';
try
    here=fileparts(mfilename('fullpath'));
    [st,o]=system(sprintf('git -C "%s" rev-parse --short HEAD',here));
    if st==0 && ~isempty(strtrim(o)), v=strtrim(o); end
catch %#ok<CTCH>
end
end

% =====================================================================
function myoStage(s,stage,detail)
if isfield(s,'stageFcn') && ~isempty(s.stageFcn)
    try, s.stageFcn(stage,detail); catch, end
else
    fprintf('[%s] %s\n',stage,detail);
end
end
function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
