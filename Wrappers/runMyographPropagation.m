%runMyographPropagation  Direction and speed of the diameter changes along the vessel
%
%   runMyographPropagation(s,fNames) asks, for every interval of every *_MYO_r.mat
%   recording in fNames, whether the constriction/dilation TRAVELS along the vessel
%   and, if so, how fast and which way.  fNames is a cell array of *_MYO_r.mat
%   paths; the wrapper iterates them itself, so there is no launcher for-loop.
%
%   The estimator is the lag at maximum cross-correlation and nothing else: each
%   location's diameter fluctuation is correlated against a robust reference built
%   from the coherent locations, and the lag is regressed on the row index (slope ->
%   speed, sign -> direction).  A speed is ALWAYS returned, together with the
%   metrics that say how much to believe it and one plain sentence summarising them,
%   so the reader decides - see getMyographPropagation.
%
%   WHICH DIAMETER IS ANALYSED.  The diameter step measures three (outer wall, wall
%   centre, lumen); s.diameterMeasures names the ones analysed here, wall centre by
%   default.  Results nest per measure, so asking for two measures answers the same
%   question twice rather than mixing them.  It is one protocol decision shared with
%   the vasomotion step.
%
%   IT READS THE FULL PER-LINE ARRAYS, not the interval's averaged trace: a
%   travelling wave is precisely the disagreement BETWEEN the lines, and an average
%   over them has none of it left.  The per-line arrays live once, in source; the
%   interval carries the frame range into them.  The calibration (pixelSize) and the
%   measured rows (rowRange) are read back from source rather than asked for again -
%   they are properties of the recording, fixed when it was measured.
%
%   INPUTS
%     s        parameter structure.  Carried through into
%              settings.runMyographPropagation.  The fields a protocol tunes:
%                .diameterMeasures  which diameters to analyse, default {'mid'}
%                .vFR               vasomotion band [lo hi] Hz - it sizes the lag
%                                   search around the dominant oscillation
%                .detrendSec        window removing slow drift before correlating, s
%                .propMinCoh        how well a location must correlate to be used
%                .propMinRows       fewest locations an estimate may rest on
%                .propNShuffle      surrogate repeats behind the p-value
%     fNames   cell array of *_MYO_r.mat paths.  Empty cells are skipped.
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail)
%     and s.cancelFcn()->tf (checked between files).
%
%   SIDE-EFFECTS (per file)
%     <name>_MYO_r.mat   results.intervals(k).propagation.<measure> - speed,
%                        speedUnit, direction, speedCI, R2, pValue, confidence,
%                        confidenceLevel, confidenceText, metrics, lagByRow, nRows,
%                        domFreq, qualityFlags
%     <name>_MYO_s.mat   settings.runMyographPropagation = s
%     The SOURCE file is not rewritten: this step reads it and adds nothing to it.
%
%   EXAMPLE
%     s.diameterMeasures = {'mid'};
%     s.vFR = [0.05 0.25];
%     D = dir(fullfile(rootFolder,'*_MYO_r.mat'));
%     runMyographPropagation(s, fullfile({D.folder}',{D.name}'));
%
%   DEPENDS ON
%     getMyographPropagation, myographMeasureIndex, myographProduct
%     (Core/Myograph) and Core/Reporting.
%
% See also: runMyographVideo, runMyographDiameter, runMyographVasomotion,
%           getMyographPropagation, myographProduct
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

% %Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - Which diameter, and which oscillation
% s.diameterMeasures={'mid'};  % 'outer' | 'mid' | 'inner', one or several
% s.vFR=[0.05 0.25];           % vasomotion band [lo hi], Hz
% %ADJUSTED IF NECESSARY - Robustness
% s.detrendSec=30;             % window removing slow drift before correlating, s
% s.propMinCoh=0.3;            % how well a location must correlate to be used
% s.propMinRows=20;            % fewest locations an estimate may rest on
% s.propNShuffle=200;          % surrogate repeats behind the p-value

function runMyographPropagation(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_MYO_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_MYO_r.mat".  Every ' ...
        'step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_MYO_r.mat'').']);
end
if ~isfield(s,'diameterMeasures') || isempty(s.diameterMeasures)
    s.diameterMeasures={'mid'};
end

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Propagation',fNames);

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end     % cooperative cancel between files
    if isempty(fNames{fidx}), continue; end
    fName=char(fNames{fidx});
    reportFile(rep,fidx,fName);
    clearvars source results settings

    [source,results,settings]=myographProduct('open',fName);
    if isempty(source.data)
        error('runMyographPropagation:noDiameter', ...
            'No diameter has been measured for %s yet - run the Diameter step first.',fName);
    end

    % the recording's own properties, not a protocol choice: they were fixed when
    % the diameter was measured
    sFile=s; sFile.fName=fName;
    sCore=sFile;
    sCore.rowRange=fieldOr(source,'rowRange',[1 Inf]);
    sCore.pixelSize=fieldOr(source,'pixelSize',[]);

    meas=cellstr(s.diameterMeasures);
    for k=1:1:numel(results.intervals)
        fr=frameRange(results.intervals(k));
        if isempty(fr), continue; end
        p=struct();
        for j=1:1:numel(meas)
            m=myographMeasureIndex(meas{j},source.measures);
            m=min(m,size(source.data,3));
            p.(meas{j})=trim(getMyographPropagation(sCore, ...
                double(source.data(fr,:,m)),source.time(fr), ...
                source.mask(fr,:),source.valid(fr)));
        end
        results.intervals(k).propagation=p;
    end

    settings.runMyographPropagation=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,[],results,settings);   % SOURCE untouched
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function prop=trim(prop)
%trim  What is KEPT of an estimate.  getMyographPropagation returns a .diag holding
%   the whole [frames x lines] fluctuation matrix it worked on - tens of megabytes
%   per interval per measure, and re-derivable from source in a second - so the
%   stored tree is the answer and the evidence for it, not the working.
if isfield(prop,'diag'), prop=rmfield(prop,'diag'); end
end

% =====================================================================
function fr=frameRange(iv)
%frameRange  The interval's own frames into source, or [] when it has none.
fr=[];
if ~isfield(iv,'frames') || numel(iv.frames)~=2, return; end
if iv.frames(2)<iv.frames(1), return; end
fr=iv.frames(1):iv.frames(2);
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
