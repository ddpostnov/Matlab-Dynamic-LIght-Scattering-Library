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
%   IT READS THE PER-LINE ARRAYS, not the interval's averaged trace: a travelling
%   wave is precisely the disagreement BETWEEN the lines, and an average over them
%   has none of it left.  They live in the window itself,
%   results.intervals(k).diameter.lines.<measure>, so this step reads exactly the
%   window it is answering about; the calibration (pixelSize) is read back from
%   results.recording rather than asked for again, because it is a property of the
%   recording and was fixed when it was measured.
%
%   THE LOCATIONS ARE THE MEASURED ROWS, and that is a change of answer as well as
%   of storage.  This step used to be handed EVERY image row of the recording,
%   including the rows outside the measured band, which carry an interpolated fill
%   and are not measurements of anything; the stored arrays hold the band and
%   nothing else, so those rows no longer reach the estimator.  On a recording whose
%   rowRange is [1 Inf] the two row sets are identical and no number moves.  On one
%   with a narrow band the speed will differ, and it differs because the old answer
%   included rows that were never measured.
%
%   A WINDOW WITH NO ARRAYS CANNOT BE ANSWERED.  runMyographDiameter's s.keepArrays
%   false stores the trace alone and there is no source copy to fall back on, so
%   this step says so and stops rather than reporting a speed derived from one
%   averaged curve.
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
% Last revision: 03-August-2026

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
    clearvars results settings

    [results,settings]=myographProduct('open',fName);
    rec=fieldOr(results,'recording',struct());
    if ~any(hasArrays(results))
        error('runMyographPropagation:noPerLineDiameter', ...
            ['No per-line diameter is stored for %s, so there is nothing to look ' ...
             'for a travelling wave in.  Run the Diameter step on it, with ' ...
             '''Keep every line''''s diameter'' switched on.'],fName);
    end

    % the recording's own properties, not a protocol choice: they were fixed when
    % the diameter was measured
    sFile=s; sFile.fName=fName;
    sCore=sFile;
    % THE STORED ROWS ARE ALREADY THE MEASURED ONES.  getMyographPropagation clips
    % the block it is handed by s.rowRange, counting from column 1 of that block,
    % and each window keeps only the rows inside the recording's band - so passing
    % the band again would apply it a second time and quietly move the speed.
    % [1 Inf] is how "these are all the locations there are" is spelled.  It is also
    % the fix the shape change makes: the wrapper used to hand over EVERY image row,
    % including rows outside the band that carry an interpolated fill and are not
    % measurements of anything.
    sCore.rowRange=[1 Inf];
    sCore.pixelSize=fieldOr(rec,'pixelSize',[]);

    meas=cellstr(s.diameterMeasures);
    for k=1:1:numel(results.intervals)
        d=fieldOr(results.intervals(k),'diameter',[]);
        if ~isstruct(d) || ~isfield(d,'lines') || ~isfield(d,'measured'), continue; end
        p=struct();
        for j=1:1:numel(meas)
            [~,nm]=myographMeasureIndex(meas{j},measuresOf(rec));
            if ~isfield(d.lines,nm), continue; end
            p.(meas{j})=trim(getMyographPropagation(sCore, ...
                double(d.lines.(nm)),d.time, ...
                logical(d.measured),logical(d.valid)));
        end
        results.intervals(k).propagation=p;
    end

    settings.runMyographPropagation=reportSettings(sFile);

    reportWriting(rep);
    myographProduct('save',fName,results,settings);
    reportSaved(rep);
end
reportClose(rep);
end

% =====================================================================
function prop=trim(prop)
%trim  What is KEPT of an estimate.  getMyographPropagation returns a .diag holding
%   the whole [frames x lines] fluctuation matrix it worked on - tens of megabytes
%   per interval per measure, and re-derivable from the window's own arrays in a
%   second - so the stored tree is the answer and the evidence for it, not the
%   working.
if isfield(prop,'diag'), prop=rmfield(prop,'diag'); end
end

% =====================================================================
function tf=hasArrays(results)
%hasArrays  Which windows carry the per-line diameter this step reads.  A window
%   written with s.keepArrays false has its trace and nothing else, and a travelling
%   wave is precisely the disagreement BETWEEN the lines - an average over them has
%   none of it left, so there is nothing here to fall back on.
tf=false;
ivs=fieldOr(results,'intervals',[]);
if isempty(ivs) || ~isstruct(ivs), return; end
tf=arrayfun(@(iv) isstruct(fieldOr(iv,'diameter',[])) && ...
    isfield(iv.diameter,'lines') && isfield(iv.diameter,'measured'),ivs);
end

% =====================================================================
function meas=measuresOf(rec)
%measuresOf  The measures, and their order.  results.recording is the one authority
%   on it, written by the diameter step in the same breath as the windows.
meas=reshape(cellstr(fieldOr(rec,'measures',{})),1,[]);
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
