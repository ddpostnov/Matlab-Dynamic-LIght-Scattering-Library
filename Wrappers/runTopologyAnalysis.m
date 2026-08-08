%runTopologyAnalysis - How much vessel is in a field, how wide it is, and how it is spread
%
%   For every *_I_r.mat product that runSegmentation has already been over, this writes
%   six published metrics of the vasculature - vessel area fraction, length density,
%   junction density, the calibre distribution, tortuosity and the distance from tissue
%   to its nearest vessel - as one table row per scope: the whole analysed area first,
%   then each region that was drawn.
%
%   IT OPENS NO DATA CUBE, and that is worth saying because every neighbouring wrapper
%   does.  Everything it needs is the segmentation seam: results.cMask, results.sMetrics
%   and results.regionsMask, plus the two labelling transients re-derived from cMask by
%   getSegmentationLabels exactly as runDynamicSegmentation re-derives them.
%   getSegmentationLabels is deterministic, so those are the same values the static step
%   consumed - AS LONG AS THE LABELLING SETTINGS MATCH THE SEGMENTATION RUN.  sMinL,
%   prchNSize, correctNodes, simR and difR must be the ones runSegmentation was given;
%   they are their own settings group for that reason, and the requirement is the same one
%   runDynamicSegmentation carries.
%
%   THE STEP IS FREE AND THE SEGMENTATION IT READS IS THE WHOLE BILL.  Measured on a
%   1300 x 1300 angiogram: getPixelCategories 195 s, getSegmentationLabels 181 s, and the
%   metrics themselves seconds.  That is the argument for reading what is on the product
%   rather than recomputing anything.
%
%   EVERY NUMBER HERE IS A WITHIN-STUDY COMPARATOR, NOT AN ABSOLUTE, and the report page
%   says so in the reader's own words.  Two reasons, both measured.
%
%   FIRST, THE SEGMENTATION SLIDERS MOVE THESE NUMBERS FURTHER THAN TWO DIFFERENT FIELDS
%   DO.  Across a 25-point sens x sSizeN grid, vessel area fraction moves 1.51x, junction
%   density 2.96x and tortuosity only 1.02x - while two healthy cortical fields from two
%   animals differ by 1.004x, 1.22x and 1.006x.  So comparing recordings on these numbers
%   requires that they were segmented the same way.  sens and sSizeN are sharedKeys of the
%   segmentation step, so the workbench CAN carry one value across a working set, and
%   nothing enforces it - therefore this step reads what each file recorded and says so on
%   the page and in the log when the batch disagrees.
%
%   SECOND, A WIDE-FIELD ANGIOGRAM IS A 2-D PROJECTION OF A 3-D VASCULATURE.  A crossing
%   at two depths is indistinguishable from a bifurcation, so the junction count is not a
%   branching count; overlapping vessels skeletonise to one line, so length density falls
%   short where the network is densest; and projection can only shrink a distance while
%   the capillary bed is not resolved at this pixel size, so the tissue-to-vessel distance
%   is not an oxygen supply distance.  Those three warnings are printed on the page, not
%   only in a settings tooltip, because a tooltip is not where the reader of a figure
%   looks.
%
%   THE SIX ARE NOT EQUALLY TRUSTWORTHY, and the page is laid out accordingly.  Ranked by
%   how much survives contact with this preparation the order is tortuosity, calibre,
%   length density, area fraction, extravascular distance, junction density - so the two
%   shape statistics sit beside the area fraction in the title rather than below it.
%
%   INPUTS
%     s        parameter structure.
%              Labelling (read by getSegmentationLabels, must match the runSegmentation
%                call): sMinL, prchNSize, correctNodes, simR, difR.
%              Calibration: pixelSize (micrometres per pixel; empty or 0 reports
%                everything in pixels).
%              What is reported: evdBeyond, calibreEdges, evdEdges, keepMaps.
%              Machine: parforSegmentationLabels.
%     fNames   cell array of *_I_r.mat paths already processed by runSegmentation (each
%              with a matching *_s.mat sibling; the *_d.mat cube is never read).
%     Optional workbench hooks in s (no-op when absent): s.stageFcn(stage,detail),
%     s.cancelFcn()->tf.
%
%   OUTPUT SIDE-EFFECTS (per file)
%     <name>_r.mat   results.topology.{metrics,calibre,evd,units} and, only when
%                    s.keepMaps, results.topology.evdMap
%     <name>_s.mat   settings.runTopologyAnalysis = s
%     <name>_rep_topology.jpg   the vessel mask on the angiogram, the distance to the
%                    nearest vessel, and the two distributions
%
%   THE CALIBRATION IS A SETTING, NOT THE PRODUCT'S RECORD.  runIntensityBolus writes
%   results.pixelSize when the operator answered it, and this step does NOT read it: the
%   workbench carries one answer here through sharedKeys, and leaving s.pixelSize empty
%   means "report in pixels", which is a live choice rather than a missing value.
%
% Syntax:
%    runTopologyAnalysis(s,fNames)
%
% Example:
%    s.libraryFolder=libraryFolder;
%    s.pixelSize=2.5;                  % or [] to report in pixels
%    fNames = getFileNamesList(rootFolder,'*_a_I_r.mat');
%    runSegmentation(s,fNames);
%    runTopologyAnalysis(s,fNames);
%
%   DEPENDS ON
%     getTopologyMetrics, getSegmentationLabels (Core) + the Core Reporting module and
%     MATLAB's Image Processing Toolbox (visboundaries).
%
% See also: getTopologyMetrics, runSegmentation, runDynamicSegmentation,
%           getSegmentationLabels, setVascularTree
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026


%%Example of s structure parametrisation
% s.libraryFolder=libraryFolder;
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - LABELLING (MUST MATCH runSegmentation)
% s.sMinL=15;         % minimum segment length
% s.prchNSize=50;     % parenchymal cell neighbourhood
% s.correctNodes=true; s.simR=0.3; s.difR=0.4;  % branch-node correction
% %ADJUSTED (OR VERIFIED) PER PROTOCOL - CALIBRATION AND REPORTING
% s.pixelSize=2.5;    % micrometres per pixel.  Leave empty or 0 for an uncalibrated
%                     % recording - results are then reported in px
% s.evdBeyond=25;     % the share of tissue further than this from the nearest vessel
%                     % is reported.  25 um is provisional: at 50 um the fraction is
%                     % 5e-5 on real fields, i.e. identically zero
% s.calibreEdges=0:5:100; % edges of the vessel-width histogram, um
% s.evdEdges=0:2:60;      % edges of the tissue-distance histogram, um
% s.keepMaps=false;   % also store the picture of the distance to the nearest vessel

%------------- BEGIN CODE --------------
function runTopologyAnalysis(s,fNames)

if ~all( cellfun(@(x) isempty(x) || contains(x,'_I_r.mat'), fNames(:)) )
    error(['One or more *non-empty* entries do not contain "_I_r.mat".  ' ...
        'Every step takes the RESULTS member of a product - list them with ' ...
        'getFileNamesList(rootFolder,''*_I_r.mat'').']);
end
if ~isfield(s,'pixelSize'), s.pixelSize=[];   end
if ~isfield(s,'keepMaps') || isempty(s.keepMaps), s.keepMaps=false; end
% Resolved here, not only in the core, so the choice is RECORDED in the saved settings
% like every other tunable - and it is the SAME switch the two segmentation steps carry,
% because all three drive getSegmentationLabels and a shared machine is one tick.
if ~isfield(s,'parforSegmentationLabels') || isempty(s.parforSegmentationLabels)
    s.parforSegmentationLabels=true;
end

% --- what each file recorded about its own segmentation, BEFORE any work -------------
% Judged from the sidecars alone: a file that was never segmented is refused by name here
% rather than on the line that reads edgeSize, and a batch segmented at two different
% settings is named once at the top and again on every page.  The sidecars are kilobytes.
prov = readSegmentationProvenance(fNames);

% reportOpen (Core/Reporting) owns the hook seam: the optional workbench callbacks
% are resolved to no-ops when absent and ride in rep.  s is never mutated, and
% reportSettings strips the hooks from the settings before saving.
rep=reportOpen(s,'Vascular density',fNames);
% An ordinary MATLAB warning and NOT a report line: the reporting module is three lines
% per recording by design and this is not a stage.  runCTTH says the same kind of thing
% the same way, so it is not a new idiom - and the page carries it for the reader who
% never sees a command window.
if prov.mixed
    warning('runTopologyAnalysis:mixedSegmentationSettings','%s',prov.warnText);
end

for fidx=1:1:numel(fNames)
    if reportCancelled(rep), break; end          % cooperative cancel between files
    if ~isempty(fNames{fidx})

    s.fName=fNames{fidx};
    reportFile(rep,fidx,s.fName);
    clearvars results settings
    load(getProductPath(s.fName,'s'),'settings');
    load(s.fName,'results');

    % --- re-derive the labelling transients from the persisted seam -------------------
    % The same two arrays runDynamicSegmentation re-derives, from the same call.  Only
    % the centre lines and the junction nodes are wanted here: the merged mask, the
    % vessel-support map and the wall-distance mask belong to steps that fit a wall.
    edgeSize = settings.runSegmentation.edgeSize;
    [~,~,sLines,~,~,~,nodes] = getSegmentationLabels(results.cMask,edgeSize,s);

    % --- the metrics ------------------------------------------------------------------
    % keepMaps is forced ON for the core call and the SETTING decides what is STORED: the
    % report page draws the distance image on every run, and computing the transform twice
    % to avoid carrying one single image would be the more expensive of the two.
    sCore = s;  sCore.keepMaps = true;
    [m,maps] = getTopologyMetrics(results.cMask,results.sMetrics,sLines,nodes, ...
                                  results.regionsMask,sCore);

    results.topology = m;
    if s.keepMaps, results.topology.evdMap = maps.evd; end

    % --- the report page ---------------------------------------------------------------
    writeTopologyPage(rep,results,m,maps.evd,prov,fidx);

    settings.runTopologyAnalysis=reportSettings(s);
    reportWriting(rep);
    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
    save(s.fName,'results','-v7.3','-nocompression');
    reportSaved(rep);
    end
end
reportClose(rep);
end

% =====================================================================
function prov = readSegmentationProvenance(fNames)
%readSegmentationProvenance  The segmentation settings each file recorded, across a batch.
%
%   THIS IS THE ONE THING THAT MAKES THE STEP MEAN ANYTHING ACROSS A STUDY.  R3 measured
%   that no metric here separates two different cortical fields by more than the
%   segmentation sliders alone can move it, so two recordings are comparable only if they
%   were segmented at the same sens and sSizeN.  Both are sharedKeys of the segmentation
%   step, so a workbench session CAN carry one value across a working set - and nothing
%   warns a user who changed them between file 1 and file 20.  This is that warning.
%
%   Reading a settings group a previous step wrote is a statement about a live CHOICE, not
%   a back-compatibility fallback: a product either has been segmented or it has not, and
%   the second is refused by name rather than worked around.
n = numel(fNames);
prov = struct('sens',nan(1,n),'sSizeN',nan(1,n),'angio',{cell(1,n)}, ...
              'mixed',false,'warnText','');
seen = false(1,n);
for i=1:1:n
    f = fNames{i};
    if isempty(f), continue, end
    sName = getProductPath(f,'s');
    if ~isfile(sName)
        error('runTopologyAnalysis:notSegmented', ...
            ['%s has no settings beside it, so it has not been through the ' ...
             'segmentation step.  Run Segmentation on it first.'], shortName(f));
    end
    L = load(sName,'settings');
    if ~isfield(L,'settings') || ~isfield(L.settings,'runSegmentation')
        error('runTopologyAnalysis:notSegmented', ...
            ['%s has not been through the segmentation step, so there is no vessel ' ...
             'mask to measure.  Run Segmentation on it first.'], shortName(f));
    end
    prov.sens(i)   = getNum(L.settings.runSegmentation,'sens');
    prov.sSizeN(i) = getNum(L.settings.runSegmentation,'sSizeN');
    % The stretch averaged into the picture, when a bolus entry step wrote one.  A
    % product from another entry step simply carries no such group, and both states are
    % current - this is not a fallback onto an older spelling.
    if isfield(L.settings,'runIntensityBolus') && ...
            isfield(L.settings.runIntensityBolus,'fAngio')
        prov.angio{i} = L.settings.runIntensityBolus.fAngio;
    end
    seen(i) = true;
end
if nnz(seen)<2, return, end
prov.mixed = numel(uniqueNum(prov.sens(seen)))>1 || numel(uniqueNum(prov.sSizeN(seen)))>1;
if prov.mixed
    prov.warnText = sprintf(['These recordings were not all found the same way: ' ...
        'sensitivity %s and small-vessel scale %s.  How much vessel a picture is ' ...
        'judged to contain moves further with those two than two different animals ' ...
        'differ, so these numbers can be read one recording at a time but not ' ...
        'compared with each other.'], listNum(prov.sens(seen)), listNum(prov.sSizeN(seen)));
end
end

% =====================================================================
function v = getNum(g,f)
%getNum  One numeric setting, or NaN when the field is absent.
v = NaN;
if isfield(g,f) && isscalar(g.(f)) && isnumeric(g.(f)), v = double(g.(f)); end
end

% =====================================================================
function u = uniqueNum(v)
%uniqueNum  Distinct values, counting NaN as one value rather than as many.
u = unique(v(~isnan(v)));
if any(isnan(v)), u = [u NaN]; end
end

% =====================================================================
function t = listNum(v)
%listNum  The distinct values of v, as text a biologist reads.  A value nothing recorded
%   is named as such rather than printed as NaN.
u = uniqueNum(v);
t = strjoin(arrayfun(@(x) numOrUnknown(x), u, 'UniformOutput',false), ' and ');
end

% =====================================================================
function t = numOrUnknown(x)
if isnan(x), t='not recorded'; else, t=num2str(x,'%g'); end
end

% =====================================================================
function writeTopologyPage(rep,results,m,evdMap,prov,fidx)
%writeTopologyPage  <name>_rep_topology.jpg - four tiles and the words that go with them.
%
%   THE TITLE PUTS CALIBRE AND TORTUOSITY BESIDE THE AREA FRACTION, NOT BELOW IT.  Vessel
%   area fraction is the number everyone asks for and it ranks FOURTH of the six: it has
%   the second-largest swing under the segmentation sliders and it took essentially the
%   same value (0.457 against 0.453) on two fields that look nothing alike.  The two
%   shape statistics are the trustworthy ones, and a page that gave area fraction the
%   top-left tile alone would say the opposite.
%
%   AND THE THREE WARNINGS ARE ON THE PAGE, not only in the settings tooltips.  A figure
%   travels; a tooltip does not.
T = m.metrics;  U = m.units;
cMask = results.cMask;

fh=reportFigure(rep,'topology');
try
    tl=tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

    % --- 1. the angiogram, with the mask the metrics were computed on outlined --------
    ax=nexttile(tl);
    img=double(results.imgI);
    imagesc(ax,img);  colormap(ax,'gray');
    lo=prctile(img(:),1); hi=prctile(img(:),99.5);
    if hi>lo, clim(ax,[lo hi]); end
    hold(ax,'on'); visboundaries(ax,cMask>3,'Color','m','LineWidth',0.5); hold(ax,'off');
    axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
    title(ax,sprintf('vessel %.1f%% of the analysed area, %d segments', ...
        100*T.areaFraction(1), round(T.nSegments(1))));

    % --- 2. how far each point of tissue is from the nearest vessel -------------------
    ax=nexttile(tl);
    imagesc(ax,double(evdMap)); axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
    colorbar(ax);
    title(ax,sprintf('distance to the nearest vessel, %s (19 in 20 within %.1f)', ...
        U.lengthUnit, T.evdP95(1)));

    % --- 3. the calibre distribution, LENGTH-weighted ---------------------------------
    ax=nexttile(tl);
    drawHist(ax,m.calibre.edges,m.calibre.counts(1,:));
    xlabel(ax,sprintf('vessel width, %s',U.lengthUnit));
    ylabel(ax,sprintf('vessel length, %s', lengthAxisUnit(U)));
    title(ax,sprintf('width, weighted by length (middle %.1f, spread %.1f %s)', ...
        T.calibreMedian(1), T.calibreIQR(1), U.lengthUnit));

    % --- 4. the tissue-distance distribution, with the threshold marked ---------------
    ax=nexttile(tl);
    drawHist(ax,m.evd.edges,m.evd.counts(1,:));
    hold(ax,'on'); xline(ax,m.evd.beyond,'--','LineWidth',1); hold(ax,'off');
    xlabel(ax,sprintf('distance to the nearest vessel, %s',U.lengthUnit));
    ylabel(ax,sprintf('tissue area, %s',U.areaUnit));
    title(ax,sprintf('%.2f%% of tissue is further than %g %s (%.0f%% of it measurable)', ...
        100*T.evdFractionBeyond(1), m.evd.beyond, U.lengthUnit, 100*T.evdCoverage(1)));

    % THE HEADLINE GOES IN THE SUBTITLE, NOT IN A TITLE.  reportSave writes the
    % recording's name across the top of every page with sgtitle, and sgtitle attaches to
    % whichever grid the figure holds - so a title set here is silently overwritten and
    % the four numbers would vanish without anything looking wrong.
    notes = wrapLines([{headline(T,U)}, pageNotes(prov,fidx)],230);
    subtitle(tl,notes,'FontSize',8,'Interpreter','none');
catch ME
    delete(fh); rethrow(ME);        % reportSave never runs, so the figure goes here
end
reportSave(rep,fh,'topology');      % reportSave titles the page and deletes fh
end

% =====================================================================
function drawHist(ax,edges,counts)
%drawHist  One histogram from bin edges and per-bin totals, drawn at the bin centres.
c = edges(1:end-1)+diff(edges)/2;
bar(ax,c,counts,1);
if numel(edges)>1, xlim(ax,[edges(1) edges(end)]); end
end

% =====================================================================
function t = headline(T,U)
%headline  The four numbers a reader takes away, and which unit they are in.
if isempty(U.pixelSize)
    cal = 'uncalibrated, values in pixels';
else
    cal = sprintf('%g um per px',U.pixelSize);
end
t = sprintf('vessel area %.3f  |  %.1f %s  |  median calibre %.1f %s  |  tortuosity %.2f  |  %s', ...
    T.areaFraction(1), T.lengthDensity(1), U.densityUnit, ...
    T.calibreMedian(1), U.lengthUnit, T.tortuosityMedian(1), cal);
end

% =====================================================================
function c = pageNotes(prov,fidx)
%pageNotes  How the field was found, and the three things the numbers do not mean.
%
%   The three warnings are word for word the settings tooltips, on purpose: a reader who
%   met one in the panel meets the same sentence under the figure.
c = {};
if prov.mixed
    c{end+1} = ['NOT COMPARABLE ACROSS THIS BATCH.  ' prov.warnText];
end
c{end+1} = provenanceLine(prov,fidx);
c{end+1} = ['Junction density counts where centre-lines meet. Two vessels crossing at ' ...
    'different depths look the same as one vessel branching, so compare fields, not animals.'];
c{end+1} = ['Extravascular distance is how far tissue is from the nearest vessel it can ' ...
    'see. Vessels below the surface are counted as if they were in the plane, and the ' ...
    'smallest ones are not resolved at all, so this is not an oxygen supply distance.'];
c{end+1} = ['Length density is total vessel length per unit area. Vessels that overlap ' ...
    'in the picture merge into one line, so this falls short where the vasculature is densest.'];
end

% =====================================================================
function t = provenanceLine(prov,fidx)
%provenanceLine  What this field was found with, and what the picture was averaged over.
%
%   The angiogram span is printed because it is a USER CHOICE on the bolus branch - the
%   stretch of the recording averaged into the picture is drawn by hand - and a shorter
%   span gives a noisier mean image, hence a different segmentation, hence different
%   densities.  A product from another entry step simply carries no such span, and the
%   line then says nothing about one.
bits = {};
if ~isnan(prov.sens(fidx))
    bits{end+1} = sprintf('found at sensitivity %g',prov.sens(fidx));
end
if ~isnan(prov.sSizeN(fidx))
    bits{end+1} = sprintf('small-vessel scale %g',prov.sSizeN(fidx));
end
a = prov.angio{fidx};
if numel(a)==2
    bits{end+1} = sprintf('picture averaged over frames %d to %d',round(a(1)),round(a(2)));
end
if isempty(bits), t=''; else, t = strjoin(bits,'   |   '); end
end

% =====================================================================
function c = wrapLines(c,n)
%wrapLines  Break each line on a word boundary so nothing runs off the page.
%   The canvas is a fixed 1400 px and a text object does not wrap itself: the
%   mixed-settings warning is the longest string on the page by far, and unwrapped it
%   loses BOTH of its ends - which is a warning nobody can read.
out = cell(1,0);
for i=1:1:numel(c)
    t = c{i};
    while numel(t)>n
        k = find(isspace(t(1:n+1)),1,'last');
        if isempty(k), k = n; end
        out{end+1} = strtrim(t(1:k-1));  %#ok<AGROW> - a handful of lines
        t = strtrim(t(k:end));
    end
    out{end+1} = t;                      %#ok<AGROW>
end
c = out;
end

% =====================================================================
function u = lengthAxisUnit(U)
%lengthAxisUnit  The unit the calibre histogram's counts are in: mm when calibrated.
if isempty(U.pixelSize), u='px'; else, u='mm'; end
end

% =====================================================================
function n=shortName(f)
%shortName  Name and extension - what the user sees on disk.
[~,b,e]=fileparts(char(f));
n=[b e];
end
%------------- END OF CODE --------------
