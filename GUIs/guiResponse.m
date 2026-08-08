%guiResponse - Look at one recording's response: the trace and the markers on it
%
% WHAT THIS TOOL DOES
%   A single-window viewer for ONE processed recording that carries a stimulus
%   response (runNVC) or a drug response (runFitVasoreactivity).  It draws the
%   measured trace and the markers the analysis measured - where the response rose,
%   how big it was, by when it had been delivered, whether it came back - each one
%   drawn where it happened rather than only listed as a number.
%
%   THERE IS NO MODEL ON THE STIMULUS PAGE.  The author's ruling of 06-Aug-2026
%   retired fitNVC, so a stimulus response is measured and not fitted, and this window
%   shows what was measured: the four windows the markers were read in, the levels,
%   the peak at the sample it actually is, and THE CUMULATIVE CURVE THE TIMES ARE
%   QUANTILES OF.  That second panel is the instrument now - every time on the screen
%   is a level crossing on it - and a viewer that hid it would leave a reader nothing
%   to check a time against.  The drug page is unchanged and still fits.
%
%   IT IS A VIEWER, NOT A SECOND ANALYSIS.  Every curve and every number comes out of
%   the same cores the pipeline used - getNVCMetrics and getNVCConfidence for a
%   stimulus response, fitVasoreactivity for a drug one - called with the settings that
%   recording was actually processed with, which are read from the SETTINGS file beside
%   it.  For one segment of one repetition the confidence numbers are therefore THE
%   PRODUCT'S OWN, to the last digit: every confidence factor is a property of one unit
%   and one repetition, so a one-segment block gives the same answer as the producer's
%   whole-field one.  Nothing here defines a marker, a window or a threshold of its own.
%
%   IT NEVER OPENS THE SOURCE.  Only the RESULTS (*_r.mat) and SETTINGS (*_s.mat)
%   members are read - a few hundred megabytes rather than the flow cube's twenty
%   gigabytes - so a recording opens in seconds and a laptop can browse one.
%
% WHY THE REPETITIONS ARE CUT HERE
%   runNVC stores no epoch-cut traces: [samples x segments x repetitions] is 18.6 GB on
%   a reference recording, so the product keeps the numbers and not the cuttings.  The
%   cut itself is cheap to repeat and completely determined by what IS stored - the
%   epoch starts in the settings and the recording clock in the results - so this window
%   repeats it, by the same nearest-frame rule, and CHECKS THE RESULT against
%   results.nvc.epochStart before it draws anything.  If the settings beside the product
%   describe a different run, that check is what says so.
%
% A REPRESENTATIVE-REPETITION PRODUCT IS A DIFFERENT SHAPE, AND IT IS NOT REFUSED
%   With s.nvcRepresentative the step replaces the recording with the average of its
%   trusted repetitions, in place: there is no epochStart, results.time IS the epoch
%   clock and every marker is [nSeg x 1].  So there is nothing to cut, no repetition to
%   choose and no check to make - one curve, one set of markers - and the title says
%   which kind of product is open.  The one number such a product cannot reproduce
%   exactly is its dropped-frame allowance: nothing per-repetition survives the collapse,
%   so this window reads it as none.
%
% HOW TO USE IT
%   1. "Open recording..." - pick a processed *_BFI_r.mat.  The window offers
%      whichever analyses the file carries.
%   2. Choose the signal (segments, vessels, vessel diameter) and whether to look at
%      the average of all of them or at one.
%   3. For a stimulus response, choose the average of the kept repetitions or a single
%      one - including a repetition that was not kept, which is often the interesting one.
%   4. Read the plot; the numbers for exactly what is drawn are in the panel below the
%      controls, ending with the check that scored lowest, in words.  "Save the plot..."
%      writes the panels to an image or a PDF.
%
% PROGRAMMATIC USE
%   getappdata(fig,'responseAPI') exposes .open .setKind .setSignal .setSegment
%   .setEpoch .setFit .setEach .state .savePlot, so the window can be driven headless
%   (guiResponse(path,'Visible','off')) exactly like guiExplore and guiExport.
%   .state() reports WHAT THE LAST REDRAW ACTUALLY DREW - .kind .signal .title .epoch
%   .collapsed .nCurves .valid .error .m .conf .confMin .trust .weakest .weakestWords
%   .numbers - including the panel text and, when a recording was refused, the sentence
%   saying why.  testRunNVC case `viewer` drives it.
%
% NO LAUNCHER IS NEEDED TO OPEN IT
%   Every time this window opens it puts the library on the MATLAB path itself, from
%   its own location (setLibraryPath), so typing guiResponse in a fresh MATLAB is the
%   whole of the setup.
%
% Syntax:
%    guiResponse                          % open the tool, then pick a recording
%    guiResponse(fName)                   % open it ON a processed *_BFI_r.mat
%    h = guiResponse(fName,'Visible','off')   % headless (programmatic drive / tests)
%
% Inputs:
%    fName - path of any member of one product triplet (*_r.mat, *_s.mat, *_d.mat);
%            the RESULTS and SETTINGS members are the two that are read.
%
% Outputs:
%    h     - the uifigure handle.
%
% Example:
%    guiResponse('C:\Data\a1\LSCI_20250318_1WT_t_BFI_r.mat');
%
% See also: runNVC, editNVCEpochs, getNVCMetrics, getNVCConfidence,
%           runFitVasoreactivity, fitVasoreactivity, guiExplore, guiExport,
%           getProductPath, setLibraryPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 06-August-2026

%------------- BEGIN CODE --------------
function h = guiResponse(varargin)

% ---- THE LIBRARY PATH, SET FROM THIS FILE'S OWN LOCATION -------------------
% A front window is usually the FIRST thing typed into a fresh MATLAB, with no
% launcher run, so it puts the library on the path itself instead of assuming
% somebody already did.  setLibraryPath keeps the .claude tooling copies OFF the
% path - they hold whole checkouts of this library at older commits.
libraryFolder = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(libraryFolder,'Utilities'));
setLibraryPath(libraryFolder);

[fileArg, vis] = parseArgs(varargin);

% ---- the palette of the report pages, so a screen and a page match --------
COL = struct('meas',[0.45 0.45 0.45], 'fit',[0.85 0.20 0.15], ...
    'mark',[0.30 0.55 0.85], 'thin',[0.35 0.35 0.35], ...
    'each',[0.78 0.80 0.86], 'map',[0.55 0.30 0.65], 'cum',[0.20 0.35 0.60]);

% ---- shared application state (seen by every nested function) --------------
app = struct();
app.fName    = '';       % the RESULTS member currently open
app.results  = [];
app.settings = [];
app.kind     = '';       % 'nvc' | 'vasoreactivity'
app.geom     = struct(); % per-signal layouts, built once and reused
app.lastRoot = pwd;
app.shown    = struct(); % what the last redraw actually drew (for the API/tests)
app.tl       = [];       % the tiled layout the panels live in - what Save writes

delete(findall(groot,'Type','figure','Tag','guiResponse'));
fig = uifigure('Name','guiResponse - the response of one recording', ...
    'Color','w','Position',[90 80 1300 780],'Visible',vis,'Tag','guiResponse');

outer = uigridlayout(fig,[1 2],'ColumnWidth',{330,'1x'}, ...
    'Padding',[6 6 6 6],'ColumnSpacing',8,'BackgroundColor','w');

left = uigridlayout(outer,[5 1],'RowHeight',{'fit','fit','fit','fit','1x'}, ...
    'Padding',[0 0 0 0],'RowSpacing',8,'BackgroundColor','w');

c = struct();

% ===================== 1 - the recording ==================================
p1 = uipanel(left,'Title','1 - Recording','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
g1 = uigridlayout(p1,[2 1],'RowHeight',{'fit','fit'},'Padding',[6 6 6 6], ...
    'RowSpacing',4,'BackgroundColor','w');
uibutton(g1,'Text','Open recording...','ButtonPushedFcn',@(~,~)onOpen(), ...
    'Tooltip','pick a processed *_BFI_r.mat that carries a response analysis');
c.fileLbl = uilabel(g1,'Text','(no recording open)','WordWrap','on', ...
    'FontColor',[0.35 0.35 0.35]);

% ===================== 2 - what to show ===================================
p2 = uipanel(left,'Title','2 - What to show','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
g2 = uigridlayout(p2,[5 2],'RowHeight',repmat({'fit'},1,5), ...
    'ColumnWidth',{88,'1x'},'Padding',[6 6 6 6],'RowSpacing',5, ...
    'ColumnSpacing',6,'BackgroundColor','w');

uilabel(g2,'Text','Response');
c.kind = uidropdown(g2,'Items',{'-'},'ItemsData',{''},'Enable','off', ...
    'ValueChangedFcn',@(~,~)onKind());

uilabel(g2,'Text','Signal');
c.signal = uidropdown(g2,'Items',{'-'},'ItemsData',{''},'Enable','off', ...
    'ValueChangedFcn',@(~,~)onSignal());

uilabel(g2,'Text','Trace');
c.scope = uidropdown(g2,'Items',{'Average of all segments','One segment'}, ...
    'ItemsData',{'all','one'},'Enable','off','ValueChangedFcn',@(~,~)onScope());

uilabel(g2,'Text','Segment');
c.segment = uispinner(g2,'Limits',[1 2],'Value',1,'Step',1,'Enable','off', ...
    'RoundFractionalValues','on','ValueChangedFcn',@(~,~)redraw());

c.epochLbl = uilabel(g2,'Text','Repetition');
c.epoch = uidropdown(g2,'Items',{'-'},'ItemsData',0,'Enable','off', ...
    'ValueChangedFcn',@(~,~)redraw());

% ===================== 3 - display ========================================
% THE FIT CHECKBOX BELONGS TO THE DRUG PAGE ALONE.  A stimulus response has no model
% behind it any more, so the control is inert there rather than offering a curve this
% window would have to invent.
p3 = uipanel(left,'Title','3 - Display','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
g3 = uigridlayout(p3,[4 1],'RowHeight',repmat({'fit'},1,4), ...
    'Padding',[6 6 6 6],'RowSpacing',3,'BackgroundColor','w');
c.showFit = uicheckbox(g3,'Text','Fitted drug response','Value',true, ...
    'Enable','off','ValueChangedFcn',@(~,~)redraw(), ...
    'Tooltip','the only slow thing here - switch it off to flick through segments');
c.showMark = uicheckbox(g3,'Text','Measured markers','Value',true, ...
    'Enable','off','ValueChangedFcn',@(~,~)redraw());
c.showEach = uicheckbox(g3,'Text','Every repetition behind it','Value',false, ...
    'Enable','off','ValueChangedFcn',@(~,~)redraw());
c.showMap = uicheckbox(g3,'Text','Arterial pressure','Value',true, ...
    'Enable','off','ValueChangedFcn',@(~,~)redraw());

% ===================== the save button ====================================
g4 = uigridlayout(left,[1 1],'RowHeight',{'fit'},'Padding',[0 0 0 0], ...
    'BackgroundColor','w');
c.saveBtn = uibutton(g4,'Text','Save the plot...','Enable','off', ...
    'ButtonPushedFcn',@(~,~)onSave(),'Tooltip','write the plot to an image or a PDF');

% ===================== 4 - the numbers ====================================
p5 = uipanel(left,'Title','4 - Numbers for the trace shown','FontWeight','bold', ...
    'BackgroundColor','w','ForegroundColor',[0.15 0.25 0.45]);
g5 = uigridlayout(p5,[1 1],'Padding',[6 6 6 6],'BackgroundColor','w');
c.numbers = uitextarea(g5,'Value',{'Open a recording to see its response.'}, ...
    'Editable','off','FontName','Courier New','FontSize',11);

% ===================== the plot ===========================================
% THE PANELS LIVE IN A TILED LAYOUT and not in a grid of uiaxes, and the reason is the
% Save button: exportgraphics writes a tiled layout whole, and refuses a uiaxes inside a
% container with 'UI components will not be included'.  A stimulus response needs two
% panels - the trace and the cumulative the times are read off - so a window that could
% only export one of them would be saving half the evidence.
c.plotP = uipanel(outer,'BackgroundColor','w','BorderType','none');

setappdata(fig,'responseAPI',struct( ...
    'open',      @(p) openFile(p), ...
    'setKind',   @(k) setKind(k), ...
    'setSignal', @(s) setSignal(s), ...
    'setSegment',@(i) setSegment(i), ...
    'setEpoch',  @(e) setEpoch(e), ...
    'setFit',    @(tf)setFit(tf), ...
    'setEach',   @(tf)setEach(tf), ...
    'savePlot',  @(p) savePlot(p), ...
    'state',     @shownState));
%   .state is a handle to a NESTED function and not @() app.shown: an anonymous
%   function captures its variables BY VALUE where it is written, so that spelling
%   would hand every caller the empty state this window started with.

if ~isempty(fileArg), openFile(fileArg); end
if nargout>0, h = fig; end

% =====================================================================
%                            THE FILE
% =====================================================================
    function onOpen()
        [f,d] = uigetfile({'*_r.mat','Processed results (*_r.mat)'}, ...
            'Open a processed recording',app.lastRoot);
        figure(fig);                       % the dialog steals focus on Windows
        if isequal(f,0), return, end
        app.lastRoot = d;
        openFile(fullfile(d,f));
    end

% ---------------------------------------------------------------------
    function ok = openFile(p)
        %openFile  Read the RESULTS and SETTINGS members and rebuild every picker.
        %   THE SETTINGS ARE NOT OPTIONAL.  They carry the windows, the stimulus
        %   geometry and the thresholds this recording was analysed with, and measuring
        %   a trace against a window guessed here would be a second definition of the
        %   analysis rather than a view of it.
        ok = false;
        rp = getProductPath(p,'r');
        sp = getProductPath(p,'s');
        if ~isfile(rp)
            say('%s is not there.',rp); return
        end
        if ~isfile(sp)
            say(['The settings written beside this recording are missing (%s), so ' ...
                 'there is nothing that says which windows it was analysed with.'],sp)
            return
        end

        busy(true);
        cleaner = onCleanup(@() busy(false));
        R = load(rp,'results');
        S = load(sp,'settings');
        if ~isfield(R,'results') || ~isfield(S,'settings')
            say('%s or its settings file holds no results/settings variable.',rp); return
        end

        [kinds,labels,missing] = availableKinds(R.results,S.settings);
        if isempty(kinds)
            say(['%s carries no response analysis.  Run the stimulus response ' ...
                 '(runNVC) or the vasoreactivity fit (runFitVasoreactivity) on it ' ...
                 'first.%s'],rp,missing); return
        end

        app.fName    = rp;
        app.results  = R.results;
        app.settings = S.settings;
        app.geom     = struct();
        [~,base,ext] = fileparts(rp);
        c.fileLbl.Text = [base ext missing];

        %Items and ItemsData are set in ONE call: assigned one at a time they are
        %briefly of different lengths, which the component rejects.
        set(c.kind,'Items',labels,'ItemsData',kinds);
        c.kind.Value  = kinds{1};
        c.kind.Enable = 'on';
        onKind();
        ok = true;
    end

% =====================================================================
%                          THE PICKERS
% =====================================================================
    function onKind()
        app.kind = c.kind.Value;
        app.geom = struct();
        [ids,labels] = signalsFor();
        if isempty(ids)
            say(['This analysis measured signals the recording no longer carries, ' ...
                 'so there is no trace to draw.']);
            set(c.signal,'Items',{'-'},'ItemsData',{''},'Enable','off');
            return
        end
        set(c.signal,'Items',labels,'ItemsData',ids);
        c.signal.Value  = ids{1};
        c.signal.Enable = 'on';
        onSignal();
    end

% ---------------------------------------------------------------------
    function onSignal()
        %onSignal  The segment count and the repetition list follow the chosen signal.
        %a spinner's range has to increase, so a one-segment signal still gets a
        %two-wide range; the draw clamps the index to the segments that exist
        n = size(app.results.(c.signal.Value),2);
        c.segment.Limits = [1 max(n,2)];
        c.segment.Value  = min(c.segment.Value,max(n,1));
        c.scope.Enable   = 'on';
        c.showMark.Enable= 'on';
        c.saveBtn.Enable = 'on';

        isNVC = strcmp(app.kind,'nvc');
        %A COLLAPSED PRODUCT HAS ONE REPETITION AND IT IS ALREADY AN AVERAGE, so there
        %is nothing to choose between and nothing to draw behind it.
        collapsed = isNVC && ~isfield(app.results.nvc,'epochStart');
        c.showFit.Enable  = onOff(~isNVC);
        c.showEach.Enable = onOff(isNVC && ~collapsed);
        c.epoch.Enable    = onOff(isNVC && ~collapsed);
        c.epochLbl.Enable = onOff(isNVC && ~collapsed);
        c.showMap.Enable  = onOff(~isNVC && isfield(app.results.vasoreactivity,'map'));

        if isNVC
            [items,data] = epochItems(collapsed);
            set(c.epoch,'Items',items,'ItemsData',data);
            if ~ismember(c.epoch.Value,data), c.epoch.Value = data(1); end
        end
        onScope();
    end

% ---------------------------------------------------------------------
    function onScope()
        c.segment.Enable = onOff(strcmp(c.scope.Value,'one'));
        redraw();
    end

% ---------------------------------------------------------------------
    function [ids,labels] = signalsFor()
        %signalsFor  The signals this analysis measured AND the recording still
        %   carries.  The names are the explorer's, so one vocabulary describes a
        %   signal everywhere the library shows one.
        switch app.kind
            case 'nvc'
                ids = fieldnames(app.results.nvc.esMetrics);
            case 'vasoreactivity'
                ids = intersect({'sData','dvsData','dvsDiameter'}, ...
                    fieldnames(app.results.vasoreactivity),'stable');
            otherwise
                ids = {};
        end
        keep = false(size(ids));
        for i=1:numel(ids)
            keep(i) = isfield(app.results,ids{i}) && ~isempty(app.results.(ids{i}));
            if keep(i) && strcmp(ids{i},'gsData')
                keep(i) = isfield(app.results,'gsTime') && ~isempty(app.results.gsTime);
            end
        end
        ids = ids(keep);
        labels = cellfun(@signalLabel,ids,'UniformOutput',false);
    end

% ---------------------------------------------------------------------
    function [items,data] = epochItems(collapsed)
        %epochItems  One entry per repetition, saying which ones the recording kept.
        %   READ STRAIGHT OFF THE TREE and never off a layout: this runs before any
        %   geometry is built, so a settings file belonging to another run must still
        %   leave a usable list of repetitions and let the DRAW be the thing that
        %   refuses.
        if collapsed
            items = {'The representative repetition'};
            data  = 1;
            return
        end
        v = logical(app.results.nvc.epochTrust(:))';
        items = {sprintf('Average of the %d repetitions kept',sum(v))};
        data  = 0;
        for k=1:numel(v)
            if v(k), items{end+1} = sprintf('Repetition %d',k);              %#ok<AGROW>
            else,    items{end+1} = sprintf('Repetition %d - not kept',k);   %#ok<AGROW>
            end
            data(end+1) = k;                                                 %#ok<AGROW>
        end
    end

% =====================================================================
%                           THE DRAWING
% =====================================================================
    function redraw()
        %redraw  Rebuild the panels and draw whatever the pickers now say.
        %   THE LAYOUT IS REBUILT rather than cleared: a stimulus page has two panels
        %   and a drug page has one with a second y-axis for the pressure covariate,
        %   and clearing an axes does not reliably take a second ruler with it.
        if isempty(app.results) || isempty(c.signal.Value), return, end
        delete(c.plotP.Children);
        busy(true);
        cleaner = onCleanup(@() busy(false));
        app.shown = struct('kind',app.kind,'signal',c.signal.Value,'title','', ...
            'epoch',0,'collapsed',false,'nCurves',0,'valid',false,'error','', ...
            'm',[],'conf',NaN,'confMin',NaN,'trust',false,'weakest','', ...
            'weakestWords','','numbers',{{}});
        try
            switch app.kind
                case 'nvc',            drawNVC();
                case 'vasoreactivity', drawVRC();
            end
        catch err
            %A VIEWER MUST SAY WHY IT IS EMPTY.  Every refusal above is a written
            %sentence about this recording, so it belongs where the numbers go rather
            %than in the console the user is not looking at.
            app.shown.error = err.message;
            c.numbers.Value = [{'Nothing could be drawn:'};{''}; ...
                cellstr(splitlines(string(err.message)))];
        end
        %THE PANEL TEXT IS PART OF WHAT WAS DRAWN, so it is read back here rather than
        %in each branch: both of them write it, the catch arm writes it too, and a
        %caller asking .state what this window is showing means the numbers as much as
        %the curves.
        app.shown.numbers = cellstr(string(c.numbers.Value));
    end

% ---------------------------------------------------------------------
    function ax = newPanels(n)
        %newPanels  The tiled layout this redraw draws into, and its axes.  One tile for
        %   a drug response, two for a stimulus one - the trace, and the cumulative the
        %   times are quantiles of.
        app.tl = tiledlayout(c.plotP,n,1,'TileSpacing','compact','Padding','compact');
        ax = gobjects(1,n);
        for i=1:n
            ax(i) = nexttile(app.tl);
            hold(ax(i),'on'); grid(ax(i),'on');
        end
    end

% ---------------------------------------------------------------------
    function drawNVC()
        %drawNVC  One repetition, or the average of the kept ones: the trace with its
        %   four windows and its markers, and the cumulative curve underneath.
        sig = c.signal.Value;
        G   = nvcGeom(sig);
        [y,E,who] = nvcSeries(G,sig);
        t   = G.t;
        L   = G.L;
        [m,cf] = nvcNumbers(G,E,y);

        ax = newPanels(2);
        hs = gobjects(0); names = {};

        %THE FOUR WINDOWS, WHERE THEY ACTUALLY ARE - drawn off the layout's own index
        %vectors, so the shading and the marker read the same window.  Same four, in the
        %same colours, as '_rep_nvcresponse.jpg'.
        %
        %EACH ONE IS LABELLED ON ITSELF AND NOT IN THE LEGEND.  Four region entries put a
        %seven-line legend inside the axes, which 'best' then parks on top of the trace -
        %and a window named where it is needs no key at all.  The labels go on after the
        %limits are fixed, because a label placed against a band that has not been sized
        %yet lands wherever the default happened to be.
        reg = {t(find(L.blIdx,1)),   t(find(L.blIdx,1,'last')),  'baseline',    COL.meas,0.08,'center'; ...
               L.tStim,              L.tStim+L.D,                'stimulus',    COL.mark,0.14,'center'; ...
               L.tStim+L.D,          t(find(L.pkIdx,1,'last')),  'peak window', COL.mark,0.06,'left'; ...
               t(find(L.finIdx,1)),  t(find(L.finIdx,1,'last')), 'end window',  COL.meas,0.08,'center'};
        for k = 1:size(reg,1)
            xregion(ax(1),reg{k,1},reg{k,2},'FaceColor',reg{k,4},'FaceAlpha',reg{k,5});
        end

        eachY = [];
        if c.showEach.Value && ~G.collapsed
            h1 = plot(ax(1),t,E,'-','Color',COL.each,'LineWidth',0.5);
            hs(end+1) = h1(1); names{end+1} = 'each repetition';
            eachY = E(:);              % the band the y-limits then have to hold
        end
        hs(end+1) = plot(ax(1),t,y,'-','Color',COL.meas,'LineWidth',1.2);
        names{end+1} = 'measured';

        if c.showMark.Value && m.valid
            [hm,nm] = nvcMarkers(ax(1),m,G,t,y);
            hs = [hs hm]; names = [names nm];
        end

        axis(ax(1),'tight');
        yl = padY(ax(1),y,eachY);
        for k = 1:size(reg,1)
            if strcmp(reg{k,6},'left'), xL = reg{k,1}; else, xL = 0.5*(reg{k,1}+reg{k,2}); end
            text(ax(1),xL,yl(2),reg{k,3},'HorizontalAlignment',reg{k,6}, ...
                'VerticalAlignment','top','FontSize',9,'Color',[0.40 0.40 0.45]);
        end
        hold(ax(1),'off');
        ylabel(ax(1),signalUnit(sig));
        title(ax(1),who);
        %THE KEY GOES OUTSIDE THE BOX, in one row.  Only the curves are in it now, so it
        %is three or four short words, and 'best' inside the axes has nowhere to put even
        %that on a trace that fills its own box - it lands on the data or on the window
        %labels every time.
        legend(ax(1),hs,names,'Location','northoutside','Orientation','horizontal', ...
            'Box','off','FontSize',9);

        nCum = drawCumulative(ax(2),m,G,t,y);

        xlabel(ax(2),'Time in the repetition, s');
        linkaxes(ax,'x');
        xlim(ax(1),[t(1) t(end)]);

        c.numbers.Value = nvcLines(m,cf,G,who);
        app.shown.title        = who;
        app.shown.epoch        = currentEpoch(G);
        app.shown.collapsed    = G.collapsed;
        app.shown.nCurves      = numel(hs)+nCum;
        app.shown.valid        = logical(m.valid);
        app.shown.m            = m;
        app.shown.conf         = cf.conf;
        app.shown.confMin      = cf.confMin;
        app.shown.trust        = cf.trust;
        app.shown.weakest      = cf.weakest;
        app.shown.weakestWords = cf.weakestWords;
    end

% ---------------------------------------------------------------------
    function [hs,names] = nvcMarkers(ax,m,G,t,y)
        %nvcMarkers  Every marker drawn where it was measured.
        %   The times in the bag are quoted from the STIMULUS MARK, which is the clock
        %   the physiology is read in; the axis is the repetition clock, which is the
        %   clock the protocol is written in.  tStim is what converts one to the other,
        %   and it is the layout's, so the mark and the metric cannot disagree.
        hs = gobjects(0); names = {};
        L  = G.L;
        Bl = double(m.Bl);

        hs(end+1) = plot(ax,[t(1) t(end)],[Bl Bl],':','Color',COL.thin,'LineWidth',1);
        names{end+1} = 'baseline';

        tf = t(L.finIdx);
        plot(ax,[tf(1) tf(end)],double(m.Fn)*[1 1],':','Color',COL.thin,'LineWidth',1.5);

        %THE PEAK IS DRAWN AT THE SAMPLE IT ACTUALLY IS.  Its VALUE is the core's - the
        %point sits at Bl + m.Pk, which is that sample's own height - and only the place
        %to put it is found here, by the core's own expression.  Reading it back any
        %other way would be this window inventing a marker.
        sg = sign(double(m.Auc)); if sg==0 || ~isfinite(sg), sg = 1; end
        r = sg*(double(y)-Bl); r(~L.pkIdx) = -Inf;
        [~,iPk] = max(r);
        vp = Bl+double(m.Pk);
        if isfinite(vp)
            plot(ax,[t(iPk) t(iPk)],[Bl vp],'-','Color',COL.fit,'LineWidth',0.75);
            hs(end+1) = plot(ax,t(iPk),vp,'o','MarkerSize',7, ...
                'MarkerFaceColor',COL.fit,'MarkerEdgeColor','none');
            names{end+1} = 'peak';
            if isfinite(m.PkRel)
                text(ax,t(iPk),vp,sprintf('  %+.1f %%',100*double(m.PkRel)), ...
                    'VerticalAlignment','bottom','FontSize',10,'Color',COL.fit);
            end
        end
        %THE TIMING LABELS HANG AT THE BOTTOM, where the crossings themselves are: the
        %top of the axes belongs to the peak and to its size, and however many timing
        %levels the protocol asked for competing for that strip is how a marker stops
        %being readable.  How many there are is READ, never listed.
        for j=1:numel(L.pcts)
            nm = sprintf('T%d',L.pcts(j));
            if ~isfield(m,nm) || ~isfinite(m.(nm)), continue, end
            xline(ax,L.tStim+double(m.(nm)),'--',sprintf('%d %%',L.pcts(j)), ...
                'Color',COL.fit,'LabelOrientation','horizontal', ...
                'LabelVerticalAlignment','bottom','FontSize',9);
        end
    end

% ---------------------------------------------------------------------
    function n = drawCumulative(ax,m,G,t,y)
        %drawCumulative  THE INSTRUMENT THE TIMES ARE READ OFF.  Every T is the moment
        %   this curve passes a percentage of its own maximum, so the panel is what lets
        %   a reader check a time instead of taking it.
        %
        %   IT IS THE CORE'S OWN CONSTRUCTION, line for line: the baseline-referenced
        %   running integral, mirrored by sign(Auc), and levels taken as fractions of its
        %   MAXIMUM rather than of Auc - the two differ wherever the trace dips below
        %   baseline late in the repetition.  The 50 % level is drawn whether or not it
        %   was requested, because it is the anchor every other level is searched out
        %   from and not a marker that happens to be asked for.
        n = 0;
        L = G.L;
        if ~m.valid
            axis(ax,'off'); return
        end
        sg = sign(double(m.Auc)); if sg==0 || ~isfinite(sg), sg = 1; end
        Cb = cumtrapz(t,sg*(double(y)-double(m.Bl)));
        cbMax = max(Cb);

        xregion(ax,L.tStim,L.tStim+L.D,'FaceColor',COL.mark,'FaceAlpha',0.14);
        plot(ax,t,Cb,'-','Color',COL.cum,'LineWidth',1.2); n = n+1;
        if cbMax>0
            yline(ax,0.5*cbMax,':','Color',COL.thin);
            for j=1:numel(L.pcts)
                yline(ax,0.01*L.pcts(j)*cbMax,'--',sprintf('%d %%',L.pcts(j)), ...
                    'Color',COL.fit,'LabelHorizontalAlignment','left','FontSize',9);
                nm = sprintf('T%d',L.pcts(j));
                if isfield(m,nm) && isfinite(m.(nm))
                    xline(ax,L.tStim+double(m.(nm)),'--','Color',COL.fit);
                end
            end
        end
        hold(ax,'off');
        axis(ax,'tight'); padY(ax,Cb);
        ylabel(ax,'Delivered so far');
        title(ax,'The accumulated change the times are read off');
    end

% ---------------------------------------------------------------------
    function drawVRC()
        %drawVRC  The whole recording, its fit, and the drug markers.
        sig = c.signal.Value;
        want = {'markers'};
        if c.showFit.Value, want = {'markers','model','reconstruction'}; end
        G = vrcGeom(want);
        [y,who] = vrcSeries(sig);

        m  = fitVasoreactivity(y,G.L,G.s);
        yd = blockDecimate(y,G.L.avgN);
        t  = G.L.time;

        ax = newPanels(1);
        hs = gobjects(0); names = {};
        hs(end+1) = xregion(ax,G.L.t1+G.L.bl(1),G.L.t1+G.L.bl(2), ...
            'FaceColor',COL.meas,'FaceAlpha',0.10);
        names{end+1} = 'baseline window';

        %THE PRESSURE IS STORED, NEVER SUBTRACTED, and it is drawn beside the trace
        %for the same reason the producer keeps it: autoregulation opposing a
        %drug-driven flow change is a confound that is visible or invisible, and a
        %covariate quietly removed would be worse than one plainly displayed.  It is
        %drawn AFTER the flow trace has its limits, because a second ruler is scaled
        %from its own data and `axis tight` only ever touches the active side.
        vr = app.results.vasoreactivity;
        mapv = [];
        if c.showMap.Value && isfield(vr,'map') && any(isfinite(vr.map))
            mapv = double(vr.map(:));
        end

        hs(end+1) = plot(ax,t,yd,'-','Color',COL.meas,'LineWidth',1.2);
        names{end+1} = 'measured';
        fd = [];
        if c.showFit.Value && isfield(m,'fData') && numel(m.fData)==numel(t)
            fd = double(m.fData);
            hs(end+1) = plot(ax,t,fd,'-','Color',COL.fit,'LineWidth',1.5);
            names{end+1} = 'fitted';
        end

        if c.showMark.Value && m.valid
            [hm,nm] = vrcMarkers(ax,m,G,t);
            hs = [hs hm]; names = [names nm];
        end

        axis(ax,'tight');
        padY(ax,yd,fd);
        if ~isempty(mapv)
            yyaxis(ax,'right');
            hs(end+1) = plot(ax,t,mapv,'-','Color',COL.map,'LineWidth',1);
            names{end+1} = 'arterial pressure';
            padY(ax,mapv);
            ylabel(ax,'Arterial pressure');
            ax.YAxis(2).Color = COL.map;
            yyaxis(ax,'left');
            ax.YAxis(1).Color = 'k';
        end
        hold(ax,'off');
        xlabel(ax,'Time, s');
        ylabel(ax,signalUnit(sig));
        title(ax,who);
        legend(ax,hs,names,'Location','best','Box','off');

        c.numbers.Value = vrcLines(m,G,who);
        app.shown.title   = who;
        app.shown.nCurves = numel(hs);
        app.shown.valid   = logical(m.valid);
        app.shown.m       = m;
    end

% ---------------------------------------------------------------------
    function [hs,names] = vrcMarkers(ax,m,G,t)
        %vrcMarkers  The injection, the baseline level, the peak and the washout.
        %   TMax and TWash are quoted from the INJECTION, so the injection mark is
        %   what puts them back on the recording clock.
        hs = gobjects(0); names = {};
        tI = G.L.tInj;
        Bl = double(m.BlMean);

        xline(ax,tI,'-','injection','Color',COL.mark,'LineWidth',1.5, ...
            'LabelOrientation','horizontal','FontSize',9);
        hs(end+1) = plot(ax,[t(1) t(end)],[Bl Bl],':','Color',COL.thin,'LineWidth',1);
        names{end+1} = 'baseline';

        if isfinite(m.TMax) && isfinite(m.Peak)
            tp = tI+double(m.TMax); vp = Bl+double(m.Peak);
            plot(ax,[tp tp],[Bl vp],'-','Color',COL.fit,'LineWidth',0.75);
            hs(end+1) = plot(ax,tp,vp,'o','MarkerSize',7, ...
                'MarkerFaceColor',COL.fit,'MarkerEdgeColor','none');
            names{end+1} = 'peak';
            if isfinite(m.PeakRel)
                text(ax,tp,vp,sprintf('  %+.1f %%',100*double(m.PeakRel)), ...
                    'VerticalAlignment','bottom','FontSize',10,'Color',COL.fit);
            end
        end
        %as on the stimulus page, the top strip belongs to the peak and its size
        if isfinite(m.TWash)
            xline(ax,tI+double(m.TWash),'--','washout','Color',COL.thin, ...
                'LabelOrientation','horizontal','LabelVerticalAlignment','bottom', ...
                'FontSize',9);
        end
        if isfield(m,'TMaxModel') && isfinite(m.TMaxModel)
            xline(ax,tI+double(m.TMaxModel),':','fitted peak','Color',COL.fit, ...
                'LabelOrientation','horizontal','LabelVerticalAlignment','bottom', ...
                'FontSize',9);
        end
    end

% =====================================================================
%                       WHAT IS BEING MEASURED
% =====================================================================
    function G = nvcGeom(sig)
        %nvcGeom  This signal's cut and its layout, built once and kept.
        %   gsData is stored at the RAW frame rate on its own clock - 100 Hz against
        %   the segment traces' 4 Hz on a reference recording - so it is cut on that
        %   clock and gets its own layout; everything else shares the primary one.
        key = matlab.lang.makeValidName(sig);
        if isfield(app.geom,key), G = app.geom.(key); return, end

        sf = app.settings.runNVC;
        N  = app.results.nvc;
        %A COLLAPSED PRODUCT IS TOLD BY WHAT IT NO LONGER HAS.  Nothing per-repetition
        %survives s.nvcRepresentative, so an absent epochStart is the fact, and results
        %.time IS the repetition clock rather than the recording's.
        G.collapsed = ~isfield(N,'epochStart');
        if strcmp(sig,'gsData'), clk = double(app.results.gsTime(:));
        else,                    clk = double(app.results.time(:));
        end
        G.clk = clk;

        if G.collapsed
            G.cut      = (1:numel(clk))';
            G.t        = clk-clk(1);
            G.trust    = true;
            %NOTHING SURVIVES THE COLLAPSE, INCLUDING THE DROPPED-FRAME ALLOWANCE the
            %average was measured with, so it is read as none.  It is the one number in
            %this window a collapsed product cannot reproduce to the last digit, and on
            %a good recording it is worth about a hundredth of a factor.
            G.timeLoss = 0;
        else
            %THE DROPPED-FRAME ALLOWANCE IS RE-DERIVED WITH THE CUT and not read off
            %results.nvc.timeLoss, and that is the more faithful of the two: the tree
            %stores it as single while the producer graded the confidence on the double
            %it computed, and the tree carries the PRIMARY signal's series only - gsData
            %runs on its own clock and lost its own frames.
            [G.cut,G.timeLoss] = cutEpochs(clk,sf.epochStartSec,sf.epochDurationSec);
            G.t     = clk(G.cut(:,1))-clk(G.cut(1,1));
            G.trust = logical(N.epochTrust(:))';

            %THE CUT IS CHECKED AGAINST THE TREE BEFORE ANYTHING IS DRAWN.  The
            %repetitions are re-cut here because the product stores none;
            %results.nvc.epochStart is where the producer cut them, so the two must land
            %on the same frames.  They cannot differ unless the settings file beside this
            %product belongs to another run, which is exactly what this says out loud
            %instead of drawing.
            if ~strcmp(sig,'gsData')
                es = double(N.epochStart(:))';
                dt = median(diff(clk));
                if size(G.cut,1)~=numel(N.time) || numel(es)~=size(G.cut,2) || ...
                        max(abs(clk(G.cut(1,:))'-es)) > 0.5*dt
                    error('guiResponse:settingsMismatch', ...
                        ['The settings beside this recording describe a different ' ...
                         'run: the repetitions they place do not land where ' ...
                         'results.nvc.epochStart says they were cut.  Re-run the ' ...
                         'stimulus response, or open the settings file that belongs ' ...
                         'to this analysis.']);
                end
            end
        end

        %ALL THREE LEVELS ARE MEASURED whatever the run stored: the levels exist to keep
        %an unwanted array out of a saved tree, nothing is saved here, and the confidence
        %needs the complete set anyway.
        sM = sf; sM.segNvcReturn = {};
        if G.collapsed, sM.rejectFirstEpoch = false; end
        G.s = sM;
        G.L = getNVCMetrics(G.t,sM);

        app.geom.(key) = G;
    end

% ---------------------------------------------------------------------
    function [y,E,who] = nvcSeries(G,sig)
        %nvcSeries  The trace on the screen, every repetition behind it, and which is
        %   which in words.
        M = app.results.(sig);
        if strcmp(c.scope.Value,'one')
            k  = min(round(c.segment.Value),size(M,2));
            y0 = double(M(:,k));
            whoSeg = sprintf('Segment %d of %d',k,size(M,2));
        else
            y0 = mean(double(M),2,'omitnan');
            whoSeg = sprintf('Average of %d segments',size(M,2));
        end
        E = y0(G.cut);

        if G.collapsed
            y = E(:,1);
            who = [whoSeg ', the representative repetition'];
            return
        end
        ep = c.epoch.Value;
        if ep==0
            keep = G.trust;
            %AN EMPTY KEPT SET IS NOT AN ERROR - runNVC says the same and writes NaN
            %columns - so the average falls back to every repetition and says so.
            if ~any(keep)
                y = mean(E,2);
                who = [whoSeg ', average of all ' num2str(size(E,2)) ...
                    ' repetitions - none was kept'];
                return
            end
            y = mean(E(:,keep),2);
            who = sprintf('%s, average of %d kept repetitions',whoSeg,sum(keep));
        else
            y = E(:,ep);
            who = sprintf('%s, repetition %d',whoSeg,ep);
            if ~G.trust(ep), who = [who ' - not kept']; end
        end
    end

% ---------------------------------------------------------------------
    function ep = currentEpoch(G)
        %currentEpoch  Which repetition is on screen, 0 for an average.
        if G.collapsed, ep = 1; else, ep = c.epoch.Value; end
    end

% ---------------------------------------------------------------------
    function [m,cf] = nvcNumbers(G,E,y)
        %nvcNumbers  The markers of what is drawn and the confidence that goes with it.
        %
        %   A SINGLE REPETITION IS MEASURED THE WAY THE PRODUCER MEASURED IT, and that
        %   means every repetition of this trace goes through the core and ONE
        %   getNVCConfidence call is made over all of them.  Two things need it: the
        %   optional across-repetition group compares a repetition with the others of the
        %   same segment, and the first-repetition rule is about column 1 of a block.
        %   Handed one column, both would answer a different question - and every factor
        %   is a property of one unit, so a one-segment block gives exactly the number
        %   the product stores for that segment.
        %
        %   AN AVERAGE IS NOT ONE OF THE REPETITIONS, so it is measured on its own, with
        %   the first-repetition rule off - leaving it on would zero the average's
        %   confidence outright.  That is runNVC>collapseToEpoch's rule, not a
        %   convenience.
        names = G.L.blockNames;
        ep = currentEpoch(G);
        if ep>=1 && ~G.collapsed
            nEp = size(E,2);
            M = nan(1,numel(names),nEp,'single');
            bags = cell(1,nEp);
            for k=1:nEp
                bags{k} = getNVCMetrics(E(:,k),G.L,G.s);
                for j=1:numel(names), M(1,j,k) = bags{k}.(names{j}); end
            end
            C = getNVCConfidence(M,names,G.timeLoss,G.L,G.s);
            m = bags{ep};
            cf = confAt(C,ep);
        else
            sAvg = G.s; sAvg.rejectFirstEpoch = false;
            m = getNVCMetrics(y,G.L,sAvg);
            M = nan(1,numel(names),1,'single');
            for j=1:numel(names), M(1,j,1) = m.(names{j}); end
            loss = 0;
            if ~G.collapsed && any(G.trust), loss = mean(G.timeLoss(G.trust)); end
            C = getNVCConfidence(M,names,loss,G.L,sAvg);
            cf = confAt(C,1);
        end
    end

% ---------------------------------------------------------------------
    function G = vrcGeom(want)
        %vrcGeom  The decimated base, the injection anchor and the windows.
        %   Keyed on the level set because switching the fit off changes the layout
        %   the core is built with - and only then.
        key = ['w' num2str(numel(want))];
        if isfield(app.geom,key), G = app.geom.(key); return, end

        s2 = app.settings.runFitVasoreactivity;
        s2.segVrcReturn = want;
        G.s = s2;
        G.L = fitVasoreactivity(double(app.results.time(:)),s2);

        %The same check the stimulus branch makes: the decimated base the settings
        %imply is the base the tree publishes, or the settings belong to another run.
        vr = app.results.vasoreactivity;
        if numel(G.L.time)~=numel(vr.time) || ...
                abs(G.L.tInj-double(vr.injection)) > G.L.dt
            error('guiResponse:settingsMismatch', ...
                ['The settings beside this recording describe a different run: the ' ...
                 'time base and injection they imply are not the ones stored in ' ...
                 'results.vasoreactivity.  Re-run the vasoreactivity fit, or open ' ...
                 'the settings file that belongs to this analysis.']);
        end
        app.geom.(key) = G;
    end

% ---------------------------------------------------------------------
    function [y,who] = vrcSeries(sig)
        %vrcSeries  The trace on the screen, at the RECORDING's rate - the core
        %   decimates it itself, so the drawn curve and the fitted one cannot be
        %   decimated by two slightly different rules.
        M = app.results.(sig);
        if strcmp(c.scope.Value,'one')
            k = min(round(c.segment.Value),size(M,2));
            y = double(M(:,k));
            who = sprintf('Segment %d of %d',k,size(M,2));
        else
            y = mean(double(M),2,'omitnan');
            who = sprintf('Average of %d segments',size(M,2));
        end
    end

% =====================================================================
%                            THE NUMBERS
% =====================================================================
    function txt = nvcLines(m,cf,G,who)
        %nvcLines  What was measured, then how much of it to believe.  Every line asks
        %   the bag whether it has the field and every timing line is named from the
        %   LAYOUT, so a protocol that asked for five levels prints five lines without a
        %   word changing here.
        txt = [{who}; {''}];
        if ~m.valid
            txt = [txt; {'This trace has gaps, so nothing can be measured on it.'}];
            return
        end
        txt = [txt; {'MEASURED'}];
        txt = [txt; num(m,'Bl',    'baseline       %.4g')];
        txt = [txt; num(m,'BlStd', 'baseline noise %.3g')];
        txt = [txt; num(m,'Fn',    'level at end   %.4g')];
        txt = [txt; num(m,'FnStd', 'noise at end   %.3g')];
        txt = [txt; num(m,'Ep',    'mean level     %.4g')];
        txt = [txt; {''}];
        txt = [txt; num(m,'St',    'response      %+.4g')];
        txt = [txt; pct(m,'StRel', 'relative      %+.1f %%')];
        txt = [txt; num(m,'Pk',    'peak          %+.4g')];
        txt = [txt; pct(m,'PkRel', 'relative      %+.1f %%')];
        txt = [txt; num(m,'Auc',   'area          %+.4g')];
        txt = [txt; num(m,'AucRel','extra flow    %+.3g s')];
        %THE TIMES ARE NAMED FROM THE SETTING and quoted from the stimulus, which is why
        %one of them can be negative - a rise that had already happened before the mark.
        txt = [txt; {''};{'DELIVERED BY, from the stimulus'}];
        for j=1:numel(G.L.pcts)
            txt = [txt; numAt(m,sprintf('T%d',G.L.pcts(j)), ...
                sprintf('%d %%',G.L.pcts(j)),'%.2f s')];                        %#ok<AGROW>
        end
        %THE CONFIDENCE IS THE PRODUCT'S OWN for a single repetition of a single segment,
        %so this block is not a second opinion about the row - it is the row.
        txt = [txt; {''};{'CONFIDENCE'}];
        txt = [txt; {sprintf('overall        %.3f',cf.conf)}];
        txt = [txt; {sprintf('weakest check  %.3f',cf.confMin)}];
        txt = [txt; {sprintf('kept           %s',yesNo(cf.trust))}];
        if ~isempty(cf.weakestWords)
            txt = [txt; {''};{'the weakest check says'}; ...
                cellstr(wrapWords(cf.weakestWords,26))];
        end
    end

% ---------------------------------------------------------------------
    function txt = vrcLines(m,G,who)
        %vrcLines  Minutes, not seconds, for every timing: a t_max of 780 s is a
        %   number to divide and 13 min is a number to read - the report page's rule.
        txt = {who;''};
        if ~m.valid
            txt = [txt; {'This trace has gaps, so nothing can be measured on it.'}];
            return
        end
        txt = [txt; {'MEASURED'}];
        txt = [txt; num(m,'BlMean', 'baseline      %.4g')];
        txt = [txt; num(m,'BlStd',  'baseline sd   %.3g')];
        txt = [txt; num(m,'Peak',   'peak          %+.4g')];
        txt = [txt; pct(m,'PeakRel','relative      %+.1f %%')];
        txt = [txt; inMinutes(m,'TMax', 'time to peak  %.1f min')];
        txt = [txt; inMinutes(m,'TWash','washout       %.1f min')];
        txt = [txt; num(m,'AUC',    'area          %.4g')];

        if isfield(m,'Gain')
            txt = [txt; {'';sprintf('FITTED  (%s)',G.L.model)}];
            txt = [txt; num(m,'Gain',     'peak change   %.4g')];
            txt = [txt; pct(m,'CVR',      'relative      %+.1f %%')];
            txt = [txt; pct2(m,'CVRSE',   '  +/- 95%%     %.2f %%')];
            txt = [txt; inMinutes(m,'Onset','onset         %.1f min')];
            txt = [txt; inMinutes(m,'TMaxModel','time to peak  %.1f min')];
            txt = [txt; num(m,'RateFast','rate, rising  %.4g 1/s')];
            txt = [txt; num(m,'RateSlow','rate, decay   %.4g 1/s')];
            txt = [txt; num(m,'Drift',   'drift         %.3g per s')];
            txt = [txt; num(m,'R2',      'R2            %.4f')];
            txt = [txt; num(m,'RMSE',    'RMSE          %.4g')];
            txt = [txt; num(m,'StartsAgree',sprintf('starts agree  %%d of %d',G.L.nStarts))];
        end
        %A FIT WITH THE STEAL FLAG RAISED IS A FIT NOT TO TRUST, so it is called out
        %rather than filed in the block - the report page says the same.
        if isfield(m,'Steal') && m.Steal>0
            txt = [txt; {'';'*** STEAL DETECTED ***'; ...
                sprintf('%.1f min below baseline',double(m.StealSec)/60); ...
                'this model cannot describe a sub-baseline'; ...
                'dip, so the peak and washout times above'; ...
                'are not to be trusted'}];
        end
    end

% =====================================================================
%                        THE PROGRAMMATIC API
% =====================================================================
    function setKind(k)
        if ~ismember(k,c.kind.ItemsData)
            error('guiResponse:noSuchKind','This recording carries no %s analysis.',k);
        end
        c.kind.Value = k; onKind();
    end
    function setSignal(sg)
        if ~ismember(sg,c.signal.ItemsData)
            error('guiResponse:noSuchSignal','This analysis did not measure %s.',sg);
        end
        c.signal.Value = sg; onSignal();
    end
    function setSegment(i)
        if isempty(i)
            c.scope.Value = 'all';
        else
            c.scope.Value = 'one'; c.segment.Value = i;
        end
        onScope();
    end
    function setEpoch(e)
        c.epoch.Value = e; redraw();
    end
    function setFit(tf)
        c.showFit.Value = logical(tf); redraw();
    end
    function setEach(tf)
        c.showEach.Value = logical(tf); redraw();
    end
    function st = shownState()
        st = app.shown;
    end
    function savePlot(p)
        %THE WHOLE LAYOUT IS WRITTEN, not the top panel: on a stimulus page the
        %cumulative curve is where every time on the screen comes from, and half the
        %evidence is not a saved plot.  Resolution is an IMAGE setting; passing it with
        %a vector format is an error rather than a no-op, so the two cases are written
        %out.
        if isempty(app.tl) || ~isgraphics(app.tl)
            error('guiResponse:nothingDrawn','There is no plot to save yet.');
        end
        [~,~,e] = fileparts(p);
        if any(strcmpi(e,{'.pdf','.eps','.emf'})), exportgraphics(app.tl,p);
        else, exportgraphics(app.tl,p,'Resolution',200);
        end
    end
    function onSave()
        [f,d] = uiputfile({'*.png';'*.pdf';'*.jpg'},'Save the plot', ...
            fullfile(app.lastRoot,'response.png'));
        figure(fig);
        if isequal(f,0), return, end
        savePlot(fullfile(d,f));
    end

% =====================================================================
%                            SMALL THINGS
% =====================================================================
    function busy(tf)
        if ~isgraphics(fig), return, end
        if tf, fig.Pointer = 'watch'; else, fig.Pointer = 'arrow'; end
        drawnow          % NOT limitrate: limitrate DROPS updates, and the one update
    end                  % that must not be dropped is the pointer coming back
    function say(fmt,varargin)
        %A REFUSAL IS THE MESSAGE, so it must never fall down the console.  Visible
        %is an on/off ENUM rather than a char, and strcmp against 'on' answers false
        %for a window that is plainly on the screen - which would send every one of
        %these sentences to a warning nobody is reading.
        msg = sprintf(fmt,varargin{:});
        if strcmp(char(string(fig.Visible)),'on'), uialert(fig,msg,'guiResponse');
        else, warning('guiResponse:message','%s',msg);
        end
    end
end

% =====================================================================
function cf = confAt(C,k)
%confAt  One (unit, repetition) pair's confidence, and THE CHECK IT SCORED LOWEST ON.
%
%   THE CHECKS ARE READ OUT OF THE CORE, never listed: there are nine of them by default
%   and fourteen with the across-repetition group on, and the timing family is named from
%   s.nvcAreaPcts.  Their sentences come from the same place, so this window has no
%   phrasebook of its own to fall out of step.
cf = struct('conf',double(C.conf(1,k)),'confMin',double(C.confMin(1,k)), ...
    'trust',logical(C.trust(1,k)),'weakest','','weakestWords','');
best = Inf;
for j=1:numel(C.factorNames)
    v = double(C.factors.(C.factorNames{j})(1,k));
    if isfinite(v) && v<best
        best = v;
        cf.weakest = C.factorNames{j};
        cf.weakestWords = C.factorWords{j};
    end
end
end

% =====================================================================
function [cut,timeLoss] = cutEpochs(clk,startSec,durSec)
%cutEpochs  Where every repetition sits in a signal's own samples, [nES x nEp], and how
%   many seconds of it the camera did not deliver.
%
%   THE SAME NEAREST-FRAME RULE runNVC cut them with, repeated because the product
%   stores the numbers and not the cuttings ([samples x segments x repetitions] is
%   18.6 GB on a reference recording).  Every repetition is the same number of samples -
%   the epoch clock and the windows are shared - and the count is the TYPICAL span rather
%   than the shortest: a recording that stops on its own final sample makes the last
%   repetition one frame short by arithmetic, and letting that shorten all twenty would
%   throw a sample away from every one of them to describe the edge of one.
%
%   AND THE SAME LOSS RULE, for the same reason: everything beyond one nominal step is
%   time the camera did not deliver.  It is re-derived here rather than read out of
%   results.nvc.timeLoss because the tree holds the PRIMARY signal's, as a single, and
%   the confidence was graded on the double of whichever signal is being drawn.
clk = double(clk(:));
nT  = numel(clk);
startSec = double(startSec(:))';
endSec   = startSec+double(durSec);

[~,startFrame] = min(abs(clk-startSec),[],1);
[~,endFrame]   = min(abs(clk-endSec),[],1);

nES = min(round(median(endFrame-startFrame)), nT-startFrame(end)+1);
if nES<4
    error('guiResponse:shortEpoch', ...
        'A repetition of %g s holds only %d samples of this signal.',durSec,nES);
end
cut = startFrame+(0:nES-1)';

gaps = [0; diff(clk)-median(diff(clk))];
gaps(gaps<0) = 0;
timeLoss = sum(gaps(cut),1);
end

% =====================================================================
function [kinds,labels,missing] = availableKinds(results,settings)
%availableKinds  Which analyses this recording carries AND can be described.
%   A tree without the settings it was made with cannot be drawn against its own
%   windows, so it is named in the label rather than silently dropped.
%
%   AND NEITHER CAN A TREE FROM A RETIRED VERSION OF THE STEP.  This is a refusal and
%   not a fallback: nothing here reads an old field name or offers an old marker, it
%   simply recognises that the branch in front of it is not the one runNVC writes and
%   says which step to run again.  Without it the shape mismatch arrives as an
%   undefined-field error out of a picker callback, where no sentence can reach the
%   person who caused it.
kinds = {}; labels = {}; missing = '';
have = {'nvc','runNVC','Stimulus response'; ...
        'vasoreactivity','runFitVasoreactivity','Vasoreactivity'};
for i = 1:size(have,1)
    if ~isfield(results,have{i,1}), continue, end
    if ~isfield(settings,have{i,2})
        missing = sprintf('%s  (the %s analysis has no settings beside it)', ...
            missing,lower(have{i,3}));
        continue
    end
    if strcmp(have{i,1},'nvc') && ~currentNvcTree(results.nvc)
        missing = sprintf(['%s  (the stimulus response in this recording was ' ...
            'written by an earlier version of the step - run it again)'],missing);
        continue
    end
    kinds{end+1}  = have{i,1};                                          %#ok<AGROW>
    labels{end+1} = have{i,3};                                          %#ok<AGROW>
end
end

% =====================================================================
function tf = currentNvcTree(N)
%currentNvcTree  Is this the branch runNVC writes now?  Two facts and no more: it holds
%   the per-signal marker sets, and it either says which repetitions the recording kept
%   or has no repetitions at all because it was collapsed into one.
tf = isstruct(N) && isfield(N,'esMetrics') && ...
    (isfield(N,'epochTrust') || ~isfield(N,'epochStart'));
end

% =====================================================================
function s = signalLabel(sig)
%signalLabel  The explorer's vocabulary, so one phrase describes a signal wherever
%   the library shows one.
switch sig
    case 'sData',       s = 'Segment signal';
    case 'dvsData',     s = 'Vessel signal';
    case 'dvsDiameter', s = 'Vessel diameter';
    case 'gsData',      s = 'Segment signal at full time resolution';
    otherwise,          s = sig;
end
end

% =====================================================================
function u = signalUnit(sig)
%signalUnit  What the y-axis is a measurement of - the wrappers' own two labels.
if strcmp(sig,'dvsDiameter'), u = 'Diameter, px'; else, u = 'Flow, BFI'; end
end

% =====================================================================
function c = num(m,name,fmt)
%num  One line, or nothing at all when this level does not emit it.
c = {};
if isfield(m,name) && isfinite(m.(name)), c = {sprintf(fmt,double(m.(name)))}; end
end

function c = numAt(m,name,label,fmt)
%numAt  One line whose LABEL is computed rather than typed - a timing level named from
%   s.nvcAreaPcts.
%
%   THE LABEL IS PLACED BY CONCATENATION AND NEVER BY A SECOND sprintf.  A label built
%   with sprintf carries a real per-cent sign, and handing that string back to sprintf as
%   a FORMAT makes the sign the start of a conversion: the number after it then has
%   nowhere to land and the line prints as a bare "10" with the time silently gone.  It
%   printed exactly that until somebody opened the window and read it.
c = {};
if isfield(m,name) && isfinite(m.(name))
    c = {[sprintf('%-15s',label) sprintf(fmt,double(m.(name)))]};
end
end

function c = pct(m,name,fmt)
%pct  The same, for a fraction printed as a percentage.
c = {};
if isfield(m,name) && isfinite(m.(name)), c = {sprintf(fmt,100*double(m.(name)))}; end
end

function c = pct2(m,name,fmt)
%pct2  A 95 % half-width from a standard error, as a percentage.
c = {};
if isfield(m,name) && isfinite(m.(name))
    c = {sprintf(fmt,100*1.96*double(m.(name)))};
end
end

function c = inMinutes(m,name,fmt)
%inMinutes  A drug timing, in the unit a reader reads it in.  NOT named `minutes`:
%   that is a MATLAB built-in, and a local function of that name would shadow it for
%   every line of this file.
c = {};
if isfield(m,name) && isfinite(m.(name)), c = {sprintf(fmt,double(m.(name))/60)}; end
end

function s = yesNo(tf)
if tf, s = 'yes'; else, s = 'no'; end
end

% =====================================================================
function lines = wrapWords(txt,width)
%wrapWords  One sentence over as many lines as a fixed-width panel needs.  The numbers
%   panel is a Courier text area and does not reflow, so a sentence handed to it whole
%   would run off the right-hand edge of the one line that explains the score above it.
lines = {};
words = strsplit(strtrim(char(txt)));
cur = '';
for i = 1:numel(words)
    if isempty(cur), trial = words{i}; else, trial = [cur ' ' words{i}]; end
    if numel(trial)>width && ~isempty(cur)
        lines{end+1} = cur;                                             %#ok<AGROW>
        cur = words{i};
    else
        cur = trial;
    end
end
if ~isempty(cur), lines{end+1} = cur; end
lines = lines(:);
end

% =====================================================================
function yl = padY(ax,varargin)
%padY  Headroom above the trace, sized from the CURVES rather than from the axes.
%
%   The peak carries its size as a text label and the timings carry theirs, and a
%   label fitted flush to the box by `axis tight` reads as part of the title rather
%   than as part of the plot.
%
%   THE RANGE IS TAKEN FROM THE DATA AND NOT READ BACK OUT OF THE AXES.  `axis tight`
%   does not resolve its limits before the next line runs, so a `ylim(ax)` query
%   immediately after it can still answer with the default [0 1] on one end - which pads
%   the trace into the top twentieth of the box.  The same rule the report pages use:
%   compute the band, then give it room.  IT HANDS THE BAND BACK for the same reason it
%   computes it: a caller placing a label at the top of the axes must use the limits that
%   were SET, not the ones a read-back would answer with.
lo = Inf; hi = -Inf;
for k = 1:numel(varargin)
    v = double(varargin{k}(:));
    v = v(isfinite(v));
    if isempty(v), continue, end
    lo = min(lo,min(v)); hi = max(hi,max(v));
end
if isfinite(lo) && isfinite(hi) && hi>lo
    yl = [lo hi]+(hi-lo)*[-0.04 0.12];
    ylim(ax,yl);
else
    yl = ylim(ax);          % a flat trace: whatever the axes settled on is the answer
end
end

% =====================================================================
function v = onOff(tf)
if tf, v = 'on'; else, v = 'off'; end
end

% =====================================================================
function [fileArg,vis] = parseArgs(a)
%parseArgs  guiResponse(path) / guiResponse(...,'Visible','off').
fileArg = ''; vis = 'on';
if ~isempty(a) && (ischar(a{1})||isstring(a{1})) && ~strcmpi(a{1},'Visible')
    fileArg = char(a{1});
    a(1) = [];
end
for i = 1:2:numel(a)-1
    if strcmpi(a{i},'Visible'), vis = char(string(a{i+1})); end
end
end
