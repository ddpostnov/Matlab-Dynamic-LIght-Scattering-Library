%meEpochs - When the animal was still, and which frames the magnifier may believe.
%
%   THE ANIMAL IS AWAKE.  It runs, whisks and grooms, and a bout displaces the whole
%   preparation by many pixels.  That is not a cosmetic problem: it leaves the range
%   the magnifier is valid over - alpha*delta < lambda/4 - and tears the picture up,
%   and on the way it drags the whole-field phase mean that meGlobal removes, so a
%   bout can corrupt frames that are themselves quiet.  Every frame therefore carries
%   a quality flag before anything else reads it.
%
%   THIS CLASSIFIES; IT DOES NOT MEASURE.  meProbe already walks the whole recording
%   with a phase correlation and produces the frame-to-frame displacement series, and
%   that walk is the expensive part - one pass over 31 000 frames.  Re-deriving it
%   here would be a second answer to a question that has one.  So the series comes in
%   and the flag goes out, and meProbe classifies its own measurement through this
%   same function rather than carrying a second copy of the rule.
%
%   A BOUT IS NOT A REJECTED CYCLE, and the two must not be added up.  A bout spans
%   many beats and is a property of the ANIMAL; a rejected cycle is a property of one
%   beat and comes out of meGate.  Reporting them separately is what lets a reader
%   tell "the animal moved for a quarter of the recording" from "the pulse was
%   irregular", which are different facts about a preparation.
%
%   THE THRESHOLD IS RELATIVE TO THE RECORDING'S OWN MEDIAN STEP, not absolute.  A
%   quiet awake animal still moves - 0.078 px between frames on the reference
%   recording - and an absolute threshold would either call all of that a bout or
%   none of a calmer animal's.  s.quietCoef multiplies the median.
%
%   THE EDGES OF A BOUT GO WITH IT.  The displacement crosses the threshold in the
%   middle of a movement, not at its start, so the flag is dilated by s.boutPad
%   seconds either way.  Without it the frames on the shoulders of every bout - the
%   ones already moving fast but not yet fast enough - are kept, and they are exactly
%   the frames whose motion is large and unflagged.
%
% Syntax:
%    ep = meEpochs(s, displacement, fs)
%
% Inputs:
%    s            - settings.  Fields read:
%                   .quietCoef  a frame belongs to a bout when it moves this many
%                               times the recording's own median step.
%                   .boutPad    seconds the flag is grown either side of a bout.
%                   .contFrames longest window the continuous mode may take.  Empty
%                               takes the whole quiet stretch.
%                   .contSpanBouts  false keeps the continuous window inside one
%                               unbroken quiet stretch; true lets it reach
%                               .contFrames across a bout, with the bad frames
%                               marked in .badInWindow for the video to burn in.
%    displacement - [frames 1] frame-to-frame displacement in pixels of the original
%                   frame, from meProbe.  Its first entry has no predecessor and is
%                   not used.
%    fs           - frames per second.
%
% Outputs:
%    ep - .bad          [frames 1] logical, true inside a bout
%         .good         its negation, kept because every caller wants one of the two
%         .median       the recording's own median step, pixels
%         .threshold    the step above which a frame is a bout, pixels
%         .quietFraction
%         .nBouts       how many separate bouts there are
%         .runs         [n 2] first and last frame of every quiet stretch
%         .longestFrames  length of the longest one
%         .longestFirst   its first frame
%         .longestStart   the same in seconds, which is what meProbe reports
%         .longestSeconds
%         .window       [first last] the frames the continuous mode runs on
%         .windowFrames how many that is
%         .badInWindow  [windowFrames 1] logical, all false unless the window was
%                       allowed to span a bout
%         .spansBout    whether it does
%
% Example:
%    probe = load('rec_ME_probe.mat').probe;
%    ep    = meEpochs(meSettings, probe.motion.displacement, probe.pass.fs);
%    fprintf('quiet %.0f%%, continuous run on frames %d to %d\n', ...
%        100*ep.quietFraction, ep.window(1), ep.window(2));
%
% Dependencies: none beyond core MATLAB.
% See also: meProbe, meGate, meCine, meEnhance, meSettings
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 07-August-2026

%------------- BEGIN CODE --------------
function ep = meEpochs(s, displacement, fs)

d  = double(displacement(:));
nT = numel(d);
if nT < 3
    error('meEpochs:frames','A displacement series of %d frames says nothing about bouts.',nT);
end

% The first frame has no predecessor, so its zero would drag the median down.
med = median(d(2:end),'omitnan');
thr = s.quietCoef*med;

bad = d > thr;
pad = max(1, round(s.boutPad*fs));
bad = movmax(double(bad), [pad pad]) > 0;

runs = quietRuns(~bad);
if isempty(runs)
    error('meEpochs:noQuiet', ...
        ['Every frame of this recording is inside a movement bout at %.2f times the ' ...
         'median step. Nothing can be averaged; look at the recording.'], s.quietCoef);
end

[longest, iLong] = max(runs(:,2)-runs(:,1)+1);
first = runs(iLong,1);

% ---- the window the continuous mode runs on --------------------------------
% THE DEFAULT IS TO STAY INSIDE ONE UNBROKEN QUIET STRETCH even when that is shorter
% than asked for.  A bout does not merely look bad in a magnified movie; it exceeds
% the bound the magnifier is valid over, and marking a torn stretch does not turn it
% back into data.  s.contSpanBouts is there because the opposite choice is defensible
% when the point of the run is to show what a bout does - and then every bad frame in
% the window is named, because a gap the viewer cannot see is worse than one they can.
want = s.contFrames;
if isempty(want), want = longest; end
want = min(want, nT);

if s.contSpanBouts
    % Centre the window on the longest quiet stretch and grow it to what was asked
    % for, clipping at the ends of the recording rather than sliding past them.
    centre = first + floor((longest-1)/2);
    w0     = min(max(centre - floor((want-1)/2), 1), nT-want+1);
    win    = [w0, w0+want-1];
else
    win = [first, first + min(want, longest) - 1];
end
inWin = win(1):win(2);

ep = struct();
ep.bad           = bad;
ep.good          = ~bad;
ep.median        = med;
ep.threshold     = thr;
ep.quietFraction = mean(~bad);
ep.nBouts        = nnz(diff([false; bad(:)])==1);
ep.runs          = runs;
ep.longestFrames = longest;
ep.longestFirst  = first;
ep.longestStart  = (first-1)/fs;
ep.longestSeconds= longest/fs;
ep.window        = win;
ep.windowFrames  = numel(inWin);
ep.badInWindow   = bad(inWin);
ep.spansBout     = any(ep.badInWindow);
end

% =====================================================================
function runs = quietRuns(tf)
%quietRuns  First and last index of every unbroken stretch of true.
tf = tf(:)';
d  = diff([false tf false]);
b  = find(d==1)';
e  = find(d==-1)'-1;
runs = [b e];
end
%------------- END OF CODE --------------
