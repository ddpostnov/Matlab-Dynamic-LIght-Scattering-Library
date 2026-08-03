%exModality - Is this a myograph recording or a speckle one?
%
%   m = exModality(results) answers with one word, and it answers from the DATA
%   rather than from the file name, the workbench session or anything a caller
%   remembered.  A myograph recording holds its analysis inside the WINDOWS it was
%   measured in; nothing else in the library does.
%
%   THE TEST IS "NON-EMPTY", AND THAT IS NOT PEDANTRY.  Both myograph shapes carry
%   BOTH field names, and one of them is always an empty placeholder:
%
%       a pressure myograph   results.intervals [1x3]   results.channel   [0x0]
%       a wire myograph       results.intervals [0x0]   results.channel   [1x4]
%
%   So `isfield(results,'intervals')` is true for a wire recording that has no
%   flat intervals at all, and would become true for a speckle file the moment
%   either name appeared anywhere in the tree.  Asking whether a window actually
%   EXISTS is the question that stays right.
%
%   The flattening is myographIntervals' job, not this function's: it already knows
%   both shapes, walks the channels in tree order, and writes .channelName onto
%   every window.  Calling it here means there is exactly one place in the library
%   that knows a wire myograph hangs its windows off the channel axis.
%
%   Syntax:
%      m = exModality(results)
%
%   INPUTS
%     results  the RESULTS member of a triplet, loaded.  There is no second answer
%              for a file too big to load: since D7 the explorer loads every
%              results file it reads, so this question is always asked of a tree.
%
%   OUTPUTS
%     m  'myograph' when the recording holds at least one analysed window,
%        'speckle' otherwise - including for anything that is not a results struct
%        at all, so a caller handed rubbish gets the harmless answer.
%
%   EXAMPLE
%      S = load('WT794 uden endothel_MYO_r.mat');
%      exModality(S.results)          % 'myograph'
%
%   DEPENDS ON
%     myographIntervals (Core/Myograph).
%
% See also: exAxes, exSchema, exScan, myographIntervals, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function m = exModality(results)

m = 'speckle';
if isempty(results) || ~isstruct(results) || ~isscalar(results), return; end

if ~isempty(myographIntervals(results)), m = 'myograph'; end
end
