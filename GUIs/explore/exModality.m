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
%   A LAZY HANDLE IS ALWAYS SPECKLE, and that is a fact about file sizes rather
%   than a shortcut.  exScan falls back to h5info only above app.sizeLimitGB; the
%   files that large are per-pixel speckle results, while a myograph result is
%   megabytes (the largest in the reference sets is 24 MB).  h5info cannot see
%   inside a struct array anyway - results.intervals appears as an opaque
%   H5T_REFERENCE dataset - so a lazy myograph could not be read even if one
%   existed, and answering 'speckle' for it would be a silent wrong answer.  The
%   handle is therefore tested for the two names directly and, on the strength of
%   the size argument, reported as speckle when neither is there.
%
%   Syntax:
%      m = exModality(results)
%
%   INPUTS
%     results  the RESULTS member of a triplet, loaded; or the lazy handle
%              struct('x_lazy_',true,'x_path_',<file>) exScan uses for a file over
%              the size limit.
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

if isfield(results,'x_lazy_')
    if hasNonEmptyDataset(results,'intervals') || hasNonEmptyDataset(results,'channel')
        m = 'myograph';
    end
    return
end

if ~isempty(myographIntervals(results)), m = 'myograph'; end
end

% =====================================================================
function tf = hasNonEmptyDataset(handle, name)
%hasNonEmptyDataset  Does the >size-limit file carry this branch with anything in
%   it?  A struct array is stored as an object-reference dataset whose dataspace is
%   the ARRAY's size, so an empty placeholder reads back as 0 elements and a real
%   one does not - which is the same non-empty test the loaded arm makes, asked of
%   the header instead of the data.
tf = false;
if ~isfield(handle,'x_path_'), return; end
try
    info = h5info(char(handle.x_path_), '/results');
catch
    return
end
if isempty(info.Datasets), return; end
k = strcmp({info.Datasets.Name}, name);
if ~any(k), return; end
sz = info.Datasets(k).Dataspace.Size;
tf = ~isempty(sz) && all(sz > 0);
end
