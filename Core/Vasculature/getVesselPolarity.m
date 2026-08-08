%getVesselPolarity - are this product's vessels bright or dark against the tissue?
%
%   [polarity,product] = getVesselPolarity(fName)
%
% DESCRIPTION
%   THE ONE PLACE THE QUESTION IS ASKED.  Every step that looks for vessels in a
%   picture has to know which way round they are, and before this function existed
%   the answer was a boolean called isK spread over five files - which conflated the
%   polarity with three other things that merely happened to be true of the same
%   files (the SLSCI edge width, the contrast trust mask, and the diffusion
%   schedule).  Those three follow the PRODUCT TOKEN and this one does not, so they
%   are asked separately now; getPixelCategories' header sets out the whole split.
%
%   THE RULE HAS TWO ARMS AND NEITHER IS A FALLBACK FOR THE OTHER.
%
%     '_K' (speckle contrast)   ->  ALWAYS 'dark'.  Fast flow decorrelates the
%           speckle and lowers the contrast, so a vessel is a dark object in a
%           contrast image as a matter of physics.
%
%     '_I' (fluorescence intensity)  ->  ALWAYS 'bright'.  A plasma label makes the
%           vessels the bright objects, and THIS VERSION OF THE LIBRARY PROCESSES
%           ONLY SUCH RECORDINGS - the author's decision of 2026-08-08, taken in
%           exchange for a fluorescence branch with no half-supported second path in
%           it.  Preparations that invert the contrast exist; running one through
%           this library produces wrong numbers that look entirely ordinary, and the
%           operator is trusted not to.  Dark polarity is a separate development, and
%           it starts from claude-docs/intensity-branch/03-regions-segmentation.md
%           and 05-bright-only.md rather than from anything left in the tree.
%
%   IT IS STILL A FUNCTION RATHER THAN A CONSTANT, and that is the point of it.  Both
%   arms are reached on every run of the two branches, and welding the answer back to
%   a file-name test inside each core is exactly the arrangement that produced
%   runDynamicSegmentation's unconditional imcomplement - a real bug on the intensity
%   path, found only because the four decisions were separated.
%
% INPUT
%   fName     the RESULTS member's path, i.e. the '*_K_r.mat' or '*_I_r.mat' the
%             step was handed.  Reading the name is this function's JOB - it is why
%             the other cores can be told never to do it.
%
% OUTPUT
%   polarity  'bright' | 'dark' - are the vessels brighter or darker than the
%             tissue.  What decides whether an image is inverted before ridges are
%             looked for in it.
%   product   'K' | 'I' - the product token.  Returned beside the polarity because
%             the two used to be one flag and every caller needs both: the token
%             decides which mean image to read, whether there is an SLSCI border to
%             measure and whether a contrast trust mask means anything.
%
% See also: getPixelCategories, setRegions, runSegmentation,
%           runDynamicSegmentation, showSegmentsPreview
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 08-August-2026

%------------- BEGIN CODE --------------
function [polarity,product] = getVesselPolarity(fName)

fName=char(fName);

if contains(fName,'_K_r.mat')
    product='K';
    polarity='dark';            % physics, not a preference - see the header
    return
end
if ~contains(fName,'_I_r.mat')
    error('getVesselPolarity:unknownProduct', ...
        ['"%s" is neither a contrast nor an intensity result, so there is no way to ' ...
         'tell which way round its vessels are.'], fName);
end

product='I';
polarity='bright';              % the branch's stated assumption - see the header
end
%------------- END OF CODE --------------
