%getProductPath - The SOURCE / RESULTS / SETTINGS member beside a data product.
%
%   Every recording this library processes yields a triplet of .mat files that
%   share one base name and differ only in a final role letter:
%
%       <base>_d.mat   SOURCE     the data cube (or the frame time base) + time
%       <base>_r.mat   RESULTS    images, masks, per-segment tables, metric trees
%       <base>_s.mat   SETTINGS   a copy of the s a step ran with
%
%   THE RESULTS FILE IS THE ONE THE LIBRARY REFERS TO.  Every step reads and
%   writes RESULTS; only some of them open the SOURCE at all, and a step that
%   was told not to keep it (runIntensity's s.saveSource) writes no SOURCE
%   whatsoever - so a file list built on the SOURCE names is a list that can be
%   missing the very recordings it is meant to describe.  fNames is therefore a
%   list of *_r.mat paths everywhere in this library, and the other two members
%   are named from it, here.
%
%   THE ROLE LETTER IS REPLACED, NOT THE FIRST THING THAT LOOKS LIKE IT.  The
%   substitution is anchored to the end of the name, so a recording whose own
%   stem ends in a role letter ('...Rat_d_t_K_r.mat') keeps its stem intact.
%   Any member of the triplet may be passed in - all three name the same
%   product - which is what lets one call site accept whichever member a caller
%   happens to be holding.
%
% Syntax:
%    p = getProductPath(fName)            % the RESULTS member (default)
%    p = getProductPath(fName, role)
%
% Inputs:
%    fName - char/string path of any member of one triplet ('..._d.mat',
%            '..._r.mat' or '..._s.mat').  A cell array is mapped elementwise
%            and keeps its size and shape; empty entries stay empty.
%    role  - 'd' (SOURCE) | 'r' (RESULTS, default) | 's' (SETTINGS).
%
% Outputs:
%    p     - the requested member's full path, or a cell array of them.
%
% Examples:
%    load(getProductPath(s.fName,'d'),'source');
%    save(getProductPath(s.fName,'s'),'settings','-v7.3','-nocompression');
%    dNames = getProductPath(fNames,'d');
%
% See also: getFileNamesList, removeProcessedFiles, myographProduct, wbFileModel
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function p = getProductPath(fName, role)

if nargin < 2 || isempty(role), role = 'r'; end
role = lower(char(role));
if ~any(strcmp(role,{'d','r','s'}))
    error('getProductPath:badRole', ...
        'Role must be ''d'', ''r'' or ''s'', not ''%s''.', role);
end

if iscell(fName)
    p = cell(size(fName));
    for idx = 1:1:numel(fName)
        if isempty(fName{idx}), p{idx} = fName{idx}; continue; end
        p{idx} = oneName(fName{idx}, role);
    end
    return
end
p = oneName(fName, role);
end

% =====================================================================
function p = oneName(fName, role)
%oneName  Swap the trailing role letter of a single name.
%   A NAME THAT CARRIES NO ROLE IS NOT SILENTLY REWRITTEN.  Handed a raw
%   recording or anything else that does not end in '_d.mat' / '_r.mat' /
%   '_s.mat', this says so, because the alternative - handing back the input
%   unchanged - turns a wrong file list into a load of the wrong file rather
%   than into an error anyone can read.
p = char(fName);
if isempty(regexp(p, '_[drs]\.mat$', 'once'))
    error('getProductPath:noRole', ...
        ['"%s" is not a member of a data triplet - a name ending in ' ...
         '"_d.mat", "_r.mat" or "_s.mat" was expected.'], p);
end
p = regexprep(p, '_[drs]\.mat$', ['_' role '.mat']);
end
