%setLibraryPath - Put the library on the MATLAB path without the tooling folders.
%
%   addpath(genpath(libraryFolder)) is the right instinct and almost the right
%   command.  The problem is that genpath walks the FILESYSTEM: it does not read
%   .gitignore, it does not read .git/info/exclude, and it does not skip folders
%   whose name begins with a dot.  So everything under .claude/ goes on the path
%   too - including .claude/worktrees/, where the coding assistant checks out
%   COMPLETE COPIES OF THIS LIBRARY at older commits.
%
%   THAT IS A SHADOW, AND IT IS THE DANGEROUS KIND.  A leftover worktree puts a
%   second getVasomotion.m, runSegmentation.m, reportSave.m and so on where MATLAB
%   can find them, and the copy can WIN: which -all then shows the worktree version
%   first and the real file marked "% Shadowed".  Nothing errors.  The run simply
%   uses code from an older commit - against real data, silently.  This has already
%   happened once (2026-07-26), where it made a green test suite go red while the
%   actual library was fine.
%
%   The fix is one rule, in one place: add everything, then drop the .claude
%   subtree.  Nothing else is removed.  Drafts/ stays on the path (it is
%   work-in-progress the launchers still reach into), 3rd party/ stays on the path
%   (Bio-Formats and friends), and the git-excluded claude-docs / claude-tests /
%   claude-literature folders stay too - they hold notes and test harnesses with
%   names of their own, and none of them can shadow a library function.
%
%   Safe to call when there is no .claude folder at all, which is the normal case
%   for anyone who checked this library out of git: it then does exactly what the
%   plain addpath(genpath(...)) did.
%
% Syntax:
%    setLibraryPath(libraryFolder)
%    dropped = setLibraryPath(libraryFolder)
%
% Inputs:
%    libraryFolder - the library's root folder (the one holding Core, Wrappers,
%                    Launchers, ...).
%
% Outputs:
%    dropped - cellstr of the folders that were kept OFF the path, so a caller or a
%              test can say which shadows were avoided.  Empty on a clean checkout.
%
% Example:
%    libraryFolder = 'C:\Data\Matlab-Dynamic-LIght-Scattering-Library';
%    setLibraryPath(libraryFolder);
%
% See also: addpath, genpath, rmpath, which
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function dropped = setLibraryPath(libraryFolder)

libraryFolder = char(libraryFolder);
if ~isfolder(libraryFolder)
    error('setLibraryPath:noFolder','Not a folder: %s', libraryFolder);
end

addpath(genpath(libraryFolder));

% The one exclusion.  rmpath complains about folders that are not on the path, so
% the list is intersected with what genpath actually returned rather than handed
% over blind.
dropped  = {};
toolRoot = fullfile(libraryFolder,'.claude');
if isfolder(toolRoot)
    dropped = strsplit(genpath(toolRoot), pathsep);
    dropped = dropped(~cellfun(@isempty,dropped));
    if ~isempty(dropped)
        rmpath(strjoin(dropped,pathsep));
    end
end
end
