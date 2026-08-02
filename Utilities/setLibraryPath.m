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
%   THE SECOND SHADOW, AND WHY DROPPING IS NOT THE ANSWER.  The claude-docs /
%   claude-tests / claude-literature trees STAY on the path - the test harnesses are
%   called by name and would stop resolving without them.  But claude-tests has held
%   files named EXACTLY after real wrappers (runSegmentation.m, runVasomotion.m,
%   runBFI.m ... , frozen "before" copies kept as refactor evidence).  Those are the
%   same hazard as a worktree, from a folder that must remain reachable.
%
%   Measured 2026-08-01: they resolved to the real wrapper, but only because MATLAB
%   searches the path in order and 'Wrappers' happens to sort before 'claude-tests' in
%   ASCII.  That is luck, not a design, and it flips the moment a folder is renamed.
%
%   So the tooling trees are DEMOTED instead: added, then explicitly moved to the END
%   of the path, below every real library folder.  A name collision then always
%   resolves to the library, deterministically.  Any collision that does exist is
%   RETURNED rather than left silent, because a test harness shadowing a wrapper is
%   worth knowing about even when the library wins.
%
%   Safe to call when there is no .claude folder at all, which is the normal case
%   for anyone who checked this library out of git: it then does exactly what the
%   plain addpath(genpath(...)) did.
%
% Syntax:
%    setLibraryPath(libraryFolder)
%    dropped = setLibraryPath(libraryFolder)
%    [dropped,shadowed] = setLibraryPath(libraryFolder)
%
% Inputs:
%    libraryFolder - the library's root folder (the one holding Core, Wrappers,
%                    Launchers, ...).
%
% Outputs:
%    dropped  - cellstr of the folders that were kept OFF the path, so a caller or a
%               test can say which shadows were avoided.  Empty on a clean checkout.
%    shadowed - cellstr of function names that exist BOTH in a library folder and in a
%               claude-* tooling tree.  The library copy wins (that is what the
%               demotion guarantees); this names the collisions so they can be fixed
%               at the source, by giving the frozen copy a distinct name.
%
% Example:
%    libraryFolder = 'C:\Data\Matlab-Dynamic-LIght-Scattering-Library';
%    [dropped,shadowed] = setLibraryPath(libraryFolder);
%
% See also: addpath, genpath, rmpath, which
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 29-July-2026

%------------- BEGIN CODE --------------
function [dropped,shadowed] = setLibraryPath(libraryFolder)

libraryFolder = char(libraryFolder);
if ~isfolder(libraryFolder)
    error('setLibraryPath:noFolder','Not a folder: %s', libraryFolder);
end

addpath(genpath(libraryFolder));

% ---- 1. the one exclusion: .claude, which can hold whole checkouts of this library ----
% rmpath complains about folders that are not on the path, so the list is intersected
% with what genpath actually returned rather than handed over blind.
dropped  = {};
toolRoot = fullfile(libraryFolder,'.claude');
if isfolder(toolRoot)
    dropped = strsplit(genpath(toolRoot), pathsep);
    dropped = dropped(~cellfun(@isempty,dropped));
    if ~isempty(dropped)
        rmpath(strjoin(dropped,pathsep));
    end
end

% ---- 2. demote the claude-* tooling trees BELOW every real library folder ----
% They must stay reachable (the test harnesses are called by name), so they cannot be
% dropped; but they must never outrank Core/ or Wrappers/.  Removing and re-adding
% with '-end' puts them at the bottom of the path in their existing order, which makes
% the library's precedence a property of this function rather than of ASCII sorting.
entries = strsplit(path, pathsep);
entries = entries(~cellfun(@isempty,entries));
isTool  = startsWith(entries, [libraryFolder filesep 'claude-']);
tooling = entries(isTool);
if ~isempty(tooling)
    rmpath(strjoin(tooling,pathsep));
    addpath(strjoin(tooling,pathsep),'-end');
end

% ---- 3. name the collisions, so a shadow is visible even though the library wins ----
shadowed = {};
if nargout>1 && ~isempty(tooling)
    libDirs  = entries(~isTool & startsWith(entries,libraryFolder));
    shadowed = collidingNames(libDirs,tooling);
end
end

% ---------------------------------------------------------------
function names = collidingNames(libDirs,toolDirs)
%COLLIDINGNAMES  function names present in BOTH a library folder and a tooling folder
names = intersect(mNames(libDirs), mNames(toolDirs));
names = names(:)';
end

% ---------------------------------------------------------------
function n = mNames(dirs)
%MNAMES  the .m base names directly inside each of DIRS (genpath already recursed)
%   fileparts, not erase: erase('.m') would strip EVERY occurrence, so a name that
%   happens to contain '.m' elsewhere would come out mangled.
n = {};
for i = 1:numel(dirs)
    f = dir(fullfile(dirs{i},'*.m'));
    b = cell(1,numel(f));
    for k = 1:numel(f)
        [~,b{k}] = fileparts(f(k).name);
    end
    n = [n, b]; %#ok<AGROW>
end
n = unique(n);
end
