%hasLabChartSDK  Is the ADInstruments SDK reachable, so a LabChart file could be read?
%
%   tf = hasLabChartSDK() answers whether the vendored ADInstruments SDK is on the
%   path.  It does not read anything and it does not throw - it is the question, not
%   the answer.
%
%   IT EXISTS SO THE QUESTION IS ASKED IN ONE PLACE.  readLabChart is the only file
%   in this library that CALLS the SDK, and that is a rule with a test behind it -
%   the SDK is third-party code that has to stay swappable.  But other steps have a
%   legitimate reason to ask whether a recording could be read AGAIN without trying:
%   the intervals step discards everything outside the analysed windows, and what it
%   tells the operator to do about a missing recording depends on whether reading it
%   back is even possible on this machine.  That is a different question from "read
%   this file", and answering it by calling the reader and catching the error would
%   mean opening a file to find out whether a folder is on the path.
%
%   IT ASKS 'which', NOT 'exist'.  exist('adi.readFile','file') answers 0 for a
%   PACKAGE function that is perfectly reachable, so the obvious test would report a
%   correctly installed SDK as missing.  setLibraryPath puts the package's PARENT on
%   the path and never the '+adi' folder itself, which is what makes 'which' the only
%   test that tracks how the library actually loads it.
%
%   Syntax:
%      tf = hasLabChartSDK()
%
%   OUTPUT
%     tf   true when a LabChart file could be read on this machine.  The SDK runs on
%          Windows with 64-bit MATLAB only, and belongs in the library's "3rd party"
%          folder as "adinstruments_sdk_matlab".
%
%   EXAMPLE
%     if ~hasLabChartSDK()
%         warning('This recording cannot be read again on this machine.');
%     end
%
% See also: readLabChart, runLabChart, setMyographIntervals, setLibraryPath
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 03-August-2026

%------------- BEGIN CODE --------------
function tf = hasLabChartSDK()

tf = ~isempty(which('adi.readFile'));
end
