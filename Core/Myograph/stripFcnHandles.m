%stripFcnHandles  Remove function-handle fields from a parameter struct before saving
%
%   sClean = stripFcnHandles(s) returns s with every field that holds a
%   function_handle removed.  Anonymous functions can capture graphics objects in
%   their closure (e.g. the workbench sets s.progressFcn/s.cancelFcn/s.stageFcn to
%   callbacks that reference its uifigure); saving such a struct serialises the whole
%   captured figure into the .mat, and loading it later spawns a stray "ghost"
%   window.  Stripping the handles keeps the saved s to plain data so it round-trips
%   cleanly.  A non-struct input is returned unchanged.
%
%   INPUT   s       parameter struct (or anything; non-structs pass through)
%   OUTPUT  sClean  the same struct with all function_handle fields removed
%
%   DEPENDS ON  base MATLAB only.  Self-contained.
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 21-July-2026

function s = stripFcnHandles(s)
if ~isstruct(s), return; end
fn=fieldnames(s);
for i=1:numel(fn)
    if isa(s.(fn{i}),'function_handle')
        s=rmfield(s,fn{i});
    end
end
end
