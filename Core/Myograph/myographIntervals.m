%myographIntervals  The analysed windows of a myograph recording, however they are stored
%
%   [iv,channels] = myographIntervals(results) is the ONE place that knows a myograph
%   result tree holds its windows in two shapes, and hands every consumer the flat
%   list it actually wants.
%
%   THE TWO SHAPES, AND WHY THERE ARE TWO
%     A pressure myograph measures ONE vessel, so a window is a stretch of that one
%     recording and the tree is flat:
%         results.intervals(k)
%     A wire myograph records SEVERAL CHAMBERS on one LabChart file, and a window may
%     belong to one of them - the drug went into chamber 2 at a different minute from
%     chamber 3 - so the CHANNEL is the outer axis and the windows hang off it:
%         results.channel(i).name
%         results.channel(i).intervals(k)
%     A window the operator left on 'all channels' appears under every channel, which
%     is what it means; its own .channels stays empty, so the editor can tell the two
%     apart again when the recording is re-opened.
%
%   WHAT COMES BACK is the windows CONCATENATED, channel by channel and in time
%   order within each, with one field added:
%       iv(k).channelName   the channel this window was analysed on, '' for a
%                           pressure myograph
%   Everything else about an element is untouched, so a consumer that already reads
%   .name / .tStart / .diameter / .vasomotion keeps working on the flat list and only
%   has to decide whether it cares about the channel.
%
%   Syntax:
%      iv              = myographIntervals(results)
%      [iv,channels]   = myographIntervals(results)
%
%   INPUTS
%     results  the RESULTS member of a _MYO triplet.
%
%   OUTPUTS
%     iv        1 x n struct array of window elements, each with .channelName.  Empty
%               (the interval prototype, plus .channelName) when the recording has no
%               windows or is not a myograph result at all.
%     channels  1 x m cellstr, the channel names in tree order; {} for a pressure
%               myograph.  It is the axis a consumer stratifies by.
%
%   DEPENDS ON
%     myographProduct (Core/Myograph) for the interval prototype.
%
% See also: cutMyographSource, setMyographIntervals, runMyographVasomotion,
%           myographProduct, exportToExcel, guiExplore
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 02-August-2026

%------------- BEGIN CODE --------------
function [iv,channels] = myographIntervals(results)

iv=emptyFlat(); channels={};
if ~isstruct(results) || ~isscalar(results), return; end

% THE WIRE MYOGRAPH FIRST, because a recording that carries results.channel is one -
% and one that has been through the interval editor carries an empty results.intervals
% beside it, so testing intervals first would answer "no windows" for a file that has
% plenty.
ch=fieldOr(results,'channel',[]);
if ~isempty(ch) && isstruct(ch)
    for i=1:1:numel(ch)
        nm=char(fieldOr(ch(i),'name',sprintf('channel%d',i)));
        channels{end+1}=nm; %#ok<AGROW>
        iv=append(iv,fieldOr(ch(i),'intervals',[]),nm);
    end
    return
end

iv=append(iv,fieldOr(results,'intervals',[]),'');
end

% =====================================================================
function iv=append(iv,src,channelName)
%append  One channel's windows onto the flat list, with the channel written on each.
%   The prototype is grown field by field rather than by concatenation, so a tree
%   written before a branch existed still reads back.
if isempty(src) || ~isstruct(src), return; end
fn=fieldnames(iv);
for k=1:1:numel(src)
    e=iv([]);                                   % one element with every field, empty
    e(1).channelName=char(channelName);
    for f=1:1:numel(fn)
        if strcmp(fn{f},'channelName'), continue; end
        if isfield(src(k),fn{f}), e(1).(fn{f})=src(k).(fn{f}); else, e(1).(fn{f})=[]; end
    end
    iv(end+1)=e; %#ok<AGROW>
end
end

% =====================================================================
function iv=emptyFlat()
%emptyFlat  The interval prototype with the channel axis written onto it.  Built by
%   naming the fields rather than by assigning one onto an empty array: a 0x0 struct
%   array has no element to assign to, so [iv(:).channelName]=deal('') would add
%   nothing and the field would appear only once a window existed.
fn=[reshape(fieldnames(myographProduct('intervals')),1,[]),{'channelName'}];
args=[fn;repmat({{}},1,numel(fn))];
iv=struct(args{:});
end

% =====================================================================
function v=fieldOr(st,f,d)
if nargin<3, d=[]; end
if isstruct(st) && isfield(st,f), v=st.(f); else, v=d; end
end
