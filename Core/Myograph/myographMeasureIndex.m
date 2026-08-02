%myographMeasureIndex  Which of the three diameter measures a setting names
%
%   [k,name] = myographMeasureIndex(edgeMode) returns the index K into the third
%   dimension of a getMyographDiameter idxL/idxR/diameter array, and the canonical
%   NAME of that measure, for the measure named by EDGEMODE.
%
%   getMyographDiameter always returns all three measures, ordered
%   {'outer','mid','inner'} (outer-wall, wall-centre and luminal diameter).
%   s.edgeMode no longer selects which one is MEASURED - it names the one that is
%   ANALYSED and plotted by default.  This function is the single place that maps
%   the setting onto that choice, including the two historical spellings:
%
%     'outer'            -> 1  outer      (widest)
%     'mid'              -> 2  mid        (default)
%     'center'           -> 2  mid        (the pre-2026-08 name of the same measure)
%     'min'              -> 2  mid        (the same measure, detected by the darkest-
%                                          point rule instead of the centroid)
%     'inner'            -> 3  inner      (the lumen, narrowest)
%
%   [k,name] = myographMeasureIndex(edgeMode,measures) matches against the
%   MEASURES cellstr carried by the interval (interval.measures) instead of the
%   built-in order, so a result loaded from disk is read by its own labels.
%
%   An unrecognised name falls back to the wall-centre measure with a warning: a
%   mistyped setting must not silently analyse a different diameter.
%
%   Syntax:
%      k        = myographMeasureIndex(s.edgeMode)
%      [k,name] = myographMeasureIndex(s.edgeMode,interval.measures)
%
% Author: Dmitry D Postnov, CFIN, Aarhus University (dpostnov@cfin.au.dk)
% Copyright 2026 Dmitry D Postnov, Aarhus University.
% Header generation and script formatting were done with Claude Code.
% Last revision: 01-August-2026

function [k,name] = myographMeasureIndex(edgeMode,measures)

if nargin<2 || isempty(measures), measures={'outer','mid','inner'}; end
if isempty(edgeMode), edgeMode='mid'; end
edgeMode=char(edgeMode);

%the two historical spellings of the wall-centre measure
if any(strcmpi(edgeMode,{'center','centre','min'})), edgeMode='mid'; end

k=find(strcmpi(measures,edgeMode),1);
if isempty(k)
    k=find(strcmpi(measures,'mid'),1);
    if isempty(k), k=min(2,numel(measures)); end
    warning('myographMeasureIndex:unknownMeasure', ...
        'Unknown diameter measure ''%s''; using ''%s''.',edgeMode,measures{k});
end
name=measures{k};
end
