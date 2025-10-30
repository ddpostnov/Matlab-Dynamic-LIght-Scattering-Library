function L = buildOverlapTable(vsMap, thresh)
% vsMap  : X×Y×N int array, 0 = background
% thresh : overlap threshold in [0,1]     (default 0.51)
% L      : R×N lookup table (R = #objects in reference slice)

if nargin < 2, thresh = 0.51; end
[X,Y,N] = size(vsMap);

%% 1. find reference slice
nObj = squeeze(sum(vsMap ~= 0, [1 2]));     % any non-zero pixel counts
nLabels = arrayfun(@(k) numel(unique(vsMap(:,:,k))) - 1, 1:N);
[~, refIdx] = max(nLabels);                 % index of reference page
refLabels  = setdiff(unique(vsMap(:,:,refIdx)), 0);
R          = numel(refLabels);

%% 2. precompute pixel lists for every object in every slice
objPix = cell(N,1);               % objPix{k} – containers.Map(label → linear indices)
for k = 1:N
    lbls = setdiff(unique(vsMap(:,:,k)), 0);
    mp   = containers.Map('KeyType','double','ValueType','any');
    for l = lbls.'
        mp(l) = find(vsMap(:,:,k) == l);
    end
    objPix{k} = mp;
end

%% 3. build lookup table
L = zeros(R,N);                   % initialise with 0 = "no match"
L(:,refIdx) = refLabels;          % trivial matches in reference slice

for r = 1:R
    refMaskIdx = objPix{refIdx}(refLabels(r));
    refArea    = numel(refMaskIdx);

    for k = setdiff(1:N, refIdx)
        mp = objPix{k};           % map for slice k
        cand = mp.keys;
        bestOv = 0; bestLbl = 0;

        for c = cand
            candIdx   = mp(c{1});
            interArea = numel(intersect(refMaskIdx, candIdx));
            ovRatio   = max(interArea/refArea, interArea/numel(candIdx));   % ≥ 51 % of candidate object
            if ovRatio >= thresh && ovRatio > bestOv
                bestOv = ovRatio; bestLbl = c{1};
            end
        end
        L(r,k) = bestLbl;         % 0 if no object reached threshold
    end
end
end