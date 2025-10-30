%In MATLAB, I have a 3D matrix named source.data and a time vector named 
% source.time. The source.data represents Blood Flow Index (BFI) obtained 
% using laser speckle contrast imaging. The time is in seconds. I want to 
% plot the image of that corresponds to the average of source.data in time. 
% The colour limits should account for the possibility of the image 
% containing some pixels with extremely high, extremely low or NaN values 
% in respect to the meaningful data. Then I want to draw multiple regions 
% of interest (ROI) interactively on the image. Finally I want all the ROI 
% signals, which are averaged over pixels for individual regions, to be 
% plotted versus time in two panels: first as a plot of actual values, 
% second as a plot of normalised values, e.g. with each signal normalised 
% from 0 to 1


%% 1. Robust mean-BFI image ------------------------------------------------
data  = source.data;                       % Y × X × T
t     = source.time(:);                    % T × 1  (s)

% Discard obvious outliers before display
data(data==0 | data==Inf) = NaN;
cleanVals = data(:);
cleanVals(~isfinite(cleanVals)) = [];   % removes NaN and Inf
pLow  = prctile(cleanVals, 1);          % 1st percentile
pHigh = prctile(cleanVals, 99);         % 99th percentile
imgBFI = mean(data,3,'omitnan');                % Y × X

figure('Color','w'), ax = axes;
imagesc(imgBFI,[pLow pHigh]), axis image off
colormap(ax,parula), colorbar
title('Mean BFI (robust colour limits)')

%% 2. Interactive ROI definition ------------------------------------------
disp('Draw ROIs; double-click to finish each.  Press Esc when done.')
ROIs = {};
while true
    h = drawpolygon(ax,'LineWidth',1,'Color','y');      % user draws ROI
    if isempty(h) || isempty(h.Position)                % <-- pressed Enter
        delete(h);                                      % tidy up empty ROI
        break                                           % exit the loop
    end

    ROIs{end+1} = poly2mask(h.Position(:,1), ...
                            h.Position(:,2), ...
                            size(imgBFI,1), size(imgBFI,2));
end

%% 3. Extract ROI time-courses --------------------------------------------
nROI=numel(ROIs);
sig      = nan(numel(t), nROI);   % each column one ROI
for k = 1:nROI
    msk = ROIs{k};
    pix = reshape(data(repmat(msk,1,1,size(data,3))),[],size(data,3));
    sig(:,k) = mean(pix,'omitnan')';           % T × 1
end

%% 4. Plot raw and 0-to-1 normalised signals ------------------------------
figure('Color','w','Position',[100 100 900 400])

% 4a Raw values
subplot(1,2,1)
plot(t, sig,'LineWidth',1.1), grid on
xlabel('Time (s)'), ylabel('BFI'), title('Raw ROI signals')
legend(compose('ROI %d',1:nROI),'Location','best')

% 4b Normalised (min–max per ROI → 0–1)
sigN = (sig - min(sig,[],1)) ./ max(sig-min(sig,[],1),[],1);
subplot(1,2,2)
plot(t, sigN,'LineWidth',1.1), grid on
xlabel('Time (s)'), ylabel('Normalised BFI (0–1)')
title('Min–max normalised ROI signals')
