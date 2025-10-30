%In MATLAB, I have a 3D matrix named source.data and a time vector named 
% source.time. The source.data represents Blood Flow Index (BFI) obtained 
% using laser speckle contrast imaging. Consider that some pixels might 
% have abnormal values – very low (e.g. 0), very high (infinitiy or simply 
% many times higher than all of the others) or NaNs. Make sure the further 
% analysis corrects for it. The data consist of baseline period of 5 
% seconds, which is followed by a stimulation and response period. Find 
% the peak response period and generate an image of response divided by 
% the baseline. Highlight the area where response is above 10%.  Side by 
% side make a plot of the average signal in the peak response area versus 
% the remaining area, both normalised by respective baseline values. 
% Explain how you found the peak of the response

%% --- INPUTS ------------------------------------------------------------
%  source.data  : Y × X × T  (BFI, arbitrary units)
%  source.time  : T × 1      (seconds, monotonic increasing)

winPeak    = 0.5;                  % seconds, width of search window
peakPadding= 0.25;                 % seconds either side of the peak
nanThresh = prctile(source.data(~isnan(source.data)), 99.9)*3;

%% --- 1. Basic cleaning -------------------------------------------------
data = source.data;                % shorthand
data(data==0 | data==Inf) = NaN;   % hard invalids
data(data > nanThresh)     = NaN;  % extreme outliers → NaN

%% --- 2. Baseline statistics -------------------------------------------
t          = source.time(:);
idxBase    =  t < 5;               % first 5 s
baseBFI    = mean(data(:,:,idxBase),3,'omitnan');      % Y×X
baseMean   = squeeze(mean(mean(data(:,:,idxBase),1,'omitnan'),2));

%% --- 3. Find the peak response period ---------------------------------
% Slide a window and compute mean ΔBFI/Baseline inside each window.
dt         = mean(diff(t));                           
wFrames    = round(winPeak/dt);                        % window size in frames
pctChange  = zeros(numel(t),1);

for k = 1:numel(t)-wFrames+1
    seg      = mean(data(:,:,k:k+wFrames-1),3,'omitnan');
    pctChange(k) = mean( (seg(:)-baseBFI(:))./baseBFI(:) ,'omitnan');
end

[~, idxPeak] = max(pctChange);                         % global maximum
tPeakStart   = t(idxPeak) - peakPadding;
tPeakEnd     = t(idxPeak+wFrames-1) + peakPadding;
idxPeakWin   = (t>=tPeakStart & t<=tPeakEnd);

%% --- 4. Percent-change image & ROI ------------------------------------
peakBFI   = mean(data(:,:,idxPeakWin),3,'omitnan');
respImg   = (peakBFI - baseBFI) ./ baseBFI;            % ΔBFI / BFI₀

roiMask   = respImg > 0.10;                            % >10 % increase
otherMask = respImg >= 0 & ~roiMask;                   % valid but ≤10 %

%% --- 5. Normalised time-courses ---------------------------------------
% --- ROI time-course ----------------------------------------------------
roiCurve = squeeze(mean(mean(data .* roiMask, 1,'omitnan'),2));  % 1×T
roiCurve = roiCurve(:);                                          % T×1

% --- Remaining-area time-course ----------------------------------------
otherCurve = squeeze(mean(mean(data .* otherMask,1,'omitnan'),2));
otherCurve = otherCurve(:);

% normalise by their own baselines
roiCurveN   = roiCurve   ./ mean(roiCurve(idxBase));
otherCurveN = otherCurve ./ mean(otherCurve(idxBase));

%% --- 6. Visualisation --------------------------------------------------
figure('Name','BFI response','Color','w','Position',[100 100 1200 450])

% ---- 6a. response/baseline image ----
subplot(1,2,1)
imagesc(respImg), axis image off
colormap(gca,parula); colorbar
hold on; visboundaries(roiMask,'Color','w','LineWidth',0.6)
title(sprintf('Peak ΔBFI/BFI₀ (%.1f–%.1f s)',tPeakStart,tPeakEnd))

% ---- 6b. time-courses ----
subplot(1,2,2)
plot(t, roiCurveN ,'r','LineWidth',1.2), hold on
plot(t, otherCurveN,'k','LineWidth',1.0)
xline(5,'--','Baseline end','LabelVerticalAlignment','bottom')
xline([tPeakStart tPeakEnd],'b:','Peak window')
xlabel('Time  (s)'), ylabel('Normalised BFI')
legend({'>10 % ROI','≤10 % area'},'Location','best'), grid on
title('Normalised mean BFI in and out of peak-response ROI')
