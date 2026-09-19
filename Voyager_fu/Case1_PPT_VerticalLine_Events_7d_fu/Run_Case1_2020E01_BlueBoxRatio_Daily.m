function result = Run_Case1_2020E01_BlueBoxRatio_Daily
%Run_Case1_2020E01_BlueBoxRatio_Daily Daily ratios without 3-day smoothing.
%   Read original CDFs through the current L1-first pipeline. Divide each
%   daily S1--S7 flux by its own July 22--26 mean. Retain daily PA and mask,
%   independent sectors and values below one. Panels a--e are unchanged.
%   Keep the 0--8 color scale for comparison with the three-day figure.
%   Modified: 2026-09-09

%% baseline and separate output folder
baselineUTC = [datetime(2020,7,22,'TimeZone','UTC'), ...
    datetime(2020,7,27,'TimeZone','UTC')];
outputRoot = ['C:\Users\Administrator\Documents\', ...
    'Recovery-Work-Voyager_betatron\2020E01_BlueBoxRatio_Daily'];

%% original daily values; no additional temporal smoothing
result = Run_Case1_2020E01_BottomVariants( ...
    'Modes','pad_ratio','BaselineUTC',baselineUTC, ...
    'PADDisplayAverageDays',1,'DifferenceNegativeAsMissing',false, ...
    'DifferenceColorLimits',[0 8],'OutputRoot',outputRoot);
fprintf('BLUE_BOX_BASELINE_DAILY_RATIO_VERIFIED\n');
end
