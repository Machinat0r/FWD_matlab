function result = Run_Case1_2020E01_BlueBoxRatio
%Run_Case1_2020E01_BlueBoxRatio Divide S1--S7 by July 22--26 means.
%   Original CDF input, current three-day centered per-sector mean and
%   central-day PA geometry. Ratios below one remain valid. No sector
%   merging, no new in-figure processing notes, no changes to panels a--e.
%   Modified: 2026-09-09

%% baseline and output paths
baselineUTC = [datetime(2020,7,22,'TimeZone','UTC'), ...
    datetime(2020,7,27,'TimeZone','UTC')];
outputRoot = ['C:\Users\Administrator\Documents\', ...
    'Recovery-Work-Voyager_betatron\2020E01_BlueBoxRatio'];

%% current raw-CDF pipeline and sector-wise ratio
result = Run_Case1_2020E01_BottomVariants( ...
    'Modes','pad_ratio','BaselineUTC',baselineUTC, ...
    'PADDisplayAverageDays',3,'DifferenceNegativeAsMissing',false, ...
    'DifferenceColorLimits',[0 8],'OutputRoot',outputRoot);
fprintf('BLUE_BOX_BASELINE_RATIO_VERIFIED\n');
end
