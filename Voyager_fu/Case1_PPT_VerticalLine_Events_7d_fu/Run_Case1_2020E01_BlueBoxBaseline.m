function result = Run_Case1_2020E01_BlueBoxBaseline
%Run_Case1_2020E01_BlueBoxBaseline User-selected July 22--26 baseline.
%   Subtract the separate five-day S1--S7 means, retain the current centered
%   three-day mean, then blank negative results per sector. Never merge
%   sectors or enlarge neighboring cells over a masked sector.
%   Source CDFs and the previously published figure are preserved.
%   Modified: 2026-09-09

%% requested time interval and current display settings
baselineUTC = [datetime(2020,7,22,'TimeZone','UTC'), ...
    datetime(2020,7,27,'TimeZone','UTC')];
outputRoot = ['C:\Users\Administrator\Documents\', ...
    'Recovery-Work-Voyager_betatron\2020E01_BlueBoxBaseline'];

%% original CDF processing and existing a--f plotting program
result = Run_Case1_2020E01_BottomVariants( ...
    'Modes','pad_difference','BaselineUTC',baselineUTC, ...
    'PADDisplayAverageDays',3,'DifferenceNegativeAsMissing',true, ...
    'DifferenceColorLimits',[0 4.8644292],'OutputRoot',outputRoot);
assert(result.NegativeDifferenceCells==0);
fprintf('BLUE_BOX_BASELINE_NEGATIVES_BLANK_VERIFIED\n');
end
