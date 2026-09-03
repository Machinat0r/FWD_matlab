function report = Case1_Test_Empty_Epoch_Selection()
%Case1_Test_Empty_Epoch_Selection Single-CDF empty exclusion-table regression.
%   Synthetic records only. No other event or satellite source is inspected.

%% one source and a window with no original Epoch
cfg = Case1_Config;
startTime = datetime(2004, 8, 8, 'TimeZone', 'UTC');
endTime = startTime+days(4);
data = struct;
data.Epoch = startTime+days([-1; 5]);
data.DeltaT = [40*86400; 86400];
data.SourceFileIndex = [1; 1];
data.SourceRecordNumber = [101; 102];
data.SourceManifest = table("synthetic.cdf", 'VariableNames', {'SourceFile'});
original = data;
[keep, audit] = Case1_Select_LECP_Epoch(data, startTime, endTime);
assert(~any(keep) && audit.InputRecords == 0 && audit.RetainedRecords == 0);
assert(height(audit.ExcludedRecords) == 0);
assert(isequal(size(audit.ExcludedRecords.SourceCDF), [0 1]));
assert(isequal(size(audit.ExcludedRecords.SourceCDFRecord), [0 1]));
assert(isequaln(data, original));

%% one original record, no rejected rows and unchanged long duration
data.Epoch = startTime+hours(1);
data.DeltaT = 40*86400;
data.SourceFileIndex = 1;
data.SourceRecordNumber = 101;
[keep, audit] = Case1_Select_LECP_Epoch(data, startTime, endTime);
assert(keep && audit.RetainedRecords == 1 && audit.NegativeDeltaTRejected == 0);
assert(isequal(size(audit.ExcludedRecords.SourceCDF), [0 1]));
assert(isequal(size(audit.ExcludedRecords.SourceCDFRecord), [0 1]));

%% a negative duration retains the exact excluded source identity
data.DeltaT = -1;
[keep, audit] = Case1_Select_LECP_Epoch(data, startTime, endTime);
assert(~keep && audit.NegativeDeltaTRejected == 1);
assert(audit.ExcludedRecords.SourceCDF == "synthetic.cdf");
assert(audit.ExcludedRecords.SourceCDFRecord == 101);
assert(audit.ExcludedRecords.EpochUTC == data.Epoch);
assert(audit.ExcludedRecords.DeltaT_s == -1);

%% numerical validation stays under the Voyager archive
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.Scope = 'Synthetic single-file empty window, no rejected rows, and negative-duration provenance.';
report.Passed = true;
report.AssertionCount = 13;
report.CodeFile = string(fullfile(cfg.CodeRoot, 'Case1_Select_LECP_Epoch.m'));
report.CodeSHA256 = string(Case1_File_SHA256(char(report.CodeFile)));
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'florinski_figure4');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
report.AuditFile = string(fullfile(folder, ['empty_epoch_selection_', stamp, '.mat']));
save(report.AuditFile, 'report');
fprintf('Empty-epoch selection regression: 13/13 passed.\n');
end
