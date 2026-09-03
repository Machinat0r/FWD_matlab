function report = Case1_Test_L1_Fallback()
%Case1_Test_L1_Fallback Synthetic L1 averaging/provenance/priority tests.
%   Only the requested 2004 rate CDF is opened, as a metadata template for
%   three explicitly synthetic tiny CDFs. All test data stay in the archive.

%% synthetic products and hourly policy
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'l1_fallback');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS'));
fixtureFolder = fullfile(folder, ['synthetic_', stamp]);
mkdir(fixtureFolder);
t0 = datetime(2004, 8, 8, 'TimeZone', 'UTC');
l2 = mockProduct(t0+minutes([10; 20; 70; 190; 230; 300]), false);
l2.FHDU_SectoredFluxes(:, 10, :) = repmat(reshape(11:18, 1, 1, 8), 6, 1, 1);
l2.FHDU_SectoredFluxes([2, 3, 4], 10, 7) = NaN;
l2.DeltaT = [7200; NaN; 3600; 0; -1; 3600];
l1 = mockProduct(t0+minutes([5; 65; 80; 90; 120; 170; 190; 245; 300; -1]), true);
amplitude = [100; 0; 2; 1e8; 3; 5; 4; 0; 1e8; 1e8];
l1.FHDU_SectoredRates(:, 10, :) = repmat(reshape(amplitude, 10, 1, 1), 1, 1, 8);
l1.FHDU_SectoredRateUncertainties(:, 10, :) = 0.2;
l1.FHDU_SectoredRateUncertainties(3, 10, 2) = NaN;
l1.FHDU_SectoredRates(7, 10, 7) = NaN;
l1.FHDU_SectoredRates([5, 6], 10, 8) = NaN;
l1.DeltaT = [4000; 5000; 90000; -1; NaN; 0; 3600; 3600; 3600; 3600];
l1.SourceFileIndex(3) = 2;
l1.SourceManifest = [l1.SourceManifest; l1.SourceManifest];
l1.SourceManifest.SourceFile(2) = "synthetic_rates_second.cdf";
l1.SourceMetadata{2} = l1.SourceMetadata{1};
originalL2 = l2; originalL1 = l1;
out = Case1_Apply_L1_Fallback(l2, l1, t0, t0+hours(5), 'hour');
nCheck = 0;
check(isequaln(l2, originalL2) && isequaln(l1, originalL1));
check(numel(out.Epoch) == 5);
check(isequaln(out.OriginalL2SourceRow, [1; 2; NaN; NaN; 4]));
check(isequal(out.Epoch, t0+minutes([10; 20; 90; 150; 190])));
check(isequal(out.SourceProduct, ["L2"; "L2"; "L1_UTC_mean"; "L1_UTC_mean"; "L2"]));
factor = 1/(0.44*(1.78-0.57));
check(max(abs(reshape(out.FHDU_SectoredFluxes(3, 10, :), 1, 8)-factor)) < 1e-12);
check(max(abs(reshape(out.FHDU_SectoredFluxes(4, 10, 1:7), 1, 7)-4*factor)) < 1e-12);
check(isnan(out.FHDU_SectoredFluxes(4, 10, 8)));
check(all(out.SectorSampleCount(3, :) == 2));
check(isequal(out.SectorSampleCount(4, :), [2, 2, 2, 2, 2, 2, 2, 0]));
check(abs(out.FHDU_SectoredFluxUncertainties(3, 10, 1)-sqrt(2)*0.2/2*factor) < 1e-12);
check(isnan(out.FHDU_SectoredFluxUncertainties(3, 10, 2)));
check(all(isnan(out.DeltaT(3:4))) && all(isnan(out.SourceRecordNumber(3:4))));
check(out.SourceToDifferentialFluxFactor(1) == 1 && out.SourceToDifferentialFluxFactor(3) == factor);
check(isequaln(out.FHDU_SectoredFluxes([1, 2, 5], :, :), l2.FHDU_SectoredFluxes([1, 2, 4], :, :)));
check(isequaln(out.FHDU_SectoredFluxUncertainties([1, 2, 5], :, :), l2.FHDU_SectoredFluxUncertainties([1, 2, 4], :, :)));
check(isequaln(out.DeltaT([1, 2, 5]), l2.DeltaT([1, 2, 4])));
check(isequaln(out.FHDU_EnergyRange([1, 2, 5], :, :), l2.FHDU_EnergyRange([1, 2, 4], :, :)));
check(isequal(out.FHDU_EnergyRange(3, :, 10), [0.57, 0.89]));
check(all(isnan(out.FHDU_SectoredFluxes(3, [1:9, 11:16], :)), 'all'));
check(isequal(out.L1SourceRecords{3}.SourceRow, [2; 3]));
check(isequal(out.L1SourceRecords{3}.SourceCDFRecord, [102; 103]));
check(isequal(out.L1SourceRecords{3}.SourceCDF, ["synthetic_rates.cdf"; "synthetic_rates_second.cdf"]));
check(isequal(out.L1SourceRecords{3}.DeltaT_s, [5000; 90000]));
check(all(out.L1SourceRecords{3}.Contributes_S1_to_S8, 'all'));
audit = out.L1FallbackAudit;
check(audit.L1AddedRecords == 2 && audit.L2ReplacedRecords == 1 && audit.L2PreservedRecords == 3);
check(audit.L1RecordSelection.NegativeDeltaTRejected == 1);
check(audit.L2RecordSelection.NegativeDeltaTRejected == 1);
check(audit.L1RecordSelection.MissingDeltaTRetained == 1);
check(audit.ReplacedL2.OriginalRows == 3);
check(isequaln(audit.ReplacedL2.FHDU_SectoredFluxes, l2.FHDU_SectoredFluxes(3, :, :)));
check(audit.Candidates.Decision(1) == "preserve_all_L2_rows_complete_L2_in_bin");
check(audit.Candidates.Decision(4) == "L1_mean_missing_or_nonpositive_sector");
check(audit.Candidates.Decision(5) == "L1_mean_missing_or_nonpositive_sector");
check(~any(cellfun(@(r) any(r == 4 | r == 9 | r == 10), audit.Candidates.L1Rows)));
check(out.BinStartUTC(3) == t0+hours(1) && out.BinEndUTC(3) == t0+hours(2));
check(isequal(out.SectorSampleCount([1, 2, 5], :), double(isfinite(reshape(l2.FHDU_SectoredFluxes([1, 2, 4], 10, :), 3, 8)))));

%% no L1, zero sigma, day center and empty window
empty = Case1_Apply_L1_Fallback(l2, l1, t0+days(9), t0+days(10), 'day');
check(isempty(empty.Epoch) && height(empty.L1FallbackAudit.Candidates) == 0);
dayL2 = mockProduct(t0-days(1), false);
dayL1 = mockProduct(t0+hours([1; 23]), true);
dayL1.FHDU_SectoredRates(:, 10, :) = 2;
dayL1.FHDU_SectoredRateUncertainties(:, 10, :) = 0;
dayOut = Case1_Apply_L1_Fallback(dayL2, dayL1, t0, t0+days(1), 'day');
check(dayOut.Epoch == t0+hours(12));
check(all(dayOut.FHDU_SectoredFluxUncertainties(1, 10, :) == 0, 'all'));
check(all(dayOut.FHDU_SectoredFluxes(1, 10, :) == 2*factor, 'all'));
noRate = mockProduct(t0-days(1), true);
unchanged = Case1_Apply_L1_Fallback(l2, noRate, t0, t0+hours(5), 'hour');
check(isequaln(unchanged.FHDU_SectoredFluxes, l2.FHDU_SectoredFluxes(1:4, :, :)));
check(any(unchanged.L1FallbackAudit.Candidates.Decision == "no_L1_Epoch_in_bin"));

%% tiny synthetic CDFs exercise the production reader directly
files = writeRateFixtures(cfg, fixtureFolder);
single = Case1_Read_LECP_Rates(char(files(1)));
check(numel(single.Epoch) == 2 && isequal(single.SourceRecordNumber, [1; 2]));
check(isnan(single.FHDU_SectoredRates(1, 10, 1)));
check(single.FHDU_SectoredRates(1, 10, 2) == 0);
check(single.SourceManifest.SHA256 == string(Case1_File_SHA256(char(files(1)))));
doubleRead = Case1_Read_LECP_Rates(files(1:2));
check(numel(doubleRead.Epoch) == 2 && doubleRead.DuplicateIdenticalRecordsRemoved == 2);
check(height(doubleRead.SourceRecordAliases{1}) == 2);
check(isequal(doubleRead.SourceRecordAliases{1}.SourceCDF, files(1:2)));
conflictCaught = false;
try
    Case1_Read_LECP_Rates(files([1, 3]));
catch ME
    conflictCaught = strcmp(ME.identifier, 'VoyagerCase1:ConflictingL1Epoch');
end
check(conflictCaught);

%% retain full test evidence and exact code hashes
report = struct('Passed', true, 'AssertionCount', nCheck, ...
    'CreatedUTC', datetime('now', 'TimeZone', 'UTC'), ...
    'SyntheticCDFs', files, 'Output', out, 'DayOutput', dayOut, ...
    'L2Input', l2, 'L1Input', l1);
codeFiles = string(fullfile(cfg.CodeRoot, ...
    {'Case1_Read_LECP_Rates.m'; 'Case1_Apply_L1_Fallback.m'; 'Case1_Test_L1_Fallback.m'}));
hash = strings(3, 1);
for ii = 1:3, hash(ii) = Case1_File_SHA256(char(codeFiles(ii))); end
report.CodeManifest = table(codeFiles, hash, 'VariableNames', {'CodeFile', 'SHA256'});
report.AuditFile = string(fullfile(folder, ['l1_fallback_synthetic_', stamp, '.mat']));
save(report.AuditFile, 'report');
fprintf('L1 fallback regression: %d/%d passed.\n', nCheck, nCheck);

    function check(value)
        nCheck = nCheck+1;
        assert(isscalar(value) && value, 'L1 fallback assertion %d failed.', nCheck);
    end
end

function p = mockProduct(epoch, isRate)
n = numel(epoch);
p = struct;
p.Epoch = epoch(:); p.DeltaT = ones(n, 1)*3600;
p.FHDU_Energy = nan(n, 16); p.FHDU_Energy(:, 10) = 0.73;
p.FHDU_EnergyRange = nan(n, 2, 16);
p.FHDU_EnergyRange(:, :, 10) = repmat([0.57, 0.89], n, 1);
p.variable_meta = struct('Synthetic', true);
p.SourceMetadata = {struct('Synthetic', true)};
p.SourceFileIndex = ones(n, 1); p.SourceRecordNumber = (101:100+n).';
p.SectorIterator = 1:8;
if isRate
    p.FHDU_SectoredRates = ones(n, 16, 8);
    p.FHDU_SectoredRateUncertainties = ones(n, 16, 8)*0.1;
    filename = "synthetic_rates.cdf"; cadence = "native";
else
    p.FHDU_SectoredFluxes = ones(n, 16, 8);
    p.FHDU_SectoredFluxUncertainties = ones(n, 16, 8)*0.1;
    filename = "synthetic_flux.cdf"; cadence = "hour";
end
p.SourceManifest = table(filename, "synthetic_SHA", n, cadence, ...
    'VariableNames', {'SourceFile', 'SHA256', 'Records', 'Cadence'});
p.ProductCadence = cadence;
end

function files = writeRateFixtures(cfg, folder)
template = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'native', 'l1', ...
    'sectored_rates', '2004', 'voyager-1_lecp_lev-1-rates_20040101_v1.1.1-01.cdf');
cdfObj = dataobj(template);
recordNames = {'Epoch', 'DeltaT', 'FHDU_SectoredRates', ...
    'FHDU_SectoredRateUncertainties', 'FHDU_Energy', 'FHDU_EnergyRange'};
names = [recordNames, {'SectorIterator', 'Hydrogen_Channels', 'Hydrogen_Channels_Label'}];
data = struct;
for ii = 1:numel(names)
    name = names{ii}; data.(name) = cdfObj.data.(name);
    if ismember(name, recordNames)
        data.(name).data = data.(name).data(1:2, :, :);
    end
end
data.FHDU_SectoredRates.data(:) = 1;
rateMeta = getv(cdfObj, 'FHDU_SectoredRates');
data.FHDU_SectoredRates.data(1, 10, 1) = rateMeta.FILLVAL;
data.FHDU_SectoredRates.data(1, 10, 2) = 0;
data.FHDU_SectoredRateUncertainties.data(:) = 0.1;
data.FHDU_Energy.data(:, 10) = 0.73;
variables = cdfObj.Variables(ismember(cdfObj.Variables(:, 1), names), :);
attrs = cdfObj.VariableAttributes;
attrNames = fieldnames(attrs);
for ii = 1:numel(attrNames)
    value = attrs.(attrNames{ii});
    attrs.(attrNames{ii}) = value(ismember(value(:, 1), names), :);
end
globalAttrs = cdfObj.GlobalAttributes;
globalAttrs.TEXT = {'SYNTHETIC regression fixture. Not observational data.'};
files = string(fullfile(folder, {'synthetic_rates_a.cdf'; ...
    'synthetic_rates_b_identical.cdf'; 'synthetic_rates_c_conflict.cdf'}));
for ii = 1:3
    if ii == 3, data.FHDU_SectoredRates.data(2, 10, 2) = 99; end
    irf.cdf.write_dataobj(char(files(ii)), globalAttrs, data, attrs, variables);
end
end
