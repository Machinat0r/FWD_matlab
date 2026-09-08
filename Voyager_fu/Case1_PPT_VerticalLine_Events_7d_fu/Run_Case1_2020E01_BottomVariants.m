function result = Run_Case1_2020E01_BottomVariants(varargin)
%Run_Case1_2020E01_BottomVariants Two requested daily overview variants.
%   Scientific inputs are original website CDF files on Z:. Result tables
%   and MAT audits are outputs only. Default window: 2020-06-22--09-08 UTC.
%   Baseline: lowest five-day mean with center near July 3 (July 1--5).
%   All original sector values remain independent in both figures.
%   Modified: 2026-09-08

%% paths and requested modes
parser = inputParser;
addParameter(parser, 'Modes', ["pad_difference", "sector"]);
addParameter(parser, 'SectorColorMode', 'normalized');
parse(parser, varargin{:});
options = parser.Results;
modes = string(options.Modes);
assert(all(ismember(modes, ["pad_difference", "sector"])));
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
root = fullfile(fileparts(cfg.OutputRoot), '2020E01_Variants');
reportFolder = fullfile(root, 'report');
if ~isfolder(reportFolder), mkdir(reportFolder); end
diary(fullfile(reportFolder, 'run.log'));
logCleanup = onCleanup(@() diary('off'));

dailyFiles = cfg.LECPNativeDailyCDFs(contains(cfg.LECPNativeDailyCDFs, ...
    [filesep, '2020', filesep]));
rateFiles = cfg.LECPLevel1CDFs(contains(cfg.LECPLevel1CDFs, ...
    [filesep, '2020', filesep]));
assert(isscalar(dailyFiles) && isscalar(rateFiles));
rawFiles = [dailyFiles; rateFiles];
for month = 6:9
    rawFiles(end+1,1) = string(fullfile(cfg.DataRoot, 'voyager1', 'coho', ...
        '1hr', 'l2', 'merged_mag_plasma', '2020', sprintf('%02d',month), ...
        sprintf('voyager1_coho1hr_merged_mag_plasma_2020%02d01_v01.cdf',month)));
end
beforeHash = strings(size(rawFiles));
for k = 1:numel(rawFiles), beforeHash(k) = Case1_File_SHA256(rawFiles(k)); end

%% select the low-flux five-day baseline from original CDFs
l2 = Case1_Read_LECP_CDFs(dailyFiles);
l1 = Case1_Read_LECP_Rates(rateFiles);
nearStart = datetime(2020,6,29,'TimeZone','UTC');
nearEnd = datetime(2020,7,8,'TimeZone','UTC');
p = Case1_Apply_L1_Fallback(l2, l1, nearStart, nearEnd, 'day', 'l1_first');
flux = reshape(p.FHDU_SectoredFluxes(:,10,1:7), numel(p.Epoch), 7);
flux(~isfinite(flux) | flux <= 0) = NaN;
candidateStart = (nearStart:days(1):datetime(2020,7,3,'TimeZone','UTC')).';
candidateEnd = candidateStart + days(5);
score = nan(size(candidateStart));
sectorMean = nan(numel(score),7);
counts = zeros(numel(score),7);
for k = 1:numel(score)
    rows = p.Epoch >= candidateStart(k) & p.Epoch < candidateEnd(k);
    values = flux(rows,:);
    counts(k,:) = sum(isfinite(values),1);
    sectorMean(k,:) = mean(values,1,'omitnan');
    score(k) = mean(sectorMean(k,:));
end
assert(all(isfinite(score)), 'The proposed local baseline candidates lack data.');
[~, best] = min(score);
baselineUTC = [candidateStart(best), candidateEnd(best)];
baselineCandidates = table(candidateStart, candidateEnd, score, sectorMean, ...
    counts, 'VariableNames', {'StartUTC','EndUTCExclusive', ...
    'MeanAcrossSevenSectors','SeparateSectorMean','SeparateSectorCount'});
writetable(baselineCandidates, fullfile(reportFolder,'baseline_candidates.csv'));
fprintf('Selected baseline: %s to %s UTC (exclusive end).\n', ...
    string(baselineUTC(1)),string(baselineUTC(2)));

%% same event catalog and current a--e renderer
catalog = table("2020-E01", "2020-07-22", "2020-08-07", "2020-08-08", ...
    'VariableNames', {'EventID','StartUTC','EndUTCInclusive','EndUTCExclusive'});
catalogFile = fullfile(reportFolder,'event_catalog.csv');
writetable(catalog,catalogFile);
common = { ...
    'Spacecraft',1, 'DataRoot',cfg.DataRoot, 'CatalogFile',catalogFile, ...
    'ContextMonths',1, 'ContextDays',[], 'ShowEventBoundaries',false, ...
    'Overwrite',true, 'Visible',false, 'ExportDPI',cfg.ExportDPI, ...
    'MAGGapHours',cfg.MAGGapHours, 'FluxGapHours',cfg.FluxGapHours, ...
    'LECPDailyAverage',false, 'CRSDisplay','lecp_p1_pitch_angle', ...
    'LECPNativeCDFs',dailyFiles, 'LECPLevel1Fallback',true, ...
    'LECPSourcePriority','l1_first', 'LECPLevel1CDFs',rateFiles, ...
    'PADCadence','day', 'AccumulationPolicy','epoch_drop_negative_deltat', ...
    'LECPSectorAverageDays',0, 'PitchAngleMethod','predicted_ck', ...
    'PredictedAttitudeApproved',true, 'NominalLECPGeometryApproved',true, ...
    'AttitudeDailyHourUTC',12, 'LECPBackgroundMode','none', ...
    'ExportPitchAngleTable',false, 'P1DisplayEnergyMeV',[0.57 1.78], ...
    'ExportPeakPAD',false, 'LineYMarginFraction',cfg.LineYMarginFraction, ...
    'ColorPercentiles',cfg.ColorPercentiles};

result = struct('BaselineUTC',baselineUTC, ...
    'BaselineCandidates',baselineCandidates, ...
    'GeneratedUTC',datetime('now','TimeZone','UTC'));
result.Reports = cell(numel(modes),1);
referenceUpper = [];
for k = 1:numel(modes)
    mode = modes(k);
    out = fullfile(root,char(mode));
    report = Voyager_Case1_Plot_Events(common{:}, ...
        'LECPBottomMode',char(mode), 'BaselineUTC',baselineUTC, ...
        'SectorColorMode',options.SectorColorMode, ...
        'OutputFolder',out, 'ReportFolder',fullfile(out,'report'), ...
        'ReportTag',char(mode), 'PitchAngleDataFolder',fullfile(out,'audit'));
    assert(string(report.Status(1)) == "ok", ...
        'Plot failed: %s',string(report.Notes(1)));
    a = load(string(report.PitchAngleAuditFile(1)),'pitchAngleTable','opts');
    T = a.pitchAngleTable;
    audit = T.Properties.UserData;
    assert(~audit.SectorMergeApplied);
    plotted = T.BottomPanelUsable;
    displayed = T{:,cellstr(compose('DisplayValue_S%d',1:7))};
    assert(all(sum(isfinite(displayed(plotted,:)),2)==7));
    if isempty(referenceUpper)
        referenceUpper = audit.UpperPanelLines;
    else
        assert(isequaln(referenceUpper,audit.UpperPanelLines), ...
            'One of the a--e line arrays or axis limits changed.');
    end
    J = T{:,cellstr(compose('Flux_S%d_1d',1:7))};
    if mode == "pad_difference"
        ref = T.EpochUTC>=baselineUTC(1) & T.EpochUTC<baselineUTC(2);
        b = mean(J(ref,:),1,'omitnan');
        assert(max(abs(b-sectorMean(best,:)))<1e-10);
        expected = J-b;
        assert(isequaln(displayed(plotted,:),expected(plotted,:)));
        result.NegativeDifferenceCells = nnz(displayed<0);
        result.BaselineSectorFlux = b;
    else
        assert(~audit.PitchAngleCalculated && ~audit.MagneticFieldRequired);
        assert(~any(startsWith(string(T.Properties.VariableNames),'PA_S')));
        if strcmp(options.SectorColorMode,'normalized')
            expected = J./max(J,[],2);
            assert(isequaln(displayed(plotted,:),expected(plotted,:)));
            assert(all(max(displayed(plotted,:),[],2)==1));
        else
            assert(isequaln(displayed(plotted,:),J(plotted,:)));
        end
    end
    result.Reports{k} = report;
    result.(char(mode)) = struct('FigureFile',string(report.FigureFile(1)), ...
        'AuditFile',string(report.PitchAngleAuditFile(1)), ...
        'UsableRecords',nnz(plotted));
end

%% source integrity and result audit
afterHash = strings(size(rawFiles));
for k = 1:numel(rawFiles), afterHash(k) = Case1_File_SHA256(rawFiles(k)); end
assert(isequal(beforeHash,afterHash),'Source bytes changed during rendering.');
result.SourceManifest = table(rawFiles,beforeHash,afterHash);
result.Modes = modes;
result.Options = options;
result.UpperPanelsIdenticalBetweenVariants = numel(modes)==2;
result.EntryPointSHA256 = Case1_File_SHA256([mfilename('fullpath'),'.m']);
save(fullfile(reportFolder,['result_',strjoin(cellstr(modes),'_'),'.mat']), ...
    'result','-v7.3');
writetable(result.SourceManifest,fullfile(reportFolder,'source_manifest.csv'));
disp(result);
clear logCleanup
end
