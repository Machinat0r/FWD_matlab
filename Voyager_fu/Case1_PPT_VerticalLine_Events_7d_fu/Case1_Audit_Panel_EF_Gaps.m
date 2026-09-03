function report = Case1_Audit_Panel_EF_Gaps()
%Case1_Audit_Panel_EF_Gaps Diagnose existing e/f gaps from original CDFs.
%   Read-only for production CDFs, plots, configuration and saved PAD tables.
%   UTC bin membership below counts available inputs; it does not rebin flux.
%   Original sector samples and the saved mean B are checked independently.

%% configuration and source readers already used by the production code
cfg = Case1_Config;
originalPath = path;
cleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
Case1_Add_IRFU_Path(cfg.IRFURoot);
catalog = Case1_Event_Catalog;
catalog = catalog(catalog.Spacecraft == 1, :);
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.TestFile = [mfilename('fullpath'), '.m'];
report.CodeSHA256 = Case1_File_SHA256(report.TestFile);
report.Policy = ['Diagnostic counts only: native Epoch membership in UTC bins; ', ...
    'no flux averaging/interpolation, new quality gate or figure change.'];
report.Records = table;
report.Bins = table;
report.Events = table;
report.Sources = cell(2, 1);
report.Peaks = table;
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
cadences = {'day', 'hour'};
folders = {'1d', '1h'};
sourceGroups = {cfg.LECPNativeDailyCDFs, cfg.LECPNativeHourlyCDFs};

%% all current V1 events, using independent original CDF reads
for ic = 1:2
    cadence = cadences{ic}; suffix = folders{ic};
    source = Case1_Read_LECP_CDFs(sourceGroups{ic});
    report.Sources{ic} = source.SourceManifest;
    auditFolder = fullfile(cfg.PitchAngleDataFolder, 'predicted_ck');
    if ic == 2
        auditFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1h', ...
            'derived', 'pitch_angle', '2013-2021', 'predicted_ck');
    end
    for ie = 1:height(catalog)
        eventID = string(catalog.EventID(ie));
        startTime = catalog.StartUTC(ie)-days(cfg.ContextDays);
        endTime = catalog.EndUTCExclusive(ie)+days(cfg.ContextDays);
        entry = dir(fullfile(auditFolder, ['V1_', char(eventID), ...
            '*_predictedCK_', suffix, '_nativeCDF_Epoch.mat']));
        assert(isscalar(entry), 'Expected one existing PAD audit per event.');
        saved = load(fullfile(entry.folder, entry.name));
        T = saved.pitchAngleTable;
        [keep, selection] = Case1_Select_LECP_Epoch(source, startTime, endTime);
        rows = find(keep); n = numel(rows);
        assert(n == height(T), 'Saved PAD/source record count mismatch.');
        C = table;
        files = string(saved.sourceMAG);
        for iFile = 1:numel(files)
            key = char(files(iFile));
            if ~isKey(cache, key)
                p = Voyager_Read_CDF_Product(key, 'coho');
                part = table(p.Epoch(:), p.BR(:), p.BT(:), p.BN(:), ...
                    p.protonFlux1_LECP(:), p.protonFlux2_LECP(:), p.protonFlux3_LECP(:), ...
                    repmat(string(key), numel(p.Epoch), 1), (1:numel(p.Epoch)).', ...
                    'VariableNames', {'EpochUTC','BR','BT','BN','P1','P2','P3', ...
                    'SourceCDF','SourceCDFRecord'});
                cache(key) = part;
            end
            C = [C; cache(key)]; %#ok<AGROW>
        end
        C = sortrows(C, 'EpochUTC');
        C = C(C.EpochUTC >= startTime & C.EpochUTC < endTime, :);
        cohoBin = dateshift(C.EpochUTC, 'start', cadence);
        sectorBin = dateshift(source.Epoch(rows), 'start', cadence);
        e1 = isfinite(C.P1) & C.P1 > 0;
        eAny = any(isfinite(C{:,{'P1','P2','P3'}}) & C{:,{'P1','P2','P3'}} > 0, 2);
        flux = reshape(source.FHDU_SectoredFluxes(rows, 10, 1:7), n, 7);
        flux(~isfinite(flux) | flux < 0) = NaN;
        fluxCount = sum(isfinite(flux) & flux > 0, 2);
        rawB = nan(n, 3); bSamples = zeros(n, 1);
        e1Present = false(n, 1); eAnyPresent = false(n, 1);
        paCount = zeros(n, 1); attitudeFound = false(n, 1);
        for ir = 1:n
            inBin = cohoBin == sectorBin(ir);
            e1Present(ir) = any(e1(inBin));
            eAnyPresent(ir) = any(eAny(inBin));
            goodB = inBin & all(isfinite(C{:,{'BR','BT','BN'}}), 2);
            bSamples(ir) = nnz(goodB);
            if any(goodB), rawB(ir, :) = mean(C{goodB,{'BR','BT','BN'}}, 1); end
        end
        if n > 0
            if ic == 1
                bFields = {'BR_daily_nT','BT_daily_nT','BN_daily_nT'};
            else
                bFields = {'BR_hourly_nT','BT_hourly_nT','BN_hourly_nT'};
            end
            actualFlux = nan(n, 7); pa = nan(n, 7);
            for is = 1:7
                actualFlux(:, is) = T.(sprintf('RawFlux_S%d_%s', is, suffix));
                pa(:, is) = T.(sprintf('PA_S%d_deg', is));
            end
            assert(isequal(T.EpochUTC, source.Epoch(rows)) && isequal(T.SourceRow, rows));
            assert(isequaln(actualFlux, flux), 'Saved flux differs from original CDF.');
            assert(isequaln(T{:, bFields}, rawB), 'B bin matching differs from original COHO.');
            assert(isequal(T.MAGVectorSampleCount, bSamples), 'B sample count mismatch.');
            paCount = sum(isfinite(pa), 2);
            attitudeFound = T.Properties.UserData.Attitude.Found;
        end
        bGood = all(isfinite(rawB), 2) & vecnorm(rawB, 2, 2) > 0;
        usable = fluxCount == 7 & paCount == 7;
        if n > 0, assert(isequal(T.PADUsable, usable)); end
        % Exclusive cause order for this table; individual flags also retained.
        cause = repmat("usable", n, 1);
        cause(paCount < 7) = "other_PA_unavailable";
        cause(~attitudeFound) = "attitude_transform_missing";
        cause(~bGood) = "magnetic_vector_missing";
        cause(fluxCount < 7) = "partial_sector_flux";
        cause(fluxCount == 0) = "no_positive_sector_flux";
        R = table(repmat(string(cadence), n, 1), repmat(eventID, n, 1), ...
            source.Epoch(rows), rows, sectorBin, fluxCount, bGood, bSamples, ...
            attitudeFound, paCount, usable, e1Present, eAnyPresent, cause, ...
            source.SourceManifest.SourceFile(source.SourceFileIndex(rows)), ...
            source.SourceRecordNumber(rows), ...
            'VariableNames', {'Cadence','EventID','EpochUTC','SourceRow','BinStartUTC', ...
            'PositiveFluxSectorCount','BValid','BSampleCount','AttitudeFound', ...
            'FinitePASectorCount','PADUsable','PanelEP1SameBin','PanelEAnySameBin', ...
            'Cause','SourceCDF','SourceCDFRecord'});
        for is = 1:7, R.(sprintf('FluxS%d', is)) = flux(:, is); end
        R.BR = rawB(:,1); R.BT = rawB(:,2); R.BN = rawB(:,3);
        report.Records = [report.Records; R]; %#ok<AGROW>

        step = days(1); if ic == 2, step = hours(1); end
        grid = (startTime:step:endTime-step).';
        nBin = numel(grid); eBin = false(nBin,1); eAnyBin = eBin;
        fCount = zeros(nBin,1); fUsable = fCount; fFluxGood = fCount;
        for ib = 1:nBin
            eBin(ib) = any(e1(cohoBin == grid(ib)));
            eAnyBin(ib) = any(eAny(cohoBin == grid(ib)));
            inBin = sectorBin == grid(ib);
            fCount(ib) = nnz(inBin);
            fUsable(ib) = nnz(inBin & usable);
            fFluxGood(ib) = nnz(inBin & fluxCount == 7);
        end
        Bins = table(repmat(string(cadence),nBin,1), repmat(eventID,nBin,1), ...
            grid, eBin, eAnyBin, fCount, fFluxGood, fUsable, ...
            'VariableNames', {'Cadence','EventID','BinStartUTC','PanelEP1', ...
            'PanelEAny','SectorRecords','CompleteFluxRecords','UsablePADRecords'});
        report.Bins = [report.Bins; Bins]; %#ok<AGROW>
        E = table(string(cadence), eventID, nnz(e1), nnz(eAny), n, ...
            nnz(fluxCount == 7), nnz(usable), nnz(fluxCount == 0), ...
            nnz(fluxCount > 0 & fluxCount < 7), nnz(fluxCount == 7 & ~bGood), ...
            nnz(~attitudeFound), nnz(eBin & fUsable == 0), nnz(eBin & fCount == 0), ...
            'VariableNames', {'Cadence','EventID','PanelEP1Records','PanelEAnyRecords', ...
            'SectorRecords','CompleteFluxRecords','UsablePADRecords','NoFluxRecords', ...
            'PartialFluxRecords','CompleteFluxButNoB','AttitudeMissingRecords', ...
            'EP1BinsWithoutPAD','EP1BinsWithoutSectorRecord'});
        report.Events = [report.Events; E]; %#ok<AGROW>
        if ic == 2
            peakFile = dir(fullfile(cfg.PeakPADDataFolder, ['V1_',char(eventID),'*_P1_peak5_hourly.mat']));
            A = load(fullfile(peakFile.folder, peakFile.name), 'peakAudit');
            if isfinite(A.peakAudit.PeakTableRow)
                report.Peaks = [report.Peaks; R(A.peakAudit.PeakTableRow,:)]; %#ok<AGROW>
            end
        end
        if ismember(eventID, ["Case1-S02-L01","Case1-S03-L04","Case1-S06-L01"])
            field = matlab.lang.makeValidName([char(eventID),'_',cadence]);
            report.Examples.(field) = struct('COHO',C,'SectorRecords',R,'Counts',E);
        end
    end
    fprintf('Panel e/f source diagnosis completed: %s, 45 events.\n', cadence);
end

%% unique source-record and UTC-bin counts avoid overlapping event windows
keys = report.Records.Cadence + "|" + string(report.Records.SourceRow);
[~, first] = unique(keys, 'stable');
report.UniqueRecords = report.Records(first,:);
keys = report.Bins.Cadence + "|" + string(posixtime(report.Bins.BinStartUTC));
[~, first] = unique(keys, 'stable');
report.UniqueBins = report.Bins(first,:);
report.CohoSourceFiles = string(cache.keys).';
report.CohoSHA256 = strings(size(report.CohoSourceFiles));
for ii = 1:numel(report.CohoSourceFiles)
    report.CohoSHA256(ii) = string(Case1_File_SHA256(char(report.CohoSourceFiles(ii))));
end
report.Summary = groupcounts(report.UniqueRecords, {'Cadence','Cause'});
report.PeakSummary = groupcounts(report.Peaks, 'Cause');
report.Passed = true;
folder = fullfile(cfg.DataRoot,'voyager1','lecp','validation','panel_ef_gaps');
if ~isfolder(folder), mkdir(folder); end
report.AuditFile = fullfile(folder, ['panel_ef_gaps_', ...
    char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS')),'.mat']);
save(report.AuditFile, 'report', '-v7.3');
disp(report.Summary); disp(report.PeakSummary);
fprintf('Read-only panel e/f audit: %s\n', report.AuditFile);
end
