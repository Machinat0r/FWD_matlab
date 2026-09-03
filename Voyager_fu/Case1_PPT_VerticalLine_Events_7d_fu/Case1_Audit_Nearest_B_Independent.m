function report = Case1_Audit_Nearest_B_Independent()
%Case1_Audit_Nearest_B_Independent Read-only nearest original COHO B timing.
%   No production PAD, flux or magnetic-field records are changed.

%% existing source audit and native CDF reader
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
auditFile = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', ...
    'panel_ef_gaps', 'panel_ef_gaps_20260903_071726_792.mat');
old = load(auditFile, 'report');
R = old.report.UniqueRecords;
Q = R(~R.BValid, :);
files = old.report.CohoSourceFiles;
C = readFiles(files);
N = matchNearest(Q, C);

%% include all on-disk months intersecting each nearest-distance search ball
% Any unexamined record that could beat the initial nearest candidate must
% lie inside one of these intervals. Reading these full months removes the
% event-window and first-source-month boundary bias without downloading.
entries = dir(fullfile(cfg.DataRoot, 'voyager1', 'coho', '1hr', 'l2', ...
    'merged_mag_plasma', '**', '*.cdf'));
extraFiles = strings(0, 1);
availableMonths = NaT(0, 1, 'TimeZone', 'UTC');
neededMonths = NaT(0, 1, 'TimeZone', 'UTC');
for iq = 1:height(N)
    lower = dateshift(N.EpochUTC(iq)-hours(N.AbsDeltaHours(iq)), 'start', 'month');
    upper = dateshift(N.EpochUTC(iq)+hours(N.AbsDeltaHours(iq)), 'start', 'month');
    neededMonths = [neededMonths; (lower:calmonths(1):upper).']; %#ok<AGROW>
end
neededMonths = unique(neededMonths);
for ii = 1:numel(entries)
    token = regexp(entries(ii).name, '_(\d{8})_v\d+\.cdf$', 'tokens', 'once');
    if isempty(token), continue; end
    month = datetime(token{1}, 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC');
    availableMonths(end+1, 1) = month; %#ok<AGROW>
    if any(neededMonths == month)
        extraFiles(end+1, 1) = string(fullfile(entries(ii).folder, entries(ii).name)); %#ok<AGROW>
    end
end
extraFiles = setdiff(extraFiles, files);
if ~isempty(extraFiles), C = [C; readFiles(extraFiles)]; end
C = sortrows(C, 'EpochUTC');
[~, first] = unique(C.EpochUTC, 'stable');
C = C(first, :);
N2 = matchNearest(Q, C);
assert(all(N2.AbsDeltaHours <= N.AbsDeltaHours));
P = old.report.Peaks;
P = matchNearest(P(~P.BValid, :), C);

%% grouped descriptive statistics, without applying a tolerance
summary = table;
for cadence = ["day", "hour"]
    for complete = [false true]
        use = N2.Cadence == cadence;
        if complete, use = use & N2.PositiveFluxSectorCount == 7; end
        delta = N2.AbsDeltaHours(use);
        q = prctile(delta, [50 90 95]);
        row = table(cadence, complete, numel(delta), min(delta), q(1), q(2), q(3), max(delta), ...
            'VariableNames', {'Cadence','RequireSevenPositiveFlux','N','MinHours', ...
            'MedianHours','P90Hours','P95Hours','MaxHours'});
        for limit = [1 3 6 12 24 48 72]
            row.(sprintf('Over%dh', limit)) = nnz(delta > limit);
        end
        summary = [summary; row]; %#ok<AGROW>
    end
end
report = struct;
report.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
report.PriorGapAudit = auditFile;
report.Policy = 'Nearest complete finite nonzero native COHO hourly BR/BT/BN vector; absolute Epoch difference; earlier tie; no tolerance applied in diagnostic.';
report.MissingBRecords = N2;
report.PeaksMissingB = P;
report.Summary = summary;
report.InitialSourceFiles = files;
report.AdditionalLocalSourceFiles = extraFiles;
report.InitialVsExpandedChanged = nnz(N.BEpochUTC ~= N2.BEpochUTC);
report.NeededMonths = neededMonths;
report.MissingLocalMonthsInsideInitialSearchIntervals = setdiff(neededMonths, availableMonths);
report.CohoData = C;
report.AllSourceFiles = [files; extraFiles];
report.SourceSHA256 = strings(numel(report.AllSourceFiles),1);
for ii = 1:numel(report.AllSourceFiles)
    report.SourceSHA256(ii) = string(Case1_File_SHA256(char(report.AllSourceFiles(ii))));
end
report.Passed = isempty(report.MissingLocalMonthsInsideInitialSearchIntervals);
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'validation', 'nearest_b');
if ~isfolder(folder), mkdir(folder); end
report.AuditFile = fullfile(folder, ['nearest_b_independent_', ...
    char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyyMMdd_HHmmss_SSS')), '.mat']);
save(report.AuditFile, 'report', '-v7.3');
disp(summary);
disp(P(:, {'EventID','EpochUTC','BEpochUTC','SignedDeltaHours','AbsDeltaHours','PositiveFluxSectorCount'}));
longest = sortrows(N2(:, {'Cadence','EventID','EpochUTC','BEpochUTC','AbsDeltaHours'}), 'AbsDeltaHours', 'descend');
disp(longest(1:min(10,height(N2)), :));
fprintf('Initial files=%d, extra files=%d, changed matches=%d, missing months=%d.\n', ...
    numel(files), numel(extraFiles), report.InitialVsExpandedChanged, ...
    numel(report.MissingLocalMonthsInsideInitialSearchIntervals));
fprintf('Independent nearest-B audit: %s\n', report.AuditFile);
end

function C = readFiles(files)
C = table;
for ii = 1:numel(files)
    p = Voyager_Read_CDF_Product(char(files(ii)), 'coho');
    part = table(p.Epoch(:), p.BR(:), p.BT(:), p.BN(:), ...
        repmat(string(files(ii)), numel(p.Epoch), 1), (1:numel(p.Epoch)).', ...
        'VariableNames', {'EpochUTC','BR','BT','BN','SourceCDF','SourceRecord'});
    C = [C; part]; %#ok<AGROW>
end
C = sortrows(C, 'EpochUTC');
end

function Q = matchNearest(Q, C)
good = all(isfinite(C{:, {'BR','BT','BN'}}), 2) & ...
    vecnorm(C{:, {'BR','BT','BN'}}, 2, 2) > 0;
G = C(good, :);
assert(~isempty(G), 'No complete original COHO vectors.');
qt = posixtime(Q.EpochUTC); bt = posixtime(G.EpochUTC);
idx = zeros(height(Q), 1);
for ii = 1:height(Q)
    [~, idx(ii)] = min(abs(bt-qt(ii)));
end
Q.BEpochUTC = G.EpochUTC(idx);
Q.SignedDeltaHours = hours(Q.BEpochUTC-Q.EpochUTC);
Q.AbsDeltaHours = abs(Q.SignedDeltaHours);
Q.NearestBR = G.BR(idx); Q.NearestBT = G.BT(idx); Q.NearestBN = G.BN(idx);
Q.BSourceCDF = G.SourceCDF(idx);
Q.BSourceRecord = G.SourceRecord(idx);
Q.EpochUTC.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
Q.BEpochUTC.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
end
