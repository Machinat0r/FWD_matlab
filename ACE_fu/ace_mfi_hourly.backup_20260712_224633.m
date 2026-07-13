function [Hourly, RunInfo, DailySelection, FileLog] = ace_mfi_hourly(dataRoot, outputDir, varargin)
%ACE_MFI_HOURLY Compute native-resolution hourly ACE/MFI mean and maximum.
%
%   [Hourly, RunInfo, DailySelection, FileLog] = ace_mfi_hourly(...)
%
% Selection and calculation rules:
%   1. Select exactly one product for each UTC day.
%   2. Try all H3 versions from newest to oldest. If a readable H3 file
%      contains at least one valid sample that day, use H3 for all 24 hours.
%   3. Only when no H3 candidate is usable, try H0 in the same way.
%   4. Never combine H3 and H0, never resample, interpolate, or fill gaps.
%   5. Calculate mean(Magnitude) and max(Magnitude) directly from the native
%      samples in each half-open UTC hour [hh:00, hh+1:00).
%   6. Valid samples satisfy Q_FLAG == 0 and 0 <= Magnitude <= 500 nT.
%
% Name-value options:
%   'StartDate'              UTC date; NaT means first database date
%   'EndDate'                UTC date; NaT means last database date
%   'ForceRebuild'           false (default) or true
%   'WriteCSV'               true (default) or false
%   'MakeFigure'             true (default) or false
%   'FigureVisible'          'off' (default) or 'on'
%   'CheckpointEveryDays'    30 (default)
%   'Verbose'                true (default) or false
%
% Example:
%   [H, R] = ace_mfi_hourly( ...
%       "Z:\SPART-WORK\Data\ACE", ...
%       "C:\ACE_MFI_hourly_output");

    if nargin < 1 || isempty(dataRoot)
        dataRoot = "Z:\SPART-WORK\Data\ACE";
    end
    if nargin < 2 || isempty(outputDir)
        outputDir = fullfile(pwd, "ACE_MFI_hourly_output");
    end

    dataRoot = string(dataRoot);
    outputDir = string(outputDir);
    if ~isscalar(dataRoot) || ~isscalar(outputDir)
        error('ACE:InvalidPath', 'dataRoot and outputDir must be scalar paths.');
    end

    opts = parseOptions(varargin{:});
    mfiRoot = locateMfiRoot(dataRoot);

    ensureDirectory(outputDir);
    cacheDir = fullfile(outputDir, "cache_by_year");
    ensureDirectory(cacheDir);

    if opts.Verbose
        fprintf('ACE/MFI source: %s\n', mfiRoot);
        fprintf('Output folder: %s\n', outputDir);
        fprintf('Building file manifest (read-only scan)...\n');
    end

    Manifest = buildManifest(mfiRoot);
    if isempty(Manifest)
        error('ACE:NoMfiFiles', 'No ACE MFI H3/H0 CDF files were found under %s.', mfiRoot);
    end

    databaseFirstDate = min(Manifest.DateUTC);
    databaseLastDate = max(Manifest.DateUTC);
    startDate = databaseFirstDate;
    endDate = databaseLastDate;
    if ~isnat(opts.StartDate)
        startDate = max(databaseFirstDate, opts.StartDate);
    end
    if ~isnat(opts.EndDate)
        endDate = min(databaseLastDate, opts.EndDate);
    end
    if startDate > endDate
        error('ACE:InvalidDateRange', ...
            'Requested range does not overlap database range %s to %s.', ...
            char(string(databaseFirstDate, 'yyyy-MM-dd')), ...
            char(string(databaseLastDate, 'yyyy-MM-dd')));
    end

    startDate = dateshift(startDate, 'start', 'day');
    endDate = dateshift(endDate, 'start', 'day');
    allDates = (startDate:caldays(1):endDate).';
    expectedRows = numel(allDates) * 24;

    if opts.Verbose
        fprintf('UTC range: %s through %s (%d days, %d hourly rows)\n', ...
            char(string(startDate, 'yyyy-MM-dd')), ...
            char(string(endDate, 'yyyy-MM-dd')), ...
            numel(allDates), expectedRows);
        fprintf('Native samples only: no H3/H0 mixing and no resampling.\n');
    end

    config = processingConfig(opts);
    yearsToRun = year(startDate):year(endDate);
    if opts.ForceRebuild
        clearRequestedYearCaches(cacheDir, yearsToRun);
    end
    yearHourly = cell(numel(yearsToRun), 1);
    yearDaily = cell(numel(yearsToRun), 1);
    yearLogs = cell(numel(yearsToRun), 1);

    for k = 1:numel(yearsToRun)
        thisYear = yearsToRun(k);
        firstOfYear = datetime(thisYear, 1, 1, 'TimeZone', 'UTC');
        lastOfYear = datetime(thisYear, 12, 31, 'TimeZone', 'UTC');
        yearStart = max(startDate, firstOfYear);
        yearEnd = min(endDate, lastOfYear);
        datesThisYear = (yearStart:caldays(1):yearEnd).';
        manifestThisYear = Manifest(year(Manifest.DateUTC) == thisYear, :);

        [yearHourly{k}, yearDaily{k}, yearLogs{k}] = processYear( ...
            thisYear, datesThisYear, manifestThisYear, cacheDir, config, opts);
    end

    Hourly = vertcat(yearHourly{:});
    DailySelection = vertcat(yearDaily{:});
    FileLog = vertcat(yearLogs{:});

    Hourly.TimeUTC.Format = 'yyyy-MM-dd''T''HH:mm:ss''Z''';
    DailySelection.DateUTC.Format = 'yyyy-MM-dd';
    Manifest.DateUTC.Format = 'yyyy-MM-dd';
    FileLog.DateUTC.Format = 'yyyy-MM-dd';

    validateFinalResult(Hourly, DailySelection, startDate, endDate, expectedRows);

    startTag = char(string(startDate, 'yyyyMMdd'));
    endTag = char(string(endDate, 'yyyyMMdd'));
    baseName = sprintf('ACE_MFI_hourly_%s_%s', startTag, endTag);
    paths = outputPaths(outputDir, baseName);

    outputStatus = struct( ...
        'HourlyCSVWritten', false, ...
        'DailyCSVWritten', false, ...
        'ManifestCSVWritten', false, ...
        'LogCSVWritten', false, ...
        'FigureWritten', false);

    if opts.WriteCSV
        try
            writetimetable(Hourly, paths.HourlyCSV);
            outputStatus.HourlyCSVWritten = true;
        catch ME
            FileLog = appendLog(FileLog, endDate, "OUTPUT", 0, paths.HourlyCSV, ...
                "WRITE_FAILED", ME.message);
            warning('ACE:CSVWriteFailed', 'Could not write hourly CSV: %s', ME.message);
        end

        try
            writetable(DailySelection, paths.DailyCSV);
            outputStatus.DailyCSVWritten = true;
        catch ME
            FileLog = appendLog(FileLog, endDate, "OUTPUT", 0, paths.DailyCSV, ...
                "WRITE_FAILED", ME.message);
            warning('ACE:CSVWriteFailed', 'Could not write daily selection CSV: %s', ME.message);
        end
    end

    try
        writetable(Manifest, paths.ManifestCSV);
        outputStatus.ManifestCSVWritten = true;
    catch ME
        FileLog = appendLog(FileLog, endDate, "OUTPUT", 0, paths.ManifestCSV, ...
            "WRITE_FAILED", ME.message);
        warning('ACE:CSVWriteFailed', 'Could not write manifest CSV: %s', ME.message);
    end

    if opts.MakeFigure
        try
            makeOverviewFigure(Hourly, paths.FigureFIG, paths.FigurePNG, opts.FigureVisible);
            outputStatus.FigureWritten = true;
        catch ME
            FileLog = appendLog(FileLog, endDate, "OUTPUT", 0, paths.FigurePNG, ...
                "PLOT_FAILED", ME.message);
            warning('ACE:PlotFailed', 'Could not create overview figure: %s', ME.message);
        end
    end

    try
        writetable(FileLog, paths.FileLogCSV);
        outputStatus.LogCSVWritten = true;
    catch ME
        warning('ACE:CSVWriteFailed', 'Could not write file log CSV: %s', ME.message);
    end

    RunInfo = struct();
    RunInfo.CreatedUTC = datetime('now', 'TimeZone', 'UTC');
    RunInfo.DataRoot = dataRoot;
    RunInfo.MfiRoot = mfiRoot;
    RunInfo.OutputDirectory = outputDir;
    RunInfo.DatabaseFirstDateUTC = databaseFirstDate;
    RunInfo.DatabaseLastDateUTC = databaseLastDate;
    RunInfo.ProcessedFirstDateUTC = startDate;
    RunInfo.ProcessedLastDateUTC = endDate;
    RunInfo.TotalDays = numel(allDates);
    RunInfo.TotalHourlyRows = height(Hourly);
    RunInfo.H3SelectedDays = nnz(string(DailySelection.Source) == "H3");
    RunInfo.H0SelectedDays = nnz(string(DailySelection.Source) == "H0");
    RunInfo.MissingDays = nnz(string(DailySelection.Source) == "MISSING");
    RunInfo.HoursWithValidData = nnz(Hourly.N_valid > 0);
    RunInfo.ManifestFileCount = height(Manifest);
    RunInfo.FileLogRows = height(FileLog);
    RunInfo.SelectionRule = "One native product per UTC day: usable H3 first, otherwise usable H0";
    RunInfo.Statistics = "Hourly arithmetic mean and maximum of native Magnitude samples";
    RunInfo.ValidSampleRule = "Q_FLAG==0, finite Magnitude, 0<=Magnitude<=500 nT";
    RunInfo.Options = opts;
    RunInfo.OutputStatus = outputStatus;
    RunInfo.OutputFiles = paths;

    atomicSaveFinal(paths.FinalMAT, Hourly, DailySelection, Manifest, FileLog, RunInfo);

    if opts.Verbose
        fprintf('\nCompleted.\n');
        fprintf('H3 days: %d; H0 days: %d; missing/unusable days: %d\n', ...
            RunInfo.H3SelectedDays, RunInfo.H0SelectedDays, RunInfo.MissingDays);
        fprintf('Hourly MAT: %s\n', paths.FinalMAT);
        if outputStatus.HourlyCSVWritten
            fprintf('Hourly CSV: %s\n', paths.HourlyCSV);
        end
        if outputStatus.FigureWritten
            fprintf('Figure PNG: %s\n', paths.FigurePNG);
        end
    end
end


function clearRequestedYearCaches(cacheDir, yearsToRun)
    for thisYear = yearsToRun
        files = dir(fullfile(cacheDir, sprintf('ACE_MFI_hourly_%04d*.mat', thisYear)));
        for i = 1:numel(files)
            if ~files(i).isdir
                delete(fullfile(files(i).folder, files(i).name));
            end
        end
    end
end


function opts = parseOptions(varargin)
    p = inputParser;
    p.FunctionName = 'ace_mfi_hourly';
    addParameter(p, 'StartDate', NaT, @(x) isDateLike(x));
    addParameter(p, 'EndDate', NaT, @(x) isDateLike(x));
    addParameter(p, 'ForceRebuild', false, @(x) isLogicalScalar(x));
    addParameter(p, 'WriteCSV', true, @(x) isLogicalScalar(x));
    addParameter(p, 'MakeFigure', true, @(x) isLogicalScalar(x));
    addParameter(p, 'FigureVisible', 'off', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(p, 'CheckpointEveryDays', 30, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && fix(x) == x);
    addParameter(p, 'Verbose', true, @(x) isLogicalScalar(x));
    parse(p, varargin{:});

    opts = p.Results;
    opts.StartDate = normalizeOptionalDate(opts.StartDate);
    opts.EndDate = normalizeOptionalDate(opts.EndDate);
    opts.ForceRebuild = logical(opts.ForceRebuild);
    opts.WriteCSV = logical(opts.WriteCSV);
    opts.MakeFigure = logical(opts.MakeFigure);
    opts.Verbose = logical(opts.Verbose);
    opts.CheckpointEveryDays = double(opts.CheckpointEveryDays);
    opts.FigureVisible = lower(string(opts.FigureVisible));
    if ~ismember(opts.FigureVisible, ["on", "off"])
        error('ACE:InvalidFigureVisible', 'FigureVisible must be ''on'' or ''off''.');
    end
end


function tf = isDateLike(x)
    tf = isempty(x) || (isdatetime(x) && isscalar(x)) || ischar(x) || ...
        (isstring(x) && isscalar(x));
end


function tf = isLogicalScalar(x)
    tf = (islogical(x) || isnumeric(x)) && isscalar(x) && isfinite(double(x)) && ...
        ismember(double(x), [0, 1]);
end


function value = normalizeOptionalDate(value)
    if isempty(value) || (isdatetime(value) && isscalar(value) && isnat(value)) || ...
            ((ischar(value) || isstring(value)) && strlength(string(value)) == 0)
        value = NaT(1, 1, 'TimeZone', 'UTC');
        return;
    end

    if ~isdatetime(value)
        value = datetime(value, 'TimeZone', 'UTC');
    elseif isempty(value.TimeZone)
        value.TimeZone = 'UTC';
    else
        value.TimeZone = 'UTC';
    end
    value = dateshift(value, 'start', 'day');
end


function mfiRoot = locateMfiRoot(dataRoot)
    candidate1 = fullfile(dataRoot, "ace", "mfi");
    candidate2 = dataRoot;
    if isfolder(fullfile(candidate1, "h3", "l2")) || ...
            isfolder(fullfile(candidate1, "h0", "l2"))
        mfiRoot = candidate1;
    elseif isfolder(fullfile(candidate2, "h3", "l2")) || ...
            isfolder(fullfile(candidate2, "h0", "l2"))
        mfiRoot = candidate2;
    else
        error('ACE:MfiRootNotFound', ...
            ['Could not find ace\mfi\h3\l2 or ace\mfi\h0\l2 under %s. ' ...
             'If Z: is not visible to MATLAB, use the corresponding UNC path.'], dataRoot);
    end
    mfiRoot = string(mfiRoot);
end


function ensureDirectory(folder)
    if ~isfolder(folder)
        [ok, message] = mkdir(folder);
        if ~ok
            error('ACE:CreateDirectoryFailed', 'Could not create %s: %s', folder, message);
        end
    end
end


function Manifest = buildManifest(mfiRoot)
    patterns = [ ...
        fullfile(mfiRoot, "h3", "l2", "**", "ac_h3_mfi_*.cdf"), ...
        fullfile(mfiRoot, "h0", "l2", "**", "ac_h0_mfi_*.cdf")];

    listings = cell(numel(patterns), 1);
    totalFiles = 0;
    for i = 1:numel(patterns)
        listings{i} = dir(patterns(i));
        listings{i} = listings{i}(~[listings{i}.isdir]);
        totalFiles = totalFiles + numel(listings{i});
    end

    dateText = strings(totalFiles, 1);
    product = strings(totalFiles, 1);
    version = zeros(totalFiles, 1, 'uint16');
    path = strings(totalFiles, 1);
    bytes = zeros(totalFiles, 1);
    modifiedDatenum = zeros(totalFiles, 1);

    row = 0;
    for i = 1:numel(listings)
        items = listings{i};
        for j = 1:numel(items)
            token = regexp(items(j).name, ...
                '^ac_(h[03])_mfi_(\d{8})_v(\d+)\.cdf$', 'tokens', 'once');
            if isempty(token)
                continue;
            end
            row = row + 1;
            product(row) = upper(string(token{1}));
            dateText(row) = string(token{2});
            version(row) = uint16(str2double(token{3}));
            path(row) = string(fullfile(items(j).folder, items(j).name));
            bytes(row) = double(items(j).bytes);
            modifiedDatenum(row) = double(items(j).datenum);
        end
    end

    dateText = dateText(1:row);
    product = product(1:row);
    version = version(1:row);
    path = path(1:row);
    bytes = bytes(1:row);
    modifiedDatenum = modifiedDatenum(1:row);
    dateUTC = datetime(dateText, 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC');

    Manifest = table(dateUTC, product, version, path, bytes, modifiedDatenum, ...
        'VariableNames', {'DateUTC', 'Product', 'Version', 'Path', 'Bytes', 'ModifiedDatenum'});
    Manifest = sortrows(Manifest, {'DateUTC', 'Product', 'Version', 'Path'}, ...
        {'ascend', 'ascend', 'descend', 'ascend'});
end


function config = processingConfig(opts)
    config.CacheFormatVersion = uint32(2);
    config.Signature = "ACE_MFI_NATIVE_HOURLY_V2_Q0_MAG_0_500_DAY_LEVEL_H3_THEN_H0";
    config.ValidMagnitudeMin = 0;
    config.ValidMagnitudeMax = 500;
    config.ValidQualityFlag = 0;
    config.CheckpointEveryDays = opts.CheckpointEveryDays;
end


function [YearHourly, YearDaily, YearLog] = processYear( ...
        thisYear, datesThisYear, manifestThisYear, cacheDir, config, opts)

    cacheFile = fullfile(cacheDir, sprintf('ACE_MFI_hourly_%04d.mat', thisYear));
    checkpointFile = fullfile(cacheDir, sprintf('ACE_MFI_hourly_%04d_checkpoint.mat', thisYear));
    expectedMeta = makeCacheMeta(datesThisYear, manifestThisYear, config);

    if opts.ForceRebuild
        if isfile(cacheFile)
            delete(cacheFile);
        end
        if isfile(checkpointFile)
            delete(checkpointFile);
        end
    end

    if ~opts.ForceRebuild && isfile(cacheFile)
        try
            loaded = load(cacheFile, 'YearHourly', 'YearDaily', 'YearLog', 'CacheMeta');
            if isfield(loaded, 'CacheMeta') && cacheMetaMatches(loaded.CacheMeta, expectedMeta)
                if hasRetryableFailures(loaded.YearLog)
                    if opts.Verbose
                        fprintf(['[%04d] Cache contains retryable I/O failures; ' ...
                            'reprocessing this year.\n'], thisYear);
                    end
                else
                    YearHourly = loaded.YearHourly;
                    YearDaily = loaded.YearDaily;
                    YearLog = loaded.YearLog;
                    if opts.Verbose
                        fprintf('[%04d] Loaded valid yearly cache.\n', thisYear);
                    end
                    return;
                end
            end
        catch ME
            warning('ACE:CacheReadFailed', 'Ignoring unreadable %04d cache: %s', thisYear, ME.message);
        end
    end

    nDays = numel(datesThisYear);
    nHours = nDays * 24;
    hourGrid = datesThisYear + hours(0:23);
    hourUTC = reshape(hourGrid.', [], 1);

    BMean = nan(nHours, 1);
    BMax = nan(nHours, 1);
    NTotal = zeros(nHours, 1, 'uint32');
    NValid = zeros(nHours, 1, 'uint32');
    Coverage = nan(nHours, 1, 'single');
    SourceCode = zeros(nHours, 1, 'uint8');
    SourceVersion = zeros(nHours, 1, 'uint16');

    DailySourceCode = zeros(nDays, 1, 'uint8');
    DailySourceVersion = zeros(nDays, 1, 'uint16');
    DailySourceFile = strings(nDays, 1);
    DailyValidSamples = zeros(nDays, 1, 'uint32');
    DailyStatus = strings(nDays, 1);
    YearLog = emptyLogTable();
    nextDayIndex = 1;

    if ~opts.ForceRebuild && isfile(checkpointFile)
        try
            loaded = load(checkpointFile, 'Checkpoint', 'CacheMeta');
            if cacheMetaMatches(loaded.CacheMeta, expectedMeta)
                C = loaded.Checkpoint;
                if checkpointIsValid(C, nHours, nDays) && ...
                        ~hasRetryableFailures(C.YearLog)
                    newBMean = C.BMean;
                    newBMax = C.BMax;
                    newNTotal = C.NTotal;
                    newNValid = C.NValid;
                    newCoverage = C.Coverage;
                    newSourceCode = C.SourceCode;
                    newSourceVersion = C.SourceVersion;
                    newDailySourceCode = C.DailySourceCode;
                    newDailySourceVersion = C.DailySourceVersion;
                    newDailySourceFile = C.DailySourceFile;
                    newDailyValidSamples = C.DailyValidSamples;
                    newDailyStatus = C.DailyStatus;
                    newYearLog = C.YearLog;
                    newNextDayIndex = C.NextDayIndex;

                    % Assign only after every checkpoint field has passed validation.
                    BMean = newBMean;
                    BMax = newBMax;
                    NTotal = newNTotal;
                    NValid = newNValid;
                    Coverage = newCoverage;
                    SourceCode = newSourceCode;
                    SourceVersion = newSourceVersion;
                    DailySourceCode = newDailySourceCode;
                    DailySourceVersion = newDailySourceVersion;
                    DailySourceFile = newDailySourceFile;
                    DailyValidSamples = newDailyValidSamples;
                    DailyStatus = newDailyStatus;
                    YearLog = newYearLog;
                    nextDayIndex = newNextDayIndex;
                    if opts.Verbose
                        fprintf('[%04d] Resuming checkpoint at day %d of %d.\n', ...
                            thisYear, nextDayIndex, nDays);
                    end
                elseif checkpointIsValid(C, nHours, nDays)
                    warning('ACE:RetryableCheckpoint', ...
                        ['Ignoring %04d checkpoint because it contains a retryable ' ...
                         'file I/O failure. The year will be recalculated.'], thisYear);
                else
                    warning('ACE:InvalidCheckpoint', ...
                        'Ignoring structurally invalid %04d checkpoint.', thisYear);
                end
            end
        catch ME
            warning('ACE:CheckpointReadFailed', ...
                'Ignoring unreadable %04d checkpoint: %s', thisYear, ME.message);
        end
    end

    for dayIndex = nextDayIndex:nDays
        dateUTC = datesThisYear(dayIndex);
        candidatesH3 = manifestThisYear( ...
            manifestThisYear.DateUTC == dateUTC & manifestThisYear.Product == "H3", :);
        candidatesH0 = manifestThisYear( ...
            manifestThisYear.DateUTC == dateUTC & manifestThisYear.Product == "H0", :);

        [dayResult, selectedRow, YearLog] = tryCandidates( ...
            dateUTC, "H3", candidatesH3, config, YearLog);
        h3HadCandidates = ~isempty(candidatesH3);

        selectedProductCode = uint8(1);
        status = "H3_SELECTED";
        if ~dayResult.Usable
            [dayResult, selectedRow, YearLog] = tryCandidates( ...
                dateUTC, "H0", candidatesH0, config, YearLog);
            selectedProductCode = uint8(2);
            if h3HadCandidates
                status = "H0_SELECTED_H3_UNUSABLE";
            else
                status = "H0_SELECTED_H3_MISSING";
            end
        end

        hourRows = (dayIndex - 1) * 24 + (1:24);
        if dayResult.Usable
            BMean(hourRows) = dayResult.BMean;
            BMax(hourRows) = dayResult.BMax;
            NTotal(hourRows) = dayResult.NTotal;
            NValid(hourRows) = dayResult.NValid;
            Coverage(hourRows) = dayResult.Coverage;
            SourceCode(hourRows) = selectedProductCode;
            SourceVersion(hourRows) = selectedRow.Version;

            DailySourceCode(dayIndex) = selectedProductCode;
            DailySourceVersion(dayIndex) = selectedRow.Version;
            DailySourceFile(dayIndex) = selectedRow.Path;
            DailyValidSamples(dayIndex) = uint32(sum(double(dayResult.NValid)));
            DailyStatus(dayIndex) = status;
        else
            if isempty(candidatesH3) && isempty(candidatesH0)
                DailyStatus(dayIndex) = "NO_FILE";
            else
                DailyStatus(dayIndex) = "NO_USABLE_FILE";
            end
            YearLog = appendLog(YearLog, dateUTC, "NONE", 0, "", ...
                "MISSING_DAY", "No usable H3 or H0 file for this UTC day.");
        end

        if opts.Verbose && (dayIndex == 1 || dayIndex == nDays || mod(dayIndex, 30) == 0)
            fprintf('[%04d] Processed day %d/%d (%s), source=%s\n', ...
                thisYear, dayIndex, nDays, char(string(dateUTC, 'yyyy-MM-dd')), ...
                char(sourceName(DailySourceCode(dayIndex))));
        end

        if dayIndex < nDays && mod(dayIndex, config.CheckpointEveryDays) == 0
            Checkpoint = struct();
            Checkpoint.NextDayIndex = dayIndex + 1;
            Checkpoint.BMean = BMean;
            Checkpoint.BMax = BMax;
            Checkpoint.NTotal = NTotal;
            Checkpoint.NValid = NValid;
            Checkpoint.Coverage = Coverage;
            Checkpoint.SourceCode = SourceCode;
            Checkpoint.SourceVersion = SourceVersion;
            Checkpoint.DailySourceCode = DailySourceCode;
            Checkpoint.DailySourceVersion = DailySourceVersion;
            Checkpoint.DailySourceFile = DailySourceFile;
            Checkpoint.DailyValidSamples = DailyValidSamples;
            Checkpoint.DailyStatus = DailyStatus;
            Checkpoint.YearLog = YearLog;
            atomicSaveCheckpoint(checkpointFile, Checkpoint, expectedMeta);
        end
    end

    Source = categorical(double(SourceCode), [0, 1, 2], {'MISSING', 'H3', 'H0'});
    YearHourly = timetable(hourUTC, BMean, BMax, Source, SourceVersion, ...
        NTotal, NValid, Coverage, ...
        'VariableNames', {'B_mean_nT', 'B_max_nT', 'Source', 'SourceVersion', ...
                          'N_total', 'N_valid', 'Coverage'});
    YearHourly.Properties.DimensionNames{1} = 'TimeUTC';

    DailySource = categorical(double(DailySourceCode), [0, 1, 2], {'MISSING', 'H3', 'H0'});
    YearDaily = table(datesThisYear, DailySource, DailySourceVersion, DailySourceFile, ...
        DailyValidSamples, DailyStatus, ...
        'VariableNames', {'DateUTC', 'Source', 'SourceVersion', 'SourceFile', ...
                          'N_valid_day', 'Status'});

    CacheMeta = expectedMeta;
    atomicSaveYear(cacheFile, YearHourly, YearDaily, YearLog, CacheMeta);
    if isfile(checkpointFile)
        delete(checkpointFile);
    end
end


function [result, selectedRow, logTable] = tryCandidates( ...
        dateUTC, product, candidates, config, logTable)

    result = emptyDayResult();
    selectedRow = emptySelectedRow();
    if isempty(candidates)
        return;
    end

    candidates = sortrows(candidates, {'Version', 'ModifiedDatenum'}, ...
        {'descend', 'descend'});
    for i = 1:height(candidates)
        candidate = candidates(i, :);
        [candidateResult, status, message] = readMfiDay( ...
            candidate.Path, dateUTC, product, config);
        if candidateResult.Usable
            result = candidateResult;
            selectedRow = struct( ...
                'Version', candidate.Version, ...
                'Path', candidate.Path);
            return;
        end
        logTable = appendLog(logTable, dateUTC, product, candidate.Version, ...
            candidate.Path, status, message);
    end
end


function [result, status, message] = readMfiDay(filePath, dateUTC, product, config)
    result = emptyDayResult();
    message = "";

    [headerOK, headerStatus, headerMessage] = quickHeaderCheck(filePath);
    if ~headerOK
        status = headerStatus;
        message = headerMessage;
        return;
    end

    try
        D = cdfread(filePath, ...
            'Variables', {'Epoch', 'Magnitude', 'Q_FLAG'}, ...
            'CombineRecords', true, ...
            'DatetimeType', 'datetime');

        if ~iscell(D) || numel(D) ~= 3
            error('ACE:UnexpectedCDFResult', ...
                'cdfread did not return the expected three variables.');
        end

        t = D{1}(:);
        b = double(D{2}(:));
        q = double(D{3}(:));
        if ~isdatetime(t)
            error('ACE:UnexpectedEpochType', 'Epoch was not returned as datetime.');
        end
        t.TimeZone = 'UTC';
        if numel(t) ~= numel(b) || numel(t) ~= numel(q)
            error('ACE:VariableLengthMismatch', ...
                'Epoch, Magnitude, and Q_FLAG lengths are inconsistent.');
        end
        if isempty(t)
            status = "NO_RECORDS";
            message = "CDF variables contain no records.";
            return;
        end

        dayEnd = dateUTC + caldays(1);
        validTime = ~isnat(t) & t >= dateUTC & t < dayEnd;
        hourIndexAll = hour(t(validTime)) + 1;
        result.NTotal = uint32(accumarray(hourIndexAll, ones(size(hourIndexAll)), ...
            [24, 1], @sum, 0));

        valid = validTime & isfinite(b) & ...
            b >= config.ValidMagnitudeMin & b <= config.ValidMagnitudeMax & ...
            q == config.ValidQualityFlag;
        if ~any(valid)
            status = "NO_VALID_SAMPLES";
            message = "File is readable but contains no valid samples for its nominal UTC day.";
            return;
        end

        validHours = hour(t(valid)) + 1;
        validB = b(valid);
        result.NValid = uint32(accumarray(validHours, ones(size(validHours)), ...
            [24, 1], @sum, 0));
        result.BMean = accumarray(validHours, validB, [24, 1], @mean, NaN);
        result.BMax = accumarray(validHours, validB, [24, 1], @max, NaN);

        if strcmpi(product, "H3")
            expectedPerHour = 3600;
        else
            expectedPerHour = 225;
        end
        result.Coverage = single(min(1, double(result.NValid) / expectedPerHour));
        result.Usable = true;
        status = "OK";
    catch ME
        status = "READ_FAILED";
        message = string(ME.identifier) + ": " + string(ME.message);
    end
end


function [ok, status, message] = quickHeaderCheck(filePath)
    ok = false;
    status = "INVALID_FILE_HEADER";
    message = "";
    try
        [fid, openMessage] = fopen(filePath, 'r');
    catch ME
        status = "OPEN_FAILED";
        message = string(ME.identifier) + ": " + string(ME.message);
        return;
    end
    if fid < 0
        status = "OPEN_FAILED";
        message = "Cannot open file: " + string(openMessage);
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    try
        header = fread(fid, 16, '*uint8');
    catch ME
        status = "READ_FAILED";
        message = string(ME.identifier) + ": " + string(ME.message);
        return;
    end
    if numel(header) < 8
        message = "File is shorter than the minimum CDF header check.";
        return;
    end
    if all(header == 0)
        message = "The first 16 bytes are all zero; file is an invalid zero-filled placeholder.";
        return;
    end
    ok = true;
    status = "OK";
end


function result = emptyDayResult()
    result = struct( ...
        'Usable', false, ...
        'BMean', nan(24, 1), ...
        'BMax', nan(24, 1), ...
        'NTotal', zeros(24, 1, 'uint32'), ...
        'NValid', zeros(24, 1, 'uint32'), ...
        'Coverage', nan(24, 1, 'single'));
end


function selected = emptySelectedRow()
    selected = struct('Version', uint16(0), 'Path', "");
end


function name = sourceName(code)
    switch double(code)
        case 1
            name = "H3";
        case 2
            name = "H0";
        otherwise
            name = "MISSING";
    end
end


function T = emptyLogTable()
    dateUTC = NaT(0, 1, 'TimeZone', 'UTC');
    product = strings(0, 1);
    version = zeros(0, 1, 'uint16');
    filePath = strings(0, 1);
    status = strings(0, 1);
    message = strings(0, 1);
    T = table(dateUTC, product, version, filePath, status, message, ...
        'VariableNames', {'DateUTC', 'Product', 'Version', 'FilePath', 'Status', 'Message'});
end


function T = appendLog(T, dateUTC, product, version, filePath, status, message)
    row = table(dateUTC, string(product), uint16(version), string(filePath), ...
        string(status), string(message), ...
        'VariableNames', T.Properties.VariableNames);
    T = [T; row];
end


function meta = makeCacheMeta(datesThisYear, manifestThisYear, config)
    meta.CacheFormatVersion = config.CacheFormatVersion;
    meta.ConfigSignature = config.Signature;
    meta.StartDateUTC = datesThisYear(1);
    meta.EndDateUTC = datesThisYear(end);
    meta.ManifestSnapshot = manifestThisYear(:, ...
        {'DateUTC', 'Product', 'Version', 'Path', 'Bytes', 'ModifiedDatenum'});
end


function tf = cacheMetaMatches(actual, expected)
    tf = false;
    required = {'CacheFormatVersion', 'ConfigSignature', 'StartDateUTC', ...
                'EndDateUTC', 'ManifestSnapshot'};
    if ~isstruct(actual) || ~all(isfield(actual, required))
        return;
    end
    try
        tf = actual.CacheFormatVersion == expected.CacheFormatVersion && ...
            string(actual.ConfigSignature) == string(expected.ConfigSignature) && ...
            actual.StartDateUTC == expected.StartDateUTC && ...
            actual.EndDateUTC == expected.EndDateUTC && ...
            isequaln(actual.ManifestSnapshot, expected.ManifestSnapshot);
    catch
        tf = false;
    end
end


function tf = hasRetryableFailures(logTable)
    tf = false;
    if ~istable(logTable) || ~ismember('Status', logTable.Properties.VariableNames)
        return;
    end
    retryableStatus = ["OPEN_FAILED", "READ_FAILED"];
    tf = any(ismember(string(logTable.Status), retryableStatus));
end


function tf = checkpointIsValid(C, nHours, nDays)
    tf = false;
    required = { ...
        'NextDayIndex', 'BMean', 'BMax', 'NTotal', 'NValid', 'Coverage', ...
        'SourceCode', 'SourceVersion', 'DailySourceCode', ...
        'DailySourceVersion', 'DailySourceFile', 'DailyValidSamples', ...
        'DailyStatus', 'YearLog'};
    if ~isstruct(C) || ~isscalar(C) || ~all(isfield(C, required))
        return;
    end

    hourlyNumeric = {'BMean', 'BMax', 'NTotal', 'NValid', 'Coverage', ...
                     'SourceCode', 'SourceVersion'};
    for i = 1:numel(hourlyNumeric)
        value = C.(hourlyNumeric{i});
        if ~isnumeric(value) || ~iscolumn(value) || numel(value) ~= nHours
            return;
        end
    end

    dailyNumeric = {'DailySourceCode', 'DailySourceVersion', 'DailyValidSamples'};
    for i = 1:numel(dailyNumeric)
        value = C.(dailyNumeric{i});
        if ~isnumeric(value) || ~iscolumn(value) || numel(value) ~= nDays
            return;
        end
    end

    if ~isstring(C.DailySourceFile) || ~iscolumn(C.DailySourceFile) || ...
            numel(C.DailySourceFile) ~= nDays
        return;
    end
    if ~isstring(C.DailyStatus) || ~iscolumn(C.DailyStatus) || ...
            numel(C.DailyStatus) ~= nDays
        return;
    end
    if ~istable(C.YearLog)
        return;
    end
    requiredLogVariables = {'DateUTC', 'Product', 'Version', 'FilePath', 'Status', 'Message'};
    if ~all(ismember(requiredLogVariables, C.YearLog.Properties.VariableNames))
        return;
    end

    nextIndex = double(C.NextDayIndex);
    if ~isscalar(nextIndex) || ~isfinite(nextIndex) || fix(nextIndex) ~= nextIndex || ...
            nextIndex < 1 || nextIndex > nDays + 1
        return;
    end
    tf = true;
end


function atomicSaveCheckpoint(targetFile, Checkpoint, CacheMeta)
    tempFile = string(targetFile) + ".tmp.mat";
    if isfile(tempFile)
        delete(tempFile);
    end
    save(tempFile, 'Checkpoint', 'CacheMeta', '-v7');
    replaceFile(tempFile, targetFile);
end


function atomicSaveYear(targetFile, YearHourly, YearDaily, YearLog, CacheMeta)
    tempFile = string(targetFile) + ".tmp.mat";
    if isfile(tempFile)
        delete(tempFile);
    end
    save(tempFile, 'YearHourly', 'YearDaily', 'YearLog', 'CacheMeta', '-v7');
    replaceFile(tempFile, targetFile);
end


function atomicSaveFinal(targetFile, Hourly, DailySelection, Manifest, FileLog, RunInfo)
    tempFile = string(targetFile) + ".tmp.mat";
    if isfile(tempFile)
        delete(tempFile);
    end
    save(tempFile, 'Hourly', 'DailySelection', 'Manifest', 'FileLog', 'RunInfo', '-v7');
    replaceFile(tempFile, targetFile);
end


function replaceFile(sourceFile, targetFile)
    [ok, message] = movefile(sourceFile, targetFile, 'f');
    if ~ok
        error('ACE:AtomicMoveFailed', 'Could not replace %s: %s', targetFile, message);
    end
end


function validateFinalResult(Hourly, DailySelection, startDate, endDate, expectedRows)
    if height(Hourly) ~= expectedRows
        error('ACE:RowCountMismatch', 'Expected %d hourly rows but got %d.', ...
            expectedRows, height(Hourly));
    end
    if height(DailySelection) * 24 ~= expectedRows
        error('ACE:DailyRowCountMismatch', 'Daily and hourly row counts are inconsistent.');
    end
    if Hourly.TimeUTC(1) ~= startDate || Hourly.TimeUTC(end) ~= endDate + hours(23)
        error('ACE:TimeRangeMismatch', 'Hourly result does not cover the requested UTC range.');
    end
    if any(seconds(diff(Hourly.TimeUTC)) ~= 3600)
        error('ACE:IrregularTimeGrid', 'Hourly timestamps are not unique consecutive 1-hour steps.');
    end

    sourceByDay = reshape(string(Hourly.Source), 24, []);
    if any(any(sourceByDay ~= sourceByDay(1, :)))
        error('ACE:MixedDailySource', 'At least one UTC day contains more than one product source.');
    end
    if any(sourceByDay(1, :).' ~= string(DailySelection.Source))
        error('ACE:DailySourceMismatch', 'Hourly and daily product selections are inconsistent.');
    end

    hasMean = isfinite(Hourly.B_mean_nT);
    hasMax = isfinite(Hourly.B_max_nT);
    if any(xor(hasMean, hasMax))
        error('ACE:StatisticPairMismatch', 'Mean and maximum availability are inconsistent.');
    end
    if any(hasMean ~= (Hourly.N_valid > 0))
        error('ACE:ValidCountMismatch', 'N_valid does not match statistic availability.');
    end
    if any(Hourly.B_max_nT(hasMean) < Hourly.B_mean_nT(hasMean))
        error('ACE:MaxBelowMean', 'At least one hourly maximum is below its mean.');
    end
    if any(Hourly.B_mean_nT(hasMean) < 0 | Hourly.B_max_nT(hasMean) > 500)
        error('ACE:StatisticOutOfRange', 'At least one saved statistic is outside 0..500 nT.');
    end
end


function paths = outputPaths(outputDir, baseName)
    paths.FinalMAT = fullfile(outputDir, baseName + ".mat");
    paths.HourlyCSV = fullfile(outputDir, baseName + ".csv");
    paths.DailyCSV = fullfile(outputDir, baseName + "_daily_selection.csv");
    paths.ManifestCSV = fullfile(outputDir, baseName + "_manifest.csv");
    paths.FileLogCSV = fullfile(outputDir, baseName + "_file_log.csv");
    paths.FigureFIG = fullfile(outputDir, baseName + ".fig");
    paths.FigurePNG = fullfile(outputDir, baseName + ".png");
end


function makeOverviewFigure(Hourly, figPath, pngPath, figureVisible)
    fig = figure('Visible', char(figureVisible), 'Color', 'w', ...
        'Position', [80, 80, 1800, 720]);
    cleanup = onCleanup(@() closeFigure(fig));
    ax = axes(fig);
    plot(ax, Hourly.TimeUTC, Hourly.B_mean_nT, ...
        'Color', [0.00, 0.35, 0.70], 'LineWidth', 0.75, ...
        'DisplayName', 'Hourly mean');
    hold(ax, 'on');
    plot(ax, Hourly.TimeUTC, Hourly.B_max_nT, ...
        'Color', [0.85, 0.25, 0.10], 'LineWidth', 0.60, ...
        'DisplayName', 'Hourly maximum');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, 'UTC time');
    ylabel(ax, '|B| (nT)');
    title(ax, sprintf('ACE/MFI hourly magnetic-field magnitude (%s to %s)', ...
        char(string(Hourly.TimeUTC(1), 'yyyy-MM-dd')), ...
        char(string(Hourly.TimeUTC(end), 'yyyy-MM-dd'))));
    legend(ax, 'Location', 'best');
    xlim(ax, [Hourly.TimeUTC(1), Hourly.TimeUTC(end)]);
    yLimits = ylim(ax);
    if yLimits(1) < 0
        ylim(ax, [0, yLimits(2)]);
    end
    set(ax, 'FontSize', 11, 'TickDir', 'out');
    savefig(fig, figPath);
    exportgraphics(fig, pngPath, 'Resolution', 300);
end


function closeFigure(fig)
    if isgraphics(fig)
        close(fig);
    end
end
