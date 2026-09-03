function audit = Case1_Audit_COHO_Archive(cfg, catalog, repair)
%Case1_Audit_COHO_Archive Compare Case 1 COHO sources with NASA raw bytes.
%   AUDIT = Case1_Audit_COHO_Archive(CFG, CATALOG, REPAIR) obtains the
%   official annual directory and original CDF for every required month.
%   SHA-256 compares complete files, without converting their contents.
%   Identical files are left untouched. REPAIR defaults to false. When true,
%   missing or different sources are installed only after CDF validation;
%   previous files are preserved in the matching COHO archive directory.
%
%   Coverage counts reuse Voyager_Read_CDF_Product and its current CDF
%   FILLVAL / VALIDMIN / VALIDMAX rules. Counts are diagnostic only. This
%   function performs no interpolation, averaging or gap replacement.
%   Only audit metadata are saved as MAT and JSON; scientific values remain
%   in the original CDF. The plot configuration and attitude are unchanged.
%
%   Author: Codex, following the manual MATLAB style in MMS_fu
%   Modified: 2026-09-03

%% input and output
if nargin < 1 || isempty(cfg), cfg = Case1_Config; end
if nargin < 2 || isempty(catalog), catalog = Case1_Event_Catalog; end
if nargin < 3 || isempty(repair), repair = false; end
validateattributes(repair, {'logical'}, {'scalar'});
Case1_Add_IRFU_Path(cfg.IRFURoot);

runID = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyyMMdd_HHmmss_SSS'));
auditRoot = fullfile(cfg.DataRoot, 'source_verification', ...
    'case1_coho', runID);
if ~isfolder(auditRoot), mkdir(auditRoot); end
auditFile = fullfile(auditRoot, 'Case1_COHO_archive_audit.mat');
jsonFile = fullfile(auditRoot, 'Case1_COHO_archive_audit.json');

audit = struct;
audit.RunID = runID;
audit.CheckedUTC = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
audit.DataRoot = cfg.DataRoot;
audit.RepairAuthorized = repair;
audit.AuditFile = auditFile;
audit.JSONFile = jsonFile;
audit.Method = ['NASA/SPDF annual directory; latest numeric version; ', ...
    'HTTPS original binary response; full SHA-256 comparison.'];
audit.CoverageRules = ['Counts only: use existing IRFU CDF reader; exact ', ...
    'FILLVAL, VALIDMIN/MAX, >=1e30 sentinel and nonfinite masking; ', ...
    'finite and strictly positive flux counts both recorded; no new ', ...
    'scientific screening, averaging, interpolation or gap replacement.'];
audit.Scope = ['COHO monthly merged products used by Case 1. Daily ', ...
    'sectored LECP CDF is outside this source audit.'];
audit.Months = struct([]);
audit.Events = struct([]);
audit.Kernels = struct([]);
audit.Complete = false;

%% unique required months from the current catalog
months = requiredMonths(cfg, catalog);
products = cell(numel(months), 1);
directoryCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

for ii = 1:numel(months)
    item = months(ii);
    record = struct;
    record.Spacecraft = item.Spacecraft;
    record.YearMonth = item.YearMonth;
    record.EventIDs = item.EventIDs;
    record.DirectoryURL = sprintf(['https://spdf.gsfc.nasa.gov/pub/', ...
        'data/voyager/voyager%d/coho1hr_magplasma/%04d/'], ...
        item.Spacecraft, item.Year);
    record.SourceURL = '';
    record.SourceFile = '';
    record.LatestVersion = NaN;
    record.ListedFiles = {};
    record.RemoteSHA256 = '';
    record.RemoteBytes = 0;
    record.RemoteLastModified = '';
    record.LocalBefore = localFiles(item.Folder);
    record.Status = 'pending';
    record.Error = '';
    record.Coverage = struct;

    try
        %% official file discovery
        if isKey(directoryCache, record.DirectoryURL)
            html = directoryCache(record.DirectoryURL);
        else
            [bytes, ~] = requestOriginal(record.DirectoryURL);
            html = native2unicode(bytes.', 'UTF-8');
            directoryCache(record.DirectoryURL) = html;
        end
        expression = sprintf(['voyager%d_coho1hr_merged_mag_plasma_', ...
            '%04d%02d\\d{2}_v\\d+\\.cdf'], ...
            item.Spacecraft, item.Year, item.Month);
        names = unique(regexp(html, expression, 'match'));
        if isempty(names)
            error('VoyagerCase1:OfficialMonthMissing', ...
                'No official monthly CDF at %s', record.DirectoryURL);
        end
        versions = cellfun(@(x) str2double(regexp(x, ...
            '(?<=_v)\d+(?=\.cdf)', 'match', 'once')), names);
        newest = names(versions == max(versions));
        nominal = sprintf('_%04d%02d01_', item.Year, item.Month);
        preferred = newest(contains(newest, nominal));
        if isscalar(preferred), newest = preferred; end
        if ~isscalar(newest)
            error('VoyagerCase1:OfficialMonthAmbiguous', ...
                'More than one latest official CDF for %s.', item.YearMonth);
        end
        fileName = newest{1};
        record.ListedFiles = names;
        record.LatestVersion = max(versions);
        record.SourceURL = [record.DirectoryURL, fileName];
        record.SourceFile = fullfile(item.Folder, fileName);

        %% complete original response; no CDF re-encoding
        [bytes, response] = requestOriginal(record.SourceURL);
        record.RemoteSHA256 = hashBytes(bytes);
        record.RemoteBytes = numel(bytes);
        header = response.getFields('Last-Modified');
        if ~isempty(header)
            record.RemoteLastModified = char(header(1).Value);
        end
        exact = isfile(record.SourceFile) && strcmp( ...
            Case1_File_SHA256(record.SourceFile), record.RemoteSHA256);
        if exact
            record.Status = 'identical_official_raw';
            products{ii} = Voyager_Read_CDF_Product(record.SourceFile, 'coho');
        else
            verificationFolder = fullfile(auditRoot, ...
                sprintf('voyager%d', item.Spacecraft), 'coho', ...
                '1hr', 'l2', 'merged_mag_plasma', ...
                sprintf('%04d', item.Year), sprintf('%02d', item.Month));
            if ~isfolder(verificationFolder), mkdir(verificationFolder); end
            verifiedFile = fullfile(verificationFolder, fileName);
            writeBytes(verifiedFile, bytes);
            cdfinfo(verifiedFile); % validate the downloaded CDF container
            products{ii} = Voyager_Read_CDF_Product(verifiedFile, 'coho');
            if ~strcmp(Case1_File_SHA256(verifiedFile), record.RemoteSHA256)
                error('VoyagerCase1:VerificationWrite', ...
                    'Downloaded bytes changed while writing %s.', verifiedFile);
            end
            if repair
                record.LocalBefore = archivePrevious( ...
                    record.LocalBefore, cfg.DataRoot, item, runID);
                if ~isfolder(item.Folder), mkdir(item.Folder); end
                [ok, message] = copyfile(verifiedFile, record.SourceFile, 'f');
                if ~ok, error('VoyagerCase1:COHOPublish', '%s', message); end
                if ~strcmp(Case1_File_SHA256(record.SourceFile), ...
                        record.RemoteSHA256)
                    error('VoyagerCase1:PublishedHash', ...
                        'Published CDF hash differs from NASA response.');
                end
                record.Status = 'installed_official_raw_previous_archived';
                products{ii}.source_file = record.SourceFile;
            else
                record.Status = 'different_or_missing_no_repair';
            end
        end
        record.Coverage = productCoverage(products{ii});
        if numel(record.LocalBefore) > 1 && exact
            record.Status = 'official_raw_present_extra_local_CDF_review';
        end
    catch exception
        record.Status = 'error';
        record.Error = exception.message;
    end

    if isempty(audit.Months), audit.Months = record;
    else, audit.Months(ii) = record;
    end
    save(auditFile, 'audit');
    fprintf('COHO audit %2d/%2d %s V%d: %s\n', ii, numel(months), ...
        item.YearMonth, item.Spacecraft, record.Status);
end

%% exact catalog-window coverage, without temporal resampling
for ii = 1:height(catalog)
    startUTC = catalog.StartUTC(ii) - days(cfg.ContextDays);
    stopUTC = catalog.EndUTCExclusive(ii) + days(cfg.ContextDays);
    event = struct('EventID', char(catalog.EventID(ii)), ...
        'Spacecraft', catalog.Spacecraft(ii), ...
        'StartUTC', utcText(startUTC), 'EndUTCExclusive', utcText(stopUTC), ...
        'Months', {{}}, 'CoverageParts', {{}});
    for jj = 1:numel(months)
        if months(jj).Spacecraft ~= event.Spacecraft || isempty(products{jj})
            continue
        end
        time = products{jj}.Epoch(:);
        inWindow = time >= startUTC & time < stopUTC;
        if ~any(inWindow), continue, end
        event.Months{end+1} = months(jj).YearMonth;
        event.CoverageParts{end+1} = productCoverage(products{jj}, inWindow);
    end
    if isempty(audit.Events), audit.Events = event;
    else, audit.Events(ii) = event;
    end
end

%% kernel provenance audit only; no kernel download or attitude changes
kernelTable = Case1_Attitude_Files(cfg);
manifestFolder = fullfile(cfg.DataRoot, 'voyager1', 'attitude', 'spice');
previous = dir(fullfile(manifestFolder, 'kernel_manifest_*.mat'));
kernelManifest = table;
if ~isempty(previous)
    [~, order] = sort({previous.name});
    priorFile = fullfile(previous(order(end)).folder, previous(order(end)).name);
    prior = load(priorFile, 'manifest');
    kernelManifest = prior.manifest;
    audit.PreviousKernelManifest = priorFile;
end
for ii = 1:height(kernelTable)
    item = struct('Role', char(kernelTable.Role(ii)), ...
        'SourceFile', char(kernelTable.SourceFile(ii)), ...
        'SourceURL', char(kernelTable.SourceURL(ii)), ...
        'Exists', isfile(kernelTable.SourceFile(ii)), ...
        'SHA256', '', 'OriginalManifestSHA256', '', ...
        'MatchesOriginalDownloadManifest', false);
    if item.Exists
        item.SHA256 = char(Case1_File_SHA256(item.SourceFile));
    end
    if ~isempty(kernelManifest)
        row = find(kernelManifest.SourceFile == string(item.SourceFile), 1);
        if ~isempty(row)
            item.OriginalManifestSHA256 = char(kernelManifest.SHA256(row));
            item.MatchesOriginalDownloadManifest = strcmp( ...
                item.SHA256, item.OriginalManifestSHA256);
        end
    end
    if isempty(audit.Kernels), audit.Kernels = item;
    else, audit.Kernels(ii) = item;
    end
end
audit.KernelAuditScope = ['Local hashes compared with the existing official ', ...
    'download manifest. Kernels were not refreshed or changed in this run.'];
audit.Complete = all(ismember({audit.Months.Status}, ...
    {'identical_official_raw', 'installed_official_raw_previous_archived'}));
save(auditFile, 'audit');
writelines(jsonencode(audit, 'PrettyPrint', true), jsonFile, 'Encoding', 'UTF-8');
fprintf('COHO source audit: %s\n', auditFile);
end

function months = requiredMonths(cfg, catalog)
months = struct([]);
keys = strings(0,1);
for ii = 1:height(catalog)
    first = dateshift(catalog.StartUTC(ii)-days(cfg.ContextDays), 'start', 'month');
    last = dateshift(catalog.EndUTCExclusive(ii)+days(cfg.ContextDays) ...
        - seconds(1), 'start', 'month');
    for time = first:calmonths(1):last
        key = sprintf('V%d-%04d-%02d', catalog.Spacecraft(ii), year(time), month(time));
        index = find(keys == key, 1);
        if isempty(index)
            index = numel(keys) + 1;
            keys(index,1) = string(key);
            months(index).Spacecraft = catalog.Spacecraft(ii);
            months(index).Year = year(time);
            months(index).Month = month(time);
            months(index).YearMonth = sprintf('%04d-%02d', year(time), month(time));
            months(index).EventIDs = {};
            months(index).Folder = fullfile(cfg.DataRoot, ...
                sprintf('voyager%d', catalog.Spacecraft(ii)), 'coho', ...
                '1hr', 'l2', 'merged_mag_plasma', ...
                sprintf('%04d', year(time)), sprintf('%02d', month(time)));
        end
        months(index).EventIDs{end+1} = char(catalog.EventID(ii));
    end
end
end

function files = localFiles(folder)
files = struct('SourceFile', {}, 'Bytes', {}, 'SHA256', {}, 'ArchiveFile', {});
entries = dir(fullfile(folder, '*.cdf'));
for ii = 1:numel(entries)
    files(ii).SourceFile = fullfile(entries(ii).folder, entries(ii).name);
    files(ii).Bytes = entries(ii).bytes;
    files(ii).SHA256 = char(Case1_File_SHA256(files(ii).SourceFile));
    files(ii).ArchiveFile = '';
end
end

function files = archivePrevious(files, dataRoot, item, runID)
archiveFolder = fullfile(dataRoot, sprintf('voyager%d', item.Spacecraft), ...
    'coho', '1hr', 'l2', 'merged_mag_plasma', 'archive', runID, ...
    sprintf('%04d', item.Year), sprintf('%02d', item.Month));
if ~isfolder(archiveFolder), mkdir(archiveFolder); end
for ii = 1:numel(files)
    [~, name, extension] = fileparts(files(ii).SourceFile);
    archiveFile = fullfile(archiveFolder, [name, extension]);
    [ok, message] = copyfile(files(ii).SourceFile, archiveFile);
    if ~ok, error('VoyagerCase1:COHOArchive', '%s', message); end
    if ~strcmp(Case1_File_SHA256(archiveFile), files(ii).SHA256)
        error('VoyagerCase1:COHOArchiveHash', 'Backup SHA-256 mismatch.');
    end
    files(ii).ArchiveFile = archiveFile;
end
% Each exact monthly source has a verified recoverable copy before removal.
% No recursive filesystem operation or broad directory target is used.
for ii = 1:numel(files)
    delete(files(ii).SourceFile);
end
end

function [bytes, response] = requestOriginal(url)
% IRFU irf.get_file uses websave; the MATLAB HTTP API additionally exposes
% original response bytes and allows direct access when a proxy truncates TLS.
options = matlab.net.http.HTTPOptions('UseProxy', false, ...
    'ConvertResponse', false, 'ConnectTimeout', 25, 'ResponseTimeout', 45);
for attempt = 1:3
    try
        response = matlab.net.http.RequestMessage('get').send(url, options);
        if response.StatusCode ~= matlab.net.http.StatusCode.OK
            error('VoyagerCase1:COHOHTTP', 'HTTP %s: %s', ...
                string(response.StatusCode), url);
        end
        bytes = response.Body.Data;
        if endsWith(url, '/') && (ischar(bytes) || isstring(bytes))
            bytes = unicode2native(char(bytes), 'UTF-8');
        end
        if ~isa(bytes, 'uint8') || isempty(bytes)
            error('VoyagerCase1:COHOHTTPBytes', 'No original binary response: %s', url);
        end
        bytes = bytes(:);
        return
    catch exception
        if attempt == 3, rethrow(exception); end
        pause(1);
    end
end
end

function value = hashBytes(bytes)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(bytes(:), 'int8'));
value = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end

function writeBytes(fileName, bytes)
fid = fopen(fileName, 'wb');
if fid < 0, error('VoyagerCase1:COHOWrite', '%s', fileName); end
cleanup = onCleanup(@() fclose(fid));
count = fwrite(fid, bytes, 'uint8');
if count ~= numel(bytes)
    error('VoyagerCase1:COHOWriteCount', 'Incomplete CDF write: %s', fileName);
end
end

function coverage = productCoverage(product, useRows)
time = product.Epoch(:);
if nargin < 2, useRows = true(size(time)); end
time = time(useRows);
coverage = struct('Records', numel(time), 'FirstUTC', '', 'LastUTC', '', ...
    'UniqueTimeRecords', 0, 'NonHourlySteps', 0, ...
    'BVectorFinite', 0, 'BVectorFiniteNonzero', 0, 'Variables', struct([]));
if isempty(time), return, end
coverage.FirstUTC = utcText(min(time));
coverage.LastUTC = utcText(max(time));
coverage.UniqueTimeRecords = numel(unique(time));
coverage.NonHourlySteps = nnz(seconds(diff(time)) ~= 3600);
if all(isfield(product, {'BR', 'BT', 'BN'}))
    field = [product.BR(:), product.BT(:), product.BN(:)];
    field = field(useRows,:);
    valid = all(isfinite(field), 2);
    coverage.BVectorFinite = nnz(valid);
    coverage.BVectorFiniteNonzero = nnz(valid & sum(field.^2,2) > 0);
end
names = fieldnames(product.variable_meta);
index = 0;
for ii = 1:numel(names)
    name = names{ii};
    if strcmp(name, 'Epoch'), continue, end
    values = product.(name);
    if ~isnumeric(values) || numel(values) ~= numel(useRows), continue, end
    values = values(useRows);
    valid = isfinite(values);
    index = index + 1;
    row = struct('Name', name, 'Finite', nnz(valid), ...
        'Positive', nnz(valid & values > 0), 'Missing', nnz(~valid), ...
        'FirstFiniteUTC', '', 'LastFiniteUTC', '', 'Units', '');
    if any(valid)
        row.FirstFiniteUTC = utcText(time(find(valid,1)));
        row.LastFiniteUTC = utcText(time(find(valid,1,'last')));
    end
    meta = product.variable_meta.(name);
    if isfield(meta.attributes, 'UNITS')
        row.Units = char(string(meta.attributes.UNITS));
    end
    if isempty(coverage.Variables), coverage.Variables = row;
    else, coverage.Variables(index) = row;
    end
end
end

function text = utcText(time)
time.Format = 'yyyy-MM-dd''T''HH:mm:ss''Z''';
text = char(time);
end
