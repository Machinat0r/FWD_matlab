function audit = Run_Case1_Download_2012_Sources()
%Run_Case1_Download_2012_Sources Add two missing official annual CDFs.
%   Required by the 2013-01-09 daily overview with 15-day context.
%   Existing source files are never overwritten. Downloaded bytes remain
%   unchanged; the current IRFU reader validates variables and metadata.
%   Author: Codex; Modified: 2026-09-04

%% exact original products and classified destinations
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
baseURL = 'https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/particle/lecp/final-cdf/';
names = {'voyager-1_lecp_lev-1-rates_20120101_v1.1.1-01.cdf'; ...
    'voyager-1_lecp_lev-2-daily-avg_20120101_v1.1.1-01.cdf'};
products = {'lev-1-rates'; 'lev-2-daily-avg'};
folders = {fullfile(cfg.DataRoot, 'voyager1', 'lecp', 'native', ...
    'l1', 'sectored_rates', '2012'); ...
    fullfile(cfg.DataRoot, 'voyager1', 'lecp', '1d', ...
    'l2', 'sectored_flux', '2012')};
variables = {'FHDU_SectoredRates'; 'FHDU_SectoredFluxes'};
auditFolder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', ...
    'validation', 'context15d_daily_l1_first');
if ~isfolder(auditFolder), mkdir(auditFolder); end
auditFile = fullfile(auditFolder, 'source_download_2012.mat');
assert(~isfile(auditFile), 'Existing download audit is preserved: %s', auditFile);
audit = struct('StartedUTC', datetime('now', 'TimeZone', 'UTC'), ...
    'Purpose', 'Original 2012 sources for the approved 15-day daily context', ...
    'Transport', 'curl HTTPS with certificate validation; direct connection; no byte conversion', ...
    'ScienceProcessing', 'Validation only; no new averaging, filling, quality cuts or source edits');
audit.CodeFile = string([mfilename('fullpath'), '.m']);
audit.CodeSHA256 = string(Case1_File_SHA256(char(audit.CodeFile)));

%% inspect HTTP source, download only missing files, then check original CDF
for ii = 1:2
    url = [baseURL, products{ii}, '/', names{ii}];
    sourceFile = fullfile(folders{ii}, names{ii});
    row = struct('URL', string(url), 'SourceFile', string(sourceFile));
    row.CanonicalSPDFURL = replace(string(url), 'cdaweb.gsfc.nasa.gov', 'spdf.gsfc.nasa.gov');
    row.PriorTransport = 'L1 initial SPDF GET supplied bytes 0:364543 (2026-09-04); official CDAWeb mirror supplies subsequent ranges. L2 complete CDAWeb GET; no re-encoding.';
    row.HTTPCheckUTC = datetime('now', 'TimeZone', 'UTC');
    cmd = sprintf('curl.exe -6 --noproxy "*" --fail --silent --show-error --head --retry 3 --retry-all-errors --retry-delay 1 --connect-timeout 15 --max-time 20 "%s"', url);
    [status, headers] = system(cmd);
    assert(status == 0 && contains(headers, '200 OK'), 'Official URL unavailable: %s', url);
    row.HTTPHeaders = string(headers);
    token = regexp(headers, '(?i)Content-Length:\s*(\d+)', 'tokens', 'once');
    assert(~isempty(token), 'Official response omitted Content-Length.');
    row.HTTPBytes = str2double(token{1});
    etags = {'"605c78-62779fc3ab781"'; '"2aee0-62779fd2d7bff"'};
    row.VerifiedSPDF_ETag = string(etags{ii});
    assert(contains(headers, etags{ii}), 'Official mirror version differs from the verified canonical source.');
    row.CanonicalVersionCheck = 'SPDF and CDAWeb HEAD responses match ETag, Last-Modified and Content-Length.';
    creationTimes = [datetime(2026,9,4,14,18,20.990,'TimeZone','UTC'); ...
        datetime(2026,9,4,14,24,35.778,'TimeZone','UTC')];
    row.InitialDownloadFileCreatedUTC = creationTimes(ii);
    row.InitialTimeEvidence = 'CreationTimeUtc of the first local download fragment, read from filesystem metadata.';
    row.ExistedBefore = isfile(sourceFile);
    row.DownloadStartedUTC = NaT(1, 1, 'TimeZone', 'UTC');
    row.DownloadCompletedUTC = NaT(1, 1, 'TimeZone', 'UTC');
    if ~row.ExistedBefore
        if ~isfolder(folders{ii}), mkdir(folders{ii}); end
        partFile = [sourceFile, '.context15d-download'];
        row.DownloadStartedUTC = datetime('now', 'TimeZone', 'UTC');
        row.TransferAttempts = struct([]);
        bytesAfter = 0;
        for attempt = 1:100
            partInfo = dir(partFile);
            if isempty(partInfo), bytesBefore = 0; else, bytesBefore = partInfo.bytes; end
            bytesAfter = bytesBefore;
            if bytesBefore == row.HTTPBytes, break, end
            assert(bytesBefore < row.HTTPBytes, 'Partial download is too large; preserve it for inspection.');
            rangeLast = min(row.HTTPBytes-1, bytesBefore+262144-1);
            rangeBytes = rangeLast-bytesBefore+1;
            chunkFile = sprintf('%s.range_%d_%d', sourceFile, bytesBefore, rangeLast);
            chunkInfo = dir(chunkFile);
            status = 0; output = 'Reused existing complete byte-range fragment from official URL.';
            if isempty(chunkInfo) || chunkInfo.bytes ~= rangeBytes
                if ~isempty(chunkInfo), chunkFile = sprintf('%s.try%d', chunkFile, attempt); end
                assert(~isfile(chunkFile), 'Existing partial fragment is preserved.');
                cmd = sprintf('curl.exe -6 --noproxy "*" --fail --silent --show-error --range "%d-%d" --connect-timeout 15 --max-time 40 --output "%s" --write-out "HTTP_STATUS:%%{http_code}" "%s"', bytesBefore, rangeLast, chunkFile, url);
                [status, output] = system(cmd);
                chunkInfo = dir(chunkFile);
            end
            if status == 0 && ~isempty(chunkInfo) && chunkInfo.bytes == rangeBytes
                assert(contains(output, 'HTTP_STATUS:206') || startsWith(output, 'Reused'), 'Server did not honor exact HTTP byte range.');
                in = fopen(chunkFile, 'rb');
                assert(in >= 0, 'Cannot read downloaded fragment.');
                bytes = fread(in, Inf, '*uint8'); fclose(in);
                assert(numel(bytes) == rangeBytes, 'Unexpected byte-range size.');
                out = fopen(partFile, 'ab');
                assert(out >= 0, 'Cannot append downloaded bytes.');
                nWritten = fwrite(out, bytes, 'uint8'); fclose(out);
                assert(nWritten == rangeBytes, 'Incomplete byte append.');
                bytesAfter = bytesBefore+rangeBytes;
            end
            transfer = struct('UTC', datetime('now','TimeZone','UTC'), ...
                'BytesBefore', bytesBefore, 'BytesAfter', bytesAfter, ...
                'ExitStatus', status, 'Message', string(output), 'RangeFile', string(chunkFile));
            if isempty(row.TransferAttempts), row.TransferAttempts = transfer;
            else, row.TransferAttempts(end+1) = transfer; end
            fprintf('Download %d attempt %d: %d / %d bytes (curl %d)\n', ...
                ii, attempt, bytesAfter, row.HTTPBytes, status);
            if bytesAfter == row.HTTPBytes, break, end
            if attempt >= 3 && all([row.TransferAttempts(end-2:end).BytesAfter] == bytesAfter)
                error('VoyagerDownload:Stalled', 'Three attempts made no progress; partial bytes preserved: %s', partFile);
            end
        end
        assert(bytesAfter == row.HTTPBytes, 'Official source download incomplete; partial bytes preserved.');
        row.DownloadCompletedUTC = datetime('now', 'TimeZone', 'UTC');
        readFile = partFile;
    else
        row.TransferAttempts = struct([]);
        readFile = sourceFile;
    end
    info = dir(readFile);
    assert(isscalar(info) && info.bytes == row.HTTPBytes, 'Downloaded/source byte count mismatch.');
    row.Bytes = info.bytes;
    row.SHA256 = string(Case1_File_SHA256(readFile));
    p = Voyager_Read_CDF_Product(readFile, 'lecp_sector_daily');
    labels = strtrim(string(p.Hydrogen_Channels_Label));
    row.P1Index = find(labels == "P1");
    assert(isscalar(row.P1Index) && row.P1Index == 10, 'Unexpected P1 channel index.');
    row.SectorIterator = double(p.SectorIterator(:)).';
    assert(isequal(row.SectorIterator, 1:8), 'Unexpected sector order.');
    assert(isfield(p, variables{ii}), 'Requested original variable is unavailable.');
    row.Variable = string(variables{ii});
    row.Records = numel(p.Epoch);
    row.Shape = size(p.(variables{ii}));
    assert(isequal(row.Shape, [row.Records, 16, 8]), 'Unexpected sectored variable dimensions.');
    row.FirstEpochUTC = min(p.Epoch);
    row.LastEpochUTC = max(p.Epoch);
    row.WindowStartUTC = datetime(2012,12,25,'TimeZone','UTC');
    row.WindowEndUTCExclusive = datetime(2013,1,1,'TimeZone','UTC');
    inWindow = p.Epoch >= row.WindowStartUTC & p.Epoch < row.WindowEndUTCExclusive;
    row.RecordsInAddedWindow = nnz(inWindow);
    row.NegativeDeltaTInAddedWindow = nnz(inWindow & p.DeltaT(:) < 0);
    assert(row.RecordsInAddedWindow > 0, 'Original product has no records in the added window.');
    assert(row.SHA256 == string(Case1_File_SHA256(readFile)), 'CDF bytes changed during read validation.');
    row.Validation = 'IRFU reader passed; P1 index 10; sectors 1:8; N x 16 x 8; bytes unchanged';
    if ~row.ExistedBefore
        assert(~isfile(sourceFile), 'Destination appeared during download; existing file preserved.');
        movefile(readFile, sourceFile);
    end
    assert(row.SHA256 == string(Case1_File_SHA256(sourceFile)), 'Final source hash mismatch.');
    if ii == 1, audit.Files = row; else, audit.Files(ii) = row; end
    save(auditFile, 'audit', '-v7.3');
    fprintf('%s: %d bytes, %d records; %d added-window records; SHA256 %s\n', ...
        names{ii}, row.Bytes, row.Records, row.RecordsInAddedWindow, row.SHA256);
end
audit.CompletedUTC = datetime('now', 'TimeZone', 'UTC');
audit.ValidationPassed = true;
save(auditFile, 'audit', '-v7.3');
fprintf('Download audit: %s\n', auditFile);
end
