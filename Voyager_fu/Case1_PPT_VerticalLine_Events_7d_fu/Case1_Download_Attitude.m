function manifest = Case1_Download_Attitude()
%Case1_Download_Attitude Download missing official Voyager attitude kernels.
%   Preserves existing files. All raw files and the source audit stay on Z:.
%   This downloads no flux CDF and makes no LECP geometry assumption.

%% configuration
cfg = Case1_Config;
manifest = Case1_Attitude_Files(cfg);
manifest.DownloadedThisRun = false(height(manifest), 1);
manifest.Bytes = zeros(height(manifest), 1);
manifest.SHA256 = strings(height(manifest), 1);
manifest.CheckedUTC = repmat(datetime('now', 'TimeZone', 'UTC'), ...
    height(manifest), 1);

%% original files; never overwrite an existing source
for ii = 1:height(manifest)
    fileName = char(manifest.SourceFile(ii));
    if ~isfile(fileName)
        folder = fileparts(fileName);
        if ~isfolder(folder), mkdir(folder); end
        temporaryFile = [tempname(folder), '.download'];
        cleanup = onCleanup(@() removePartial(temporaryFile));
        fprintf('Download %s\n', manifest.SourceURL(ii));
        websave(temporaryFile, char(manifest.SourceURL(ii)), ...
            weboptions('Timeout', 180));
        checkHeader(temporaryFile, manifest.Role(ii));
        [ok, message] = movefile(temporaryFile, fileName);
        if ~ok, error('VoyagerCase1:KernelPublish', '%s', message); end
        clear cleanup
        manifest.DownloadedThisRun(ii) = true;
    end
    checkHeader(fileName, manifest.Role(ii));
    info = dir(fileName);
    manifest.Bytes(ii) = info.bytes;
    manifest.SHA256(ii) = Case1_File_SHA256(fileName);
end

%% audit only; no intermediate CSV of attitude values
folder = fullfile(cfg.DataRoot, 'voyager1', 'attitude', 'spice');
auditFile = fullfile(folder, ['kernel_manifest_', ...
    char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyyMMdd_HHmmss_SSS')), '.mat']);
save(auditFile, 'manifest');
fprintf('Kernel source audit: %s\n', auditFile);
end

function checkHeader(fileName, role)
fid = fopen(fileName, 'rb');
if fid < 0, error('VoyagerCase1:KernelOpen', '%s', fileName); end
cleanup = onCleanup(@() fclose(fid));
header = char(fread(fid, 512, '*uint8').');
switch role
    case 'CK', valid = startsWith(header, 'DAF/CK');
    case 'SPK', valid = startsWith(header, 'DAF/SPK');
    case 'SCLK', valid = startsWith(header, 'KPL/SCLK');
    case 'FK', valid = startsWith(header, 'KPL/FK');
    case 'PCK', valid = startsWith(header, 'KPL/PCK');
    case 'LSK', valid = startsWith(header, 'KPL/LSK');
    otherwise, valid = contains(header, 'Voyager-1');
end
if ~valid
    error('VoyagerCase1:KernelHeader', 'Unexpected %s header: %s', ...
        role, fileName);
end
end

function removePartial(fileName)
if isfile(fileName), delete(fileName); end
end
