function result = Case1_Validate_2012_Download_Integrity()
%Case1_Validate_2012_Download_Integrity Compare ranged and complete responses.
%   The separate complete official L1 HTTP response is byte-identical to
%   the original source reconstructed from HTTP byte ranges. No data edits.
%   Author: Codex; Modified: 2026-09-04

%% existing source audit and independently downloaded complete HTTP body
cfg = Case1_Config;
folder = fullfile(cfg.DataRoot, 'voyager1', 'lecp', ...
    'validation', 'context15d_daily_l1_first');
auditFile = fullfile(folder, 'source_download_2012.mat');
saved = load(auditFile, 'audit');
audit = saved.audit;
assert(audit.ValidationPassed && numel(audit.Files) == 2);
checkFile = fullfile(folder, ...
    'voyager-1_lecp_lev-1-rates_20120101_v1.1.1-01.cdf.full-response-check');
info = dir(checkFile);
assert(isscalar(info) && info.bytes == audit.Files(1).Bytes);
sha = string(Case1_File_SHA256(checkFile));
assert(sha == audit.Files(1).SHA256, 'Complete HTTP response and joined ranges differ.');
for ii = 1:2
    assert(string(Case1_File_SHA256(char(audit.Files(ii).SourceFile))) == audit.Files(ii).SHA256);
end

%% retain the independent retrieval, exact hash and validation provenance
result = struct('SourceURL', audit.Files(1).URL, 'SourceFile', string(checkFile), ...
    'Retrieval', 'Separate complete HTTPS GET through curl -6, with certificate validation', ...
    'Bytes', info.bytes, 'SHA256', sha, ...
    'ValidatedUTC', datetime('now', 'TimeZone', 'UTC'), ...
    'CompleteResponseMatchesRangedSource', true);
result.CodeFile = string([mfilename('fullpath'), '.m']);
result.CodeSHA256 = string(Case1_File_SHA256(char(result.CodeFile)));
audit.IndependentFullResponseCheck = result;
save(auditFile, 'audit', '-v7.3');
disp(result);
end
