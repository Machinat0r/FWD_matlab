function audit = Case1_Audit_L1_Peak_Example()
%Case1_Audit_L1_Peak_Example Trace one large replot peak to original rates.
%   Diagnostic only. No new quality flag rule or outlier removal is applied.
cfg = Case1_Config; Case1_Add_IRFU_Path(cfg.IRFURoot);
file = fullfile(cfg.PeakPADDataFolder, ...
    'V1_Case1-S03-L02_20190527_P1_peak5_hourly.mat');
saved = load(file,'peakAudit'); p = saved.peakAudit;
selected = p.SelectedRows{3}; records = selected.L1SourceRecords{1};
assert(numel(unique(records.SourceCDF)) == 1);
cdf = dataobj(char(records.SourceCDF(1)));
quality = getv(cdf,'FHDU_SectoredQuality');
rate = getv(cdf,'FHDU_SectoredRates');
audit = struct('EventID',p.EventID,'PeakAuditFile',string(file), ...
    'PeakRow',selected,'SourceRecords',records, ...
    'P1Quality',reshape(quality.data(records.SourceCDFRecord,10,:),height(records),8), ...
    'QualityMetadata',rmfield(quality,'data'), ...
    'RateMetadata',rmfield(rate,'data'), ...
    'SourceSHA256',string(Case1_File_SHA256(char(records.SourceCDF(1)))), ...
    'Note',"Original large rate values retained under existing CDF bounds and approved rules. No physical interpretation or extra quality/outlier filter is inferred.");
folder = fullfile(cfg.DataRoot,'voyager1','lecp','validation','l1_fallback');
audit.File = string(fullfile(folder,'large_peak_Case1-S03-L02.mat'));
save(audit.File,'audit','-v7.3');
fprintf('Diagnostic saved: %s\n',audit.File);
end
