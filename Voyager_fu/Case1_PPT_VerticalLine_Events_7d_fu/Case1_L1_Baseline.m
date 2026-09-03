function baseline = Case1_L1_Baseline()
%Case1_L1_Baseline Preserve only the requested V1 figures before overwrite.
%   Numerical provenance stays in the classified Voyager data archive.
cfg = Case1_Config;
catalog = Case1_Event_Catalog;
catalog = catalog(catalog.Spacecraft == 1,:);
baseline = struct('CreatedUTC',datetime('now','TimeZone','UTC'), ...
    'Catalog',catalog,'Items',struct([]));
cadences = {'1d','1h'};
kk = 0;
for cc = 1:2
    folder = fullfile(cfg.DataRoot,'voyager1','lecp',cadences{cc}, ...
        'derived','pitch_angle','2013-2021','predicted_ck');
    for ii = 1:height(catalog)
        entries = dir(fullfile(folder,sprintf('V1_%s_*_%s_nativeCDF_Epoch.mat', ...
            string(catalog.EventID(ii)),cadences{cc})));
        assert(numel(entries) == 1,'Expected exactly one existing event audit.');
        file = fullfile(entries.folder,entries.name);
        saved = load(file);
        kk = kk+1;
        baseline.Items(kk).EventID = string(catalog.EventID(ii));
        baseline.Items(kk).Cadence = string(cadences{cc});
        baseline.Items(kk).AuditFile = string(file);
        baseline.Items(kk).AuditSHA256 = string(Case1_File_SHA256(file));
        baseline.Items(kk).Saved = saved;
        baseline.Items(kk).Rows = height(saved.pitchAngleTable);
        baseline.Items(kk).Usable = 0;
        if ~isempty(saved.pitchAngleTable)
            baseline.Items(kk).Usable = nnz(saved.pitchAngleTable.PADUsable);
        end
    end
end
folder = fullfile(cfg.DataRoot,'voyager1','lecp','validation','l1_fallback');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
baseline.File = string(fullfile(folder,['baseline_',stamp,'.mat']));
save(baseline.File,'baseline','-v7.3');
fprintf('Saved %d original event audits: %s\n',numel(baseline.Items),baseline.File);
end
