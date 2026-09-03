function result = Run_Case1_L1_Replots(baselineFile,cadence)
%Run_Case1_L1_Replots Overwrite requested V1 windows with missing PAD bins.
%   Baseline is immutable. Source selection and L1 fallback are in the
%   existing production runner; this entry selects only previously requested
%   events having at least one unrepresented valid PAD time bin.
if nargin < 2, cadence = 'day'; end
cadence = validatestring(cadence,{'day','hour'});
saved = load(baselineFile,'baseline'); baseline = saved.baseline;
suffix = "1d";
if strcmp(cadence,'hour'), suffix = "1h"; end
items = baseline.Items([baseline.Items.Cadence] == suffix);
selectedIDs = strings(0,1);
for ii = 1:numel(items)
    row = baseline.Catalog(string(baseline.Catalog.EventID) == items(ii).EventID,:);
    begin = row.StartUTC-days(7); finish = row.EndUTCExclusive+days(7);
    if strcmp(cadence,'day'), expected = days(finish-begin);
    else, expected = hours(finish-begin); end
    t = items(ii).Saved.pitchAngleTable;
    validBins = 0;
    if ~isempty(t)
        validBins = numel(unique(dateshift(t.EpochUTC(t.PADUsable), ...
            'start',cadence)));
    end
    if validBins < expected
        selectedIDs(end+1,1) = items(ii).EventID; %#ok<AGROW>
    end
end
fprintf('L1 retry: %s, %d requested events with missing PAD bins.\n',cadence,numel(selectedIDs));
result = struct('BaselineFile',string(baselineFile),'Cadence',string(cadence), ...
    'EventIDs',selectedIDs,'CreatedUTC',datetime('now','TimeZone','UTC'));
if ~isempty(selectedIDs)
    result.Run = Run_Case1_PPT_VerticalLine_Events_7d('RunPlots',true, ...
        'Overwrite',true,'Visible',false,'PADCadence',cadence, ...
        'EventIDs',selectedIDs,'ExportPeakPAD',strcmp(cadence,'hour'), ...
        'ContextDays',7,'LECPSourcePriority','l2_first');
end
cfg = Case1_Config;
folder = fullfile(cfg.DataRoot,'voyager1','lecp','validation','l1_fallback');
if ~isfolder(folder), mkdir(folder); end
stamp = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
result.File = string(fullfile(folder,['replot_',cadence,'_',stamp,'.mat']));
save(result.File,'result','-v7.3');
fprintf('Replot audit: %s\n',result.File);
end
