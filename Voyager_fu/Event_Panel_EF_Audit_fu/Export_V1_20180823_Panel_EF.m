function report = Export_V1_20180823_Panel_EF
%Export_V1_20180823_Panel_EF Export sources of the existing hourly event figure.
% Uses IRFU dataobj and the existing Voyager_Read_CDF_Product reader.
% Selects original Epoch in [2018-08-20,2018-08-27) UTC. No new scientific
% averaging, filtering, PA calculation, calibration or plotting is performed.
% Raw CSV values retain CDF fill values. Existing plot values are separate.
% Full original CDF copies preserve native time encoding and all variables.

%% paths and original saved plot audit
addpath(genpath('C:/Users/Administrator/Documents/irfu-matlab-master'));
addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Case1_PPT_VerticalLine_Events_7d_fu');
root = 'Z:/SPART-WORK/Data/Voyager';
out = fullfile(root,'voyager1/derived/panel_ef_audit/1h/2018/20180823');
if ~isfolder(out), mkdir(out); end
auditFile = fullfile(root,['voyager1/lecp/1h/derived/pitch_angle/2013-2021/predicted_ck/' ...
    'V1_Case1-S02-L06_20180823_20180823_COHO1h_raw_LECP_P1_pitch_angle_predictedCK_1h_nativeCDF_Epoch.mat']);
original = load(auditFile);
startUTC = datetime(2018,8,20,'TimeZone','UTC');
endUTC = datetime(2018,8,27,'TimeZone','UTC');
assert(original.opts.ContextDays == 3);
assert(strcmp(original.opts.LECPSourcePriority,'l1_first'));
files = {fullfile(root,'voyager1/coho/1hr/l2/merged_mag_plasma/2018/08/voyager1_coho1hr_merged_mag_plasma_20180801_v01.cdf'); ...
    fullfile(root,'voyager1/lecp/native/l1/sectored_rates/2018/voyager-1_lecp_lev-1-rates_20180101_v1.1.1-01.cdf'); ...
    fullfile(root,'voyager1/lecp/1h/l2/sectored_flux/2018/voyager-1_lecp_lev-2-hourly-avg_20180101_v1.1.1-01.cdf')};
productNames = {'coho_l2_1h','lecp_l1_native','lecp_l2_1h'};
raw = struct;
metadata = struct;
for i = 1:3
    obj = dataobj(files{i});
    epoch = getv(obj,'Epoch');
    if contains(lower(string(epoch.type)),"tt2000")
        epochUnix = EpochTT(int64(epoch.data(:))).epochUnix;
        epochMeaning = 'Native CDF TT2000 int64 values retained in Epoch; EpochUTC uses IRFU EpochTT conversion.';
    else
        epochUnix = double(epoch.data(:));
        epochMeaning = 'IRFU dataobj converts CDF_EPOCH to POSIX seconds; native encoding retained in original CDF copy.';
    end
    time = datetime(epochUnix,'ConvertFrom','posixtime','TimeZone','UTC');
    rows = find(time >= startUTC & time < endUTC);
    vars = {'Epoch','DeltaT','FHDU_SectoredRates','FHDU_SectoredRateUncertainties', ...
        'FHDU_SectoredFluxes','FHDU_SectoredFluxUncertainties','FHDU_SectoredQuality', ...
        'FHDU_Energy','FHDU_EnergyRange','Hydrogen_Channels','Hydrogen_Channels_Label', ...
        'SectorIterator','SectorIterator_Label','F','ABS_B','BR','BT','BN', ...
        'elevAngle','azimuthAngle','protonFlux1_LECP','protonFlux2_LECP','protonFlux3_LECP'};
    available = string(obj.Variables(:,1));
    item = struct('SourceCDF',files{i},'SourceCDFRecord',rows,'EpochUTC',time(rows));
    meta = struct('GlobalAttributes',obj.GlobalAttributes,'Variables',struct);
    for k = 1:numel(vars)
        name = vars{k};
        if ~any(available == string(name)), continue, end
        v = getv(obj,name);
        meta.Variables.(name) = rmfield(v,'data');
        value = v.data;
        if v.nrec == numel(time)
            value = value(rows,:,:);
        end
        item.(name) = value;
    end
    item.EpochDataMeaning = epochMeaning;
    raw.(productNames{i}) = item;
    metadata.(productNames{i}) = meta;
    copyDir = fullfile(out,'source_cdf',productNames{i});
    if ~isfolder(copyDir), mkdir(copyDir); end
    copyfile(files{i},copyDir);
end
copyfile(auditFile,fullfile(out,'original_figure_audit.mat'));

%% panel e: original records plus existing reader values
r = raw.coho_l2_1h;
panelE = table(r.SourceCDFRecord,r.EpochUTC,'VariableNames',{'SourceCDFRecord','EpochUTC'});
clean = Voyager_Read_CDF_Product(files{1},'coho');
for i = 1:3
    name = sprintf('protonFlux%d_LECP',i);
    panelE.([name '_raw']) = double(r.(name));
    panelE.([name '_plot']) = clean.(name)(r.SourceCDFRecord);
end
for name = {'F','ABS_B','BR','BT','BN','elevAngle','azimuthAngle'}
    if isfield(r,name{1}), panelE.([name{1} '_raw']) = double(r.(name{1})); end
end
writeCSV(panelE,fullfile(out,'panel_e_COHO_original_records.csv'));

%% panel f: raw P1 (CDF channel index 10), including S8 for provenance
panelF_L1 = rawP1Table(raw.lecp_l1_native,'Rate');
panelF_L2 = rawP1Table(raw.lecp_l2_1h,'Flux');
writeCSV(panelF_L1,fullfile(out,'panel_f_L1_P1_original_records.csv'));
writeCSV(panelF_L2,fullfile(out,'panel_f_L2_P1_original_records.csv'));

%% exact cached plotting input, without recomputing PA or flux
panelF = original.pitchAngleTable;
flatF = removevars(panelF,'L1SourceRecords');
flatF.L1SourceCDFRecords = strings(height(panelF),1);
for i = 1:height(panelF)
    source = panelF.L1SourceRecords{i};
    if istable(source) && ~isempty(source)
        flatF.L1SourceCDFRecords(i) = join(string(source.SourceCDFRecord(:)'),';');
    end
end
writeCSV(flatF,fullfile(out,'panel_f_actual_plot_inputs.csv'));
candidates = original.l1FallbackAudit.Candidates;
flatCandidates = removevars(candidates,{'L1Rows','L2Rows','SourceRecords'});
writeCSV(flatCandidates,fullfile(out,'panel_f_existing_L1_hour_means.csv'));
mapping = table;
for i = 1:height(candidates)
    source = candidates.SourceRecords{i};
    if isempty(source), continue, end
    if ismember('IdenticalRecordAliases',source.Properties.VariableNames)
        source = removevars(source,'IdenticalRecordAliases');
    end
    source.CandidateIndex = repmat(i,height(source),1);
    source.CandidateApplied = repmat(candidates.Applied(i),height(source),1);
    source.CandidateBinStartUTC = repmat(candidates.BinStartUTC(i),height(source),1);
    mapping = [mapping;source]; %#ok<AGROW>
end
writeCSV(mapping,fullfile(out,'panel_f_L1_contribution_mapping.csv'));

%% comparison aligned only by the existing UTC-hour assignment
comparison = table(panelF.EpochUTC,panelF.MAGBinStartUTC,panelF.SourceProduct,panelF.PADUsable, ...
    'VariableNames',{'PanelF_EpochUTC','UTC_HourStart','PanelF_SourceProduct','PanelF_PADUsable'});
[found,eIndex] = ismember(comparison.UTC_HourStart,panelE.EpochUTC);
comparison.PanelE_P1 = nan(height(comparison),1);
comparison.PanelE_P1(found) = panelE.protonFlux1_LECP_plot(eIndex(found));
fNames = arrayfun(@(s)sprintf('Flux_S%d_1h',s),1:7,'UniformOutput',false);
comparison.PanelF_S1_to_S7 = panelF{:,fNames};
comparison.PanelF_S4_L1_MeanRate = panelF.RawRate_S4_1h;
comparison.PanelF_ConversionFactor = panelF.SourceToDifferentialFluxFactor;
[found,l2Index] = ismember(comparison.UTC_HourStart,dateshift(panelF_L2.EpochUTC,'start','hour'));
comparison.Original_L2_S4 = nan(height(comparison),1);
comparison.Original_L2_CDFRecord = nan(height(comparison),1);
comparison.Original_L2_S4(found) = panelF_L2.Flux_S4_raw(l2Index(found));
comparison.Original_L2_CDFRecord(found) = panelF_L2.SourceCDFRecord(l2Index(found));
writeCSV(comparison,fullfile(out,'panel_e_f_same_UTC_hour_comparison.csv'));

%% lightweight verification of export content against the exact inputs
assert(height(panelE) == 168 && height(panelF_L1) == 110 && height(panelF_L2) == 46);
assert(height(panelF) == 46 && nnz(panelF.PADUsable) == 17);
assert(all(panelF.SourceProduct(panelF.PADUsable) == "L1_UTC_mean"));
back = readtable(fullfile(out,'panel_e_COHO_original_records.csv'),'VariableNamingRule','preserve');
assert(all(abs(back.protonFlux2_LECP_raw - panelE.protonFlux2_LECP_raw) <= ...
    1e-12 * max(1,abs(panelE.protonFlux2_LECP_raw))));
report = struct('OutputFolder',out,'WindowStartUTC',char(string(startUTC)), ...
    'WindowEndUTCExclusive',char(string(endUTC)),'PanelERecords',height(panelE), ...
    'PanelEValidPerChannel',sum(isfinite(panelE{:,{'protonFlux1_LECP_plot','protonFlux2_LECP_plot','protonFlux3_LECP_plot'}})), ...
    'L1OriginalRecords',height(panelF_L1),'L2OriginalRecords',height(panelF_L2), ...
    'PanelFInputRecords',height(panelF),'PanelFVisiblePADRecords',nnz(panelF.PADUsable), ...
    'L1AppliedHourCount',nnz(candidates.Applied),'P1ChannelIndex',10);
[peak,index] = max(panelE.protonFlux1_LECP_plot);
report.PanelE_P1_Peak = peak;
report.PanelE_P1_PeakUTC = char(string(panelE.EpochUTC(index)));
report.MetadataNote = 'LECP P1 source band 0.57-0.89 MeV; displayed band and historical L1 conversion use 0.57-1.78 MeV. Equivalence unresolved.';
save(fullfile(out,'event_raw_and_plot_data.mat'),'raw','metadata','panelE','panelF_L1','panelF_L2','panelF', ...
    'candidates','mapping','comparison','report','startUTC','endUTC','-v7.3');
writeJSON(metadata,fullfile(out,'source_variable_metadata.json'));
writeJSON(report,fullfile(out,'export_summary.json'));
disp(report);
disp(comparison(ismember(comparison.UTC_HourStart, ...
    [datetime(2018,8,23,10,0,0,'TimeZone','UTC');datetime(2018,8,24,8,0,0,'TimeZone','UTC')]),:));
end

function result = rawP1Table(r,kind)
labels = strtrim(string(r.Hydrogen_Channels_Label));
assert(labels(10) == "P1");
if strcmp(kind,'Rate')
    value = r.FHDU_SectoredRates;
    sigma = r.FHDU_SectoredRateUncertainties;
else
    value = r.FHDU_SectoredFluxes;
    sigma = r.FHDU_SectoredFluxUncertainties;
end
result = table(r.SourceCDFRecord,r.EpochUTC,double(r.DeltaT(:)), ...
    'VariableNames',{'SourceCDFRecord','EpochUTC','DeltaT_s_raw'});
result.P1_LowerEnergy_MeV_raw = double(r.FHDU_EnergyRange(:,1,10));
result.P1_UpperEnergy_MeV_raw = double(r.FHDU_EnergyRange(:,2,10));
for s = 1:8
    result.(sprintf('%s_S%d_raw',kind,s)) = double(value(:,10,s));
    result.(sprintf('%sSigma_S%d_raw',kind,s)) = double(sigma(:,10,s));
    result.(sprintf('SectoredQuality_S%d_raw',s)) = r.FHDU_SectoredQuality(:,10,s);
end
end

function writeCSV(t,file)
% Preserve UTC explicitly, expand vector-valued columns with sector suffixes.
names = t.Properties.VariableNames;
out = table;
for i = 1:numel(names)
    value = t.(names{i});
    if isdatetime(value)
        value.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''';
        value = string(value);
    end
    if (isnumeric(value) || islogical(value)) && size(value,2) > 1
        for s = 1:size(value,2)
            out.(sprintf('%s_S%d',names{i},s)) = value(:,s);
        end
    else
        out.(names{i}) = value;
    end
end
writetable(out,file,'Encoding','UTF-8');
end

function writeJSON(value,file)
fid = fopen(file,'w','n','UTF-8');
assert(fid >= 0);
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end

