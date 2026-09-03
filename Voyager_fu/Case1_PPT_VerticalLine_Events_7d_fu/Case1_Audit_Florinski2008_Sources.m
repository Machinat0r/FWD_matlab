function report = Case1_Audit_Florinski2008_Sources(DataRoot)
%Case1_Audit_Florinski2008_Sources Check only the Figure 4 four-day event.
%   Read original NASA/SPDF annual LECP and August COHO CDFs using the
%   current IRFU readers. Report raw Epoch, P1 sector values, DeltaT and
%   same-UTC-day/hour magnetic coverage for 2004 DOY 221 through 224.
%   No source is modified. Flux and uncertainty remain at original Epoch.
%   Only DeltaT < 0 is rejected; all other source rules follow current code.
%   Magnetic vector means reproduce the approved existing bin policy.
%   No nearest-field replacement, particle interpolation or new threshold.

%% directories and original source files
if nargin < 1, DataRoot = 'Z:/SPART-WORK/Data/Voyager'; end
IRFRoot = 'C:/Users/Administrator/Documents/irfu-matlab-master';
Case1_Add_IRFU_Path(IRFRoot);
StartUTC = datetime(2004,8,8,'TimeZone','UTC');
EndUTC = datetime(2004,8,12,'TimeZone','UTC');
RunID = char(datetime('now','TimeZone','UTC','Format','yyyyMMdd_HHmmss_SSS'));
AuditFolder = fullfile(DataRoot,'voyager1','lecp','validation', ...
    'florinski2008_figure4','2004');
if ~isfolder(AuditFolder), mkdir(AuditFolder); end
report = struct;
report.Scope = 'Voyager 1; 2004 DOY 221-224; [2004-08-08,2004-08-12) UTC';
report.StartUTC = StartUTC;
report.EndUTC = EndUTC;
report.RunID = RunID;
report.CodeFile = [mfilename('fullpath'),'.m'];
report.CodeSHA256 = Case1_File_SHA256(report.CodeFile);
report.SourceManifest = table;
report.Summary = table;
report.AuditFile = fullfile(AuditFolder,['source_audit_',RunID,'.mat']);
CohoFile = fullfile(DataRoot,'voyager1','coho','1hr','l2', ...
    'merged_mag_plasma','2004','08', ...
    'voyager1_coho1hr_merged_mag_plasma_20040801_v01.cdf');
[Mag,MagManifest] = Case1_Read_MAG_CDFs(string(CohoFile));
report.COHOManifest = MagManifest;
report.COHOURL = ['https://spdf.gsfc.nasa.gov/pub/data/voyager/', ...
    'voyager1/coho1hr_magplasma/2004/', ...
    'voyager1_coho1hr_merged_mag_plasma_20040801_v01.cdf'];
report.MagRecords = Mag(Mag.EpochUTC >= StartUTC & Mag.EpochUTC < EndUTC,:);
Coho = Voyager_Read_CDF_Product(CohoFile,'coho');
report.COHOFourDayCoverage = table;
for iday = 0:3
    DayStart = StartUTC+days(iday);
    OnDay = Coho.Epoch >= DayStart & Coho.Epoch < DayStart+days(1);
    B = [Coho.BR(OnDay),Coho.BT(OnDay),Coho.BN(OnDay)];
    Flux = Coho.protonFlux1_LECP(OnDay);
    Coverage = table(DayStart,nnz(OnDay), ...
        nnz(all(isfinite(B),2) & vecnorm(B,2,2)>0), ...
        nnz(isfinite(Flux) & Flux>0),'VariableNames', ...
        {'DayUTC','OriginalHourlyRecords','ValidBRecords','PositiveLECPP1Records'});
    report.COHOFourDayCoverage = [report.COHOFourDayCoverage;Coverage]; %#ok<AGROW>
end
report.COHOSourceStatus = 'Existing original CDF; not replaced in this task.';

%% read the two official products without averaging their particles
Cadences = ["day","hour"];
Folders = ["1d","1h"];
Products = ["daily","hourly"];
for ic = 1:2
    FileName = sprintf('voyager-1_lecp_lev-2-%s-avg_20040101_v1.1.1-01.cdf', ...
        Products(ic));
    SourceFile = fullfile(DataRoot,'voyager1','lecp',Folders(ic), ...
        'l2','sectored_flux','2004',FileName);
    SourceURL = sprintf(['https://spdf.gsfc.nasa.gov/pub/data/voyager/', ...
        'voyager1/particle/lecp/final-cdf/lev-2-%s-avg/%s'],Products(ic),FileName);
    P = Case1_Read_LECP_CDFs(SourceFile);
    % Pass only time arrays: the existing helper's optional provenance
    % branch does not handle an empty source-index vector consistently.
    TimeInput = struct('Epoch',P.Epoch,'DeltaT',P.DeltaT);
    [Keep,Selection] = Case1_Select_LECP_Epoch(TimeInput,StartUTC,EndUTC);
    NativeEpoch = spdfcdfread(char(SourceFile),'Variables',{'Epoch'}, ...
        'ConvertEpochToDatenum',true,'Validate',true);
    if iscell(NativeEpoch), NativeEpoch = NativeEpoch{1}; end
    NativeEpoch = datetime(NativeEpoch(:),'ConvertFrom','datenum','TimeZone','UTC');
    NativeInWindow = NativeEpoch >= StartUTC & NativeEpoch < EndUTC;
    assert(numel(NativeEpoch)==numel(P.Epoch), ...
        'Native NASA and current reader record counts differ.');
    Difference_s = abs(seconds(sort(NativeEpoch)-P.Epoch));
    assert(all(Difference_s < 0.001),'Independent Epoch conversion differs.');
    NativeAudit = struct('Reader','NASA spdfcdfread; native datenum conversion', ...
        'TotalRecords',numel(NativeEpoch),'WindowRecords',nnz(NativeInWindow), ...
        'EpochAgreementMaxSeconds',max(Difference_s), ...
        'FirstEpochUTC',min(NativeEpoch),'LastEpochUTC',max(NativeEpoch));
    [~,NearRows] = sort(abs(NativeEpoch-StartUTC));
    NearRows = sort(NearRows(1:min(12,numel(NearRows))));
    NativeAudit.NearestRecords = table(NearRows,NativeEpoch(NearRows), ...
        P.DeltaT(NearRows),'VariableNames',{'SourceCDFRecord','EpochUTC','DeltaT_s'});
    Before = find(NativeEpoch < StartUTC,1,'last');
    After = find(NativeEpoch >= EndUTC,1,'first');
    BracketRows = [Before;After];
    NativeAudit.BracketingRecords = table(BracketRows,NativeEpoch(BracketRows), ...
        P.DeltaT(BracketRows),'VariableNames',{'SourceCDFRecord','EpochUTC','DeltaT_s'});
    report.([char(Products(ic)),'NativeEpoch']) = NativeAudit;
    InWindow = P.Epoch >= StartUTC & P.Epoch < EndUTC;
    Rows = find(InWindow);
    N = numel(Rows);
    F = reshape(P.FHDU_SectoredFluxes(Rows,10,1:7),N,7);
    Sigma = reshape(P.FHDU_SectoredFluxUncertainties(Rows,10,1:7),N,7);
    BM = nan(N,3);
    Count = zeros(N,1);
    BinStart = dateshift(P.Epoch(Rows),'start',char(Cadences(ic)));
    for ii = 1:N
        if ic == 1, BinEnd = BinStart(ii)+days(1);
        else, BinEnd = BinStart(ii)+hours(1); end
        InBin = Mag.EpochUTC >= BinStart(ii) & Mag.EpochUTC < BinEnd;
        V = Mag{InBin,{'BR','BT','BN'}};
        V = V(all(isfinite(V),2),:);
        Count(ii) = size(V,1);
        if Count(ii) > 0, BM(ii,:) = mean(V,1); end
    end
    All7Positive = all(isfinite(F) & F > 0,2);
    BAvailable = all(isfinite(BM),2) & vecnorm(BM,2,2) > 0;
    T = table(P.Epoch(Rows),Rows,P.SourceRecordNumber(Rows), ...
        P.DeltaT(Rows),Keep(Rows),F,Sigma,All7Positive,BinStart,Count, ...
        BM,BAvailable,Keep(Rows) & All7Positive & BAvailable, ...
        'VariableNames',{'EpochUTC','SourceRow','SourceCDFRecord','DeltaT_s', ...
        'TimeRetained','P1FluxS1ToS7','P1SigmaS1ToS7','All7Positive', ...
        'MAGBinStartUTC','MAGVectorSampleCount','BMeanRTN_nT', ...
        'BAvailable','FluxAndBEligible'});
    T.SourceCDF = repmat(string(SourceFile),N,1);
    T.Properties.UserData = struct('Selection',Selection, ...
        'SourceURL',SourceURL,'ProductCadence',Cadences(ic), ...
        'P1EnergyRange',P.FHDU_EnergyRange(Rows,:,10), ...
        'SourceMetadata',{P.SourceMetadata});
    report.(char(Products(ic))) = T;
    Manifest = P.SourceManifest;
    Manifest.SourceURL = string(SourceURL);
    Manifest.FirstEpochUTC = min(P.Epoch);
    Manifest.LastEpochUTC = max(P.Epoch);
    FileInfo = dir(SourceFile);
    Manifest.Bytes = FileInfo.bytes;
    Manifest.SourceStatus = "Downloaded official raw CDF in this task; no re-encoding";
    report.SourceManifest = [report.SourceManifest;Manifest]; %#ok<AGROW>
    for iday = 0:3
        DayStart = StartUTC+days(iday);
        OnDay = T.EpochUTC >= DayStart & T.EpochUTC < DayStart+days(1);
        Summary = table(Cadences(ic),DayStart,nnz(OnDay), ...
            nnz(OnDay & T.TimeRetained),nnz(OnDay & T.TimeRetained & T.All7Positive), ...
            nnz(OnDay & T.TimeRetained & T.BAvailable), ...
            nnz(OnDay & T.FluxAndBEligible), ...
            'VariableNames',{'Cadence','DayUTC','OriginalRecords', ...
            'RetainedRecords','All7PositiveRecords','BAvailableRecords', ...
            'FluxAndBEligibleRecords'});
        report.Summary = [report.Summary;Summary]; %#ok<AGROW>
    end
end
report.Method = ['Current IRFU readers; CDF fill/valid-range rules; ', ...
    'retain original Epoch/flux/sigma; discard DeltaT<0; ', ...
    'same-UTC-day/hour finite complete vector mean; S1-S7 positive; ', ...
    'S8 omitted; no nearest field, no flux interpolation. ', ...
    'FluxAndBEligible does not independently certify attitude.'];
report.Passed = true;
save(report.AuditFile,'report');
disp(report.dailyNativeEpoch);
disp(report.hourlyNativeEpoch);
disp(report.Summary);
disp(report.daily(:,{'EpochUTC','DeltaT_s','All7Positive', ...
    'MAGVectorSampleCount','BAvailable','FluxAndBEligible'}));
fprintf('Florinski source audit: %s\n',report.AuditFile);
end
