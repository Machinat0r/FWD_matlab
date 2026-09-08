function [result,raw] = Run_V1_Interstellar_Daily_Overview(varargin)
% Four daily panels from original NASA CDFs; 2012-08-25 onward.
% Run_V1_Interstellar_Daily_Overview('SyncArchive',false) uses local sources.
% Author: Codex; 2026-09-07. See README_四面板日统计.md.

%% 路径与运行参数
root = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(root), 'Case1_PPT_VerticalLine_Events_7d_fu'));
cfg = Case1_Config;
Case1_Add_IRFU_Path(cfg.IRFURoot);
p = inputParser;
addParameter(p,'SyncArchive',true,@islogical);
addParameter(p,'Visible',true,@islogical);
addParameter(p,'MakePlot',true,@islogical);
addParameter(p,'DataRoot',cfg.DataRoot,@(x)ischar(x)||isstring(x));
addParameter(p,'OutputRoot','C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview',@(x)ischar(x)||isstring(x));
parse(p,varargin{:});
cfg.DataRoot = char(p.Results.DataRoot);
startUTC = datetime(2012,8,25,'TimeZone','UTC');
base = fullfile(cfg.DataRoot,'voyager1');
out = char(p.Results.OutputRoot);
if ~isfolder(out)
    mkdir(out);
end
if p.Results.SyncArchive
    archive = V1_Sync_Overview_Sources(cfg.DataRoot);
else
    archive = struct('VerifiedOnline',false);
    cached = fullfile(cfg.DataRoot,'source_verification','V1_daily_overview','official_archive_inventory.mat');
    if isfile(cached)
        a = load(cached,'audit'); % 仅目录核查信息，通量仍读取 CDF
        archive = a.audit;
        archive.FromCachedInventory = true;
    end
end
coho = latestFiles(fullfile(base,'coho','1hr','l2','merged_mag_plasma'), ...
    'voyager1_coho1hr_merged_mag_plasma_*_v*.cdf',2012);
assert(~isempty(coho),'No original COHO CDFs found.');

%% 读取原始 CDF：小时总磁场和 P1 通量
raw = table; manifest = table;
for k = 1:numel(coho)
    fprintf('COHO %d/%d: %s\n',k,numel(coho),coho(k));
    q = Voyager_Read_CDF_Product(coho(k),'coho');
    b = nan(numel(q.Epoch),1);
    field = "missing";
    if isfield(q,'ABS_B') && any(isfinite(q.ABS_B))
        b = q.ABS_B(:);
        field = "ABS_B";
    elseif isfield(q,'F') && any(isfinite(q.F))
        b = q.F(:);
        field = "F";
    end
    j = nan(size(b));
    if isfield(q,'protonFlux1_LECP')
        j = q.protonFlux1_LECP(:);
    end
    n = numel(b);
    part = table(q.Epoch(:),b,j,repmat(k,n,1),(1:n).', ...
        'VariableNames',{'EpochUTC','B_nT','P1','FileIndex','CDFRecord'});
    raw = [raw;part(part.EpochUTC>=startUTC,:)]; %#ok<AGROW>
    row = table(coho(k),string(Case1_File_SHA256(coho(k))),field, ...
        {q.variable_meta},{q.global_attributes}, ...
        'VariableNames',{'SourceFile','SHA256','MagnitudeVariable','VariableMetadata','GlobalAttributes'});
    manifest = [manifest;row]; %#ok<AGROW>
end
raw = sortrows(raw,'EpochUTC');
assert(numel(unique(raw.EpochUTC))==height(raw), ...
    'Overlapping original COHO Epochs require inspection; no implicit deduplication.');
assert(~isempty(raw),'No data after heliopause crossing.');
stopUTC = dateshift(max(raw.EpochUTC),'start','day')+caldays(1);
day = (startUTC:caldays(1):stopUTC-caldays(1)).';
nDay = numel(day);
[meanB,~,nB] = dailyStats(raw.EpochUTC,raw.B_nT,day);
[meanJ,medianJ,nJ] = dailyStats(raw.EpochUTC,raw.P1,day);
sectorMean = nan(nDay,8);
sectorCount = zeros(nDay,8);
sectorSource = repmat("missing",nDay,1);
sectorAuditFiles = strings(0,1);

%% 读取原始 CDF：按既有 L1 优先规则计算扇区日均
l1 = latestFiles(fullfile(base,'lecp','native','l1','sectored_rates'), ...
    'voyager-1_lecp_lev-1-rates_*_v*.cdf',2012);
l2 = latestFiles(fullfile(base,'lecp','1d','l2','sectored_flux'), ...
    'voyager-1_lecp_lev-2-daily-avg_*_v*.cdf',2012);
sourceYears = unique([fileYears(l1);fileYears(l2)]);
for yy = sourceYears.'
    f1 = l1(fileYears(l1)==yy);
    f2 = l2(fileYears(l2)==yy);
    assert(numel(f1)==1 && numel(f2)==1, ...
        'Incomplete L1/L2 source pair for %d; synchronize official archive.',yy);
    fprintf('LECP daily sectors: %d\n',yy);
    rates = Case1_Read_LECP_Rates(f1);
    flux = Case1_Read_LECP_CDFs(f2);
    begin = max(startUTC,datetime(yy,1,1,'TimeZone','UTC'));
    finish = min(stopUTC,datetime(yy+1,1,1,'TimeZone','UTC'));
    if finish<=begin
        continue
    end
    merged = Case1_Apply_L1_Fallback(flux,rates,begin,finish,'day','l1_first');
    values = reshape(merged.FHDU_SectoredFluxes(:,10,:),[],8);
    rows = find(day>=begin & day<finish);
    for ss = 1:8
        [sectorMean(rows,ss),~,sectorCount(rows,ss)] = ...
            dailyStats(merged.Epoch,values(:,ss),day(rows));
    end
    for ii = rows.'
        use = merged.Epoch>=day(ii) & merged.Epoch<day(ii)+caldays(1);
        if any(use)
            sectorSource(ii) = join(unique(merged.SourceProduct(use)),';');
        end
    end
    auditFile = fullfile(out,sprintf('sector_source_audit_%d.mat',yy));
    save(auditFile,'merged','f1','f2','-v7.3');
    sectorAuditFiles(end+1,1) = string(auditFile); %#ok<AGROW>
end
selectedSectors = [1 2 3 5 6 7];
sumSix = sum(sectorMean(:,selectedSectors),2); % NaN propagates: all six required.
validSix = all(isfinite(sectorMean(:,selectedSectors)),2);
sumSix(~validSix) = NaN;
daily = table(day+hours(12),meanB,meanJ,medianJ,sumSix,nB,nJ, ...
    sectorMean,sectorCount,sectorSource,validSix, ...
    'VariableNames',{'EpochUTC','BMean_nT','P1Mean','P1Median','P1SixSectorSum', ...
    'MAGSampleCount','P1SampleCount','SectorDailyMean','SectorRecordCount', ...
    'SectorSource','SixSectorsComplete'});

%% 保存统计、源文件信息和方法记录
method = struct;
method.StartUTC = startUTC;
method.EndUTCExclusive = stopUTC;
method.TimeBins = '[UTC 00:00, next UTC 00:00); points at 12:00 UTC';
method.PanelA = 'Arithmetic mean of original COHO hourly scalar ABS_B (F fallback).';
method.PanelBC = 'Arithmetic mean and median of the same finite original COHO protonFlux1_LECP hourly values; no log transform before aggregation.';
method.PanelD = 'Sum of daily S1,S2,S3,S5,S6,S7 differential fluxes. L1-first seven-sector-complete source selection reused unchanged. Multiple retained L2 rows in one UTC day are averaged per sector. No PA bins, solid-angle weighting or normalization.';
method.Missing = 'CDF fill/valid-range screening inherited from IRFU reader; finite-only daily statistics, zero retained if source-valid; no minimum-count threshold, interpolation, gap filling or background subtraction. Missing required sector leaves sum NaN.';
method.L1 = 'Drop negative DeltaT, anchor at original Epoch; per-sector nonnegative rates arithmetic mean / (0.44*(1.78-0.57)). Historical conversion not established equivalent to calibrated L2. Complete S1-S7 day required for L1 substitution; original L2 otherwise.';
method.Energy = 'Display 0.57-1.78 MeV; preserve source P1 EnergyRange (previously 0.57-0.89 MeV) and all source metadata unchanged.';
method.Attitude = 'Sector-flux sum requires no pitch-angle calculation; no new attitude approximation or PAD usability mask applied.';
method.Uncertainty = 'No error bars requested. Original sector sigma and approved L1 propagated sigma retained in annual merged audits.';
method.BoundarySource = 'https://www.nasa.gov/news-release/nasa-spacecraft-embarks-on-historic-journey-into-interstellar-space/';
coverage = table(["a";"b";"c";"d"],NaT(4,1,'TimeZone','UTC'), ...
    NaT(4,1,'TimeZone','UTC'),zeros(4,1), ...
    'VariableNames',{'Panel','FirstValidUTC','LastValidUTC','ValidDays'});
v = [meanB meanJ medianJ sumSix];
for k=1:4
    good = isfinite(v(:,k));
    coverage.ValidDays(k) = nnz(good);
    if any(good)
        coverage.FirstValidUTC(k)=day(find(good,1));
        coverage.LastValidUTC(k)=day(find(good,1,'last'));
    end
end
codeFiles = string(fullfile(root,{'Run_V1_Interstellar_Daily_Overview.m', ...
    'V1_Sync_Overview_Sources.m','Plot_V1_Interstellar_Daily_Overview.m', ...
    'Overview_V1_Interstellar_Daily.m'})).';
codeHash = strings(size(codeFiles));
for k = 1:numel(codeFiles)
    codeHash(k) = Case1_File_SHA256(codeFiles(k));
end
result = struct('Daily',daily,'Coverage',coverage,'Method',method,'Archive',archive, ...
    'COHOSources',manifest,'SectorAuditFiles',sectorAuditFiles,'CodeFiles',codeFiles, ...
    'CodeSHA256',codeHash,'OutputFolder',string(out));

writetable(daily,fullfile(out,'V1_daily_overview.csv'));
writetable(coverage,fullfile(out,'V1_daily_coverage.csv'));

%% 四面板绘图
if p.Results.MakePlot
    result.DisplayAudit = Plot_V1_Interstellar_Daily_Overview(result,p.Results.Visible);
else
    result.DisplayAudit = struct('FigureExported',false);
end
save(fullfile(out,'V1_daily_overview_audit.mat'),'result','raw','-v7.3');
disp(coverage);
fprintf('Saved all outputs: %s\n',out);
end

function [avg,med,count]=dailyStats(time,value,day)
avg = nan(numel(day),1);
med = avg;
count = zeros(numel(day),1);
for k=1:numel(day)
    x=value(time>=day(k) & time<day(k)+caldays(1) & isfinite(value));
    count(k)=numel(x);
    if ~isempty(x)
        avg(k) = mean(x);
        med(k) = median(x);
    end
end
end

function yy=fileYears(files)
yy=zeros(numel(files),1);
for k=1:numel(files)
    t=regexp(char(files(k)),'_(\d{4})\d{4}_v','tokens','once');
    yy(k)=str2double(t{1});
end
end

function files=latestFiles(folder,pattern,minimumYear)
e=dir(fullfile(folder,'**',pattern));
files=sort(string(fullfile({e.folder},{e.name})).');
if isempty(files)
    return
end
files=files(fileYears(files)>=minimumYear);
keys=strings(size(files));
for k=1:numel(files)
    t=regexp(char(files(k)),'_(\d{8})_v','tokens','once'); keys(k)=string(t{1});
end
[~,keep] = unique(keys,'last');
files = files(keep);
end





