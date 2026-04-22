%------ MP crossing search (SDCFilenames-based) --------------------------
% 条件保持不变：
% 2) 读位置：若 x>0 或 R>8Re 才继续
% 3) 磁层：Bz 连续>=3s >20nT 且同时 Ne<2 cm^-3
%    磁鞘：|Vi| 连续>=3s >100 km/s 且同时 Ne>8 cm^-3
% 4) 同一磁场文件时间段内同时出现磁层与磁鞘 -> 磁层顶穿越；记录时间并 SDCPlot 画整段磁场文件

clear; clc;
global ParentDir OutputDir
ParentDir   = 'Z:\SPART-WORK\Data\MMS/';
DownloadDir = 'C:\MMS/';
TempRoot    = [DownloadDir,'temp/']; if ~isfolder(TempRoot), mkdir(TempRoot); end

OutputDir = [ParentDir,'MPBoundarySearch_20190101_20191231/'];
if ~isfolder(OutputDir), mkdir(OutputDir); end
if ~isfolder([OutputDir,'OverviewFig/']), mkdir([OutputDir,'OverviewFig/']); end

units = irf_units;

% -------------------- 参数区（按需修改） --------------------------------
% 检索时间范围（包含 2023-12-31 当天）
DateStart = '2019-06-09';
DateEnd   = '2020-01-01';
DateEndPlus1 = datestr(datenum(DateEnd,'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [DateStart,'/',DateEndPlus1];

% 卫星设置
ic  = 1;
iic = 1; 

% 判据阈值
MIN_DUR_S  = 3.0;
BZ_TH_NT   = 20.0;
NE_MS_MAX  = 1.0;
VI_TH_KMS  = 100.0;
NE_SH_MIN  = 3.0;
R_MIN_RE   = 6.0;
X_MIN_RE   = 0.0;

% 用于定义“抓取一个 brst 文件的时间窗口”（通常 brst 文件 < 10min）
WIN_SEC_FOR_ONE_FILE = 12*60;  % 12min，保证覆盖文件

% 下载参数
DL_THREADS = 16;
DL_CHECKSIZE = 0;

caseLog  = fullfile(OutputDir,'MP_caselist.txt');
errorLog = fullfile(OutputDir,'MP_errorlog.txt');

% -------------------- 1) 用 SDCFilenames 获取文件名清单 -----------------
% 注：这里先“取清单”，不等于全部下载；后面会按每个磁场文件逐个过滤并下载
try
    % FGM：用于定义“逐个磁场文件”的循环对象（只用 ic）
    filenames_fgm = SDCFilenames(Date, ic,  'inst','fgm','drm','brst');

    % FPI moments：判据需要 Ne 与 Vi（只用 ic）
    filenames_fpi = SDCFilenames(Date, ic,  'inst','fpi','drm','brst','dpt','des-moms,dis-moms');

    % 画图可能用到：SCM/EDP（按你旧脚本习惯）
    filenames_scm = SDCFilenames(Date, ic,  'inst','scm','drm','brst','dpt','scb');
    filenames_edp = SDCFilenames(Date, iic, 'inst','edp','drm','brst','dpt','dce');

    % （可选但建议）位置数据可能依赖 MEC；如果你本地已有可注释掉
    % 某些环境下 mms.get_data('R_gsm',...) 会自动找 mec srvy
    % filenames_mec = SDCFilenames(Date, ic, 'inst','mec','drm','srvy','dpt','epht89d');

    filenames_all = [filenames_fgm, filenames_fpi, filenames_scm, filenames_edp];
catch ME
    error('SDCFilenames 获取清单失败：%s', ME.message);
end

if isempty(filenames_fgm)
    error('在 %s 范围内未找到 fgm brst 文件名。', Date);
end

% 从 FGM 文件名中解析每个文件的起始时刻（14位时间戳），作为循环对象
tStr14_list = regexp(filenames_fgm, '_(\d{14})_v', 'tokens');
tStr14_list = cellfun(@(c) c{1}{1}, tStr14_list, 'UniformOutput', false);
tStr14_list = unique(tStr14_list, 'stable');

fprintf('FGM brst 文件数（ic=%d）：%d\n', ic, numel(tStr14_list));

% 初始化本地数据库（下载/移动后会重复 init）
mms.db_init('local_file_db', ParentDir);

% -------------------- 主循环：逐个磁场文件检索 --------------------------
for k = 1:numel(tStr14_list)

    tStr14 = tStr14_list{k};
    isoStart = tstr14_to_iso8601z(tStr14);
    isoEnd   = add_seconds_iso8601z(isoStart, WIN_SEC_FOR_ONE_FILE);
    TT = [isoStart,'/',isoEnd];

    clc; fprintf('(%d/%d) 当前检索文件起点：%s  (sc=%d)\n', k, numel(tStr14_list), isoStart, ic);

    % 过滤出与该 TT 相交的文件名，并拿到 des-moms 文件（供 SDCPlot 或你内部使用）
    try
        [fn_sel, desmoms1, ~] = findFilenames(TT, filenames_all, 'brst', ic);
    catch ME
        append_error(errorLog, sprintf('%s  findFilenames失败：%s', isoStart, ME.message));
        continue
    end

    if isempty(fn_sel)
        continue
    end

    % 为该次检索创建独立临时目录
    tempDir = TempRoot;
    if ~isfolder(tempDir), mkdir(tempDir); end

    % 先下载“最小集合”用于判据：fgm + des-moms + dis-moms + (可选)mec
    fn_min = fn_sel( contains(fn_sel,'fgm') | contains(fn_sel,'des-moms') | contains(fn_sel,'dis-moms') | contains(fn_sel,'mec') );

    try
        if ~isempty(fn_min)
            SDCFilesDownload_NAS(fn_min, tempDir, 'CheckSize', DL_CHECKSIZE, 'Threads', DL_THREADS);
            SDCDataMove(tempDir, ParentDir);
            mms.db_init('local_file_db', ParentDir);
        end
    catch
        % 下载失败不一定致命（可能本地已有），继续尝试读数据
        warning('本轮最小集合下载/移动失败或无文件下载，继续尝试本地读取。');
    end

    % -------------- 读取磁场，确定“该磁场文件真实时间段” -----------------
    try
        tint0 = irf.tint(TT);
        B_ts = mms.get_data('B_gsm_brst', tint0, ic);
        if isempty(B_ts)
            cleanup_tempdir(tempDir, errorLog, tStr14);
            continue
        end
        Bmat = irf.ts2mat(B_ts); % [epoch, Bx, By, Bz]
        tint = irf.tint(B_ts.time.epoch(1), B_ts.time.epoch(end)); % 真实文件时间段
        tEpoch0 = Bmat(1,1); tEpoch1 = Bmat(end,1);
    catch ME
        append_error(errorLog, sprintf('%s  磁场导入失败：%s', isoStart, ME.message));
        cleanup_tempdir(tempDir, errorLog, tStr14);
        continue
    end

    % -------------------- 2) 读取位置并筛选 ------------------------------
    try
        Pos = mms.get_data('R_gsm', tint);
        posField = ['gsmR', num2str(ic)];
        if ~isfield(Pos, posField)
            cleanup_tempdir(tempDir, errorLog, tStr14);
            continue
        end
        Posmat = Pos.(posField); % [epoch, x, y, z] km
        if isempty(Posmat)
            cleanup_tempdir(tempDir, errorLog, tStr14);
            continue
        end

        if size(Posmat, 2) == 3
            Posmat = irf_abs(Posmat);
        end

        Re_km = units.RE/1e3;
        x_Re = Posmat(1,2)/Re_km;
        R_Re = sqrt(Posmat(1,2)^2 + Posmat(1,3)^2 + Posmat(1,4)^2)/Re_km;
        if ~(x_Re > X_MIN_RE || R_Re > R_MIN_RE)
            cleanup_tempdir(tempDir, errorLog, tStr14);
            continue
        end
    catch ME
        append_error(errorLog, sprintf('%s  位置导入失败：%s', isoStart, ME.message));
        cleanup_tempdir(tempDir, errorLog, tStr14);
        continue
    end

    % -------------------- 3) 判据：磁层/磁鞘（保持不变） ------------------
    try
        Ne_ts = mms.db_get_ts( ...
            sprintf('mms%d_fpi_brst_l2_des-moms', ic), ...
            sprintf('mms%d_des_numberdensity_brst', ic), ...
            tint);
        Vi_ts = mms.db_get_ts( ...
            sprintf('mms%d_fpi_brst_l2_dis-moms', ic), ...
            sprintf('mms%d_dis_bulkv_gse_brst', ic), ...
            tint);

        if isempty(Ne_ts) || isempty(Vi_ts)
            cleanup_tempdir(tempDir, errorLog, tStr14);
            continue
        end

        Ne = irf.ts2mat(Ne_ts);          % [epoch, Ne] cm^-3
        ViAbs = irf.ts2mat(Vi_ts.abs);   % [epoch, |Vi|] km/s

        % “与此同时”：重采样到 Ne 时间轴
        Bmat = irf_abs(Bmat);
        Babs = [Bmat(:,1), Bmat(:,5)];
        Babs_res = irf_resamp(Babs, Ne);
        Vi_res = irf_resamp(ViAbs, Ne);

        tNe = Ne(:,1);
        bz  = Babs_res(:,2);
        ne  = Ne(:,2);
        vi  = Vi_res(:,2);

        condMS = (bz > BZ_TH_NT)  & (ne < NE_MS_MAX);
        condSH = (vi > VI_TH_KMS) & (ne > NE_SH_MIN);

        intMS = find_intervals(tNe, condMS, MIN_DUR_S);
        intSH = find_intervals(tNe, condSH, MIN_DUR_S);

        if isempty(intMS) || isempty(intSH)
            cleanup_tempdir(tempDir, errorLog, tStr14);
            continue
        end

        [crossEpoch, pairInfo] = pick_best_crossing(intMS, intSH);
        crossUTC = irf_time(crossEpoch,'epoch>utc');

    catch ME
        append_error(errorLog, sprintf('%s  判据计算失败：%s', isoStart, ME.message));
        cleanup_tempdir(tempDir, errorLog, tStr14);
        continue
    end

    % -------------------- 4) 命中：下载全量并出图（保持不变） -------------
    try
        % 记录事件
        append_case(caseLog, ic, tStr14, x_Re, R_Re, crossUTC, pairInfo, intMS, intSH);

        % 命中后下载该 TT 内的“全量文件”（你示例的做法）
        try
            SDCFilesDownload_NAS(fn_sel, tempDir, 'CheckSize', DL_CHECKSIZE, 'Threads', DL_THREADS);
            SDCDataMove(tempDir, ParentDir);
            mms.db_init('local_file_db', ParentDir);
        catch
            warning('命中后全量下载失败或部分缺失，将尝试直接出图。');
        end

        % 出图：整段磁场文件时间段
        PlotTint = irf_time([tEpoch0, tEpoch1], 'epoch>epochTT');
        NameTag = [tStr14]; % 与旧脚本风格一致
        SDCPlot(PlotTint, desmoms1, desmoms1, ic, NameTag, crossEpoch);

    catch ME
        append_error(errorLog, sprintf('%s  记录/出图失败：%s', crossUTC, ME.message));
    end

    % 清理临时目录
    % cleanup_tempdir(tempDir, errorLog, tStr14);

end

fprintf('\n完成：结果写入 %s\n', caseLog);

% ======================= local functions ================================

function iso = tstr14_to_iso8601z(tStr14)
% 'yyyymmddHHMMSS' -> 'YYYY-mm-ddTHH:MM:SS.000Z'
dn = datenum(tStr14, 'yyyymmddHHMMSS');
iso = [datestr(dn,'yyyy-mm-dd'), 'T', datestr(dn,'HH:MM:SS'), '.000Z'];
end

function iso2 = add_seconds_iso8601z(iso1, sec)
% iso1: 'YYYY-mm-ddTHH:MM:SS.000Z'
dn = datenum(iso1(1:19), 'yyyy-mm-ddTHH:MM:SS');
dn2 = dn + sec/86400;
iso2 = [datestr(dn2,'yyyy-mm-dd'), 'T', datestr(dn2,'HH:MM:SS'), '.000Z'];
end

function intervals = find_intervals(tEpoch, cond, minDurSec)
tEpoch = tEpoch(:); cond = cond(:);
good = ~isnan(tEpoch);
tEpoch = tEpoch(good); cond = cond(good);
cond(~isfinite(cond)) = false;

if numel(tEpoch) < 2
    intervals = [];
    return
end

dt = diff(tEpoch);
dtMed = median(dt(isfinite(dt) & dt>0));
if isempty(dtMed) || ~isfinite(dtMed), dtMed = 0; end

edges = diff([false; cond; false]);
iStart = find(edges == 1);
iEnd   = find(edges == -1) - 1;

intervals = [];
for i = 1:numel(iStart)
    t0 = tEpoch(iStart(i));
    t1 = tEpoch(iEnd(i)) + dtMed;
    if (t1 - t0) >= minDurSec
        intervals(end+1,:) = [t0, t1]; %#ok<AGROW>
    end
end
end

function [crossEpoch, info] = pick_best_crossing(intMS, intSH)
bestGap = inf;
crossEpoch = NaN;
info = [];

for i = 1:size(intMS,1)
    s1 = intMS(i,1); e1 = intMS(i,2);
    for j = 1:size(intSH,1)
        s2 = intSH(j,1); e2 = intSH(j,2);

        if e1 <= s2
            gap = s2 - e1; cross = (e1 + s2)/2; order = 1; % MS->SH
        elseif e2 <= s1
            gap = s1 - e2; cross = (e2 + s1)/2; order = 2; % SH->MS
        else
            gap = 0; cross = max(s1,s2); order = 0;        % overlap
        end

        if gap < bestGap
            bestGap = gap;
            crossEpoch = cross;
            info = [i, j, order, gap];
        end
    end
end
end

function append_case(caseLog, ic, tStr14, x_Re, R_Re, crossUTC, pairInfo, intMS, intSH)
fid = fopen(caseLog, 'a');
if fid < 0, return; end
fprintf(fid, 'UTC_cross=%s | sc=%d | fileTag=%s | x=%.2fRe | R=%.2fRe | nMS=%d | nSH=%d\n', ...
    crossUTC, ic, tStr14, x_Re, R_Re, size(intMS,1), size(intSH,1));
if ~isempty(pairInfo)
    i = pairInfo(1); j = pairInfo(2); order = pairInfo(3); gap = pairInfo(4);
    fprintf(fid, '   MS[%d]=[%s, %s], SH[%d]=[%s, %s], order=%d, gap=%.3fs\n', ...
        i, irf_time(intMS(i,1),'epoch>utc'), irf_time(intMS(i,2),'epoch>utc'), ...
        j, irf_time(intSH(j,1),'epoch>utc'), irf_time(intSH(j,2),'epoch>utc'), ...
        order, gap);
end
fclose(fid);
end

function append_error(errorLog, msg)
fid = fopen(errorLog, 'a');
if fid < 0, return; end
fprintf(fid, '%s\n', msg);
fclose(fid);
end

function cleanup_tempdir(tempDir, errorLog, tag)
try
    if isfolder(tempDir)
        rmdir(tempDir, 's');
    end
catch
    % append_error(errorLog, sprintf('删除临时目录失败：%s  (%s)', tempDir, tag));
end
end
