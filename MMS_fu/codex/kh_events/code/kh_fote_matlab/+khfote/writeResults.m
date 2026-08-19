function result = writeResults(fig,rawAligned,aligned,mag,vel, ...
    fluxField,agreement,summary,smoothingInfo,opt)
%KH FOTE 输出模块：统一保存图、时间序列、MAT数据和摘要。
%
% 输出目录结构
%   pdf/              单事件矢量PDF
%   png/              与PDF同布局的预览图
%   data/timeseries/  每个共同时间点的零点、误差、连续掩码和一致性
%   data/mat/         完整MATLAB结构体，便于后续统计
%   data/summaries/   JSON运行摘要
%
% 距离字段写入CSV/MAT作为诊断信息，不用于good筛选。

%% 1. 建立分类输出目录与稳定文件名
pdfDir = fullfile(opt.OutputRoot,'pdf');
pngDir = fullfile(opt.OutputRoot,'png');
csvDir = fullfile(opt.OutputRoot,'data','timeseries');
matDir = fullfile(opt.OutputRoot,'data','mat');
summaryDir = fullfile(opt.OutputRoot,'data','summaries');
folders = {pdfDir,pngDir,csvDir,matDir,summaryDir};
for k = 1:numel(folders)
    if ~exist(folders{k},'dir'), mkdir(folders{k}); end
end

smoothTag = numberTag(opt.SmoothSeconds);
qualityTag = numberTag(opt.QualityPercent);
eigenTag = numberTag(opt.EigenvalueStabilityThreshold);
if opt.MinQualityDurationSeconds>0
    durationTag = ['run' numberTag(opt.MinQualityDurationSeconds) 's'];
else
    durationTag = 'pointwise';
end
% 使用短文件名，避免Windows传统260字符路径限制。
base = sprintf('%s_%s_%s_%s_sm%ss_q%s_eig%s_%s', ...
    safeName(opt.eventID),compactUtc(opt.startUtc), ...
    compactUtc(opt.endUtc),upper(opt.VelocityField),smoothTag,qualityTag, ...
    eigenTag,durationTag);
pdfPath = fullfile(pdfDir,[base '.pdf']);
pngPath = fullfile(pngDir,[base '.png']);
csvPath = fullfile(csvDir,[base '_timeseries.csv']);
matPath = fullfile(matDir,[base '.mat']);
jsonPath = fullfile(summaryDir,[base '_summary.json']);

%% 2. 导出图像
exportgraphics(fig,pdfPath,'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,pngPath,'Resolution',140,'BackgroundColor','white');
if ~isempty(opt.CombinedPdf)
    appendFlag = isfile(opt.CombinedPdf);
    exportgraphics(fig,opt.CombinedPdf,'ContentType','vector', ...
        'BackgroundColor','white','Append',appendFlag);
end
if strcmpi(opt.FigureVisible,'off'), close(fig); end

%% 3. 写出逐时刻时间序列表
timeUtc = datetime(aligned.t,'ConvertFrom','posixtime','TimeZone','UTC');
timeUtc.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
outTable = table(timeUtc,mag.distanceKm,mag.distanceL,string(mag.type), ...
    mag.eta,mag.xi,mag.good,mag.eigenvalueStability, ...
    mag.eigenvalueStable,mag.spiralQualified,mag.sustainedGood,mag.near, ...
    vel.distanceKm,vel.distanceL,string(vel.type),vel.alpha,vel.good, ...
    vel.eigenvalueStability,vel.eigenvalueStable,vel.spiralQualified, ...
    vel.sustainedGood,vel.near,agreement,'VariableNames', ...
    {'TimeUTC','B_NullDistance_km','B_NullDistance_L','B_Type', ...
    'B_Eta_percent','B_Xi_percent','B_Good40','B_EigenvalueStability', ...
    'B_EigenvalueStable','B_StableSpiral','B_AcceptedSpiral','B_Near', ...
    'V_NullDistance_km','V_NullDistance_L','V_Type','V_Alpha_percent', ...
    'V_Good40','V_EigenvalueStability','V_EigenvalueStable', ...
    'V_StableSpiral','V_AcceptedSpiral','V_Near','ScrewAgreement'});
writetable(outTable,csvPath);

%% 4. 保存完整MAT数据
save(matPath,'rawAligned','aligned','mag','vel','fluxField', ...
    'agreement','summary','smoothingInfo','opt','-v7.3');

%% 5. 生成机器可读运行摘要
result = struct('EventID',opt.eventID,'StartUTC',opt.startUtc, ...
    'EndUTC',opt.endUtc, ...
    'Engine','FOTE_Taylor_Expansion + IRFU c_4_grad', ...
    'SmoothSeconds',opt.SmoothSeconds, ...
    'EigenvalueStabilityThreshold',opt.EigenvalueStabilityThreshold, ...
    'MinQualityDurationSeconds',opt.MinQualityDurationSeconds, ...
    'RawCommonSamples',numel(rawAligned.t), ...
    'CommonSamples',numel(aligned.t), ...
    'FigurePdf',pdfPath,'FigurePng',pngPath, ...
    'TimeSeries',csvPath,'MatFile',matPath,'Summary',summary);
fid = fopen(jsonPath,'w');
if fid>=0
    cleaner = onCleanup(@()fclose(fid));
    fwrite(fid,jsonencode(result,'PrettyPrint',true),'char');
    clear cleaner
end

fprintf('MATLAB/IRFU result: %s\n',pdfPath);
fprintf(['B quality40=%d, V quality40=%d; M_lambda >=%g: B=%d, V=%d; ' ...
    'accepted (run >=%g s): B=%d, V=%d.\n'], ...
    summary.BQuality40,summary.VQuality40, ...
    opt.EigenvalueStabilityThreshold,summary.BStableSpiral, ...
    summary.VStableSpiral,opt.MinQualityDurationSeconds, ...
    summary.BAcceptedSpiral,summary.VAcceptedSpiral);
end


function out = safeName(in)
out = regexprep(char(string(in)),'[^A-Za-z0-9_-]','_');
end


function out = compactUtc(in)
out = regexprep(char(string(in)),'[^0-9]','');
if numel(out)>14, out = out(1:14); end
end


function out = numberTag(in)
out = strrep(sprintf('%g',in),'.','p');
end
