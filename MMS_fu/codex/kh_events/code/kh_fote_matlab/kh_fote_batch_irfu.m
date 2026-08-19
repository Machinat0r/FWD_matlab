function batchReport = kh_fote_batch_irfu(varargin)
%KH_FOTE_BATCH_IRFU 批处理模块：逐行运行KH事件表并合并PDF。
%
% 模块职责：读取EventID/PlotStartUTC/PlotEndUTC，逐事件调用统一单事件
% 入口，记录成功或失败状态，并持续写出batch_status.csv。
%
% Example:
%   report = kh_fote_batch_irfu('SmoothSeconds',5, ...
%       'EigenvalueStabilityThreshold',0.5, ...
%       'MinQualityDurationSeconds',5);

p = inputParser;
p.addParameter('Catalog', ...
    'C:\Users\Administrator\Documents\KH\selected_burst_windows_all.csv', ...
    @(x)ischar(x) || isstring(x));
p.addParameter('OutputRoot',fullfile( ...
    'C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\codex\kh_events', ...
    'all_events','matlab_irfu_q40_eig0p5_run5s'), ...
    @(x)ischar(x) || isstring(x));
p.addParameter('SmoothSeconds',5,@(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('QualityPercent',40,@(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('EigenvalueStabilityThreshold',0.5, ...
    @(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('MinQualityDurationSeconds',5, ...
    @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('VelocityField','Vi',@(x)any(strcmpi(string(x),["Vi","Ve"])));
p.addParameter('FigureVisible','off',@(x)any(strcmpi(string(x),["on","off"])));
p.addParameter('ContinueOnError',true,@(x)islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
catalogPath = char(string(opt.Catalog));
outputRoot = char(string(opt.OutputRoot));
if ~isfile(catalogPath), error('Catalog not found: %s',catalogPath); end
if ~exist(outputRoot,'dir'), mkdir(outputRoot); end

catalog = readtable(catalogPath,'TextType','string');
required = ["EventID","PlotStartUTC","PlotEndUTC"];
if ~all(ismember(required,string(catalog.Properties.VariableNames)))
    error('Catalog must contain EventID, PlotStartUTC, and PlotEndUTC.');
end

pdfDir = fullfile(outputRoot,'pdf');
dataDir = fullfile(outputRoot,'data');
if ~exist(pdfDir,'dir'), mkdir(pdfDir); end
if ~exist(dataDir,'dir'), mkdir(dataDir); end
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
combinedPdf = fullfile(pdfDir,sprintf( ...
    'KH_all_%s_sm%gs_q%g_eig%g_run%gs_%s.pdf', ...
    upper(char(string(opt.VelocityField))),opt.SmoothSeconds, ...
    opt.QualityPercent,opt.EigenvalueStabilityThreshold, ...
    opt.MinQualityDurationSeconds,stamp));

n = height(catalog);
eventID = strings(n,1); startUtc = strings(n,1); endUtc = strings(n,1);
success = false(n,1); bGood = nan(n,1); vGood = nan(n,1);
bContinuous = nan(n,1); vContinuous = nan(n,1);
figurePdf = strings(n,1);
message = strings(n,1);
for k = 1:n
    eventID(k) = string(catalog.EventID(k));
    startUtc(k) = normalizeUtc(catalog.PlotStartUTC(k));
    endUtc(k) = normalizeUtc(catalog.PlotEndUTC(k));
    fprintf('BATCH %d/%d %s\n',k,n,eventID(k));
    try
        one = kh_fote_event(char(eventID(k)),char(startUtc(k)),char(endUtc(k)), ...
            'SmoothSeconds',opt.SmoothSeconds, ...
            'QualityPercent',opt.QualityPercent, ...
            'EigenvalueStabilityThreshold', ...
            opt.EigenvalueStabilityThreshold, ...
            'MinQualityDurationSeconds',opt.MinQualityDurationSeconds, ...
            'VelocityField',char(string(opt.VelocityField)), ...
            'FigureVisible',char(string(opt.FigureVisible)), ...
            'OutputRoot',fullfile(outputRoot,char(eventID(k))), ...
            'CombinedPdf',combinedPdf);
        success(k) = true;
        bGood(k) = one.Summary.BQuality40;
        vGood(k) = one.Summary.VQuality40;
        bContinuous(k) = one.Summary.BQuality40Continuous;
        vContinuous(k) = one.Summary.VQuality40Continuous;
        figurePdf(k) = string(one.FigurePdf);
        message(k) = "OK";
    catch ME
        message(k) = string(ME.message);
        warning('KH batch: %s failed: %s',eventID(k),ME.message);
        if ~opt.ContinueOnError, rethrow(ME); end
    end
    status = table(eventID,startUtc,endUtc,success,bGood,vGood, ...
        bContinuous,vContinuous, ...
        figurePdf,message,'VariableNames',{'EventID','StartUTC','EndUTC', ...
        'Success','BQuality40','VQuality40','BQuality40Continuous', ...
        'VQuality40Continuous', ...
        'FigurePdf','Message'});
    writetable(status(1:k,:),fullfile(dataDir,'batch_status.csv'));
end

batchReport = struct('Catalog',catalogPath,'OutputRoot',outputRoot, ...
    'CombinedPdf',combinedPdf,'TotalEvents',n,'Succeeded',nnz(success), ...
    'Failed',nnz(~success),'StatusCsv',fullfile(dataDir,'batch_status.csv'));
save(fullfile(dataDir,'batch_report.mat'),'batchReport','status','opt');
disp(batchReport);
end


function out = normalizeUtc(value)
out = strip(string(value));
out = replace(out,' ','T');
if ~endsWith(out,'Z'), out = out+"Z"; end
end
