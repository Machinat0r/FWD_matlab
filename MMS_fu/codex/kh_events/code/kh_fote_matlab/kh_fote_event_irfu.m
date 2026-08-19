function result = kh_fote_event_irfu(eventID,startUtc,endUtc,varargin)
%KH_FOTE_EVENT_IRFU MMS KH事件的MATLAB/IRFU FOTE与FOTE-V主入口。
%
% 本文件只负责组织各模块的执行顺序。具体算法位于+khfote package：
%   parseOptions                 参数与默认路径
%   setupEnvironment             IRFU/FOTE环境及MMS数据库
%   loadMmsData                  四星原生数据读取
%   prepareData                  物理时间平均、FOTE前irf_resamp和缺口保留
%   runOriginalFote              原FOTE函数及零点类型
%   computeFoteVError            nV连续性误差alpha
%   applyQualityAndSummarize     40%、特征值稳定度及持续时间筛选
%   plotEvent                    7-panel绘图
%   writeResults                 PDF/PNG/CSV/MAT/JSON输出
%
% 示例
%   result = kh_fote_event_irfu('KH005', ...
%       '2015-10-01T18:01:24Z','2015-10-01T18:09:00Z', ...
%       'SmoothSeconds',5,'QualityPercent',40, ...
%       'EigenvalueStabilityThreshold',0.5, ...
%       'MinQualityDurationSeconds',5);

%% 模块1：参数配置与运行环境
opt = khfote.parseOptions(eventID,startUtc,endUtc,varargin{:});
[tStart,tStop,tintLoad] = khfote.setupEnvironment(opt);

%% 模块2：读取四星原生数据
fprintf('Loading %s (%s to %s), velocity=%s ...\n', ...
    opt.eventID,opt.startUtc,opt.endUtc,opt.VelocityField);
raw = khfote.loadMmsData(tintLoad,opt.VelocityField);

%% 模块3：同一物理时间窗平均及FOTE前irf_resamp对齐
[rawAligned,aligned,smoothingInfo] = khfote.prepareData( ...
    raw,tStart,tStop,opt.SmoothSeconds);

%% 模块4：调用原FOTE函数识别磁场和流场零点
fprintf('Calling original FOTE_Taylor_Expansion for B ...\n');
mag = khfote.runOriginalFote(aligned.t,aligned.R,aligned.B);
fprintf('Calling original FOTE_Taylor_Expansion for %s ...\n', ...
    opt.VelocityField);
vel = khfote.runOriginalFote(aligned.t,aligned.R,aligned.V);

%% 模块5：FOTE-V误差
[vel.alpha,fluxField] = khfote.computeFoteVError( ...
    aligned.t,aligned.R,aligned.V,aligned.N);

%% 模块6：误差、特征值稳定度及持续时间筛选
[mag,vel,agreement,summary] = ...
    khfote.applyQualityAndSummarize(mag,vel,aligned.t,opt);

%% 模块7：绘图
fig = khfote.plotEvent(rawAligned,aligned,mag,vel,opt,tStart,tStop);

%% 模块8：分类保存全部结果
result = khfote.writeResults(fig,rawAligned,aligned,mag,vel, ...
    fluxField,agreement,summary,smoothingInfo,opt);
end
