function opt = parseOptions(eventID,startUtc,endUtc,varargin)
%KH FOTE 配置模块：解析并标准化单事件分析参数。
%
% 模块职责
%   1. 检查事件号、UTC时间和可选参数；
%   2. 集中维护数据、IRFU、FOTE及输出目录的默认值；
%   3. 将 string/char 输入统一转换为 char，方便旧版MATLAB函数使用。
%
% 输出
%   opt.eventID          事件号，例如 KH005
%   opt.startUtc/endUtc  ISO-8601 UTC字符串
%   opt.SmoothSeconds    物理时间平均窗，单位 s
%   opt.QualityPercent   FOTE/FOTE-V误差阈值，单位 %
%   opt.EigenvalueStabilityThreshold  As/Bs唯一实特征值稳定度阈值
%   opt.MinQualityDurationSeconds 误差判据必须连续成立的最短时间，单位 s
%   opt.ReliableDistanceL 仅供诊断的零点距离阈值，不参与质量筛选

p = inputParser;
p.FunctionName = 'kh_fote_event_irfu';
p.addRequired('eventID',@(x)ischar(x) || isstring(x));
p.addRequired('startUtc',@(x)ischar(x) || isstring(x));
p.addRequired('endUtc',@(x)ischar(x) || isstring(x));
p.addParameter('DataRoot','Z:\SPART-WORK\Data\MMS', ...
    @(x)ischar(x) || isstring(x));
p.addParameter('OutputRoot',fullfile(pwd,'outputs'), ...
    @(x)ischar(x) || isstring(x));
p.addParameter('IrfuRoot', ...
    'C:\Users\Administrator\Documents\irfu-matlab-master', ...
    @(x)ischar(x) || isstring(x));
p.addParameter('FoteRoot','C:\Users\Administrator\Documents\FWD_matlab\FOTE', ...
    @(x)ischar(x) || isstring(x));
p.addParameter('VelocityField','Vi', ...
    @(x)any(strcmpi(string(x),["Vi","Ve"])));
p.addParameter('SmoothSeconds',5, ...
    @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('QualityPercent',40, ...
    @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('EigenvalueStabilityThreshold',0.5, ...
    @(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('MinQualityDurationSeconds',5, ...
    @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('ReliableDistanceL',2, ...
    @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('FigureVisible','off', ...
    @(x)any(strcmpi(string(x),["on","off"])));
p.addParameter('CombinedPdf','',@(x)ischar(x) || isstring(x));
p.parse(eventID,startUtc,endUtc,varargin{:});

opt = p.Results;
opt.eventID = char(string(eventID));
opt.startUtc = char(string(startUtc));
opt.endUtc = char(string(endUtc));

% 统一文本类型，避免部分IRFU旧函数对string类型处理不一致。
textFields = {'DataRoot','OutputRoot','IrfuRoot','FoteRoot', ...
    'VelocityField','FigureVisible','CombinedPdf'};
for k = 1:numel(textFields)
    opt.(textFields{k}) = char(string(opt.(textFields{k})));
end
end
