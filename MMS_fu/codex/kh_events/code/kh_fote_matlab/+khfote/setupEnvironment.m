function [tStart,tStop,tintLoad] = setupEnvironment(opt)
%KH FOTE 环境模块：配置IRFU/FOTE路径、MMS数据库和分析时间。
%
% 模块职责
%   1. 添加完成本任务所需的最小MATLAB路径；
%   2. 验证原FOTE函数和IRFU函数均可调用；
%   3. 注册本地MMS数据库；
%   4. 将UTC字符串转换为Unix秒，并为时间平均增加读取边界。
%
% 时间输出均采用Unix秒；tintLoad为IRFU GenericTimeArray时间区间。

addRequiredPaths(opt.IrfuRoot,opt.FoteRoot);
registerMmsDatabase(opt.DataRoot);

tStart = iso2epoch(opt.startUtc);
tStop = iso2epoch(opt.endUtc);
if ~isfinite(tStart) || ~isfinite(tStop) || tStop<=tStart
    error('Invalid UTC interval: %s / %s',opt.startUtc,opt.endUtc);
end

% 读取区间向两端扩展，保证中心平均窗在事件边界处有完整数据。
loadMargin = max(10,opt.SmoothSeconds/2+2);
tintLoad = irf.tint(sprintf('%s/%s',epoch2iso(tStart-loadMargin), ...
    epoch2iso(tStop+loadMargin)));
end


function addRequiredPaths(irfuRoot,foteRoot)
% 仅加入当前流程用到的目录，减少IRFU启动时的额外路径检查。
dirs = {irfuRoot,fullfile(irfuRoot,'irf'),fullfile(irfuRoot,'plots'), ...
    fullfile(irfuRoot,'plots','mms'), ...
    fullfile(irfuRoot,'mission','cluster'), ...
    fullfile(irfuRoot,'mission','mms'), ...
    fullfile(irfuRoot,'contrib','nasa_cdf_patch'), ...
    fullfile(foteRoot,'FOTE_function')};

required = {'mms.get_data','irf.tint','irf_resamp','c_4_grad', ...
    'FOTE_Taylor_Expansion'};
available = cellfun(@(name)~isempty(which(name)),required);
if ~all(available)
    % 已由startup或MATLABPATH配置的会话无需重复addpath。这个判断也能
    % 避免大型用户路径在每个批处理事件前重复刷新。
    for k = 1:numel(dirs)
        if exist(dirs{k},'dir'), addpath(dirs{k}); end
    end
end

leapFile = fullfile(irfuRoot,'contrib','nasa_cdf_patch', ...
    'CDFLeapSeconds.txt');
if exist(leapFile,'file'), setenv('CDF_LEAPSECONDSTABLE',leapFile); end

for k = 1:numel(required)
    if isempty(which(required{k}))
        error('Required MATLAB function not found: %s',required{k});
    end
end
end


function registerMmsDatabase(dataRoot)
% IRFU通过全局MMS_DB访问本地CDF。此全局变量是IRFU现有接口要求。
global MMS_DB;
MMS_DB = mms_db;
localDb = mms_local_file_db([char(string(dataRoot)) filesep]);
MMS_DB.add_db(localDb);
end
