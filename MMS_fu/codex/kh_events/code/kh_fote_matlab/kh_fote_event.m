function result = kh_fote_event(eventID,startUtc,endUtc,varargin)
%KH_FOTE_EVENT 对外入口模块：转入MATLAB/IRFU模块化单事件流程。
%
% 保留这个短入口可兼容此前的调用方式。所有计算模块位于+khfote目录，
% 流程控制位于kh_fote_event_irfu.m。
result = kh_fote_event_irfu(eventID,startUtc,endUtc,varargin{:});
end
