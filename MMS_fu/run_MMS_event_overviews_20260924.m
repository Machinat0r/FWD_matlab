function run_MMS_event_overviews_20260924(eventNumbers,spacecraft)
% 重现2026-09-24 PPT中的MATLAB图。默认重画有CDF的15个事件、MMS1-4。
% 用法：run_MMS_event_overviews_20260924; 或 run_MMS_event_overviews_20260924(4,1)
% 路径在MMS_event_overview_20260924.m第1节中设置；不下载新CDF。
if nargin<1,eventNumbers=1:15;end
if nargin<2,spacecraft=1:4;end
addpath(fileparts(mfilename('fullpath')));
MMS_event_overview_20260924(eventNumbers,spacecraft,false);
end
