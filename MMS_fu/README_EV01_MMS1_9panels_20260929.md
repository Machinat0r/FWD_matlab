# EV01 MMS1 九个 panel（2026-09-29）

本次按用户要求重画2026-07-21事件的MMS1，删除全部FEEPS和HPCA panel。保留B、Vi、Ve、E、N、Ti、Te、离子能谱、电子能谱。时间为2026-07-20 23:50至2026-07-21 04:10 UTC；矢量GSM；FPI能谱范围20–40000 eV。数据读取、burst优先和缺口处理沿用20260927脚本。

运行程序：Overview_EV01_MMS1_9panels_20260929.m。将事件JSON与脚本放在C:\Users\Administrator\Documents\FWD_matlab\MMS_fu，在MATLAB运行：

    EventList=1; Spacecraft=1;
    Overview_EV01_MMS1_9panels_20260929

沿用本机IRFU包及Z:\SPART-WORK\Data\MMS中的原始CDF。路径集中在脚本开头。输出到KH\MMS_EV01_MMS1_20260929；记录到Z盘derived\EV01_MMS1_9panels_20260929。没有新增MATLAB function；修改差异见overview_remove_feeps_hpca.diff。