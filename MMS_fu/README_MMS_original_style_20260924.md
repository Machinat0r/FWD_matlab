# 20个事件重新跑图（2026-09-24）

本版按用户原 Overview_download.m、Overview_download_mms4.m 和 MMS_orbit.m 的直接脚本形式重写。两个新 MATLAB 文件均无 function 定义，不调用上一版自行封装的科学读取、拼接、绘图或轨道函数。

## 实际运行的程序

- Overview_events_original_style_20260924.m：按 %% load data、%% B plot、%% Vi plot、%% Ve plot、%% N plot、%% Ti plot、%% Te plot、能谱与保存等分区执行。
- MMS_orbit_original_style_20260924.m：直接调用 mms.mms4_pl_conf(tint) 与 mms.mms4_pl_conf('gsm')。
- MMS_event_overview_20260924_events.json：20个事件的原始时间与候选解释。
- references：原程序副本，仅供对照。用户原文件未覆盖。

在 MATLAB 中执行：
~~~
cd('C:\Users\Administrator\Documents\FWD_matlab\MMS_fu')
EventList=1:20; Spacecraft=1:4;
Overview_events_original_style_20260924
EventList=1:20;
MMS_orbit_original_style_20260924
~~~

单独重画某一事件/卫星：
~~~
EventList=4; Spacecraft=1;
Overview_events_original_style_20260924
~~~

每次都会重画选定范围。两脚本开头的路径集中列出，可直接修改。

## 使用的现成接口

| 步骤 | 已有程序或函数 |
| --- | --- |
| 初始化路径和数据库 | irf('check_path'), mms.db_init |
| 读取原CDF | mms.get_data, mms.db_get_ts, mms.db_get_variable |
| CDF填充值与时间定义 | mms.variable2ts |
| 时间/矩阵/坐标转换 | irf.tint, irf_time, EpochTT, EpochUnix, irf.ts2mat, irf_gse2gsm |
| overview | irf_subplot, irf_plot, irf_spectrogram, irf_legend, irf_zoom, irf_plot_axis_align, irf_pl_mark |
| 密集曲线显示 | irf_plot(...,'reduce'), IRFU包自带 reduce_to_width |
| 四星位置与构型 | mms.mms4_pl_conf（内部继续使用其已有坐标转换、插值与边界模型） |

## 绘图设置

每个事件两端各扩展10分钟。9个panel为B、Vi、Ve、E、N、Ti、Te、Ei、Ee；y轴标签沿用原overview，panel间隙0.002归一化图高，没有底部模式时间条。温度沿用(Tparallel+2*Tperpendicular)/3及平行/垂直分量。

先读取burst与survey/fast，绘制时在有效burst连续段内隐藏低时间分辨率曲线，其余时段保留survey/fast。测量值不平滑、不插值。脚本只插入用于断开曲线的NaN，避免跨真实缺口连线。模式信息与样本数保存在Z盘对应JSON中。

B、Vi、Ve使用GSM。MMS1–3 E使用GSM；MMS4使用现有DSL XY电场，缺失电子数据留白。FPI时间中心和CDF填充值沿用mms.variable2ts的原逻辑，包括该函数给出的时间修正提示。

burst电场按600秒分段读取，使用IRFU现有2GB有界内存缓存。显示压缩调用包内reduce_to_width，该函数输出显示时间箱中的极小/极大值。显示时间箱仅用于绘图，不作为新的物理数据输出。

## 轨道图

事件单时刻直接使用；区间事件取中点，UTC，GSM。EV01–EV15使用归档MEC L2 CDF。EV16–EV20使用已归档NASA SSCWeb原始GSE星历；直接解包JSON字段并将四星相同的原生时间轴传给mms.mms4_pl_conf。未自行重采样SSC数据。

本机IRFU的mms.mms4_pl_conf中，日期格式 HH:MM:SS.mmm 与现有时间格式解析器冲突。轨道脚本保留原文件备份后，仅把两处显示格式改为 HH:MM:SS。图形计算没有改动。备份位于：
Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924\mms4_pl_conf_before_date_format_fix.m

图上的磁层顶/弓激波为原程序的模型边界，模型参数及来源在图内显示。

## 数据与成果位置

- 原始CDF及SSCWeb JSON：Z:\SPART-WORK\Data\MMS，沿用原分类。
- 本版小型数据记录：Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924。
- 图和PPT：C:\Users\Administrator\Documents\KH\MMS_events_original_style_20260924。
- 新程序：C:\Users\Administrator\Documents\FWD_matlab\MMS_fu。

本次重跑读取现有归档数据。8月11日现有所需L2数据仅能重画B；8月28–29日和9月16日缺少所需L2科学数据，PPT保留此前已提供的官方quicklook参考页。重跑不会将这些缺测变为已获得的数据。FPI Vi/Ve使用现有fast数据，现有数据库没有这些事件的对应burst moments。

PPT中的科学图为嵌入的高清图片，曲线内容通过所附MATLAB脚本修改重画。PPT封面、目录及说明文字可编辑。
