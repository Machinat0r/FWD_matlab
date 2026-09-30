# 20个事件重新跑图（2026-09-24，v2）

本版按用户原 Overview_download.m、Overview_download_mms4.m 和 MMS_orbit.m 的直接脚本形式重写。两个新 MATLAB 文件均无 function 定义，不调用上一版自行封装的科学读取、拼接、绘图或轨道函数。

## 实际运行的程序

- Overview_events_original_style_20260924_v2.m：按 %% load data、%% B plot、%% Vi plot、%% Ve plot、%% N plot、%% Ti plot、%% Te plot、能谱与保存等分区执行。
- MMS_orbit_original_style_20260924_v2.m：直接调用 mms.mms4_pl_conf(tint) 与 mms.mms4_pl_conf('gsm')。
- MMS_event_overview_20260924_events.json：20个事件的原始时间与候选解释。
- references：原程序副本，仅供对照。用户原文件未覆盖。

在 MATLAB 中执行：
~~~
cd('C:\Users\Administrator\Documents\FWD_matlab\MMS_fu')
EventList=1:20; Spacecraft=1:4;
Overview_events_original_style_20260924_v2
EventList=1:20;
MMS_orbit_original_style_20260924_v2
~~~

单独重画某一事件/卫星：
~~~
EventList=4; Spacecraft=1;
Overview_events_original_style_20260924_v2
~~~

每次都会重画选定范围。两脚本开头的路径集中列出，可直接修改。

## 使用的现成接口

| 步骤 | 已有程序或函数 |
| --- | --- |
| 初始化路径和数据库 | irf('check_path'), mms.db_init |
| 读取原CDF | mms.get_data, mms.db_get_ts, mms.db_get_variable |
| CDF填充值与时间定义 | mms.variable2ts |
| 时间/矩阵/坐标转换 | irf.tint, irf_time, EpochTT, EpochUnix, irf.ts2mat, irf_gse2gsm |
| overview | irf_subplot, irf_plot, irf_spectrogram, irf_legend, irf_zoom, irf_plot_axis_align |
| 密集曲线显示 | irf_plot(...,'reduce'), IRFU包自带 reduce_to_width |
| 四星位置与构型 | mms.mms4_pl_conf（内部继续使用其已有坐标转换、插值与边界模型） |

## 绘图设置

每个事件两端各扩展10分钟。9个panel为B、Vi、Ve、E、N、Ti、Te、Ei、Ee；y轴标签沿用原overview，panel间隙0.0015归一化图高，没有底部模式时间条，也没有事件竖线；使用原程序的set(gca,"XTickLabelRotation",0)将底部时间标签设为水平。温度沿用(Tparallel+2*Tperpendicular)/3及平行/垂直分量。

先读取burst与survey/fast，绘制时在有效burst连续段内隐藏低时间分辨率曲线，其余时段保留survey/fast。测量值不平滑、不插值。脚本只插入用于断开曲线的NaN，避免跨真实缺口连线。模式信息与样本数保存在Z盘对应JSON中。

本次重画的B、Vi、Ve和E均为GSM。电场只读取GSE产品并通过irf_gse2gsm转换，只有DSL时留空；标题不再标坐标系。缺失电子数据留白。FPI时间中心和CDF填充值沿用mms.variable2ts的原逻辑，包括该函数给出的时间修正提示。

burst电场按600秒分段读取，使用IRFU现有2GB有界内存缓存。显示压缩调用包内reduce_to_width，该函数输出显示时间箱中的极小/极大值。显示时间箱仅用于绘图，不作为新的物理数据输出。

## 轨道图

事件单时刻直接使用；区间事件取中点，UTC，GSM。EV01–EV15使用归档MEC L2 CDF。EV16–EV20使用已归档NASA SSCWeb原始GSE星历；直接解包JSON字段并将四星相同的原生时间轴传给mms.mms4_pl_conf。未自行重采样SSC数据。

本机IRFU的mms.mms4_pl_conf中，日期格式 HH:MM:SS.mmm 与现有时间格式解析器冲突。轨道脚本保留原文件备份后，仅把两处显示格式改为 HH:MM:SS。图形计算没有改动。备份位于：
Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924\mms4_pl_conf_before_date_format_fix.m

图上的磁层顶/弓激波为原程序的模型边界，模型参数及来源在图内显示。投影顺序为X–Y、X–Z、Y–Z；下排相对构型采用相同顺序。X–Y/X–Z的纵轴范围关于零对称，地球居中；X方向向地球和卫星两侧至少扩展10 RE，Y/Z方向及相对构型保留至少约14%的边缘空间。关闭轨道图网格。

## 数据与成果位置

- 原始CDF及SSCWeb JSON：Z:\SPART-WORK\Data\MMS，沿用原分类。
- 本版小型数据记录：Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924_v2。
- 图和PPT：C:\Users\Administrator\Documents\KH\MMS_events_original_style_20260924_v2。
- 新程序：C:\Users\Administrator\Documents\FWD_matlab\MMS_fu。

原有仪器复用归档数据；新增FEEPS、HPCA原始CDF按官方产品目录补下载至同一Z盘数据库。8月11日原有仪器仅取得B，本轮另加入可用FEEPS、HPCA；8月28–29日和9月16日缺少所需L2科学数据，PPT保留此前已提供的官方quicklook参考页。重跑不会将这些缺测变为已获得的数据。FPI Vi/Ve使用现有fast数据，现有数据库没有这些事件的对应burst moments。

PPT中的科学图为嵌入的高清图片，曲线内容通过所附MATLAB脚本修改重画。PPT封面、目录及说明文字可编辑。

用户确认：仅修改本次重画的图，44页官方quicklook保持原图、原坐标系和原布局。

## FEEPS、HPCA（本轮追加要求）

有数据时，在同一张overview下方增加FEEPS电子/离子全向能谱、HPCA H+/He+/He++/O+密度和各离子能谱，最多7个panel。全向通量和能谱结构直接调用mms.get_data、PDist.specrec，不另写积分或通量算法。只增加含有效数据的panel。burst优先，survey补充。FEEPS通量单位为1/(cm^2 s sr keV)，HPCA为1/(cm^2 s sr eV)，按CDF和现成接口返回的单位标注。能量轴沿用eV。

IRFU现有get_data的HPCA He++ survey白名单有一处拼写错误：Omnifluxheplusplus_srvy_brst_l2。将其改为Omnifluxheplusplus_hpca_srvy_l2，使该名称进入原有HPCA读取分支，计算逻辑保持原样。备份在Z:\SPART-WORK\Data\MMS\derived\events_original_style_20260924_v2\get_data_before_hpca_name_fix.m。

官方文件目录用用户原SDCFilenames查询。新增原CDF仍存Z盘官方分类目录，查询、候选覆盖和下载记录位于上述derived目录下的particles。下载复用之前已有MMS_event_overview_20260923_download.py中的原CDF下载/校验代码，按卫星分开执行；全部科学处理仍在MATLAB中。

## FEEPS survey占位值兼容修正

实测2026-07-28 MMS1 survey未启用探头的强度和能量值为-2147483648，与CDF声明的FILLVAL=-1e31不同。原mms.get_data会把该值参加探头平均，使电子全向通量和能量为负。本版MMS_IRFU_particle_compatibility_20260924_v2.m仅在现有读取器中将该精确占位值置NaN，能量仍按原有探头算术平均（忽略缺失），通量仍使用原有mean(...,3,'omitnan')和增益系数。未设置新的通量阈值，也未补写测量值。所有科学处理继续由MATLAB与现成IRFU类/函数完成。

修改前文件get_data_before_feeps_pad_fix.m和逐行差异get_data_feeps_pad_fix.diff保存在本版Z盘derived目录。附带compatibility/get_data_feeps_pad_fix.diff供核对。该修正后的样例survey电子能量范围为约47.7–575.8 keV；正常burst样例的通量和能量范围保持一致。此处理仅解决已核实的占位值问题，不额外加入未验证的探头质量筛选。
