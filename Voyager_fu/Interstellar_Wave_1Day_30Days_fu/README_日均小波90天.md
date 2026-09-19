# 日均磁场 Morlet 小波：2–90 天，延长至最新公开磁场日期

## 运行与结果

2026-09-19 用户要求把所发图从 2012-08-25 延长到 Voyager 1 最后有磁场数据的日期。
本次直接核查 NASA CDAWeb 官方目录和原始 CDF，确认当前最后有效磁场日为 **2025-06-30**。

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Wave_1Day_30Days_fu');
    result = Run_V1_Wave_1Day_90Days;

默认每次先在线核查官方 COHO/VIM 目录，并从原始 CDF 的有限磁场值判定终点，再读取原始 CDF 重新计算。
仅在需要复现旧时间范围时使用：
    
    result = Run_V1_Wave_1Day_90Days('UseLatestAvailable',false,'OutputRoot','指定独立结果目录');

正式程序：C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Wave_1Day_30Days_fu
结果目录：C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_90Days
图件：V1_daily_B_Morlet_2day_90days.png / .pdf / .fig
先前截至2021-12-16的图件、输入和分析结果保存在 before_extension_2021-12-16/。

图形沿用用户本次提供的参考：周期轴2–90天（顶部2天），log10功率色标0–3，两个上下对齐面板，无新增方法脚注。

## 终点核查（2026-09-19）

- 官方 COHOWeb 可用时间表：Voyager 1 磁场到 2025/06/30。
- COHO 最新月文件：voyager1_coho1hr_merged_mag_plasma_20250601_v01.cdf。
- 原始 ABS_B 最后有限记录：2025-06-30 00:00:00 UTC，0.466 nT。
- 独立核对最新审阅后的48秒产品 voyager1_48s_mag-vim_20250101_v01.cdf，其 F1 最后有限记录为 2025-06-30 00:01:28 UTC。
- 图中最后一天包含该日真实记录。计算区间为 [2012-08-25 00:00:00, 2025-07-01 00:00:00) UTC。
- 日期代表目前公开数据的可用终点，不表示航天器或磁强计在该日停止运行。
- 末日的有效样本数完整保留，继续沿用有限小时样本的算术日均定义，不人为补齐末日其余小时。
- 不向第一个或最后一个真实观测日之外外推。

[官方 COHOWeb 可用时间表](https://omniweb.gsfc.nasa.gov/coho/html/cw_data.html)
[官方小时 CDF 目录](https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/coho1hr_magplasma/)
[官方48秒 CDF 目录](https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/magnetic_fields_cdaweb/vim_48secmag/)

目录 HTML、下载文件清单、核查时间与 SHA256 均归档在 Z:/SPART-WORK/Data/Voyager/source_verification/V1_daily_morlet_latest/ 对应运行时间子目录。
结果中另存 official_availability_audit.mat，完整分析中的 result.Availability 和 result.Config.AvailabilityAudit 也保留核查依据。
本次所需原始文件已在归档中，无须用 CSV 或旧 MAT 替代原始科学输入。

## 原始数据、日均与缺测

使用152个原始 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA 月 CDF，沿用 Epoch 和 ABS_B（nT），按 UTC 日对有效小时标量磁场做算术平均，代表时刻12:00。
不从平均矢量重新计算标量模值，不改变 CDF 的填充值/有效范围规则，不增加覆盖率或样本数阈值。

官方目录未发布2023年10、11、12月的该产品 CDF。只有经过在线目录确认的这三个月允许作为原始数据缺口；本地遗漏其他已发布 CDF 会报错或先下载，不会被静默视为科学缺测。

本次共有4693个日格，4299个真实日均值、394个原始缺测日。
沿用前次用户明确要求的线性缺测插值，对延长后的图使用相同处理。
最长连续缺口为2023-09-28至2024-05-29，共245天；该段的小波结果依赖插值，不能提供独立的实测波动证据。

每个缺测日由前后最近真实日均值线性连接，alpha=(t-t_left)/(t_right-t_left)。
显式调用 irf_resamp(...,'linear')，关闭自动平均分支。只赋值到原始 NaN 日，已有观测日的数值和样本数保持不变。
没有新增长缺口截断阈值，不外推。末端依照真实磁场终点截取后才插值。
原始NaN保存在 result.ObservedDaily 和 Daily.BMeanObserved_nT；Daily.IsInterpolated 标记394个估计日，MAGSampleCount在这些日子仍为0。
V1_daily_interpolation_audit.csv 逐日记录两端真实日期/数值、线性权重、插值结果和缺口长度。

## 小波与显示保持的设置

继续复用本地原始 irf_wavelet Morlet：
- 100个对数频点，频率为1/(90*86400)=1.286008230453e-7至1/(2*86400)=5.787037037037e-6 Hz。
- wavelet_width=5.36，returnpower=1，cutedge=1。
- 原归一化 P=2*pi*abs(W)^2/f，单位 nT²/Hz。
- 完整插值序列统一变换，不在原缺测日处分段。
- 原生边界屏蔽保留，长周期在整个时段两端的空白更宽。
- 原函数奇数样本末点舍弃行为保留：本次4693个输入日格，4692点进入变换，原始末日值仍保留在数据与上面板。
- 不新增去均值、去趋势、预滤波、平滑、背景扣除、姿态或仪器几何假设。
- 波谱不作时间压缩平均，按原始日格平色显示。仅纵轴从频率换算成周期，功率值及单位仍按每Hz。
- 超过色标的数值只作颜色饱和，完整数值保留在分析 MAT 中。

输出还包括 source_manifest.csv、V1_daily_magnetic_input.csv、V1_daily_morlet_segments.csv、V1_daily_morlet_analysis.mat、运行日志及延长区间验证记录。

