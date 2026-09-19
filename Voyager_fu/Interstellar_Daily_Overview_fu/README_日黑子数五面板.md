# 日黑子数与 Voyager 五面板图

入口：`Run_V1_Daily_WithS4_Sunspots(false)`。正式程序位于 FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu。

本图对应用户提供的日统计、包含 S4 的四面板图，上方新增每日太阳黑子数，panel 顺次编号 a–e。时间为 [2012-08-25 00:00 UTC, 2021-12-17 00:00 UTC)。没有对太阳黑子数作太阳风传播时延校正。

太阳数据实际目录为 Z:/SPART-WORK/Data/Solar_Indices。直接读取 raw/sunspot/SN_d_tot_V2.0.csv，来源 WDC-SILSO, Royal Observatory of Belgium，Daily total sunspot number V2.0，https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv ，DOI 10.24414/qnza-ac80。CSV 是该源的原始发布格式。第 1–3 列为日期，第 5 列为日黑子数；-1 留空，零保留；不平滑、不插值、不重采样。时间点标在当天 UTC 中午以对齐已有日统计。源行号、全部源列和 SHA256 保存在审计中。

Voyager 使用既有入口直接重读原始 CDF。原来的磁场日平均、P1 日平均、P1 日中位数保持定义；最下方为 S1–S7 日均通量之和，包含 S4、排除 S8，七扇区任何缺测则该日总和留空。沿用既有 L1 优先和历史换算规则，详细方法与源 CDF 清单见 source_recompute/V1_daily_overview_audit.mat。图的源能段标签差异沿用项目约定，不重新标定数据。

太阳黑子数和磁场使用线性轴，三个粒子面板使用对数轴；0.4 点宽细线连接，NaN 断线；对数轴隐藏非正值但审计保留原值。不在图内新增方法注释。

输出：V1_daily_5panels_with_S4_sunspots.png/.pdf/.fig、daily_5panels_values.csv、daily_5panels_audit.mat；日志位于上一级 sunspots_run_20260909.log。CSV/MAT 为结果交付和审计，不作为默认科学输入。

已检索 IRFU irf_get_data_omni 的 ssn 通路；它读取 OMNI 产品，本任务使用用户指定的本地 SILSO 原始文件，故用 MATLAB readmatrix 直接读取，并继续复用现有 IRFU CDF 读取通路。
