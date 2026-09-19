# 1990年起日图：终端激波与日球层顶穿越标记

入口：Run_Voyager_Daily_Boundaries(false)。正式源码位于Interstellar_Daily_Overview_fu；结果独立存入C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Boundary_Markers_1990/V1或V2，文件名后缀_TS_HP。false跳过网站目录同步；所有科学曲线仍直接读取Z:/SPART-WORK/Data/Voyager内原始CDF重算。太阳黑子数仍读SILSO原始发布文件。

| 航天器 | 终端激波TS | 日球层顶HP |
|---|---|---|
| Voyager 1 | 2004-12-16 | 2012-08-25 |
| Voyager 2 | 2007-08-30 | 2018-11-05 |

日期来源：NASA Voyager CRS Mission，https://voyager.gsfc.nasa.gov/mission.html，核查日期2026-09-11。V2日球层顶日期也见NASA发布的https://www.nasa.gov/news-release/nasas-voyager-2-probe-enters-interstellar-space/。V2终端激波采用任务年表列出的2007-08-30作为所要求的单条事件标记。该标记没有展开多次边界往返过程。

本图时间精度为日；虚线按相应UTC日期00:00定位，仅作日期锚点，不表示精确到小时的实测穿越时刻。每个panel中添加相同两条竖虚线：暗红色为Termination shock，深绿色为Heliopause；仅顶部panel显示名称和日期。图内没有额外方法说明或脚注。

所有原日均统计、日中位数、含S4的S1-S7求和、S8排除、缺测以及V2实测MAG补充沿用原规则。逐日表和原图曲线/坐标范围作一致性检查；旧MAT审计只用于验证，未作为科学计算输入。新图源文件、变量、SHA256、样本数与处理规则记录在新V1/V2子目录中。旧图及其统计审计不覆盖；所有原有36个PNG/PDF/FIG文件用SHA256验证保持一致。

新图提供PNG、PDF及MATLAB FIG。验证记录为validation_audit.mat和original_figures_unchanged.csv，运行日志为run_daily_boundaries.log。