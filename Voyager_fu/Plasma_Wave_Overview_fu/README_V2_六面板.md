> 后续更新：默认入口现已加入官网全部 PWS 密度点并延长至 2025-09-01，新结果位于 With_Official_PWS_Density。本文保留原 PLS 六面板说明；复现此旧版使用 Run_V2_Plasma_Overview(true,false)。详见 README_V2_官网密度补充.md。

# Voyager 2 全时段日均六面板

## 本次交付

按用户已确认的“全部时段从 1990 年开始”约定，横轴为 1990-01-01 至 2024-12-31。数据按完整 UTC 日统计，点放在当日 12:00 UTC，以细线相连，缺测打断连线。最新官网 COHO 年目录为 2024。本次新增文件位于 V2_Plasma_Overview，原 Voyager_Extended_Overviews/V2 图件与统计表保持不变。

| panel | 参数 | 单位 | 日统计和显示 |
|---|---|---|---|
| a | SILSO 太阳黑子数 | 无量纲 | 原始日值；线性轴 |
| b | 磁场总强度 | nT | 原始小时标量的有限值算术均值；线性轴 |
| c | LECP P1 质子通量，0.52–1.45 MeV | cm⁻² s⁻¹ sr⁻¹ MeV⁻¹ | 原始小时通量的有限值算术均值；对数轴 |
| d | PLS 质子数密度 | cm⁻³ | 原始小时参数的有限值算术均值；对数轴 |
| e | PLS 质子温度 | K | 原始小时温度的有限值算术均值；对数轴 |
| f | PLS 总流速 | km/s | 原始小时标量流速的有限值算术均值；线性轴 |

数据表共 12,784 个 UTC 日。磁场、P1、密度、温度和速度分别保留有效样本数。太阳黑子数保留源文件行号。各参数独立按有效值统计，不要求所有参数同时有效。

## V2 的 P1 能段核对

用户于 2026-09-15 要求按照已发表文章标注能段。新图使用 **0.52–1.45 MeV**：

1. Decker, R. B., Roelof, E. C., and Krimigis, S. M. (1999), *Solar Energetic Particles from the April 1998 Activity: Observations from 1 to 72 AU*, Proceedings of the 26th International Cosmic Ray Conference, vol. 6, p. 328, SH.1.6.10。图 1 在同一张图中分别标注 Voyager 2 为 0.52–1.45 MeV、Voyager 1 为 0.57–1.78 MeV。
   https://galprop.stanford.edu/elibrary/icrc/1999/proceedings/root/vol6/s1_6_10.pdf
2. Rice, W. K. M., Zank, G. P., Richardson, J. D., and Decker, R. B. (2000), *Ion injection and shock acceleration in the outer heliosphere*, Geophysical Research Letters, 27, 509–512。论文摘要明确使用 Voyager 2 的 0.52–1.45 MeV。
   https://www.research.ed.ac.uk/en/publications/ion-injection-and-shock-acceleration-in-the-outer-heliosphere/
3. NASA 官方 V2 LECP 小时产品和 COHO 页面亦使用 0.52–1.45 MeV：
   https://omniweb.gsfc.nasa.gov/ftpbrowser/v2_lecp_flux_1h.html
   https://omniweb.gsfc.nasa.gov/coho/form/voyager2.html

本次遍历 420 个 COHO CDF，逐文件检查 protonFlux1_LECP 的 FIELDNAM 与 UNITS。继续使用原 P1 通量数值，不因标签变更而缩放通量、变更校准或合并能道。本次结论针对 V2 COHO 的 protonFlux1_LECP；其他扇区产品的能段元数据应单独核实。原 V2 图的 0.57–1.78 MeV 标签未在此次任务中覆盖修改。

## 数据来源、覆盖及数据库核查

| 参数 | 本图第一个有效日 | 最后一个有效日 | 有效日数 |
|---|---|---|---|
| 太阳黑子数 | 1990-01-01 | 2024-12-31 | 12,784 |
| 总磁场 | 1990-01-03 | 2024-12-31 | 12,079 |
| P1 | 1990-01-01 | 2023-08-30 | 7,548 |
| PLS 质子密度 | 1990-01-01 | 2018-11-04 | 10,170 |
| PLS 质子温度 | 1990-01-01 | 2018-11-04 | 10,170 |
| PLS 总流速 | 1990-01-01 | 2018-11-04 | 10,170 |

默认绘图入口直接读取网站原始 CDF：

- COHO: Z:/SPART-WORK/Data/Voyager/voyager2/coho/1hr/l2/merged_mag_plasma/YYYY/MM/
- MAG 补充: Z:/SPART-WORK/Data/Voyager/voyager2/mag/48s/reviewed_vim/YYYY/
- 太阳黑子数: Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv。该产品原生为 SILSO 日值文本；使用既有原始文件。

2026-09-15 比对官方目录后，本任务范围内的 464 个 CDF 全部已在数据库中，新增 CDF 下载数为 0：

- 420 个 COHO 月文件（1990–2024）
- 18 个常规太阳风 PLS 高分辨率年文件（1990–2007）
- 12 个日鞘 PLS 高分辨率年文件（2007–2018）
- 14 个 reviewed MAG 48 秒年文件（2009–2022）

下载核查选择官网同日期最高版本，检查 CDF 文件头，保存各文件 SHA-256、大小、URL 和核查时间。本地存在文件未被覆盖。原始目录页面及产品说明快照、清单保存在：
Z:/SPART-WORK/Data/Voyager/voyager2/source_verification/plasma_overview_1990_2024/

官方目录入口：
- https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/coho1hr_magplasma/
- https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma_cdaweb/hires_plasma/
- https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma_cdaweb/hires_plasma/heliosheath/
- https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/magnetic_fields_cdaweb/vim_48secmag/
- https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma/hour/
- https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma/daily/

高分辨率 PLS 文件用于确认数据库完整性。本图科学输入采用 COHO 中同一来源、同一小时分辨率的 PLS 参数，不将高分辨率记录混入小时日均。另查常规小时、日值、高分辨率及 ions 产品目录，未找到可补入本图的 2018-11-04 之后常规 PLS 三参数序列。PLS 仪器持续运行到 2024 年不代表该常规参数产品持续发布至 2024；本说明不把 2018 年后的仪器电流观测等同于缺测，也不声称不存在专门分析或论文数据。仪器运行背景：https://www.space.mit.edu/news/an-interstellar-instrument-takes-a-final-bow/

## 处理规则及限制

- 所有粒子和等离子体参数沿用原始 CDF 的 FILLVAL、VALIDMIN、VALIDMAX 及既有 CDF 读取器规则。无新质量阈值、插值、平滑、拟合、背景扣除、零填充或人为删除尖峰。
- 若 Epoch 重复，仅逐变量完全相同的记录保留一次，重复记录单独存入审计；冲突会停止运行。本次重复情况见 MAT。
- 太阳黑子数 -1 表示缺测，0 是有效值。对数面板中的非正值仅在显示副本中隐藏，统计表保留原有效值；本次三个对数面板均无非正日值。
- 温度直接取 protonTemp，源单位为 K。产品温度与热速的关系为 T=60.5 Vth²；不先平均热速再平方，不将 K 误认为 eV。流速取 V，保留标量含义，不由日均 RTN 分量计算速度模长。
- 密度为 PLS 质子数密度。本图未以 PWS 电子密度替代，也未采用电子密度与质子密度相等的假设。
- 沿用用户先前批准的磁场补齐：仅在 COHO 整日缺测时，使用 reviewed MAG F1 的 48 秒有效记录先求 UTC 小时均值，再求日均值。已有效的 COHO 日值保持原样。没有插值磁场。共补 340 日，仍有 705 日缺测。
- 原 a/b/c 数值同以前交付的 all_daily_values.csv 核对。该旧 CSV 只用于验证，未参与科学计算。残余差异限于 CSV 文本保存精度；记录最大绝对差和相对差。
- 本图未新增方法脚注。源文件清单、假设和限制在本说明及 MAT 审计中保留。

## 运行

正式程序目录：
C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Plasma_Wave_Overview_fu/

在 MATLAB 中：
~~~matlab
addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Plasma_Wave_Overview_fu');
Run_V2_Plasma_Overview;
Validate_V2_Plasma_Overview;
~~~

默认入口先运行 Sync_V2_Plasma_Overview_Archive.ps1 核对官方目录并补下载缺失 CDF，然后从原始 CDF 计算。若官网核查已经完成，可用 Run_V2_Plasma_Overview(false) 复用源文件清单；此时仍重新读取原始 CDF，不用 MAT/CSV 代替输入。若官网最新年份改变，程序停止并提示显式更新日期范围。

依赖现有 IRFU dataobj 和项目 Voyager_Read_CDF_Product、Voyager_Supplement_V2_MAG 等函数；路径和代码 SHA-256 保存在审计中。

## 输出与验证

输出目录：
C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V2_Plasma_Overview/

- V2_1990_20241231_daily_PLS_6panels.png / .pdf / .fig：新图。
- V2_daily_six_panel_values.csv：全部 12,784 日参数、有效样本数与磁场来源。
- coverage.csv：各参数覆盖。
- missing_date_intervals.csv：逐参数连续缺测日期段。
- abc_preservation.csv：原 a/b/c 数值核对。
- V2_plasma_overview_audit.mat：原始选中 COHO 记录、CDF 来源/变量元数据/全局属性/哈希、MAG 补充记录、日值、覆盖、假设及旧文件保护哈希。
- source_recompute/：MAG 补充通路的独立审计。
- final_validation.json：验证结果。校验六条图形曲线、时间范围、对数轴、缺测掩码、样本数守恒和旧图哈希，并重新读取 11 个跨时期原始 CDF 抽查日统计（包括 2007 年终端激波前后、2018 年最后有效日及随后空白日）。
- matlab_final_run_20260915.log：最终绘图运行日志，保留图形对象审计调试过程。matlab_validation_20260915.log 和 final_validation.json：修正验证报告字段后单独运行的最终验证结果。早期调试日志继续保留。
