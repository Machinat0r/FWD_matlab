# Voyager 2 六面板：加入官网全部 PWS 密度点

用户要求：只要数据网站直接提供参数就绘制，未提供的参数留空，图内不增加额外说明标注。本次结果单独保存在 With_Official_PWS_Density 子目录，先前六面板和五面板图件及其统计文件保留。

## 结果

时间范围：1990-01-01 至 2025-09-01。沿用太阳黑子数、磁场、LECP P1 通量、密度、温度、总流速六面板。

密度面板保留全部 PLS 质子密度日均值，并加入 NASA PDS 当前发布的 Voyager 2 星际空间电子密度产品的全部 5 条原始记录：

| 原始 UTC 时刻 | 官网 N_e（cm^-3） |
|---|---:|
| 2019-02-18 11:22:00 | 0.039 |
| 2020-06-26 21:16:00 | 0.087 |
| 2021-05-23 23:26:00 | 0.120 |
| 2022-11-10 04:45:00 | 0.120 |
| 2025-09-01 13:29:00 | 0.120 |

5 条记录的原始时刻和 N_e 数值全部原样保留。最新一条位于 2025 年，因此将图的数据时间范围延长到 2025-09-01。太阳黑子数从已有原始 SILSO 文件读到同一天；其他参数只在其实际有效范围内绘制。

温度、总流速和 PLS 质子密度最后有效日仍为 2018-11-04；MAG 最后有效日为 2024-12-31，P1 最后有效日为 2023-08-30。官网提供的 PLS L/M 模式电流属于仪器观测量，未提供可直接使用的连续后续温度、流速参数。本次未自行开展电流反演。

## 图内显示与科学处理

- 图中只保留正常标题、坐标轴、单位和 panel 编号，没有增加来源图例、方法注释或脚注。
- 密度轴改为通用 Density（cm^-3）。PLS 质子密度与 PWS 电子密度分别保留在数据及审计中，不将电子密度改写为质子密度，不将其填入 PLS 数组，也不假定两者数值相等。
- 原日均曲线仍以细线和点显示。PWS 五点使用官网原始时刻，跨越没有观测的时间段保留空白，不生成日均、插值、拟合或桥接曲线。
- PWS 源文件的 N_e、频率及上下限全部保留；本图直接使用已发布 N_e，不自行通过频率公式重算或修正。
- P1、密度和温度使用对数轴；其他面板为线性轴。密度轴上方和时间轴右侧留出纯显示余量，确保最新及最高的点完整可见。这些余量不产生新数据。
- 原 COHO 日均、CDF 填充值规则、MAG 48 秒→小时→日均的既有补充处理保持不变。对 1990–2024 原 a/b/c 的比较通过。
- 新时间范围包含 13,028 个完整 UTC 日。派生日表存放原 COHO/PLS 参数；PWS 原始记录另表交付，便于区分两种密度。

## 原始数据及网站核查

密度官方产品：
https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/vg2-vlism-density-2019-2025.csv

PDS 标签：
https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/vg2-vlism-density-2019-2025.lblx

官方数据目录：
https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/

该密度产品原生格式为 CSV，目录中没有对应密度 CDF。按用户“网站提供就画”的要求直接读取官方 CSV，未转成自制 CDF；其余航天器参数继续直接读取网站原始 CDF。标签声明 5 条、22 列，MD5 为 76287d3eb8388b6816848b357ce927e8，下载文件与官方校验和一致。文件名含 2019–2025；实际时刻按记录和标签确认，不由文件名推定。

原始文件保存于：
Z:/SPART-WORK/Data/Voyager/voyager2/pws/derived/electron_density/native/PDS_release_20260910/

原始 URL、SHA-256、大小、核查时间和官方目录快照：
Z:/SPART-WORK/Data/Voyager/voyager2/source_verification/plasma_overview_density_extension/official_sources.json

同时核查：
- MIT 参数与电流产品说明：https://web.mit.edu/space/www/voyager/voyager_data/voyager_data.html
- SPDF PLS 日鞘：https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma_cdaweb/hires_plasma/heliosheath/
- SPDF PLS 小时：https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma/hour/
- SPDF PLS 日值：https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/plasma/daily/
- SPDF PWS：https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager2/wave_spectra_pws/

## 程序与运行

正式代码目录：
C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Plasma_Wave_Overview_fu/

~~~matlab
addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Plasma_Wave_Overview_fu');
Run_V2_Plasma_Overview;
Validate_V2_PWS_Density_Extension;
~~~

默认入口先核查并下载官方文件，再直接读取原始科学数据。已在当前任务核查完成时，可运行 Run_V2_Plasma_Overview(false) 复用源文件清单；科学数据仍从原始 CDF 和官方原生 CSV 读取。

新增配套程序：
- Sync_V2_PWS_Density_Archive.ps1：官网下载、版本保留、源清单与校验。
- V2_Read_Official_PWS_Density.m：读取标签和原生 CSV，每条记录保留。
- Validate_V2_PWS_Density_Extension.m：验证交付图与所有原始记录。

第二个参数控制是否纳入 PWS：Run_V2_Plasma_Overview(true,false) 可重新生成以前截止 2024 年的 PLS 六面板。旧验证器 Validate_V2_Plasma_Overview 对应那一版本。

## 输出和验证

结果目录：
C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V2_Plasma_Overview/With_Official_PWS_Density/

- V2_1990_20250901_daily_6panels_PWS.png / .pdf / .fig：新图。
- V2_daily_six_panel_values.csv：原太阳黑子数和 COHO/PLS 日值及样本数；PWS 未写入质子密度列。
- V2_PWS_density_native_records.csv：完整 5 个 PWS 原始观测点。
- coverage.csv：七种源参数的覆盖；两个密度源共用一个显示面板。
- abc_preservation.csv：与旧 a/b/c 的数值核对。
- missing_date_intervals.csv：原 COHO/PLS 各列的缺测段；PWS 的稀疏覆盖由原始记录表表示。
- V2_plasma_overview_audit.mat：源 CDF 路径/变量/元数据/哈希、原始 COHO 记录、MAG 补充审计、官方 PWS 全表及标签、原生点、处理规则、代码哈希和旧文件保护。
- final_validation.json：源文件 MD5、5/5 点逐点同值同时刻、所有新增点位于轴内、未填造 PLS 参数、旧图及旧统计哈希不变。
- matlab_final_20260915.log：最终绘图与验证日志。

数值审计与图内曲线一一对照，同时保留所有来源和方法说明，不在图中添加小字。
