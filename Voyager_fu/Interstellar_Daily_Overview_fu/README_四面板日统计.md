# Voyager 1 日球层顶外四面板

## 目录与运行方式

- 程序：`C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu/`
- 结果：`C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/`
- 原始 CDF：`Z:/SPART-WORK/Data/Voyager/` 下按航天器、仪器、产品和年份分类。
- 官方源文件目录核查：`Z:/SPART-WORK/Data/Voyager/source_verification/V1_daily_overview/official_archive_inventory.mat`

在 MATLAB 打开并运行 `Overview_V1_Interstellar_Daily.m`。此脚本参照 MMS_fu 中 `Overview_lmn.m` 和 `SDCFilenames.m` 的人工写法，集中设置 CodeDir、ParentDir、OutputDir，按路径、参数、读取计算、绘图分节。统计和已核验的 CDF 接口保留为函数，方便复用与检查。

```matlab
run('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu/Overview_V1_Interstellar_Daily.m');
```

主函数 `Run_V1_Interstellar_Daily_Overview` 和独立无参调用的 `Plot_V1_Interstellar_Daily_Overview` 都会重新读取原始 CDF。`SyncArchive=false` 表示读取本地原始 CDF；`true` 会先查询 NASA 目录并下载缺失文件。数值输入不依赖已导出的 CSV 或 MAT。源目录清单 MAT 只提供下载核查信息，不提供磁场或通量样本。

## 科学定义

a：参考图同源 COHO 一小时总磁场强度 ABS_B 的 UTC 日算术平均；如文件提供 F 则沿用正式代码的 F 回退。使用源标量强度，不计算日均矢量的模。

b/c：同一组原始 `protonFlux1_LECP` 小时样本，分别求 UTC 日算术平均和日中位数。在线性通量上计算；偶数个样本时中位数取中间两数均值。

d：P1（官方氢通道索引 10）的日均通量 S1+S2+S3+S5+S6+S7。先逐扇区日均，再求和；任何所需扇区缺测则总和 NaN。没有角度权重、立体角积分、背景扣除、归一化、PA 分箱或对数值求和。求和只需扇区通量，不应用 PADUsable 姿态掩码。

复用 `Voyager_Read_CDF_Product`，通过 IRFU `dataobj` 直接读取原始 CDF，沿用 FILLVAL、VALIDMIN、VALIDMAX 和极大填充值筛选。日窗 `[D 00:00,D+1 00:00)`，图上每天中午一个值。有限源有效零值保留，缺测不作为零；没有最低样本数或覆盖率门限、插值、补值、平滑或额外异常值阈值。

d 复用 `Case1_Apply_L1_Fallback(...,'day','l1_first')`：舍弃 L1 DeltaT<0 记录，按原始 Epoch 所属日逐扇区求有限非负率（含零）的算术平均。S1–S7 日均率全部正且完整时，整日使用 L1，通过既有公式 `J=mean(R)/(0.44*(1.78-0.57))` 换算；否则保留原始 L2 记录，不逐扇区混合级别。同一天有多个保留 L2 记录时按本次用户要求逐扇区日均。完整参与记录、DeltaT、被替换 L2 payload、源文件索引和 Sigma 在年度审计中保留。该历史转换尚未验证等价于官方 L2 校准。

显示能段沿用用户指定的 `0.57–1.78 MeV`；源 P1 EnergyRange（此前核验为 0.57–0.89 MeV）和其它源元数据保持原样。没有因标签进行额外重标定。d 与 b/c 的源产品定义和转换通路不同。

本任务不需要新的姿态计算，没有引入新的扫描面假设。图中未要求误差棒；原始扇区 Sigma 及已批准的 L1 独立误差传播结果保存在审计中。

## 当前显示规则

仅显示 `[2012-08-25 00:00,2021-12-17 00:00)`，即包含 2021-12-16 全日。a 线性轴，b/c/d 对数轴；点以 0.4 pt 细线连接相邻有效日，缺测处断开。

对数图显示数组中的非正值置为 NaN；数值结果保留原值。本次 d 有 9 个非正日值无法显示。没有用正小数替换零值，也没有跨这些点插值。完整日统计保留至公开 COHO 源末日 2025-06-30，图上日期范围按用户要求截取。

## 输出与验证

结果目录保存 PNG、PDF、FIG、每日统计 CSV、覆盖 CSV、主审计 MAT、年度扇区来源 MAT、显示审计 MAT、验证结果和 logs。原始 CDF 不复制到结果目录；CSV/MAT 仅用于交付、审计和验证。

主审计保留原始小时输入、CDF 路径和记录索引、SHA256、源变量/全局元数据、统计定义、覆盖和代码 SHA256。年度扇区审计保留完整 L1/L2 来源与替代决策。

`Validate_V1_Daily_Overview` 用独立 accumarray 分组复核日均和日中位数，并与参考图 Case1-S01-L01 的 29 条日均扇区记录逐项比较。

2026-09-07 目录修正：此前误放于 Z 盘 daily_overview 的 18 个结果文件已迁入结果目录，源目录清单归入 source_verification；每个迁移文件的 SHA256 已核对。`path_relocation_20260907.csv` 保存原路径、新路径和哈希。两份历史重复 MATLAB 代码经哈希一致性检查后清理，保留 Voyager_fu 中的正式版本，记录见 `duplicate_code_cleanup_20260907.csv`。

官方目录核对为 159 个 COHO 月 CDF、10 个 L1 年 CDF、10 个 L2 日均年 CDF。日统计共 4693 日；有效日数 a 为 4299，b/c 各 3035，d 为 3272。末次有效日分别为 2025-06-30、2023-09-02、2021-12-16。

## 官方来源

- COHO：<https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/coho1hr_magplasma/>
- LECP：<https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/particle/lecp/final-cdf/>
- 起始日依据：<https://www.nasa.gov/news-release/nasa-spacecraft-embarks-on-historic-journey-into-interstellar-space/>

## 与磁场的相关系数

入口 `Correlate_V1_Interstellar_Daily` 从原始 CDF 重读并沿用上述日统计，计算当前图时段 2012-08-25 至 2021-12-16 的 a-b、a-c、a-d 三组 Pearson 相关系数。每组分别按同日有限值配对，使用原始日值并保留源有效零值；不取对数、不去趋势、不平滑、不移时、不增加异常值筛选。

结果保存在本任务结果目录的 `correlation` 子目录，包含系数、各对有效天数、全部参与日期、未参与日期及来源审计。MATLAB corrcoef 的输出另外用中心化向量公式核对。这里只报告描述性相关，不把相邻天当独立样本进行显著性推断。

## 包含 S4 的独立版本

运行 `Run_V1_Interstellar_Daily_WithS4`，从原始 CDF 重新计算，另存至结果目录的 `with_S4` 子目录。a/b/c 与原图保持一致；d 改为 S1–S7 日均通量之和，排除 S8，七个扇区均有限才给出日值。新图保持 2012-08-25 至 2021-12-16 的显示范围、b/c/d 对数纵轴和 0.4 pt 细线，缺测处断开。原六扇区图保留。

包含 S4 图的相关系数使用 `Correlate_V1_Interstellar_Daily(true)`；结果独立保存在 `with_S4/correlation`。仍使用 2012-08-25 至 2021-12-16、同日有限值分别配对的原始值 Pearson 相关；d 使用 S1–S7，S8 排除，源有效零值保留。无参数调用继续计算原六扇区版本。

## 六扇区版本的三日平均和月平均

入口 `Run_V1_ThreeDay_Monthly_Overview` 从原始CDF重算日统计，结果另存于 `averaged_no_S4`。d始终使用S1、S2、S3、S5、S6、S7；各面板分别对已有日值做等权算术平均，c为日中位数的平均。

三日平均采用从2012-08-25 UTC开始的不重叠三日窗口；月平均按UTC自然月。范围截至2021-12-16全日，首尾不足整窗口的部分保留并记录。缺测不参与均值，不赋零，不插值；不设最低有效天数门限。相关系数为原始窗口均值上的同窗口配对Pearson，仍计算a-b、a-c、a-d。

## 三日 panel c 最新定义

用户随后要求三日图 c 改为三日中位数：在原来的不重叠三日窗口内，对原始 COHO `protonFlux1_LECP` 小时样本直接求中位数；有限源有效零值保留，缺测不参与，没有最低样本数门限。偶数个样本取排序后中间两个值的平均。

仅更新三日图可运行 `Run_V1_ThreeDay_Monthly_Overview(false,'three_day')`。该入口仍从原始 CDF 读取小时和日统计；月图结果仅作原样保留，未作为三日计算输入。a/b/d 继续是三日算术平均，d继续排除S4、S8。月图c仍为日中位数的月平均。

三日表中的 `ThreeDayP1Median` 为新 c 值，`P1HourlySampleCount_C` 为参与的有效小时记录数；旧的日中位数均值保留为 `PreviousMeanOfDailyP1Medians` 供审计。三日a-c相关系数按新c值重算，仍采用未取对数的同窗口Pearson。
