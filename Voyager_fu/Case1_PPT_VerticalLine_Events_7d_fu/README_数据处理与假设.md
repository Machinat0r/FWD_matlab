# Voyager Case 1 七天窗口图：数据、处理与假设

## 最新配置：前后各3天、优先L1

用户最新要求将45个V1事件全部改为前后各3天（共7个UTC日），三组图全部覆盖重画，优先采用L1。当前 `ContextDays=3`、`LECPSourcePriority='l1_first'`；完整L1格优先，L1不完整时保留L2。以 [三天窗口及L1优先说明](README_三天窗口_L1优先.md) 为准。下文涉及±7天、L2优先和旧计数的段落均为之前版本的历史记录。

## 2026-09-03 最新：恢复 Level-1 扇区计数率补充通路

用户已要求接入旧 Figure 4 使用的 L1 通路，并覆盖重画指定事件中缺失 PAD 的图。最新规则见 [L1 补充通路说明](README_Level1补充通路.md)。本节优先于下文仅针对 L2 的历史运行描述与统计。

完整 L2 时间格保留原 Epoch、通量及误差；缺少完整 L2 的 UTC 日/小时格，可由官方 L1 原始计数率逐扇区算术平均补充。L1 派生记录在时间格中点计算现有三维预测姿态与投掷角，采用旧程序 `J=mean(R)/[0.44*(1.78-0.57)]`，无背景扣除。该换算与官方 L2 校准的一致性没有在本轮确立，源能段元数据保持不变。每个输出行分别标记 `L2` 或 `L1_UTC_mean`，禁止单扇区混用产品；源 CDF、原始记录、样本数、原始 DeltaT 与被替换 L2 行均保留在 MAT 审计。

以下既有运行计数描述的是本次 L1 补充前的基线；本次重画统计另存输出目录的 L1 重画结果说明。

## 本轮新增要求：显示能段与五时刻 PAD

2026-09-03 用户明确要求新图 panel e/f 对应 P1 能段显示为 `0.57–1.78 MeV`，显示上限按其指定的 Mosley 参考采用。源 CDF 的 `FHDU_EnergyRange`、通量和不确定度原样保留。原始元数据中的 P1 上限为 `0.89 MeV`，当前核查的 V1/V2 官方相关元数据取值相同，其与显示能段的差异原因尚未解释。显示字段和源元数据须分开记录，不能把标签调整解释成通量重标定、能宽修正或能道合并，也不能据此确定真实仪器能段。

用户同时将第 3.3 组明确为每个 V1 事件五个并排小时 PAD 黑色点/误差棒图，具体峰值、邻点与归一化规则见第 5.5 节。原先事件中央日及 Florinski 精确原图版式的要求以本次明确版式为准。现已重绘 90 张 V1 概览图并新增 45 张五时刻图：10 个事件五格完整，2 个事件峰值可绘但邻点不足，15 个事件峰值缺少完整 PA，17 个窗口无官方 L2 扇区记录，1 个窗口没有正扇区通量。峰值缺 PA 的 15 个事件仍绘出有效邻点，其中 3 个事件还缺少部分邻点。V2 四个事件未套用未经核验的 V1 扇区几何。

最终回归：概览图 90/90 事件、1283/1283 检查；五时刻图 45/45 事件、717/717 检查。审计分别为 `Z:\SPART-WORK\Data\Voyager\voyager1\lecp\validation\epoch_anchor\epoch_plot_audits_regression_20260903_062022_662.mat` 和 `Z:\SPART-WORK\Data\Voyager\voyager1\lecp\validation\peak_pad\peak_plot_audits_regression_20260903_062157_077.mat`。

## 2026-09-03 最新时间规则与数据审计状态

已逐字节核验当前 54 个必要 COHO 月 CDF，均与 NASA/SPDF 当前原始文件一致。
新下载 22 个 LECP 原始年 CDF：V1 2013–2021、V2 2019–2020 各自日/小时产品，保留官网文件名并分类存放在 Z 盘。
并已补齐这两艘航天器相同年份的11份Level-1原始扇区计数率CDF；49个事件窗口均有下层扇区记录。尚未将计数率校准成通量，也未用它们替换L2数值。
原来的 V1 长时间段 CDF 来自 NASA 服务器生成的子集下载，保留用于比对；新入口已改为读取官网原始年文件。

LECP `Epoch` 是累计开始，`DeltaT` 是停止减开始。2016 日产品有一条七扇区有效通量记录，其元数据区间从 2016-01-07 延至 2016-12-29；该区间不能直接解释为同样时长的连续实际曝光。
用户最新明确批准：舍去 `DeltaT < 0` 的记录，其余直接按原始 `Epoch` 锚定。
当前 `AccumulationPolicy='epoch_drop_negative_deltat'` 已解除此前 `review_required` 的统一出图暂停条件；日图读取官方日产品，小时图读取官方小时产品，通量及不确定度保持原始 CDF 的有效数值，不再执行本地 UTC 日/小时通量平均。
不因跨整点、跨午夜或正累计时长较长而剔除，不拆分区间，不依据 `DeltaT` 填满时间轴，也不插值磁场或通量。零或缺失 `DeltaT` 不额外排除；既有 CDF 填充值与有效范围处理继续适用。
预测姿态在每条原始 `Epoch` 求值，磁场按该 `Epoch` 所属 UTC 日/小时匹配现有有效矢量均值。所有实际缺测仍留空，不恢复已删除的诊断阈值。
此前 Epoch 时间规则的 90 张 V1 日/小时概览图已完成，旧验证记录保留。本轮能段显示更新与五时刻图也已通过独立验证，记录见本文开头。Florinski 2008 Figure 4 出版商原图读取受限保留为历史事项，本轮直接使用用户所附参考图的明确版式。
底部仍取原来的 P1 第 10 个氢能道。源 `FHDU_EnergyRange` 记录 `0.57–0.89 MeV`，新图按用户要求显示 `0.57–1.78 MeV`；两项信息分别保留，不修改原 CDF、通道或原始数值。

## 1. 程序用途与当前状态

这套程序对应 Case1_PPT_VerticalLine_Events_7d 中的 49 个事件图，包含 Voyager 1 的 45 个事件和 Voyager 2 的 4 个事件。每个事件来自 Case1.pptx 中一条竖线，其日期已在旧程序中固化，并按最近的 UTC 日记录。现有记录没有保留 PPT 横坐标、年度轴范围和逐线换算残差，因此事件时刻只有整日精度。

每幅图覆盖 15 个完整 UTC 日：

~~~text
[事件日 - 7 天 00:00, 事件日 + 8 天 00:00)
~~~

文件夹名称中的 VerticalLine 描述事件来源。旧图配置 ShowEventBoundaries=false，图内没有额外绘制事件边界竖线。重构程序保留这一设置。

默认入口仅执行预检并生成可复现清单，不覆盖现有 PNG。第 5 节记录当前已批准的时间处理、保留的磁场均值与显示假设；新增科学处理仍须先取得用户同意。

## 2. 原始数据

### 2.1 COHO 一小时合并 CDF

正式数据位于：

~~~text
Z:\SPART-WORK\Data\Voyager\voyager1\coho\1hr\l2\merged_mag_plasma\YYYY\MM\*.cdf
Z:\SPART-WORK\Data\Voyager\voyager2\coho\1hr\l2\merged_mag_plasma\YYYY\MM\*.cdf
~~~

49 个精确半开时间窗共需 54 个唯一月份文件。旧下载清单还包含 2017-06 和 2020-10 两个恰好落在窗口右端点的月份，因此旧归档总数为 56 份。程序使用 IRFU-MATLAB dataobj 直接读取原始 CDF，不生成 CSV、BIN 或 JSON 数值缓存。

图中使用的变量：

- Epoch：CDF 记录时间。
- ABS_B：优先用于磁场模；没有有效 ABS_B 时回退到 F。
- BR、BT、BN：RTN 坐标系磁场分量，单位 nT。
- V、protonDensity、protonTemp：程序会检查这些直接观测变量。49 个图窗中它们都没有有效值，所以现有图没有速度、密度或温度面板。
- protonFlux*_LECP：三个 LECP 质子通量通道。Voyager 1 的标称能段为 0.57–1.78、3.40–17.6、22.0–31.0 MeV；Voyager 2 为 0.52–1.45、3.04–17.3、22.0–30.0 MeV。
- protonFlux*_CRS：Voyager 2 图的质子能谱。共有 18 个约 3–599 MeV 能道，具体上下边界从各 CDF 的变量标签读取。

每个 15 天窗口有 360 个规则的一小时 Epoch。磁场有效样本为每事件 51–242 个。三条 LECP 线通量在全部 49 个窗口内都有有效数据；4 个 Voyager 2 窗口均有 CRS 数据。

### 2.2 Voyager 1 LECP 官方日/小时扇区 CDF

正式数据位于：

~~~text
Z:\SPART-WORK\Data\Voyager\voyager1\lecp\1d\l2\sectored_flux\YYYY\voyager-1_lecp_lev-2-daily-avg_*.cdf
Z:\SPART-WORK\Data\Voyager\voyager1\lecp\1h\l2\sectored_flux\YYYY\voyager-1_lecp_lev-2-hourly-avg_*.cdf
~~~

当前入口直接读取官网原始年 CDF，V1 2013–2021 共 9 份日产品和 9 份小时产品，数据集为 `VOYAGER-1_LECP_LEV-2-DAILY-AVG` 与 `VOYAGER-1_LECP_LEV-2-HOURLY-AVG`。小时图不由日产品拆分或插值生成。主要变量：

- FHDU_SectoredFluxes
- FHDU_SectoredFluxUncertainties
- FHDU_Energy
- FHDU_EnergyRange
- Epoch
- DeltaT

V1 全部 9 年原始日产品共有 1870 条记录，小时产品共有 14874 条；本轮实核日产品有 3 条负 `DeltaT`，按新规则保留 1867 条，小时产品未发现负值。之后仍按各事件窗口的原始 Epoch 选择记录并计算 PAD。负值筛选不修改源文件，剔除明细保留来源与原始记录行号。

旧 NASA 生成日子集 `Z:\SPART-WORK\Data\Voyager\voyager1\lecp\1d\l2\sectored_flux\2013-2021\voyager1_lecp_p1_sectored_daily_case1_20130102_20210525.cdf` 保留作历史对照。它含 1688 个 Epoch、16 个氢能道和 8 个扇区，在原请求区间内与官网年文件逐值一致。P1 选取能量中心最接近 0.73 MeV 的第 10 个氢能道，源 CDF `FHDU_EnergyRange` 字段为 0.57–0.89 MeV。当前新图按用户指定的 Mosley 上限显示 0.57–1.78 MeV；源字段与显示值分别保存，差异尚未解释。通量单位为 Counts/(MeV-cm^2-s-sr)，程序不做几何因子率转换或因显示能段不同而重算通量。

以下为历史日产品的窗口覆盖统计，不能用于推断官方小时产品的覆盖：45 个 Voyager 1 窗口中，14 个有完整 15 天记录，14 个有部分记录，17 个没有任何该日产品的 Epoch。完全无记录的事件为：

~~~text
Case1-S03-L01/L02/L03/L05/L06
Case1-S05-L03/L04/L05/L06
Case1-S06-L01/L02
Case1-S08-L03/L04/L05/L06
Case1-S09-L02/L03
~~~

Case1-S03-L04 有一个 Epoch，但 P1 开放扇区没有有效样本。缺测区保持空白。

### 2.3 LECP 扇区指向

2026-09-03 用户明确批准使用官方预测姿态、将约 7 天间隔的姿态插值到每天，并批准使用官方手册的标称 LECP 安装轴和扇区中心。用户要求图内无需标注这些方法；本说明、代码和审计仍完整记录来源与限制。

当前默认 `PitchAngleMethod='predicted_ck'`，MATLAB 通过 IRFU 自带 MICE 直接读取原始姿态 CK，主文件为：

~~~text
Z:\SPART-WORK\Data\Voyager\voyager1\attitude\spice\ck\7d_predicted\1977-2027\vgr1_super.bc
~~~

该 CK 在 1990-04-16 以后的部分使用预测姿态，当前事件都位于此段。原始约 7 天节点之间采用官方 CK type-3 的旋转插值；依据用户最新 Epoch 规则，在每条保留记录的原始 Epoch 直接求值，不平移到每天 12:00 或每小时 30 分。`AttitudeDailyHourUTC=12` 仅为旧接口兼容保留字段，当前 Epoch 模式不使用它。零时间容差，不外推、不越过 CK 定义的缺口。所得方向为该 Epoch 的预测方向，不能称为日均或小时实测姿态。预测数据不能恢复实际滚转、姿态机动和未记录的姿态误差。

旧 CSV 只是预留接口，未对应到已生成文件。现在默认流程不需要该 CSV，不把原始 CDF 转为 CSV。具体扇区角度、方向符号、RTN 变换及来源见 [预测姿态与扇区几何说明](README_预测姿态与扇区几何.md)。

## 3. 图件面板

Voyager 1：

1. 磁场模 |B|：黑色一小时点/线，灰色 UTC 日均。
2. BR：黑色一小时值，灰色 UTC 日均。
3. BT：黑色一小时值，灰色 UTC 日均。
4. BN：黑色一小时值，灰色 UTC 日均。
5. 三个 COHO LECP 能带的原始一小时通量。
6. P1 投掷角—通量图；使用官方预测姿态和官方标称扇区中心。磁场、通量或可用姿态不足时继续留空。

本轮 V1 概览图 panel e 对应的低能通量标识与 panel f 的 P1 能段均按用户要求显示 0.57–1.78 MeV，其他 panel 和原始数值处理维持现有规则。

新增第 3.3 组为独立五格图，左至右为峰值前的第二近可用 PAD、前一条可用 PAD、峰值、后一条可用 PAD、后的第二近可用 PAD。每格横轴为 PA 0–180°，纵轴为本格最大扇区通量归一化后的值，绘制 S1–S7 的黑色点与误差棒，不做拟合或角度合并。时间间隔由实际可用源记录决定，可不等于一小时。

Voyager 2：

1–5. 与上述磁场和 LECP 线图一致。
6. COHO CRS 质子通量能谱。

## 4. CDF 读取与缺测处理

程序首先搜索了 C:\Users\Administrator\Documents\irfu-matlab-master，随后复用 dataobj 和 getv 直接读取 CDF。

每个数值变量执行以下质量处理：

1. 与 CDF FILLVAL 精确相等的值改为 NaN。
2. 低于 VALIDMIN 或高于 VALIDMAX 的值改为 NaN。
3. 绝对值大于等于 1e30 的残留哨兵值改为 NaN。
4. NaN 不按零参与平均。
5. 多个文件的数据按 Epoch 排序；LECP 完全相同的重复记录保留第一条，同一 Epoch 存在冲突数值时报告错误。
6. LECP 在事件半开窗口内舍去 `DeltaT < 0` 的整条记录，保留其来源和行号审计。其余记录不因 DeltaT 再设筛选，不改写源 CDF。

磁场和通量没有执行插值、补缺、重采样、平滑、外推或背景扣除。唯一批准的插值为官方预测姿态的旋转插值。第 3.3 组另外采用用户批准的逐记录显示归一化，原始通量与不确定度仍保留。没有调用会隐式平均/插值科学时间序列的 irf_resamp，也没有调用会把 NaN 改为零的 irf_waverage。

## 5. 平均、筛选和显示假设

以下操作会影响科学解释，均在 Case1_Config.m 中显式配置。

### 5.1 磁场均值与 PAD 时间匹配

- UTC 自然日分箱。
- 每个分量只平均该分量的有限一小时样本。
- 灰色日均点放在当天 12:00 UTC。
- 没有设定磁场日均的最低覆盖率，一个有效样本也可以形成该变量的日均。该规则继承自旧程序。
- 从粒子到达方向矢量计算投掷角时，日均磁场方向只接收 BR/BT/BN 同时有效的记录。指向文件直接提供 PA_S1_deg...PA_S8_deg 时，不要求本地日均磁场有效。
- 日内方向离散度按各有效单位方向与日均单位方向的角距 RMS 计算。
- 日 PAD 用原始 LECP Epoch 所属 UTC 日的有效磁场矢量均值；小时 PAD 用该 Epoch 所属 UTC 小时的有效磁场矢量均值，均只接收 BR/BT/BN 同时有效的源记录。
- 上述所属日/小时只用于匹配磁场，不移动 LECP 时间，也不重新平均扇区通量。对应时间格没有有效三分量磁场或均值模为零时，无法求夹角，PA 保留 NaN。

### 5.2 一小时线图

- 磁场和 LECP 线图保留原始有效值。
- 相邻有效记录间隔超过 1.5 小时时，在显示序列中插入 NaN，使曲线断开。
- 断线处理只改变连线方式，不产生新数据点。
- 纵轴范围在有效数据范围外增加 8% 显示边距。

### 5.3 Voyager 1 LECP 扇区与原始 Epoch

- 日 PAD 读取官方日平均 CDF；小时 PAD 读取官方小时平均 CDF。每条记录的时间直接使用原始 Epoch，通量和不确定度直接取该源记录，不执行新的时间平均或独立误差传播。
- 时间选择仅为事件半开窗口内的原始 Epoch 及 `~(DeltaT < 0)`。不设正时长上限，不排除跨整点或跨午夜记录，不将长记录扩展成多个时间点。DeltaT 为零或缺失时也不新增排除条件。
- S1–S7 作为开放视场。
- S8 是受遮挡/校准相关扇区，只保存为背景诊断量。
- LECPBackgroundMode='none'，不从 S1–S7 中减去 S8。
- 样本数和扇区有效性仅作为诊断量导出，不参与门限筛选。源记录按扇区分别保留，不以缺失误差排除有效通量；不确定度缺失时保留 NaN。
- 当前默认用原始 Epoch 的预测扇区中心方向与该 Epoch 所属 UTC 日/小时的有效磁场均值计算 alpha=acos(B_hat dot u_particle)。这是官方平均通量在 Epoch 处配对出的角度，不等同于逐曝光投掷角的时间平均；也未按各扇区实际曝光时间加权。
- 官方标称 S1 中心相对 HGA 望向为 22.5°，随后按已核验的右手方向每扇区增加 45°。本体安装平面经完整姿态转换到 RTN，未施加 RTN 的 u_N=0。
- 当前 P1 属于低能前孔径，u_particle=-u_look。扇区中心线用于代表有限视场，未积分仪器角响应，也未引入未知实际安装偏差的修正。
- 通过时间选择后，投掷角图唯一的记录有效性条件为：S1–S7 每个扇区均具有有限正通量和有限投掷角。
- 投掷角跨度、磁场方向 RMS 和中位相对不确定度仅作为诊断量导出，不参与筛选。
- 日/小时概览时间色图中，投掷角中心相差不超过 2° 的点按算术平均合并。该合并仅用于色图显示；原扇区值保留在可选派生表中。新增五时刻黑色点图不执行该合并。
- 时间着色块以原始 Epoch 为中心，日产品使用 1 天、小时产品使用 1 小时的标称显示宽度。该宽度仅承担图形显示，不代表该记录的实际曝光区间，也不按 DeltaT 拉伸或补齐缺口。

2° 投掷角中心合并仅为概览时间色图的显示假设。样本数、覆盖率、投掷角跨度、磁场方向 RMS 和相对不确定度只保留原始诊断值，不设置阈值。

### 5.4 对数色标

- 只使用有限且大于零的通量，颜色值为 log10(flux)。
- 自动色限取第 1 和第 99 百分位，向外舍入到 0.25 dex。
- 色限最小跨度为 1 dex。
- Voyager 1 投掷角图按事件计算色限。
- Voyager 2 CRS 在当前所选事件集合上计算共同色限。重叠窗口中的同一原始样本可能被重复计入百分位。
- Voyager 2 CRS 记录映射到一小时时间格；CDF 标签表明这些通量可能代表六小时产品，所以着色块的显示宽度需谨慎解释。

### 5.5 第 3.3 组：峰值与相邻五时刻小时 PAD

本节为用户最新明确批准的处理，限当前已经核验几何的 V1 P1。数据直接读取官方小时 CDF；负 DeltaT 剔除、原 Epoch 锚定、原通量与不确定度读取、预测姿态及磁场匹配均沿用上述规则。

1. 时间范围固定为当前事件的 `[事件日 - 7 天 00:00, 事件日 + 8 天 00:00)`。不得搜索窗口外的记录，也不把官方日产品插值为小时值。
2. 对窗口内每条保留记录计算峰值指标 `M(t)=max(J_S1(t),...,J_S7(t))`，参与指标的通量须为有限正原始值，S8 不参与。取整个窗口 `M(t)` 最大的记录为中心峰值；并列取最早原始 Epoch。峰值搜索不预先限定 `PADUsable`，避免因 PA 缺测而改选较弱峰值。
3. 用户明确同意邻点跳过缺测：从峰值之前选择时间上最近的两条 `PADUsable` 记录，从之后选择最近的两条 `PADUsable` 记录。`PADUsable` 沿用 S1–S7 全部有限正通量和有限 PA 的唯一记录门槛。按时间顺序排为五格，不增加间隔、覆盖率、角跨度、磁场 RMS 或不确定度门限。
4. 峰值自身无有效 PAD 时中格留空，并保留已识别峰值的 Epoch 和状态；不换用较弱峰值。任一侧不足两条可用邻点时，对应位置留空。窗口内完全没有可选峰值时五格留空。缺测不填零，不补值。
5. 每条可绘记录独立取其 S1–S7 原始通量最大值 `Jmax(t)`，显示 `J_s(t)/Jmax(t)`。五个时刻各用自身分母，不共用窗口峰值作分母；原始七扇区通量及每格分母写入审计。
6. 误差棒显示 `sigma_s(t)/Jmax(t)`，分母按固定值处理，未传播最大值估计的误差及其与各扇区的相关性。因此它只表示按固定显示尺度缩放后的源误差，不代表完整比值不确定度。源 sigma 缺失时该误差棒缺失，不因此剔除有效通量点。
7. 每格显示七个扇区的黑色点及已有误差棒，横轴 PA 固定 0–180°。不拟合，不连入推测曲线，不把接近的角度合并；没有背景扣除、有限视场反卷积或新科学质量筛选。
8. 五格使用实际原始 Epoch，所选邻点允许不等时间间隔。审计记录源 CDF 与源行号、时间窗、峰值指标与并列规则、峰值状态、各位置源时间、缺失位置、原始 Flux/Sigma/PA、归一化分母及结果，以及源能段元数据和独立显示能段。
9. 同一张图的五格共用纵坐标范围，并自动包含已有误差棒。较大的对称误差棒可能延伸到零以下，这不表示测得负通量，也不触发数据筛选。图中文字时间显示到秒；MAT 审计保留原始 Epoch 精度。

## 6. 旧 Voyager 1 第六面板的限制

现有目标目录中的旧 Voyager 1 PNG 使用了固定的 8 个 45° 扇区中心、扫描面平行 RT、u_N=0 和固定 clock angle。该方法没有逐时 SEDR/AACS/NAVMAG 姿态，BN 只通过磁场模进入分母。旧第六面板只能视作名义二维 RT 几何展示，不能用于真实三维投掷角结论。

旧代码生成的 28 份投掷角 CSV 共 420 个日记录，352 日通过旧门限；17 个事件没有 LECP Epoch。旧 manifest 仍把全部运行标为 ok，且给部分无数据事件填写了实际不存在的 CSV 路径。重构版把派生表默认关闭，并在报告中记录完整 CDF 路径。

## 7. 可复现输出

默认预检会在图件目录新增：

- Case1_event_catalog_refactored.csv
- Case1_event_catalog_V1_refactored.csv
- Case1_event_catalog_V2_refactored.csv
- Case1_source_manifest_refactored.csv
- Case1_data_processing_assumptions_refactored.txt
- README_数据处理与假设_refactored.md

完整源路径、产品层级、时间分辨率、文件大小、使用该文件的 EventID、缺测政策、有效性条件和科学设置均写入这些文件。

新图文件名包含 `predictedCK_1d` 或 `predictedCK_1h`，并用 `nativeCDF_Epoch` 标记本次时间方法，与旧图分开保存，默认不覆盖既有图。原生 MATLAB 审计 MAT 按 `Z:\SPART-WORK\Data\Voyager\voyager1\lecp\1d` 或 `1h` 下的 `derived\pitch_angle` 分类存放。审计保留原始 Epoch、DeltaT、源 CDF 与源记录行号、负 DeltaT 剔除清单、所配磁场时间格、各扇区通量/不确定度/PA/完整 RTN 方向，以及旋转矩阵、几何定义和实际内核 SHA256。缺少 LECP 记录或 PA 所需输入的图仍留空，不伪造扇区数据。旧日中点验证图与旧审计仅作历史结果保留。

入口按 `PADCadence` 将新图和运行清单分别放入 `cfg.OutputRoot` 下的 `Epoch_daily` 或 `Epoch_hourly` 子目录。默认图件根目录为 `C:\Users\Administrator\Documents\Recovery-Work-Voyager_betatron\Case1_PPT_VerticalLine_Events_7d`；旧根目录已有 PNG 保留。

五时刻图位于 `C:\Users\Administrator\Documents\Recovery-Work-Voyager_betatron\Case1_PPT_VerticalLine_Events_7d\Peak_PAD_5times_hourly`。其派生数据和审计按 `Z:\SPART-WORK\Data\Voyager\voyager1\lecp\1h\derived\peak_pitch_angle\2013-2021` 分类保存。小时组运行清单中的 `PeakPADFigureFile`、`PeakPADAuditFile`、`PeakPADStatus`、`PeakPADEpochUTC`、`PeakPADSelectedCount` 和 `PeakPADPlottedCount` 记录每个事件的实际结果。输出目录中的 `README_五时刻图_结果.md` 汇总本轮结果。

## 8. 运行方式

只做预检：

~~~matlab
cd('C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Case1_PPT_VerticalLine_Events_7d_fu')
result = Run_Case1_PPT_VerticalLine_Events_7d;
~~~

首次在其他机器运行前，可执行 `Case1_Download_Attitude` 获取缺少的原始官方内核；已有源文件不会覆盖。`Case1_Test_Predicted_Attitude` 执行全时段和历史几何数值回归，并把测试报告保存到 Z 盘 attitude/spice/validation。MICE 读取需空的内核池；建议使用独立 MATLAB 进程，程序不会清空调用者已有内核或程序变量。

使用当前已批准的第 5 节处理生成日产品图：

~~~matlab
result = Run_Case1_PPT_VerticalLine_Events_7d( ...
    'RunPlots', true, ...
    'PADCadence', 'day', ...
    'Overwrite', false, ...
    'Visible', false);
~~~

小时产品图直接选择官方小时 CDF：

小时入口默认 `ExportPeakPAD=true`，同时生成第 3.3 组五时刻图；设为 false 可只生成小时概览。日入口不生成五时刻小时图。以下命令仅选择目前已核验扇区几何的 V1 事件：

~~~matlab
catalog = Case1_Event_Catalog;
v1IDs = catalog.EventID(catalog.Spacecraft == 1);
result = Run_Case1_PPT_VerticalLine_Events_7d( ...
    'RunPlots', true, ...
    'PADCadence', 'hour', ...
    'ExportPeakPAD', true, ...
    'EventIDs', v1IDs, ...
    'Overwrite', false, ...
    'Visible', false);
~~~

将 Overwrite 设为 true 会允许覆盖同名旧图。建议先选择一个 EventID 做单图验证：

~~~matlab
result = Run_Case1_PPT_VerticalLine_Events_7d( ...
    'RunPlots', true, ...
    'PADCadence', 'day', ...
    'Overwrite', false, ...
    'EventIDs', "Case1-S01-L01");
~~~
