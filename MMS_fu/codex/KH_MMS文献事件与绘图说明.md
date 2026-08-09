# MMS 磁层顶 Kelvin–Helmholtz 事件清单与 overview 绘图说明

## 1. 范围与结论

本项目把已发表文献中的 MMS 磁层顶 Kelvin–Helmholtz instability/wave/vortex（KHI/KHW/KHV）时段整理为一个去重清单，并针对本地 MMS burst 数据为 MMS1–MMS4 分别制作九面板 overview。

当前 kh_event_catalog.m 中共有 **89 个互不重复的时段**：

- Rice et al. (2022) 的 45 个 clear KHI encounters 构成基础清单，CatalogID 为 R01–R45。
- 另有 44 个不与基础清单重复的已发表 MMS 时段，CatalogID 为 X01–X44，主要来自 Kavosi et al. (2023) 的补充事件表以及若干专门的案例研究。
- 后续论文若对 Rice 清单中的同一时段给出了更强的卷起涡证据，代码会在原 R 行中追加 Source/DOI、提升 Tier 或增加 PreferredStartUTC/PreferredEndUTC，而不是再加入重复行。
- 所有行按 StartUTC 排序后重新分配 EventID（KH001、KH002……KH089）；CatalogID 保留原始 R/X 身份。因此 EventID 是绘图和输出主键，CatalogID 是清单溯源键。

其中 **A_rolled_up 共 18 个优先时段**：7 个由 Rice 基础行提升，另有 11 个来自额外案例研究。priority 模式和 smoke 模式只选择 Tier 以 A_ 开头的行；all 模式先处理这 18 行，再处理其余行。2015-09-08（KH001）保留为明确 KHI，但根据 Wilder et al. (2023) 的事件级判据降为 B_clear_KHI：该文把它列入 confirmed KHI，却指出它是 12 例中唯一没有 low-density high-speed/nonlinear roll-up 特征的一例，因而不能当作成熟卷起涡。

这里的“89 个时段”不等于“89 个都已被证明是成熟卷起的 KH 涡”。KHI/KHW 目录项、协调观测表中的 KHV 标签、线性波动和论文明确确认的非线性卷起结构，证据强度并不相同，必须结合 Tier 与 Confidence 使用。

## 2. 主要文献来源

核心清单和高置信度案例的主要来源如下：

| 文献 | 在清单中的作用 | DOI / 链接 |
|---|---|---|
| Rice et al. (2022) | 45 个 clear KHI encounters 的基础清单 | [10.1029/2021JA029685](https://doi.org/10.1029/2021JA029685) |
| Kavosi et al. (2023), Nature Communications | 补充事件表中的 MMS KHW 时段；与 Rice 重叠者用于合并证据 | [10.1038/s41467-023-37485-x](https://www.nature.com/articles/s41467-023-37485-x) |
| Wilder et al. (2023) | 以 low-density high-speed（LDHS）特征加强若干 rolled-up KHI 案例的证据，并提供额外强案例 | [10.1029/2023JA031583](https://doi.org/10.1029/2023JA031583) |
| Hwang et al. (2022) | 磁层顶 KHV 案例和协调观测表时段 | [10.3389/fspas.2022.895514](https://doi.org/10.3389/fspas.2022.895514) |
| Kieokaew et al. (2020) | 2017-05-05 晨侧非线性卷起 KH 涡 | [10.1029/2019JA027527](https://doi.org/10.1029/2019JA027527) |
| Hwang et al. (2020) | 2017-05-05 案例的补充物理证据 | [10.1029/2019JA027665](https://doi.org/10.1029/2019JA027665) |
| Blasl et al. (2022) | 2017-09-23 南向 IMF 条件下的 KH 波与卷起涡 | [10.1063/5.0067370](https://doi.org/10.1063/5.0067370) |
| Li et al. (2023) | 2020-06-26 非线性 KHW、重复磁层顶穿越及重联案例 | [10.1029/2023GL105539](https://doi.org/10.1029/2023GL105539) |
| Gurram et al. (2026) | 2022-04-14 磁暴期间早期非线性/卷起 KH 波 | [10.3389/fspas.2026.1788081](https://doi.org/10.3389/fspas.2026.1788081) |

清单还引用 Nykyri et al. (2021)、Tang et al. (2018)、Lu et al. (2019)、Hou et al. (2026)、Guo et al. (2026) 和 Settino et al. (2026) 等论文。逐事件的完整作者、DOI 和证据描述均保存在 Source、DOI 与 Confidence 字段中；应以导出的 CSV 为最终逐行溯源入口。

## 3. 清单字段与证据层级

MMS_KH_published_event_catalog.csv 包含以下字段：

| 字段 | 含义 |
|---|---|
| EventID | 按开始时间排序后的唯一绘图编号 KH001–KH089 |
| CatalogID | 原始基础/扩展编号 R01–R45 或 X01–X44 |
| StartUTC、EndUTC | 文献事件的宽时段边界，UTC |
| DurationMinutes | EndUTC − StartUTC，单位 min |
| PreferredStartUTC、PreferredEndUTC | 文献给出或人工整理的优先绘图子时段；没有时为 NaT |
| Side | dawn、dusk 或 unknown |
| Source | 支持该行的论文；同一时段可合并多篇来源 |
| DOI | 与 Source 对应的 DOI；多篇用分号分隔 |
| Tier | 证据层级 |
| Confidence | 对证据性质、非线性程度或替代解释的文字说明 |

代码中的层级及行数为：

| Tier | 行数 | 解释 |
|---|---:|---|
| A_rolled_up | 18 | 已发表案例研究明确报告非线性/卷起结构，或以多星、LDHS 等特征强力支持 rolled-up KHI/KHV |
| B_clear_KHI | 36 | Rice 基础清单中的 clear KHI encounter，但并非每一行都确立了成熟卷起状态 |
| B_clear_KHW | 27 | 已发表目录中的 MMS KHW 时段，主要来自 Kavosi 补充表 |
| B_KHV_table | 5 | 协调观测表标注为 magnetopause KHV；部分为围绕论文列出时刻设置的工程绘图窗口 |
| C_probable | 1 | KHI 一致的边界层案例，但论文讨论了压力扰动等替代驱动 |
| C_linear | 2 | 报告为线性或未充分卷起的 KH 结构 |

这些 Tier 是证据分类，不是统计概率。尤其不能把 B 或 C 行在没有进一步事件级分析时写成“已确认成熟涡”。

## 4. 代表性四星共同 burst 窗口

kh_find_burst_window.m 不读取 CDF 内容，而是索引本地文件名中的 burst 文件起始时刻。它针对四颗卫星的四类必要产品建立 16 个“卫星 × 产品”集合：

1. FGM brst L2；
2. EDP brst L2 dce；
3. FPI brst L2 des-moms；
4. FPI brst L2 dis-moms。

窗口选择逻辑如下：

1. 检索事件开始前 3 min 至事件结束之间的本地文件，并优先寻找 16 个集合都具有相同文件起始时刻的严格交集；成功时 Mode 为 strict_4sc_all_products，CoverageOf16 为 16。
2. 若严格交集为空，则合并所有候选时刻，选择覆盖“卫星 × 产品”组合数最多的时刻，Mode 记录为 best_effort_Nof16，CoverageOf16 记录 N。
3. 把相邻起始时刻间隔不超过 180 s 的候选连成连续 run，并把最后一个文件起始时刻向后估计 120 s 作为可用尾端。
4. 有 PreferredStartUTC/PreferredEndUTC 时优先保留与该子时段的重叠；否则以完整事件时段为目标。run 的评分先强烈奖励与目标时段有重叠，再综合重叠长度、可用持续时间和距目标中心的距离。
5. 最终 PlotStartUTC–PlotEndUTC 最长 20 min。若本地文件索引为空，则回退到优先时段，或从事件起点开始截取最多 20 min；扫描器异常时 run_kh_overviews.m 也执行同样的 fallback。

因此，图中时间轴是“适合四星共同产品比较的代表性 burst 子窗口”，不一定覆盖论文给出的整个 KH 事件。每次实际选择会写入 selected_burst_windows_MODE.csv，其中同时保留事件宽时段、绘图子时段、Mode、CoverageOf16 和 EmptyProducts，便于审计。窗口选择依据文件名起始时刻而非逐条读取 CDF Epoch；最后一个文件向后 120 s 只是估计，所以即使状态为 ok，窗口尾部也可能保留真实空白。ok 的含义是九类面板在窗口中都至少有有限数据，不代表整个横轴连续覆盖。

同一绘图窗口可能由多个 burst 文件拼接。折线绘制前会检测相邻采样的时间断点；若间隔超过 1 s，或超过本底采样间隔的 10 倍，则在断点处插入 NaN，使真实数据缺口保持空白，而不是用直线跨接缺口两端。

## 5. 九个绘图面板、坐标系与公式

每个事件按 sc = 1:4 循环，目标是分别输出 MMS1、MMS2、MMS3、MMS4 图。九个面板为：

1. **磁场**：B_x、B_y、B_z 与 |B|，GSM，单位 nT。
2. **电场**：E_x、E_y、E_z；从 EDP 的 GSE 数据转换为 GSM，单位 mV/m。
3. **数密度**：电子 n_e 与离子 n_i，单位 cm⁻³。
4. **离子流速**：V_ix、V_iy、V_iz；从 GSE 转换为 GSM，单位 km/s。
5. **电子流速**：V_ex、V_ey、V_ez；从 GSE 转换为 GSM，单位 km/s。
6. **粒子矩电流密度**：J_x、J_y、J_z，单位 nA/m²。
7. **能量转换**：J·E 与 J·(E + V_e × B)，单位 pW/m³。
8. **电子全向能量通量谱**：FPI DES moments 中的 energyspectr_omni；纵轴电子能量，单位 eV。
9. **离子全向能量通量谱**：FPI DIS moments 中的 energyspectr_omni；纵轴离子能量，单位 eV。

粒子矩电流密度采用

J = e (n_i V_i − n_e V_e)。

代码先把 n_i、V_i、V_e 重采样到 n_e 的时间轴，再使用离子密度 n_i 计算离子电流；这避免了原 Overview_download.m 中用 n_e 代替 n_i 的隐含准中性近似。这里没有直接使用会外推的 irf_resamp，而是按每个原始产品的连续时间段分别插值：段外不外推、跨 burst 缺口不插值，只有 n_e、n_i、V_i、V_e 都有共同有限支撑时才保留 J。E、B 与 V_e 到 J 时间轴的重采样采用同一规则，J·E 与 J·E′ 也只在所有输入共同有效处保留。由于输入密度为 cm⁻³、流速为 km/s，代码使用 e × 10¹⁸ 的换算因子直接得到 nA/m²。

电子参考系电场为

E′ = E + V_e × B。

在代码单位中，V_e[km/s] × B[nT] 乘 10⁻³ 后与 E[mV/m] 相加。随后绘制

- J·E；
- J·E′ = J·(E + V_e × B)。

nA/m² 与 mV/m 的乘积正好对应 pW/m³。这里绘制的是这两个点积量；它们不应在没有额外项与物理检验时被自动等同于更完整的耗散判据。

电子/离子谱的纵轴能量取自 CDF 的 DEPEND_1，颜色数据取相应 CDF 变量的原生 omni energy-flux 数值。当前函数没有在图上另写颜色条的物理单位，因此定量引用通量时应回查 MMS/FPI 产品元数据。

## 6. 运行入口

三个入口脚本都会添加 IRFU、MMS_fu 和 codex 路径，并初始化本地数据库 Z:\SPART-WORK\Data\MMS：

| 入口 | 实际调用 | 用途 |
|---|---|---|
| run_kh_smoke.m | run_kh_overviews('smoke',1) | 只重画按时间排序后的第一个 A_rolled_up 事件，共尝试四星四图；用于验证代码和版式 |
| run_kh_priority.m | run_kh_overviews('priority',inf) | 处理全部 18 个 A_rolled_up 优先时段 |
| run_kh_all.m | run_kh_overviews('all',inf) | 处理完整 89 时段，A_rolled_up 优先，随后处理其余层级 |

priority 和 all 模式具有续跑机制：只有文件名与本轮 PlotStartUTC 精确对应、PNG 大于 50 kB、且生成时间不早于当前绘图函数时才标为 existing；旧窗口、旧代码或疑似半写入文件都会重画。smoke 模式会主动重画四幅测试图。也可直接调用 run_kh_overviews(mode,maxEvents) 限制本轮事件数。

## 7. 最终运行结果与输出文件

最终 all 模式已完成 **89 个事件 × 4 星 = 356 张 PNG**，没有绘图级失败：

| 状态 | 图数 | 含义 |
|---|---:|---|
| ok | 285 | 九类面板在窗口中都至少有有限数据；不保证整个横轴无时间缺口 |
| partial | 31 | 成功保存，且仍有部分科学面板可用；缺失面板明确显示 unavailable，Warnings 记录缺项 |
| no_data | 40 | 选定已发表时段内四类 burst 源产品均无目标窗数据；保存审计图但不计为可用 overview |
| failed | 0 | 无绘图或保存失败 |

按事件计，49 个事件四星全部为 ok，30 个事件四星均有可用图但其中一或两星为 partial，另 10 个事件（KH017、KH018、KH027、KH030、KH036、KH056、KH059、KH060、KH062、KH063）四星均为 no_data。也就是说，**79/89 个事件取得了四颗卫星各一张可用 overview，共 316 张可用图**；其余 40 张是明确的无源数据审计图。

对全部 31 个 partial 与 40 个 no_data 行都用 MMS SDC public `file_names` API 做了目标窗口核对。优先批次与三轮扩展审计均表明：相关产品在公共源目标窗口也无文件（部分项目连续三日文件数为 0，另一些是当天有 burst 但目标窗口为空），**可补下载项目为 0**。因此没有把相邻时段文件误当作目标事件下载，也没有用跨缺口插值补造物理量。

所有结果写到 C:\Users\Administrator\Documents\KH：

| 输出 | 内容 |
|---|---|
| MMS_KH_published_event_catalog.csv | 完整 89 时段清单；每次运行入口都会写出 |
| figures\KHxxx_yyyymmdd_HHmmss\KHxxx_PLOTSTART_MMSn.png | 按事件目录保存的四星 overview PNG，200 dpi |
| selected_burst_windows_MODE.csv | 每个事件实际采用的代表性 burst 子窗口及 16 组合覆盖信息 |
| run_status_MODE.csv | 每个 EventID × Spacecraft 的 PlotStartUTC、PlotEndUTC、OutputFile、Status 与 Warnings |
| summary_MODE.csv | 本模式的请求数、可用数、complete/partial/no_data/reused/failed 汇总 |
| progress_MODE.txt | 当前/最近处理到的事件、卫星、状态和 UTC 更新时间；它是覆盖更新的进度快照，不是完整历史 |
| MATLAB_run_MODE.log | MATLAB diary 运行日志 |
| missing_products_public_audit_priority.csv | 优先批次缺失 FPI 产品在 MMS SDC 公共 API 的逐项审计；区分整日无文件与事件窗口附近无文件 |
| MMS_KH_event_result_inventory.csv | 89 个事件逐行的四星完整/部分/无数据计数与图目录 |
| MMS_SDC_public_missing_audit_all.csv | all 模式全部 71 个异常图的缺项、公共源结论与下载结论 |
| audit\SDC_audit_round1.md～round3.md | 公共 SDC 查询的精确窗口、相邻文件时刻、三日文件数与逐轮结论 |

MODE 为 smoke、priority 或 all。状态表按每颗卫星逐行更新，即使某颗卫星失败也会继续处理其余卫星和事件。九类面板齐全记为 ok；成功导出但部分变量缺失记为 partial；九类都无有限数据记为 no_data；绘图级异常才记为 failed。最终可用图片数量应以对应的 summary_MODE.csv 和 run_status_MODE.csv 为准，而不是按“事件数 × 4”预先假定。

## 8. 限制与使用注意

- 89 行是文献时段的人工整理与去重结果，并非自动事件识别算法的输出；新增论文或修订文献时间边界时应同步更新清单及溯源字段。
- 宽事件时段内可能包含多个 KH 波周期、磁层顶穿越或不同演化阶段；当前每图最多 20 min，只代表所选 burst 子窗口。
- strict_4sc_all_products 只说明本地文件名索引在四星四产品上共同覆盖，不等于所有 CDF 变量均无缺失或数据质量完全一致。
- 当前批处理只读取本地数据库；另提供 kh_download_missing_product.m，使用 MMS_fu 中已有的 SDCFilenames、SDCFilesDownload 和 SDCDataMove 定向补齐公开源端存在的 FGM、EDP、FPI DES/DIS burst L2 文件。若公共 API 在目标窗口本身没有文件，则无法用下载消除空面板，应以缺失审计、窗口表、状态表及日志为准。
- 某些 B_KHV_table 行是围绕论文表中列出时刻设置的工程窗口；C_probable/C_linear 行也不能与成熟卷起结构混为一谈。
- J 来自粒子矩，而不是四星 curlometer；J·E 与 J·E′ 的解释受时间重采样、矩数据质量、坐标转换和局地结构演化影响。
- 原始 Overview_download.m 保持不变；plot_mms_kh_overview.m 是参数化派生实现，修正了原脚本中 E + V_e × B 被 E 覆盖以及离子电流使用 n_e 的问题。

## 9. 建议的结果核查顺序

1. 先查看 summary_MODE.csv，确认可用和失败数量。
2. 在 run_status_MODE.csv 中筛选 failed、Warnings 非空和 existing 项。
3. 对照 selected_burst_windows_MODE.csv 检查 Mode、CoverageOf16、EmptyProducts 与文献优先时段。
4. 逐图确认九个面板的时间对齐、坐标系、谱图能量范围和缺失面板。
5. 对准备写入论文的具体涡事件，再回到 Source/DOI 阅读原论文，确认其非线性阶段、卷起证据及物理解释。
