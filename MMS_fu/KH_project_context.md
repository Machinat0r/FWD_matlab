# MMS 磁层顶 KH 统计：项目背景与历史记录

整理日期：2026-09-03。配套项目规则：C:\Users\Administrator\Documents\KH\AGENTS.md。

## 来源与阅读范围

用户指定的历史任务标题：
C:\Users\Administrator\Documents\KH\figures 这个文件夹中有很多kh事件，都…

任务 ID：019fc834-87ee-7ce2-8489-d784c54c6154。
完整本地原始记录：
C:\Users\Administrator\.codex\sessions\2026\08\03\rollout-2026-08-03T23-18-33-019fc834-87ee-7ce2-8489-d784c54c6154.jsonl

本次从原始记录按时间顺序读完了全部用户与助手可见聊天文本：30 条用户消息、26 条最终回复、123 条过程说明。已与日志中的 30 条 user_message 事件逐一核对，无额外遗漏的用户消息。环境描述和插件清单不计入聊天文本数；内部推理与逐条工具日志不作为用户确认的方法。也查看了用户最后圈出红/蓝色块的 KH005 附图。

阅读范围为 2026-08-03 至 2026-08-12。下表日期使用原日志 UTC，MMS 事件时间也使用 UTC。本文件保存可承接研究的完整主题和决策脉络；精确原话、代码历史和工具输出仍可回到原 JSONL 查询，无需在 KH 复制庞大的过程记录。

当前项目 KH-FOTE 已指向 C:\Users\Administrator\Documents\KH，项目 ID 为 f06b60ba-3cf0-4491-b940-97a392019e8e。本次沿用该项目保存规则，没有另建重复任务。

## 当前研究目标

利用 MMS 磁层顶 Kelvin–Helmholtz 事件库，比较同一事件中磁场 FOTE 与流场 FOTE-V 的 screw-in/screw-out 性质，进一步研究类型对应关系与能量转换的联系，目标期刊包括 GRL / JGR: Space Physics。

原对话最初提出的 B–Vi–Ve 联合拓扑、冻结/解耦解释、De、J·E、压力–应变及分层统计均属科学研究建议。当前实际已运行的主线为 B 与 Vi 的类型和误差图。电子流、能量转换统计与论文结论尚未完成，不得把设想写成已得到的发现。

最后的重点已从逐采样点分类转向单个涡旋的整体识别。用户希望短时间内的一次磁场扰动或单涡旋具有一个最终标签，避免在同一结构内反复显示 As/Bs。

## 全部 30 条用户消息的顺序与承接状态

| 序号 | UTC 日期 | 用户请求或反馈 | 后续状态 |
| --- | --- | --- | --- |
| 1 | 08-03 | 用 KH\figures 中 MMS 事件做 FOTE/FOTE-V、screw-in/out 与能量转换统计，评估 GRL/JGR 潜力 | 建立研究方向；论文可发表性只是当时评估 |
| 2 | 08-06 | 学习 FWD_matlab\FOTE；参考 Sundry\njq 的 FOTE-V；先出一张再经确认批量跑；不要学习过度平滑 | 以 KH005 为预览，先前一度因 MATLAB 故障使用 Python |
| 3 | 08-06 | 适度平滑 5–10 点、删原 panel 5、标题只留事件号、保存 PDF、跑全部事件 | 生成 7 点平滑的 73 页旧版 |
| 4 | 08-06 | 继续 | 完成批处理收尾 |
| 5 | 08-06 | 继续 | 完成代码同步、说明及交付 |
| 6 | 08-06 | 询问 panel 3/4 的 d/L | 解释最近卫星距离和四面体六边平均尺度 |
| 7 | 08-06 | 搜索整理 Cluster 已发表磁层顶 KH 事件库 | 交付独立 Excel；属于支线 |
| 8 | 08-06 | 搜索晨昏侧 KH 涡旋方向一致性的研究 | 汇总观测/理论线索，区分宏观旋转与 screw 类型 |
| 9 | 08-06 | 全部 MMS 事件不平滑再画一份新 PDF | 生成 73 页无平滑旧版 |
| 10 | 08-07 | “-02” | 含义未澄清，不应推定参数或事件修改 |
| 11 | 08-07 | 询问减背景是否能让零点接近卫星 | 仅讨论物理含义；没有批准用减背景挽救距离筛选 |
| 12 | 08-08 | 取消远距离零点筛选，改用论文 40% 误差，删距离 panel、保留误差 panel | 生成 73 页 40% 无距离筛选、无平滑旧版 |
| 13 | 08-10 | 询问 PDF 文件夹 | 告知旧 Codex 工作目录路径 |
| 14 | 08-10 | 此任务图片以后放 MMS_fu\codex\kh_events，可分子目录 | 这是历史约定，已被 09-03 的 KH 图件目录覆盖 |
| 15 | 08-10 | B/V 用相同物理时间窗平滑至 10 s，只跑 KH005 | 原生时间轴平均后对齐，保持 40% |
| 16 | 08-10 | 程序也放到上述文件夹 | 形成旧代码归档；09-03 代码目录更新为 MMS_fu |
| 17 | 08-10 | 改为 5 s 平滑，上方加原始 B/V 两个 panel | 形成 7-panel 基础版式 |
| 18 | 08-10 | 询问同一涡旋为何出现多种零点类型 | 讨论局部梯度、分类边界、噪声、多尺度及采样相关性 |
| 19 | 08-10 | 从 IRFU 加 Poincaré index 识别及 panel | 曾实现，之后由用户要求完全删除 |
| 20 | 08-11 | 程序改为 MATLAB，尽量调用原 FOTE/FOTE-V 函数 | 完成 MATLAB 单事件实跑 |
| 21 | 08-11 | 模块化分区，并为每个模块加备注 | 建立旧 +khfote 模块；新代码风格遵循 09-03 的原函数参考要求 |
| 22 | 08-11 | prepareData 不插值，缺口直接显示空缺 | 曾改为最近实际样本；后续允许计算前 irf_resamp 对齐 |
| 23 | 08-11 | FOTE 前可用 irf_resamp “对其” | 紧接下一条更正 |
| 24 | 08-11 | 更正为 FOTE 前用 irf_resamp 对齐 | 最终按连续块线性对齐，保留缺口和原始 panel |
| 25 | 08-11 | 去掉 Poincaré 判据；底部只显示连续 >=5 s 满足 40% 的点 | 加连续区段筛选 |
| 26 | 08-11 | Poincaré 曲线也删掉 | 从计算、输出、图中删除，回到 7-panel |
| 27 | 08-11 | 用现在程序跑一个事件 | KH005 跑通，当时只要求误差持续合格 |
| 28 | 08-11 | 根据红色/蓝色圈出的整体特征识别单涡旋 As/Bs，避免短时反复变类 | 助手提出加权投票、纯度/覆盖率方案，尚未成为完整已验证算法 |
| 29 | 08-11 | 加特征值稳定度 0.5，画保留/取消持续 5 s 两版 | KH005 两版完成并留有结果 |
| 30 | 08-12 | 认可加权投票和单涡旋一个标签，但单涡旋识别并框出很难 | 最终未解决问题；最后分段方案仅由助手提出 |

2026-09-03 的新要求：数据统一 Z:\SPART-WORK\Data\MMS；代码尽量 MATLAB、先查 IRFU、模仿 MMS_fu 非 codex 原函数并放 MMS_fu；FOTE/FOTE-V 不得使用未审核方法；尽量直接读取 CDF；图和指定交付放 KH、减少过程文件。以上新要求优先于表内旧保存位置和中间实现习惯。

## 事件库与历史完成情况

原对话盘点为 89 个目录事件、356 张四星 overview PNG、76 个独立观测日期，覆盖 2015–2022。89 个目录不是 89 个独立单涡旋，原图存在也不证明四星 burst 产品齐全。

根目录已有事件索引：
- C:\Users\Administrator\Documents\KH\MMS_KH_published_event_catalog.csv
- C:\Users\Administrator\Documents\KH\MMS_KH_event_result_inventory.csv
- C:\Users\Administrator\Documents\KH\selected_burst_windows_all.csv
- C:\Users\Administrator\Documents\KH\KH_MMS文献事件与绘图说明.md

旧全库 FOTE/FOTE-Vi 批处理：73 个事件成功，10 个当时无本地数据索引，6 个缺至少一星离子矩。6 个为 KH007、KH019、KH033、KH034、KH035、KH075。这是旧批处理状态；新任务应检查当时所选时窗及当前 CDF，不能据此断言网站永久无数据。

旧全库 PDF 根路径：
C:\Users\Administrator\Documents\Codex\2026-08-03\c-users-administrator-documents-kh-figures

其下：
- output\pdf\KH_all_FOTE_FOTEV_Vi_7pt.pdf
- output\no_smoothing\pdf\KH_all_FOTE_FOTEV_Vi_1pt.pdf
- output\quality40_no_distance_no_smoothing\pdf\KH_all_FOTE_FOTEV_Vi_quality40_no_distance_1pt.pdf

这些 73 页旧版采用早期流程。之后的 MATLAB、5 s 原生时间平均、irf_resamp、特征值稳定度等改动主要在 KH005 验证，尚不能声称最终设置已覆盖全库。不要把旧 Python 图集当成当前 MATLAB 流程的全库产物。

## 最后 MATLAB 实现与已核对的原函数

历史代码目录：
C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\codex\kh_events\code\kh_fote_matlab

入口及模块：
- kh_fote_event_irfu.m / kh_fote_event.m：主流程与兼容入口；
- kh_fote_batch_irfu.m：旧批处理；
- run_kh005_matlab_irfu.m：单事件；
- run_kh005_eigenstability_compare.m：持续 5 s / 逐点对照；
- +khfote\parseOptions.m、setupEnvironment.m、loadMmsData.m、prepareData.m；
- +khfote\runOriginalFote.m、computeFoteVError.m；
- +khfote\applyQualityAndSummarize.m、plotEvent.m、writeResults.m、Contents.m。

以上用于历史结果追溯。未来新代码放在 MMS_fu，风格从 codex 外原函数学习，不能照搬历史 codex 架构作为用户偏好的样板。运行旧入口之前必须落实最新目录和最少输出要求：目前旧 parseOptions 的 OutputRoot 仍为 fullfile(pwd,'outputs')，旧 writeResults 默认写 PDF/PNG/CSV/MAT/JSON，直接运行会违背新的整理偏好。本次仅建立项目规范，没有修改这些历史程序。

原函数与风格参考已经检查：
- FWD_matlab\FOTE\FOTE_function\FOTE_Taylor_Expansion.m；
- MMS_fu\FOTE_V_paras.m、FOTE_V_overview.m 的分区、命名与 IRFU 调用；
- MMS_fu\SDCDataMove.m 的现有 MMS 分类方式；
- Sundry\njq 中 FOTE-V 程序的相关段落。

FOTE_Taylor_Expansion 返回距离结构、类型标记结构和误差结构；8 参数形式跳过内部 smooth。
历史包装读取 err1 为 eta、err2 为 xi。原函数还返回另一个实部形式的 err3，不能自行把不同定义互换。
原函数中还存在类型简化（O/X 等）和 marker/face-color 编码。未来必须保持并核对这些逻辑；单纯从特征值再造分类可能与原函数不同。

IRFU 已有：
- mission\mms\+mms\get_data.m、db_init.m，以及 mms.db_get_ts；
- mission\cluster\c_4_grad.m、c_4_j.m；
- irf\irf_resamp.m；
- c_4_poincare_index / irf.solidangle 曾使用，当前流程已取消 PI。

当前原始数据已按标准目录分类，例如：
- Z:\SPART-WORK\Data\MMS\mms1\fgm\brst\l2\2015\10\01
- Z:\SPART-WORK\Data\MMS\mms1\fpi\brst\l2\dis-moms\2015\10\01

同一 CDF 按产品归档共享，无需逐 KH 事件复制一套。

## 处理方法及版本差异

### 数据与时间

KH005：2015-10-01T18:01:24Z / 2015-10-01T18:09:00Z。
最后实际比较为 GSE 磁场、离子速度 Vi、同种粒子密度 Ni 和四星位置。Ve 只是保留的接口，不能写成已完成 B–Vi–Ve 全统计。

5 s 平滑指相同物理时间长度。旧实测 FGM 约 127.998 Hz，对应约 640 点；DIS 矩约 6.667 Hz，对应约 33 点。数值仅用于 KH005，其他数据必须从时间戳确定窗口。
历史演变依次有无平滑、7 点移动中值、无平滑、10 s 时间平均、5 s 时间平均，后续不能混用标签和结果。

用户先明确要求不插值、缺口留空，随后允许 FOTE 前 irf_resamp。因此最终规则为：原始 panel 用观测样本；分析场在原生时间轴做平均后按连续数据块 irf_resamp(...,'linear') 对齐到 MMS1 离子时刻，屏蔽外推和跨缺口插值，保留完整时间轴中的 NaN。不能恢复全局插值补齐。

### 误差和类型

历史最后误差为百分数：
- eta_B = 100 |div B| / |curl B|；
- xi_B = 100 |lambda1+lambda2+lambda3| / max_j|lambda_j|，使用原函数 err2；
- alpha_V = 100 |div(nV)| / |curl(nV)|，使用对应离子或电子密度。

磁场要求 eta_B、xi_B 均 <=40%，流场要求 alpha_V <=40%。必须保留原函数公式与单位。此处记录现有实现，不替代对原 FOTE/FOTE-V 方法的审核。

d_min 为重构零点到最近一颗卫星的距离；L 为四面体六条边平均长度。早期 d_min/L<=2 入选规则已被用户取消，距离 panel 也删除。取消距离限制后类型应谨慎解释为局地梯度推断拓扑，不能自动宣称存在近场真实零点。正式筛选不得擅自新增几何、条件数或距离门槛。

历史编码：As 对应 screw-in，Bs 对应 screw-out。宏观顺/逆时针、涡量符号、磁场方向和 As/Bs 应分别说明，不能互相替代。

用户要求的稳定度试验：
M_lambda = |lambda_r| / max_j|lambda_j| >=0.5，
lambda_r 为 As/Bs 唯一实特征值。该参数获得用户明确试验授权；其物理可靠性和全库最优性尚未验证。

持续 5 s 曾有两个不同版本：
1. 只要求 40% 误差连续合格，区段内类型可变化；
2. 加稳定度后，最后实现分别要求 As 或 Bs 连续合格，类型切换会中断区段。
最新对照展示版本 2 与完全逐点筛选。不能将两次统计差异全部归因于稳定度阈值。

### 图

最后 7-panel 顺序：
1. 原始 B；
2. 原始 V；
3. 5 s 平均 B；
4. 5 s 平均 V；
5. FOTE 误差和类型；
6. FOTE-V 误差和类型；
7. screw-in/out。

标题只保留 KH 事件号，主要输出 PDF，真实数据缺口要可见。PI 和距离图均已删除。检查 PDF 版面产生的预览页应放临时目录，避免进入 KH 正式交付目录。

## KH005 最后结果：已核对保存的 JSON 摘要

历史结果目录：
C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\codex\kh_events\KH005\eigenvalue_stability_0p5_compare

PDF：
- pdf\KH005_20151001180124_20151001180900_VI_sm5s_q40_eig0p5_pointwise.pdf
- pdf\KH005_20151001180124_20151001180900_VI_sm5s_q40_eig0p5_run5s.pdf

相应摘要位于 data\summaries，PNG 在 png，旧 MAT/CSV 在 data 内。2026-09-03 已直接读回两份 JSON 核对下列数据，本次没有重新运行 MATLAB：

| 参数/计数 | 逐点 | 同类型连续 >=5 s |
| --- | ---: | ---: |
| 完整公共时间轴行数 | 2772 | 2772 |
| 平滑物理时间窗 | 5 s | 5 s |
| 稳定度阈值 | 0.5 | 0.5 |
| 磁场 40% 合格点 | 2235 | 2235 |
| 流场 40% 合格点 | 1474 | 1474 |
| 最终 B As | 266 | 40 |
| 最终 B Bs | 230 | 35 |
| 最终 V As | 164 | 0 |
| 最终 V Bs | 254 | 0 |
| 同时同向 screw 样本点 | 67 | 0 |
| 同时反向 screw 样本点 | 42 | 0 |

2772 是完整时间轴行数，包含缺口位置的表示，不能称作 2772 个完整有效四星样本。
上述计数是相关的时间采样点，不能等同独立涡旋数。
稳定度单独对 B As/Bs 未额外剔除，对 V 从 As=172、Bs=256 降到 As=164、Bs=254，共减少 10 点。持续同类型 5 s 的要求造成了更明显的收缩。

此前仅 40% 误差连续 >=5 s 的版本：B 保留 1682 点、V 保留 221 点，同向 41、反向 4。这属于旧持续规则，不应与最新同类型连续规则混算。

## 单涡旋识别：认可的方向与未定内容

用户附图在误差 panel 中圈出红色 As 团块和蓝色 Bs 团块，希望程序识别整体特征，避免同一个短扰动被反复赋予不同最终类型。图中左右有多个团块，B 和 V 在右侧某些时段呈相应的红/蓝聚集。

用户明确认可：
- 用加权投票综合局部类型证据；
- 对一个单涡旋给一个最终标签。

仍未解决：
- 如何客观识别、框出单涡旋；
- 窗口怎样区分同一涡旋的内部结构与两个相邻涡旋；
- 选择什么权重、覆盖率、置信度和不确定类别规则；
- 同一涡旋 B 与 V 怎样配对以及怎样验证。

助手曾建议两组尚未定稿的方案：
1. 误差权重 [max(0,1-e/40)]^2，磁场 e=max(eta,xi)、流场 e=alpha；5 s 内加权投票，主导纯度 0.75、螺旋覆盖率 0.30，再用持续时间与事件整体标签。小空隙合并、短反类型小岛处理等也仅为建议。
2. 用 dB/dt、dVi/dt、涡量、swirling strength 的稳健归一化构造活动度；按峰和低谷寻找中心/窗口，周期约束候选长度；最终使用权重差的置信度 C，示例阈值 0.6。

两组方案的阈值和定义并不完全一致；用户没有逐项批准。活动度峰也尚未验证对应真实涡旋中心，不能把这个提议称为成熟识别方法。不要直接将这些算法植入正式 FOTE/FOTE-V 统计。后续先依据原程序与可审阅事件实例提出具体方案，再按用户“未审核方法不得使用”的要求处理。

旧助手用 KH005 色块给过纯度示例（B 约 18:05:37–18:05:53 偏 As，B 约 18:07:38–18:07:51 偏 Bs，V 有相应区段）；这些手选窗口和示例统计尚不足以证明自动涡旋识别已成立。

## 科学讨论和文献支线

以下为历史讨论线索，本次没有重新做文献查新，也没有认定此前助手的所有结论均正确。

- 均匀的四星公共背景相减会改变线性外推零点位置；固定时刻的局地梯度不变。四星瞬时均值相减会将中心场强制置零，不应作为提高真实零点可信度的方法。
- 流场参考系可改变停滞点位置，公共平移速度不改变速度梯度。旧回复曾笼统称类型随参考系变化；后续讨论已说明梯度类型在均匀平移下不变。密度不均匀时 nV 的 alpha 可能随参考速度改变，不能沿用“所有误差不变”的笼统说法。
- Poincaré 曾在 KH005 结果中全为零；后来用户要求删除。PI=0 的解释涉及拓扑度抵消，不能简单当作绝对无零点的证明。
- 初始研究建议涉及 B–Vi–Ve 拓扑对应、能量转换、MMS 四面体构型与尺度、非独立采样、事件/日期分层统计、避免用同一梯度同时构建解释量和能量指标。均未完成正式统计或构成新方法授权。
- 旧回复将 JGR 作为较完整统计工作的潜在目标，GRL 取决于是否产生集中稳健的机制结论。这是选题判断，没有发表保证。

代表性文献线索（来自原对话，正式引用前核对）：
- FOTE：Fu et al.，10.1002/2015JA021082；
- FOTE-V：10.3847/1538-4365/ab95a0；
- MMS KH 库：Rice et al. 2022，10.1029/2021JA029685；
- KH 重联统计：Wilder 2023，10.1029/2023JA031583；
- 磁鞘梯度拓扑：10.1029/2022JA031064；
- 磁尾拓扑/能量转换：10.3847/1538-4357/acf847；
- KH 能量转换：10.26464/epp2026042、10.1029/2025JA034952；
- 晨昏旋转讨论：Hwang et al. 2022，10.3389/fspas.2022.895514；原对话特别指出其摘要和正文的方向表述可能矛盾，需核对坐标和观察方向；
- 晨昏比较：Nishino et al.，10.1016/j.pss.2010.03.011；Grimmich et al.，10.1016/j.pss.2025.106182；
- 四星曲率/涡量：10.1029/2019JA026484。

Cluster 支线：
原对话整理出 14 个日期、15 个区间和 28 篇论文，含 9 个 A1、3 个 A2、3 个 B 级区间及额外候选记录。事件级证据和四星离子矩可用性需要分别核对，不能把卫星拥有四星磁场等同于具备标准四星 FOTE-Vi 条件。
已生成的工作簿仍位于：
C:\Users\Administrator\Documents\Codex\2026-08-03\c-users-administrator-documents-kh-figures\outputs\cluster_kh_catalog\Cluster_magnetopause_KH_published_event_catalog.xlsx
本次已确认文件存在；未扩大或重跑这项支线。

## 环境历史与后续承接

旧任务中 MATLAB R2026a 多次在 Home Session / MathWorks Service Host 附近崩溃，有些轮次只有静态检查、没有完整实跑。最后两套稳定度对照已经成功生成结果。不能把早期“环境阻断”当作当前仍存在的事实，也不能把某轮静态检查声称为该轮完整回归。
后续启动问题应先检查实际日志，不照搬旧环境临时方案，不因故障自行改用未经审核的科学算法或 Python 替代流程。

本次仅新增项目规则和本背景说明，未下载新 MMS 数据、未重跑事件、未迁移旧代码/图、未删除现有文件。未来执行实际任务时按用户最新目录规则落实输出，不再复制旧目录习惯。

AGENTS.md 的用途依照 [OpenAI 官方项目说明文档](https://learn.chatgpt.com/docs/agent-configuration/agents-md)：在项目任务开始时提供持久的项目指令。本项目规则放在 KH 根目录，详细历史放在本文件并由规则引用。


## 2026-09-04：KH事件、文献与四星数据汇总交付

本次最终工作簿：C:\Users\Administrator\Documents\KH\MMS_KH事件_文献与四星数据_20260904.xlsx。
共90个主表事件、30篇文献、115条事件与论文对应关系、1440条逐星产品核验记录。保留原KH001–KH089；原CSV仍为89条，未擅自改写。以后读取旧CSV时须注意KH090只在本次工作簿及90事件核查输入中。

- 新增确证事件KH090：2018-11-06 13:00–14:10 UTC，Petrinec et al. (2022)，10.3389/fspas.2022.827612。采用原文图4的70分钟展示窗口，不能等同于已逐波精确划定的KH边界。
- MATLAB直接读取原CDF，90个事件中80个具有四星B有效记录，74个具有四星B及Vi记录，52个具有四星Ve记录，80个具有四星E记录；49个的四星四类产品均读到有效记录。判据为事件内抽查记录的时间、有限值和填充值检查，未证明完整/共同连续覆盖、质量标记通过或FOTE适用。
- KH090的B、Vi、E均确认四星；Ve仅确认MMS1–3。新下载原CDF保存在Z盘标准MMS层级，未将科学数据转换为CSV。
- CAND01（2015-10-02 08:20–11:00 UTC，10.1029/2026JA035233）为DMC/KH调制候选，未加入90事件主表。
- Radhakrishnan 2024（10.1029/2024JA032869）和2025（10.1029/2025GL116901）的统计事件表/补充表尚未逐项取得；43事件不能直接当成新增，仍需去重核验。
- KH003、006、017、018、030保留原目录工程窗口标记；KH048有Rice与Wilder不同时间窗，工作簿已注明。
- 核验证据：Z:\SPART-WORK\Data\MMS\derived\KH\catalog_audit_20260903\cdf_presence.json 和 public_availability.json。后者为官方目录证据，burst只统计窗口内起始文件；survey/fast只到UTC日期级，不能据此断言事件内覆盖。目录返回截断标志时，文件数仅表示已返回部分。
- 本次仅整理目录、核查CDF和检索论文；未重画原事件图，未改变FOTE/FOTE-V方法。

