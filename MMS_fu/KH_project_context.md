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


## 2026-09-19：SMILE_MMS 官网图片下载（后台运行中）

用户要求下载MMS官网burst和quicklook页面中2026-07-20（含）之后所有卫星、所有时段的全部图片，并明确包含各仪器的单独图。
目标目录：C:\Users\Administrator\Documents\KH\SMILE_MMS。
用户随后要求按仪器分类，不按时间建立文件夹。最终结构：仪器（综合图、ASPOC、EIS、FEEPS、FIELDS、FPI、HPCA） / MMS1–MMS4或多星综合 / burst或quicklook / 官网原文件名。日期和时段仅保留在原文件名中。

2026-09-19读取官网完整索引：burst 2252张（563时段×4星，最新2026-08-15），quicklook 39235张（41种图类型，最新2026-09-16），总计41487张。全部记录起始日期>=2026-07-20；包含官网目录中的fpi_mms_summ多星图。
官网连续返回HTTP429，程序自动等待、降低请求频率并断点续传；完整任务可能跨夜或更久。不得把“已启动后台下载”说成已全部完成。

代码：
- C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\SMILE_MMS_download_plots_20260919.ps1
- C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\SMILE_MMS_download_controller_20260919.ps1

后台控制程序2026-09-19 15:18（北京时间）启动，首次PID 76460；后续不要只凭PID判断进程，应核查程序命令行或状态文件。控制程序连续下载、核验PNG首尾和完整清单，最多六轮补下缺失项，完成时写download_complete.json。
过程文件在 $env:TEMP\MMS_plots_20260919：controller_status.json（complete/phase）、all_progress.json（当轮进度）、all_download.jsonl（逐图记录/哈希）、expected_manifest.json、missing_images.json、controller_stdout.log、controller_stderr.log、controller_pid.txt。
目标文件夹只放图片。旧日期层级已清除，已有图片已移入仪器目录。
查询完成情况时必须读取这些状态文件并核对实际文件数量；截至本记录写入时任务仍在下载中。


## 2026-09-21：SMILE_MMS 下载加速

用户要求“能否想办法加速”。已检查实际吞吐量并修改同一下载程序，保留全部41487张清单及按仪器/卫星/类型分类。先前程序的请求间隔只会增加，最终停留在6000 ms；最近两小时实际平均约7.1秒/张。

已改为全局协调的自适应请求间隔：下限2000 ms，收到HTTP429/503时共同暂停并优先遵守Retry-After，间隔增加25%；连续成功40次且距上次调整至少2分钟后，间隔减少15%，最低恢复到2000 ms。排队线程在实际发出请求前重新检查全局暂停，避免在限流后仍按旧排期发请求。保留连接复用、已下载文件验证、原子写入和最终全清单核验。

20张短测均成功；随后连续下载67张，2026-09-21 08:36:29–08:40:10 UTC，约220秒，平均3.34秒/张，包含官网限流等待，未包含程序重启和既有文件核验耗时。这是短时实测，不能保证未来始终保持该速度或按固定时间完成。

最终程序北京时间16:40重新启动，PID 109648（查询时仍应核查进程/状态，不仅依赖PID）；已确认后台新增图片且错误日志为空。过程目录仍为 $env:TEMP\MMS_plots_20260919。最新输出日志为 controller_stdout_accelerated_v2.log / controller_stderr_accelerated_v2.log；controller_pid.txt、controller_status.json、all_progress.json、all_download.jsonl 仍为原有入口。加速记录另见 acceleration_restart.json、acceleration_measurement.json。throttle.json version=2保存当前间隔和暂停时刻，后续重启遵守未到期的暂停。

截至北京时间16:41，实际已下载约24202张，仍有约17285张；尚未完成。进度文件中的success是当前轮已处理数量，刚重启扫描时低于磁盘文件总数；controller_status.json在下载阶段的valid/missing只是该轮起始值，不能用于实时计数。统计当前完成量时应核对目标目录和完整清单。完成后仍由原控制程序生成 download_complete.json。


## 2026-09-21 21:08：优先下载120分钟Quicklook

用户明确要求“优先下载120min长度的quicklook”。已将工作程序排序的最高优先级改为kind=quicklook且原文件名以_0120.png结尾，覆盖全部仪器、全部卫星；随后继续其他时长和原有完整清单。目标目录分类保持仪器/卫星/模式，不建立日期目录。

重排时120分钟Quicklook共24780张，已有15086张，余9694张。原下载进程已替换为PID 86736，仍使用原控制程序、相同完整清单和限流状态；未重置尚未到期的官网等待。最新日志为TEMP\MMS_plots_20260919下的controller_stdout_priority120.log和controller_stderr_priority120.log；启动记录priority120_start.json，controller_pid.txt及其他通用状态文件继续更新。2026-09-21 13:09:30 UTC已确认新下载的文件是120分钟Quicklook，HTTP200；错误日志为空。尚未完成全部120分钟图片或全任务，后续查询仍须核对磁盘与完整清单。


## 2026-09-21：120分钟Quicklook已下载图片打包并邮件交付

用户要求将目前已下载的120分钟Quicklook图片打包并通过邮件发给本人。已按北京时间2026-09-21 21:12:41快照打包15115张，均来自原清单的quicklook且文件名以_0120.png结尾；保留仪器/卫星/quicklook目录和官网原文件名。逐张核对PNG首尾，并对压缩包解压流与源图进行SHA256比对，15115张全部通过。

本地完整ZIP：C:\Users\Administrator\Documents\KH\MMS_Quicklook_120min_20260921_211241.zip，969860346字节，SHA256 AF747E1D6768C7910DFF4CD294BEBE86A8A24E02632C5A8DAAA867830D6A65C9。此文件是本次用户要求的交付物，不能作为无用过程文件清理。

上传通道限制单文件512MiB，约500MB分包上传又遇到60秒超时，最终使用10个约100MB、可独立解压的ZIP交付；全部分包内容合起来恰为同一批15115张，无重复，逐张内容校验通过。云端读取核对10个文件的名称、大小和下载能力一致；云端文件夹保持本人私有。
下载文件夹：https://drive.google.com/drive/folders/1vo3qgHTm_18PD2qEMnmBur0yq3EdhU8_
收件人已从连接账号确认：fuwending429@gmail.com。已通过Gmail发送给本人，主题“MMS 120分钟 Quicklook 图片压缩包（15,115张，2026-09-21）”，消息/线程ID 1a0c42380b489389，返回标签含SENT和INBOX。邮件包含下载链接、全部分包解压到同一目录的说明及本地完整ZIP位置。

本次交付仅为已下载图片快照，不代表全部120分钟图已下载完成。原后台优先下载任务保持运行。过程快照、分包和回执仍位于TEMP\MMS_plots_20260919，回执package120_email_receipt.json。

## 2026-09-22：120分钟Quicklook新增图片第二次邮件交付

用户要求检查是否下载完并把剩下的发送过来，沿用上次120分钟Quicklook和本人Gmail范围。2026-09-22 23:30:55北京时间快照：120分钟图片已下载19170/24780张，尚缺5610张；其余时长Quicklook 14455张和Burst 2252张按原清单均已有文件。控制程序PID 86736仍正常运行，处于第2轮补下载；不得声称全部完成。

本次新增4055张，严格排除首批15115张，逐张校验PNG首尾及压缩包内解压后的SHA256，均通过。完整新增ZIP：
C:\Users\Administrator\Documents\KH\MMS_Quicklook_120min_additional_20260922_233055.zip
大小279387272字节，SHA256 346C3AB06135A333E012F0AE1165CBD0259016156C5CB3B091C48DEF56B6489D。
云端使用3个约100MB可独立解压ZIP，文件名、字节大小及下载能力均已读回核对。新增子文件夹位于上次交付文件夹下，仍为用户本人私有：
https://drive.google.com/drive/folders/1aTwgwvnT-ijFE-soe67Cqeqz6fV8_-bo

已发送至fuwending429@gmail.com，主题“MMS 120分钟Quicklook新增图片（4,055张，2026-09-22）”，消息/线程ID 1a0c9c19b43de4c4，回执标签含SENT和INBOX。邮件明确标注打包时间与尚未全部下载完成。

过程目录仍为TEMP\MMS_plots_20260919。本次快照、结果、回执分别为MMS_Quicklook_120min_additional_20260922_233055_snapshot.json、MMS_Quicklook_120min_additional_20260922_233055_result.json、MMS_Quicklook_120min_additional_20260922_233055_email_receipt.json。新增累计已发送清单package120_sent_paths.json，共19170条：后续增量发送必须先读取该累计清单和成功发送回执，排除已发送图片，不能仅排除最初package120_snapshot.json。打包期间和发送后新下载的图片不自动计入本批。

## 2026-09-23：原下载清单全部完成，120分钟图片全部邮件交付

用户再次询问是否下载完。下载控制程序于北京时间2026-09-23 02:39:28完成第2轮，download_complete.json和controller_status.json均确认complete。今日15:18再次遍历原清单逐文件检查PNG首8字节和IEND尾12字节：41487张全部有效存在，无缺失、无首尾损坏；后台进程已正常退出。

分类完成数：120分钟Quicklook 24780张；其他时长Quicklook 14455张；Burst 2252张。总计2689226533字节（约2.69 GB），图片仍在C:\Users\Administrator\Documents\KH\SMILE_MMS。完成范围是2026-07-20（含）之后、2026-09-19获取的官网清单，不代表已同步后续网站新增图片。

沿用用户此前“把剩下的给我发过来”的授权，已向本人fuwending429@gmail.com补发最后5610张120分钟Quicklook，排除累计已经发送的19170张。最终三批15115+4055+5610=24780张，覆盖完整120分钟清单，不重复。

本批本地完整ZIP：C:\Users\Administrator\Documents\KH\MMS_Quicklook_120min_final_20260923_151937.zip
大小323882184字节；SHA256 93E92F002ECCDA368BF144D97483739663B1246680C59616FFF90A54731CC248。逐张解压流SHA256与源图核对通过。云端使用4个可独立解压的ZIP，已读取核对名称、文件大小及下载能力，目录仍为本人私有。
本批下载文件夹：https://drive.google.com/drive/folders/1JbkBY369DFfFbG8OlEQ9i83PBjDuE22E
完整交付总目录：https://drive.google.com/drive/folders/1vo3qgHTm_18PD2qEMnmBur0yq3EdhU8_

邮件主题“MMS图片已全部下载完成：补发最后5,610张120分钟Quicklook”，消息/线程ID 1a0cd25bf1b1d174，发送回执标签包含SENT和INBOX。TEMP\MMS_plots_20260919中的MMS_Quicklook_120min_final_20260923_151937_email_receipt.json记录本批完整回执；package120_sent_paths.json现为24780条并逐项对应完整120分钟清单。此后不要重复发送这些图片；如要求下载新日期，先明确这是扩展原网站清单的后续任务。

## 2026-09-24：20 个事件 overview 改为 PPT，恢复原程序图形风格

用户要求把上一版107页PDF改为PPT，删除底部模式时间条，压紧panel，y轴标注使用原程序，并提供实际跑图程序。60张原始CDF科学图已由MATLAB重新绘制；107页PPT保持原事件/卫星顺序和44页官方参考图。9个物理panel依次为B、Vi、Ve、E、N、Ti、Te、Ei、Ee；原程序轴标注与温度三条曲线保留，panel间隙为0.002归一化图高。原有Overview_download.m / Overview_download_mms4.m未覆盖。

交付目录：C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924
PPT：MMS_20_events_MMS1-4_overviews_20260924.pptx（107页，32818779字节）
程序包：MMS_overview_MATLAB_programs_20260924.zip（40132字节，含实际.m、事件/原CDF清单和README）
本机代码位于MMS_fu；入口run_MMS_event_overviews_20260924，读取主函数MMS_event_overview_20260924，绘图函数Overview_download_events_20260924。默认15个有CDF事件、MMS1–4；完成记录存在会跳过，第三参数true可预览重画。原CDF沿用Z盘归档；本次未新下载。

数据来源、时间范围、坐标系、模式覆盖和样本数与上一版核对一致，60图验证无错误。PPT文件验证、107页PowerPoint实际渲染和逐页视觉检查均完成。科学图以高清图片嵌入，曲线需用MATLAB修改重画。官网参考图的原始坐标系、标签与顶部状态条保留。8月11日仍仅取得所需L2磁场；8月28–29日及9月16日所需CDF在上一版查询中未取得，继续保留官网参考页，未宣称数据补齐。

PPT Drive ID：1ZV6UAytT6BVPwfbihTkqL-ziSoAPsl2D
程序ZIP Drive ID：1u7yXt9-ElDCs5Cafwh6oR8zNwsayBzXt
两者在原20事件交付文件夹1tOsVUrPLXCjnQhRKCGtUtqpkI3MvCLtP内，权限仅本人。
已发往fuwending429@gmail.com，主题“MMS overview 修订版：107页PPT与MATLAB跑图程序”，PPT用下载链接，MATLAB ZIP作为附件。邮件ID/线程ID：1a0d2916dc69a791，读回确认SENT和INBOX。
验证与完整交付回执：Z:\SPART-WORK\Data\MMS\derived\event_overview_20260924。



## 2026-09-24：PPT加入20个事件的四星位置图

用户要求使用原MMS_orbit程序把各事件卫星位置放入PPT。已完成127页新版：每个事件overview前新增1页MMS1–4地心位置/相对构型图，其余107页在PowerPoint渲染后逐像素与上一版一致。20页新增页均已单独检查。使用GSM，单一时刻直接使用，区间取中点，UTC；RE=6372 km沿用原函数。

原稿参考：FWD_matlab\新建文件夹\MMS_orbit.m。调用的mms.mms4_pl_conf因本机毫秒时间格式错误，使用兼容副本MMS_orbit_pl_conf_20260924.m：仅改名、自调用/回调名以及UTC显示格式，原文件未覆盖。布局函数将原7个物理panel及图例排为两行，扩展X范围，保留原坐标轴物理量和符号。入口MMS_orbit_events_20260924(eventNumbers,force)，默认1:20、force=false，已完成项跳过；true强制重画。实际代码位于MMS_fu。

位置数据全部归档Z:\SPART-WORK\Data\MMS。EV01–EV15新增24个MEC CDF（66216155字节），在原mms1–4/mec/srvy/l2/epht89d/2026/07或08目录。EV16–EV20暂未取得MEC，采用NASA SSCWeb官方GSE星历，存于ancillary\sscweb\2026，原生60秒，样本包围范围内调用irf_resamp对齐、原IRFU函数转换到GSM。EV15额外SSC数据用于交叉核对，绝对位置差最大27.505km，星间相对矢量差最大0.165km，此核对仅适用于这个交集样例。下载的20个DEFEPH原始文件保留在mms1–4/ancillary/defeph/2026，但未用于最终图形，避免旧IRFU转换未完整处理J2000的问题。

边界线为原程序模型；EV06、07、08、09、11、12、19、20采用默认P=2nPa/Bz=0nT，均在图上标注assumed model。其余事件用原函数取得的OMNI参数。不能把模型线当成实际边界观测；卫星坐标全部来自星历。

交付目录：C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924
PPT：MMS_20_events_MMS1-4_overviews_with_orbits_20260924.pptx，127页，35814481字节，
SHA256 8502016850D847F1E98FF55FEA3FE0C1AA1ECC3B74C9AD4047CC22F637A9AC33。
程序：MMS_orbit_MATLAB_programs_20260924.zip，15561字节，含实际3个MATLAB函数、事件清单、原稿参考及说明。
orbits子目录含20个PNG和20个MATLAB FIG。
审计、来源、每页坐标/时刻、PPT页索引及交付回执在Z:\SPART-WORK\Data\MMS\derived\event_orbits_20260924。

Drive PPT ID：1HiMteiLX4vfi7I5HKuhLZELoV-dHmelU
Drive ZIP ID：19i18ODYRZh6XcaoT7rDpYvnyuXpVclRC
两者仍在本人私有文件夹1tOsVUrPLXCjnQhRKCGtUtqpkI3MvCLtP，未改变分享权限。
已发送fuwending429@gmail.com，主题“MMS 20事件PPT已加入四星位置（127页）及轨道绘图程序”。
邮件ID/线程ID1a0d2bc940b15039，读回确认SENT和INBOX。PPT以链接、轨道程序ZIP以附件交付。

## 2026-09-24：遵照原 overview 直接脚本格式重写并重新跑完

用户再次明确要求尽可能直接调用 IRFU/用户已有程序，避免新写 function，按原 overview 的 %% 分区和变量、标注方式写代码。此要求已补充到 KH/AGENTS.md。

两份实际新脚本（均无 function 定义，替代本任务上一版自定义函数入口）：
- MMS_fu/Overview_events_original_style_20260924.m：EventList=1:20; Spacecraft=1:4; 直接读取 Z 盘原CDF、调用IRFU读取和绘图。Burst优先，其余survey/fast；不平滑、不插值科学数据，缺口保留。E按600秒块读取，调用库中reduce_to_width做绘图min/max降采样，显示bins不作为原生时间戳。B/Vi/Ve使用irf_plot的reduce选项。9 panel原标签和颜色，间距0.002，无底部模式条。
- MMS_fu/MMS_orbit_original_style_20260924.m：直接按原MMS_orbit.m调用mms.mms4_pl_conf；EV01–15用MEC，EV16–20直接读取已归档SSCWeb JSON，把原生60秒时间轴/GSE位置交给现成IRFU函数。无另行30秒预重采样。

原IRFU mms.mms4_pl_conf.m仅两处日期显示格式由HH:MM:SS.mmm改为HH:MM:SS，修复本机GenericTimeArray格式兼容问题；完整原文件备份在Z盘derived/events_original_style_20260924/mms4_pl_conf_before_date_format_fix.m。科学算法未改。用户原Overview_download.m、Overview_download_mms4.m、MMS_orbit.m保持原样。

全部80个事件/卫星任务完成：60张原始数据overview（EV01–15），20个无L2记录（EV16–20），20张四星位置图。EV13–15仍只有磁场L2，MMS4电子缺测保持空白；本地归档无FPI Vi/Ve burst，Vi/Ve使用可用fast。本次未补下载科学CDF，沿用前次数据可用性查询日期，未宣称新的在线查询。
B/Vi/Ve/E共480项原生采样点数/模式比较无差异。谱图使用mms.variable2ts居中时间及原能量矩阵，裁到请求区间。

交付目录：C:/Users/Administrator/Documents/KH/MMS_events_original_style_20260924
PPT：MMS_20_events_original_style_20260924.pptx，127页，30930172字节，SHA256 ab4f736a641a8cb93fe99bfd13e2b725bddadbe91a7330c8a67e9ce77b564813。
MATLAB包：MMS_original_style_MATLAB_20260924.zip，36618字节，SHA256 44505f44f7d0ccbb6e5210a90212fc6d9f6b7c20269b010b75ef2723433af22b。
包内含实际两份脚本、20事件JSON、README、3份用户原稿参考，已逐字节核对与本地实际运行代码一致。
PPT构建代码：MMS_fu/MMS_original_style_20260924_ppt.mjs。127页均通过结构/布局检查及PowerPoint实际渲染；81个修改页逐页目视检查，46个原官方参考/索引页与原PPT渲染逐像素一致；80张科学/轨道嵌图逐字节与新输出一致。

Drive PPT ID 1QK5MU9tkPOKJH92Ufej-m0GzTvZOgpXL；
Drive ZIP ID 1YsWdcq7kjZJD_9INI55-Ir7Yo1Efxbjc。
仍在本人私有文件夹1tOsVUrPLXCjnQhRKCGtUtqpkI3MvCLtP，权限核对仅fuwending429@gmail.com本人。
已发往fuwending429@gmail.com，主题“MMS 20个事件：按原overview格式重写并重跑的PPT及MATLAB程序（2026-09-24）”，PPT使用下载链接，MATLAB ZIP附在邮件中。
邮件ID/线程ID 1a0d2f7c73be4b81，回读确认SENT、收件人及主题。
完整运行状态、每页来源/坐标、包哈希、最终验证与邮件交付回执：Z:/SPART-WORK/Data/MMS/derived/events_original_style_20260924。

## 2026-09-24：v2 完成坐标/布局修订，并加入 FEEPS、HPCA

按用户确认，在同一张 overview 下方增加 FEEPS 电子/离子全向能谱、HPCA H+/He+/He++/O+ 密度及四种离子能谱（7个新panel）。60张原始CDF overview全部重跑为16个panel，20张轨道图重画，127页PPT检查完成。44页官方quicklook按用户要求保持原样；与原版PowerPoint渲染逐像素一致。两页索引也未改。

重画图移除事件竖线、模式时间条，时间标签水平，panel紧排，原物理量标签保留。B、Vi、Ve、E均为GSM；E只读取GSE并由irf_gse2gsm转换，只有DSL则留空；标题不标坐标系。轨道XY/XZ相邻，地球纵轴居中，卫星位置留足边缘空间。科学数据未平滑/插值，burst有效段优先、其余survey/fast，真实缺口保留。

实际科学代码：
- MMS_fu/Overview_events_original_style_20260924_v2.m
- MMS_fu/MMS_orbit_original_style_20260924_v2.m
- MMS_fu/MMS_IRFU_particle_compatibility_20260924_v2.m

以上为直接脚本，无新增 MATLAB function 定义，沿用用户原%%分区、IRFU读取绘图和mms.mms4_pl_conf。原用户程序未覆盖。IRFU get_data.m做两项已核实的兼容修正：HPCA He++ survey入口白名单拼写；FEEPS未启用探头精确占位值-2147483648（CDF声明FILLVAL不同）先置NaN再由原有探头均值处理。正常burst样例不变，修正前文件与逐行diff均保存在本版Z盘derived目录，程序包附说明及diff。不加入新探头质量门槛/物理过滤。

本轮2094个FEEPS/HPCA原始CDF、23612476209字节全部保存在Z:/SPART-WORK/Data/MMS官方产品分类中，保留文件名/版本；无全量CSV/MAT转换。下载校验CDF签名和长度，最终归档校验通过。查询清单/候选覆盖/下载记录在derived/events_original_style_20260924_v2/particles。8月11日原有仪器仍只有B，新增可用FEEPS/HPCA；EV16–20无所需L2、仅保留官方参考页与位置图。FPI Vi/Ve仍使用现有fast，未因新增FEEPS/HPCA burst而宣称FPI moments已有burst。

所有80个事件/卫星记录完成（60张图+20个无L2记录），20个轨道完成，无错误、无过期输出。对照v1核对原生点数、时间区间、GSM、缺测与位置；除授权移除的MMS4 DSL电场外，原有原生点数一致。127页PowerPoint实渲染；60张overview、20张轨道、封面目视检查，其余46页与原版逐像素一致，80张嵌图逐字节与输出PNG一致。源代码ZIP逐字节核对通过。

交付目录：C:/Users/Administrator/Documents/KH/MMS_events_original_style_20260924_v2
PPT：MMS_20_events_original_style_20260924_v2.pptx，127页，65326695字节。
SHA256 f2fe892d72f4d9694077e015beacb3ca162327ae7368e802343905bcb4bef6ee
程序ZIP：MMS_original_style_MATLAB_20260924_v2.zip，54221字节。
SHA256 e8e329403f59a2750eae89c1b5100d8f0a99e3069892e66cbea894e5a63b8f82
包含实际脚本、事件JSON、运行说明、原程序参考、现成下载器及兼容差异文件。

Drive PPT ID：1hPi-q6jg0N-2Hr6yM-OJysYK37PiadsN
Drive ZIP ID：1wtFqblFhzLL51eb95cu_FqpG8UhZFjv2
仍在本人私有文件夹1tOsVUrPLXCjnQhRKCGtUtqpkI3MvCLtP，读回大小一致，权限仅fuwending429@gmail.com本人。
邮件主题“【更新 v2】20 个事件 MMS1–4 overview：加入 FEEPS/HPCA，PPT 与实际跑图程序”。
已于北京时间2026-09-24 23:08:41发往fuwending429@gmail.com；邮件/线程ID 1a0d3f6315b28119，回读确认SENT、收件人与54221字节ZIP附件。PPT通过私有下载链接交付。
完整科学验证、最终逐页核对、源文件包哈希和邮件/上传回执在Z:/SPART-WORK/Data/MMS/derived/events_original_style_20260924_v2。

## 2026-09-27：能谱下限按用户要求调整为20 eV，已交付

用户要求：所有原下限低于10 eV的能谱panel下限改为20 eV。仅调整自行绘制的图，保留各能谱原上限；FEEPS原下限不低于10 eV的能谱保持原范围。官方quicklook沿用用户此前要求保持原样。

实际主程序：MMS_fu/Overview_events_original_style_20260927.m，从20260924_v2直接脚本复制，仅新增能谱YLim条件设置、轴范围记录及新输出路径，无新MATLAB function。全部80个事件/卫星任务完成，60张overview重画、20个无L2记录；本次复用Z盘CDF，未新下载或查询可用性。数据点数、模式、时间、缺测、坐标、上限与v2一致。修改diff附程序包。

交付目录：C:/Users/Administrator/Documents/KH/MMS_events_original_style_20260927
PPT：MMS_20_events_original_style_20260927.pptx，127页，64718876字节，SHA256 74dd57b3a6f97b2cd20cfa4be0f58462c7e1afd1cf777ff3189f1690a25f1d0e。
程序ZIP：MMS_original_style_MATLAB_20260927.zip，55863字节，SHA256 3aafad62730db4c83bf3395766738a1f41176a5057e00eb2e0aaa36e8b1fefaa。
127页PowerPoint实际渲染；60张重画图逐页目视核对，60张嵌图逐字节与新PNG一致；其余67页与v2逐像素一致（含20轨道和44官方quicklook）。结构、布局、导入检查均通过。

私有Drive PPT ID：1Gw4RyNagcaFG71Et4OSV0SM0lPG26eEt；ZIP ID：1s67NmQEtZC9YPy6nvhEVrqSmrKn1M_3h。读回大小一致，权限仅fuwending429@gmail.com本人。
已发送至fuwending429@gmail.com，主题“MMS overview 更新：能谱下限调整为 20 eV（20260927）”，PPT私有链接、ZIP附件。邮件/线程ID 1a0e38c47a66bc20，回读确认SENT、收件人和55863字节附件。
科学/逐页验证、包校验和交付回执保存在Z:/SPART-WORK/Data/MMS/derived/events_original_style_20260927。
## 2026-09-29：EV01 MMS1去掉FEEPS/HPCA，已完成本地重画

用户提供EV01 MMS2截图，并明确要求同一事件重画MMS1、去掉下方全部FEEPS/HPCA。已从Overview_events_original_style_20260927.m复制为Overview_EV01_MMS1_9panels_20260929.m，移除这些仪器的读取/绘图及无关兼容段，默认EventList=1、Spacecraft=1。保留前9个panel，原load data及B至FPI能谱绘图段逐字相同，无新MATLAB function。

时间2026-07-20 23:50至2026-07-21 04:10 UTC；矢量GSM；两张FPI能谱范围20–40000 eV；标题MMS1、时间标签水平、无事件竖线。MATLAB实跑完成，九个panel均有数据；时间、原生点数、模式、电场分段与20260927版完全一致。PNG目视检查通过。

输出：C:/Users/Administrator/Documents/KH/MMS_EV01_MMS1_20260929/EV01_MMS1_overview.png（579676字节，SHA256 7f53ea4b3572b67ba4232558a2947a207d3bbc380bd631f549b38cf7f5b31f8c）；同目录MMS_EV01_MMS1_MATLAB_20260929.zip（11811字节，含实际脚本、事件JSON、说明及修改diff，逐字节核对通过）。科学记录/校验在Z:/SPART-WORK/Data/MMS/derived/EV01_MMS1_9panels_20260929。原整批PPT保持原样。

尝试按先前授权发送图片与程序至fuwending429@gmail.com，但Gmail工具返回user rejected MCP tool call。此次邮件未发送，没有重试；成果本地交付。
## 2026-09-29：Case_for_SMILE.pptx中5个时间段的AE图已完成

用户要求按C:/Users/Administrator/Documents/Recovery-Work_SMILE-MMS/Case_for_SMILE.pptx中给出的事件时间分别画AE。已提取并目视核对全部8页，以第2、3、4、5、8页“事件时间”文字为准：Case1 2026-07-21 00:00–04:00；Case2 2026-07-24 08:00至07-25 05:00；Case3 07-28 05:10–06:00；Case4 08-04 01:50–03:00；Case5 08-11 06:50–08:30，全部UTC。

数据为NASA CDAWeb OMNI_HRO_1MIN的AE_INDEX（nT，1分钟，WDC Kyoto quicklook）。原CDF omni_hro_1min_20260701_v01.cdf、omni_hro_1min_20260801_v01.cdf已保存到Z:/SPART-WORK/Data/MMS/ancillary/omni/hro_1min/2026，每个8772893字节。来源说明 https://omniweb.gsfc.nasa.gov/html/omni_min_data.html 。

新直接MATLAB脚本 MMS_fu/AE_Case_for_SMILE_20260929.m 使用dataobj/get_ts/irf.ts2mat/irf_tlim/irf_plot/irf_zoom/irf_subplot，无新function。AE为整数CDF，按原始FILLVAL=99999掩码在double绘图数组中保留NaN并逐点核对原变量；不平滑、不插值、不重新平均。沿用irf_tlim的[start,end)选择，坐标轴与PPT范围一致。5个事件分别240/1260/50/70/100个有效分钟样本，没有缺测。跨日Case2标出两个日期。

交付目录C:/Users/Administrator/Documents/KH/Case_SMILE_AE_20260929，含Case1–5独立PNG及Case_for_SMILE_AE.pdf（5页矢量图，每个事件一页）。PDF5页目视检查通过，最终Case2日期标注检查通过，其余4页与已检查渲染逐像素一致。来源、哈希、运行与验证记录在Z盘derived/Case_SMILE_AE_20260929。代码说明README_AE_Case_for_SMILE_20260929.md在MMS_fu。原PPT未修改。

用户已明确不需要邮件发送，当前及本次后续交付直接保存本地。
## 2026-09-29：MMS1 2026-07-25 01:00–04:00 UTC九panel overview

按用户给出的MMS1九panel参考图及精确时间范围，直接复制Overview_EV01_MMS1_9panels_20260929.m为Overview_MMS1_20260725_0100_0400.m。仅修改事件/时间（本次不另加10分钟）、标题和输出路径；原数据读取、GSM转换、9个panel及20 eV能谱下限不变，无新MATLAB function。

复用已有MMS_event_overview_20260923_download.py的查询/下载函数，通过MMS1_20260725_0100_0400_download.py补齐6个CDF、164503698字节，全部按官方产品层级存入Z:/SPART-WORK/Data/MMS。当天FGM/EDP burst唯一候选从04:46:33开始，在所选区间外；FPI burst moments未查询到。本图使用FGM survey及FPI/EDP fast。

MATLAB已实跑完成，9个panel均有数据，B 172798点、FPI各2400点、E 345596点；矢量GSM、无事件竖线、时间标签水平、能谱20–40000 eV、无FEEPS/HPCA。原始缺口处理沿用前脚本。PNG目视检查及记录/原CDF长度签名检查通过。

输出 C:/Users/Administrator/Documents/KH/MMS1_20260725_0100_0400/20260725_0100_0400_MMS1_overview.png。查询、下载、运行、代码diff和校验记录在Z盘derived/MMS1_20260725_0100_0400。仅本地保存，无邮件。
## 2026-09-29：MMS1磁尾2026-07-20—08-15四小时分段图已完成

用户要求X<0时每4h两张图，分别为B/Vi/AE三panel，以及B/Vi/AE/Ni/离子能谱五panel；仅MMS1。用户确认从每次进入X<0开始，每4h分段，末尾不足4h保留。日期包括8月15日全天，UTC/GSM。官方MEC epht89d原生30s轨道通过相邻点interp1估计X=0穿越；四小时网格以实际进入时刻为起点，再与日期范围取交集。首个进入发生在7月19日，因此首段截取为7月20日00:00—01:44:41.649，保留原进入时刻的网格。

得到8次磁尾经过、162段，磁尾累计620.716小时，实际生成324张PNG，全部完成无绘图错误。输出C:/Users/Administrator/Documents/KH/MMS1_tail_20260720_0815，子目录B_Vi_AE及B_Vi_AE_Ni_Ei；index.html可逐段浏览两张图，MMS1_tail_windows.csv为时间和中点位置清单，MMS1_tail_MATLAB_programs.zip为实际程序包。

MATLAB直接脚本MMS_fu/MMS1_tail_20260720_0815_intervals.m及Overview_MMS1_tail_20260720_0815.m，沿用原overview的%%结构及B/Vi/Ni/FPI能谱绘图段，调用现成IRFU函数，无新增MATLAB function。标题为时段中点GSM位置（RE=6372km）；紧凑panel、水平UTC时间标签、无事件标记竖线；Ni按本次截图使用红色，离子能谱20—40000eV。未平滑/填补观测缺口；FGM burst优先，survey补缺。

原始MMS文件800个（MEC 30个+科学产品770个）均存Z:/SPART-WORK/Data/MMS官方层级；复用185个、补下载615个。科学产品共3348316093字节，MEC共80984666字节。47段有FGM burst。SDC本批FPI dis-moms burst无文件，fast到8月8日；单独复查8月9/11/15日仍空。CDAWeb orig_data接口查询8月9—15日fast及全时段burst也返回空清单，7月25日fast返回7文件作为接口有效对照。本批75段完全没有离子L2数据，3个离子panel明确标注No available L2 data；其他段中的局部缺口保留。

AE复用Z:/SPART-WORK/Data/MMS/ancillary/omni/hro_1min/2026下7/8月OMNI原始CDF。AE_INDEX为1分钟Kyoto quicklook，按整数CDF的FILLVAL掩码在double中保留NaN并逐点核对原变量；本批所有AE分钟记录有效，没有缺测。各段时长、AE计数、图片解码/哈希、CDF长度/文件头、程序包逐字节验证通过；162张五panel图已通过14张缩略检查页审阅，并检查了完整尺寸的3/5panel、缺测和短尾段样例。

查询、分段、来源及完整验证/交付记录位于Z:/SPART-WORK/Data/MMS/derived/MMS1_tail_20260720_0815。过程日志/检查页保存在TEMP/MMS1_tail_20260720_0815，未混入KH交付目录。本次仅保存本地，无邮件。
