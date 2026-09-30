# Cases for SMILE：聊天承接与当前工作规则

整理日期：2026-09-30，北京时间。当前结果根目录：C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS。

## 来源与完整阅读范围

- 用户指定聊天：Cases for SMILE，任务ID为01a0c9c5-4b4c-72b1-987c-d0e41c2f40e6，hostId为local。
- 已调用read_thread，逐页读取到hasMore=false。该接口对七个较早回合返回空正文，因此同时读取指定聊天的本地原始记录补全。
- 原始记录：C:\Users\Administrator\.codex\sessions\2026\09\22\rollout-2026-09-22T23-39-05-01a0c9c5-4b4c-72b1-987c-d0e41c2f40e6.jsonl。
- 按原始记录核对全部23条用户可见消息、15条助手最终回复、168条过程说明，并查看三张用户参考图。用户消息计数排除自动注入的环境、插件清单及AGENTS文本，包含四次用户选择回复。
- 阅读时间范围为2026-09-22至2026-09-29。下文聊天日期使用北京时间；事件观测时间使用UTC。
- 原聊天的工具日志、已有脚本和KH_project_context.md用于核对实现、文件路径和历史结果。助手推断、经验阈值、方法建议及当时的数据可用性均按各自来源理解。

## 2026-09-30用户明确变更

其余要求保持不变，今后本任务的图件、PDF、PPT、浏览索引、时间清单及用户要求的交付包统一保存到C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS，可按任务建立子目录。

代码仍保存到C:\Users\Administrator\Documents\FWD_matlab\MMS_fu；原始MMS数据仍按原产品层级保存到Z:\SPART-WORK\Data\MMS。OMNI、SSCWeb等辅助数据及必要衍生记录沿用该Z盘根目录的ancillary、derived分类。运行日志、临时检查页和过程文件放临时目录。

当前项目规则保存在C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\AGENTS.md。原KH/FOTE研究背景仍参考MMS_fu\KH_project_context.md，其历史KH路径用于追溯。本次未移动、删除、重画旧结果或补下载科学数据。

## 全部23条用户消息与承接状态

| 序号 | 日期 | 用户内容 | 承接要求或状态 |
| --- | --- | --- | --- |
| 1 | 09-22 | 猜测用户身份 | 助手猜测未经用户确认，不作为身份事实 |
| 2 | 09-23 | 列出20个事件，逐事件画overview，前后各至少10min，有burst优先，无则survey/fast | 保留事件顺序、时间范围与逐仪器实际覆盖 |
| 3 | 09-23 | 确认每个事件画MMS1–4 | 原20事件批次按四星绘制 |
| 4 | 09-23 | 尽量复制并修改用户已有程序，少自行写代码 | 直接复用原overview与IRFU接口 |
| 5 | 09-23 | 询问下载数据是否已存数据库 | 当时1117个CDF已归档Z盘，约21.2GB |
| 6 | 09-23 | 通过邮件交付结果 | 历史曾发给本人；09-29已改为直接本地保存 |
| 7 | 09-23 | 询问Vi、Ve是否仍无burst | 当时重新查询SDC，未取得对应FPI burst L2 |
| 8 | 09-23 | 要求同时核对CDAWeb | 当时8个DIS/DES burst矩数据集均未覆盖所选日期 |
| 9 | 09-24 | PDF改PPT；去模式时间条；压紧panel；y轴沿用原程序；给实际代码 | 107页PPT完成；后续版本继续沿用这些格式要求 |
| 10 | 09-24 | 每个事件加卫星位置，用MMS_orbit | 每个事件overview前新增一页四星位置，共127页 |
| 11 | 09-24 | 必须模仿原overview格式，尽量现成函数，重写重跑 | 新代码采用直接脚本与%%分区，减少自写处理及封装 |
| 12 | 09-24 | 询问一段处理代码在做什么 | 助手解释了burst重叠屏蔽、NaN断线与自定3*dt阈值；这些实现不能据此成为已审核科学标准 |
| 13 | 09-24 | 删竖线；DSL-only电场不画；水平时间标签；矢量GSM；位置投影相邻且地球居中 | 使用set(gca,"XTickLabelRotation",0)，标题省略坐标系；X–Y、X–Z相邻并留边距 |
| 14 | 09-24 | 确认只改自行画图，官方quicklook保持原样 | 原始quicklook不转换、重标或重绘 |
| 15 | 09-24 | 有FEEPS、HPCA数据也画 | 使用现成读取/绘图接口，按有效数据增加仪器panel |
| 16 | 09-24 | 确认新增仪器加到同一overview下方 | 原批次最多16个panel |
| 17 | 09-27 | 所有下限低于10eV的能谱改为20eV | 保留原上限；下限不低于10eV的保持原范围 |
| 18 | 09-29 | EV01改画MMS1，去掉FEEPS、HPCA | 只保留九panel，时间仍为07-20 23:50至07-21 04:10 |
| 19 | 09-29 | 不需要邮件，直接保存 | 当前默认本地交付，不自动沿用旧邮件发送流程 |
| 20 | 09-29 | 按Case_for_SMILE.pptx事件时间分别画AE | 五个精确时窗，独立PNG与每事件一页PDF已完成 |
| 21 | 09-29 | MMS1在07-25 01:00–04:00画同类图 | 精确时窗，九panel，不再另加10min |
| 22 | 09-29 | 07-20至08-15，X<0时每4h画两张MMS1图，AE用CDAWeb OMNI 1min | 三panel为B/Vi/AE；五panel增加Ni及离子能谱 |
| 23 | 09-29 | 确认从每次进入X<0起分4h，保留不足4h尾段 | 保留进入时刻的分段网格，再与请求日期范围取交集 |

## 代码、数据与方法要求

1. 新代码尽量使用MATLAB。优先检索irfu-matlab-master、MMS_fu和用户指定原程序，阅读接口后直接调用。代码采用原overview的%% load data、%% B plot、%% Vi plot等顺序分区、熟悉变量及中文备注。
2. 原程序参考包括Overview_download.m、Overview_download_mms4.m、Overview_srvydownload.m及FWD_matlab\新建文件夹\MMS_orbit.m。用户允许复制后修改，保留原文件。
3. 没有新增function并不足以证明满足要求；自写分段、过滤、重叠屏蔽等处理同样需说明依据。历史3*dt缺口阈值来自助手经验，用户未把它确立为正式科学判据。本次目录变更未改变这些历史处理段。
4. 原始CDF直接由MATLAB/IRFU读取，保留文件名与版本。优先复用本地Z盘CDF，缺失时可从官方网站补下载。不同事件共享归档文件，避免默认全量CSV/MAT转换。
5. 模式优先级按仪器分别确认。FGM、EDP或FEEPS的burst存在，以及官网burst预览存在，均不能用于证明FPI burst速度L2已发布。网站缺文件只说明当时公开可取得状态。
6. 缺测保持空白/NaN，不补造数据。矢量转换使用现成IRFU接口，当前自行画图默认GSM；只有DSL电场时留空。
7. FEEPS能量/通量以原数据单位为准，HPCA同样保留原单位。历史IRFU兼容修改包括HPCA He++ survey入口拼写、FEEPS精确占位值-2147483648及位置函数日期显示格式；有备份、diff和记录，不引入新质量门槛。
8. KH/FOTE背景约束继续沿用AGENTS.md与KH_project_context.md。KH005的GSE、5s平均及FOTE筛选属于该科学分析版本；不能自动套入本任务原始overview。
9. 图与PPT必须核对实际时间范围、卫星、单位、坐标、数据模式及缺口；交付程序对应实际源代码。过程检查页不混入正式结果目录。

## 原20事件清单

全部为UTC。以下物理描述沿用用户的初步判断，问号和待核查状态保留。

| 事件 | 日期 | 原指定时间 | 用户描述 |
| --- | --- | --- | --- |
| EV01 | 2026-07-21 | 00:00–04:00 | 磁层顶KH |
| EV02 | 2026-07-24 | 07:00 | 重联或接近等离子体片中心？ |
| EV03 | 2026-07-24 | 09:20 | 偶极化伴流涡旋 |
| EV04 | 2026-07-28 | 05:20 | 偶极化伴流涡旋，用户标有burst |
| EV05 | 2026-07-28 | 06:15 | 连续偶极化及连续流爆发？ |
| EV06 | 2026-07-28 | 21:30 | 偶极化伴流涡旋 |
| EV07 | 2026-07-28 | 23:10 | 重联 |
| EV08 | 2026-07-31 | 11:10 | 重联锋面伴电子加速 |
| EV09 | 2026-07-31 | 12:00–13:00 | 湍流重联 |
| EV10 | 2026-07-31 | 18:40 | 重联伴锋面 |
| EV11 | 2026-08-04 | 01:10–01:30 | 连续偶极化及连续流爆发？用户标有burst |
| EV12 | 2026-08-04 | 07:00 | 重联？ |
| EV13 | 2026-08-11 | 08:00–09:00 | 连续偶极化及连续流爆发？用户标有burst |
| EV14 | 2026-08-11 | 15:00±5min | 偶极化伴流涡旋，用户标有burst |
| EV15 | 2026-08-11 | 17:00–17:10 | 偶极化伴流涡旋 |
| EV16 | 2026-08-28 | 23:05 | 可能为连续偶极化锋面，待核查 |
| EV17 | 2026-08-29 | 04:00–05:00 | 可能为连续偶极化锋面，待核查 |
| EV18 | 2026-08-29 | 20:40–21:00 | 可能为PSBL的KH涡旋？ |
| EV19 | 2026-09-16 | 00:00–00:10 | 偶极化伴流涡旋 |
| EV20 | 2026-09-16 | 01:40–01:50 | 连续偶极化及连续流爆发？ |

实际事件JSON为MMS_fu\MMS_event_overview_20260924_events.json。EV14原区间为14:55–15:05，overview再延长为14:45–15:15；跨日窗口正常保留。

## 最后绘图格式与历史交付状态

- 九个基本panel为B、Vi、Ve、E、N、Ti、Te、Ei、Ee；保留原y轴名称、单位、温度标量/平行/垂直曲线与蓝绿红分量配色、黑色|B|和原能谱配色。panel间距接近零；没有底部数据模式时间条及事件标记竖线，底部UTC标签水平。
- 原20事件整批在2026-09-27最终版得到60张自行绘制overview、20张位置图、44页官方quicklook；PPT共127页。FEEPS/HPCA在同图下方，原批次为16panel。
- 当时EV01–EV15可取得部分L2；08-11原有仪器主要只有B，新增可用FEEPS/HPCA。EV16–EV20当时未取得所需科学L2，保留官方参考图。MMS4电子缺测及DSL-only电场留空。这些是历史查询结果，后续新任务应重新核对发布状态。
- 位置图调用mms.mms4_pl_conf；优先官方MEC，较新事件历史用SSCWeb原生GSE星历交给现成IRFU转换，存档在Z盘ancillary。原DEFEPH文件保留，但历史最终图未采用其J2000转换结果。
- 单时刻位置直接使用指定时刻，区间取中点；GSM，RE=6372km沿用原程序。X–Y和X–Z相邻，纵轴对称使地球居中，卫星留边距。边界线为模型，默认P=2nPa、Bz=0nT的样例已标assumed model，不能当作实际边界观测。
- 用户最后确认修改范围只限自行绘图，原44页quicklook保持其原坐标、标签和状态条。
- 09-29单独EV01 MMS1与07-25图均采用九panel、去掉FEEPS/HPCA。EV01范围23:50–04:10，07-25范围01:00–04:00；后者当时未取得覆盖时窗的burst，采用survey/fast。

## Case_for_SMILE.pptx的五个AE事件

源PPT在当前结果根目录，第2、3、4、5、8页的“事件时间”为依据；全部UTC，不另扩展。

| Case | 开始 | 结束 | 历史有效分钟数 |
| --- | --- | --- | ---: |
| 1 | 2026-07-21 00:00 | 2026-07-21 04:00 | 240 |
| 2 | 2026-07-24 08:00 | 2026-07-25 05:00 | 1260 |
| 3 | 2026-07-28 05:10 | 2026-07-28 06:00 | 50 |
| 4 | 2026-08-04 01:50 | 2026-08-04 03:00 | 70 |
| 5 | 2026-08-11 06:50 | 2026-08-11 08:30 | 100 |

AE来源为NASA CDAWeb OMNI_HRO_1MIN的AE_INDEX，nT，原生1min，WDC Kyoto quicklook。原始CDF位于Z:\SPART-WORK\Data\MMS\ancillary\omni\hro_1min\2026，文件为omni_hro_1min_20260701_v01.cdf及omni_hro_1min_20260801_v01.cdf。

整数FILLVAL=99999在double绘图数组中转NaN并核对原变量，不平滑、不插值、不重新平均；历史脚本采用irf_tlim的[start,end)截取。五段当时均无缺测，跨日Case2标出两天。已交付5张独立PNG和每事件一页的5页PDF。

## MMS1磁尾四小时分段批次

- 日期为2026-07-20 00:00至2026-08-16 00:00 UTC，包含08-15全天。仅MMS1，GSM X<0。
- 官方MEC epht89d原生30s；历史程序用相邻轨道点线性估计X=0进入/离开时刻，从每次实际进入起建立4h网格，保留短尾段，再与日期范围取交集。首个进入在07-19，因此07-20首段为00:00–01:44:41.649，保留原进入时刻的网格。
- 每段三panel为B/Vi/AE，五panel为B/Vi/AE/Ni/离子能谱；标题为时段中点GSM位置，RE=6372km，Ni按该参考图使用红色，离子能谱20–40000eV。
- 09-29历史结果：8次磁尾经过、162段、620.716h、324张PNG。47段有FGM burst；75段完全无可用FPI离子L2，对应panel标No available L2 data，其余局部缺口保留。AE分钟记录全部有效。
- 当时本批MMS原始CDF共800个，MEC30个及科学产品770个；复用185个、补下载615个，均归档Z盘。FPI burst未取得，fast当时到08-08，SDC及CDAWeb均作过空清单和正对照核查。
- 结果子目录MMS1_tail_20260720_0815内含B_Vi_AE、B_Vi_AE_Ni_Ei、index.html、MMS1_tail_windows.csv及MMS1_tail_MATLAB_programs.zip。查询、窗口及验证记录在Z盘derived\MMS1_tail_20260720_0815。

## 当前入口与结果路径更新

2026-09-30已将以下入口的结果根目录从KH改为Recovery-Work_SMILE-MMS，原子目录名沿用；CodeDir、ParentDir、DataDir、RecordDir保持原约定。

| 用途 | MMS_fu内入口 |
| --- | --- |
| 原20事件最终overview | Overview_events_original_style_20260927.m |
| 原20事件最终位置图 | MMS_orbit_original_style_20260924_v2.m |
| EV01 MMS1九panel | Overview_EV01_MMS1_9panels_20260929.m |
| PPT五事件AE | AE_Case_for_SMILE_20260929.m |
| 07-25精确时窗九panel | Overview_MMS1_20260725_0100_0400.m |
| 磁尾进入与4h分段 | MMS1_tail_20260720_0815_intervals.m |
| 磁尾三/五panel | Overview_MMS1_tail_20260720_0815.m |
| 最终PPT输入引用及输出 | MMS_original_style_20260927_ppt.mjs |
| 最终批处理完成记录检查 | MMS_overview_final_driver_20260927.py |

磁尾SkipExisting及最终批处理复用检查增加了当前结果根目录核对，防止Z盘旧完成记录引用KH图片而跳过新目录输出。这属于文件路径/复用判断变更，科学读取、计算与绘图段保持原状。本次未执行重画或生成新PPT。

较早版本、历史ZIP中的代码副本、旧README及TEMP中的打包程序保留作为历史资料。今后如运行它们，先修改实际输出路径与输入结果引用，再核对完成记录；不能只依赖本说明。新交付包应收录届时实际运行的代码。

本地默认交付依据09-29用户原话：“不需要邮箱发送，直接保存即可”。历史私有Drive文件夹及邮件仅为当时交付记录，后续不自动上传或发送。

## 2026-09-30：B_Vi_AE图片汇总为纯图片PPT

用户要求将Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\B_Vi_AE内图片做成PPT，每页一图，禁止新增文字或备注。已按W001至W162顺序嵌入162张原PNG，保持比例、没有裁切、未重绘原图，不添加封面、文字框、页码、备注或评论。

交付：C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1_B_Vi_AE_20260720_0815.pptx。实际构建代码为MMS_fu\MMS1_B_Vi_AE_PPT_20260930.mjs，使用现成Artifact Tool。验证了162页原图哈希与顺序、每页仅一图、无文字和备注、比例与边界；PowerPoint渲染全部162页并逐页与源图作完整尺寸像素比较通过。过程文件与渲染保存在TEMP\MMS1_B_Vi_AE_PPT_20260930，本次仅本地交付。

### 2026-09-30：PPT图片缩小至75%

- 用户要求把每页PPT图片缩小至大约3/4；按当前原PPT图片的宽和高各乘0.75，保留比例并在页面居中。
- 原有162页、图片顺序、图像内容及无新增文字/备注要求保持，原PPT保留。
- 交付：C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1_B_Vi_AE_20260720_0815_75percent.pptx。
- 编辑代码：C:\Users\Administrator\Documents\FWD_matlab\MMS_fu\MMS1_B_Vi_AE_PPT_75percent_20260930.mjs；构建、验证及渲染记录在TEMP\MMS1_B_Vi_AE_PPT_75percent_20260930。
- 已核对全部162页图片字节不变、宽高75%、居中、无新增文字或用户备注，已用PowerPoint渲染全部页面并核对显示。
