# Cases for SMILE：聊天承接与当前工作规则

整理日期：2026-09-30，北京时间。当前结果根目录：C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS。

当前代码目录（2026-10-03更新）：SMILE 卫星仪器程序位于 C:\Users\Administrator\Documents\FWD_matlab\SMILE_fu；MMS 与辅助分析程序位于 MMS_fu。下文2026-09-30的统一代码目录规定已由此细化。

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

### 2026-09-30：离子流速完整性核查与CDAWeb补齐

- 用户要求检查 B_Vi_AE 目录的流速空白是否由网站有数据但本地漏下载引起。本次只核查、补齐原CDF及输出核查表，未重画PNG或更新PPT。
- 2026-09-30实际新查SDC：07-19至08-17的MMS1 FPI fast L2 dis-moms仍155个、burst为0；file_info与file_names一致，原155个文件本地版本/大小均匹配。
- 同时实际新查CDAWeb：请求期07-20至08-16有204个fast文件、burst为0。扩展边界范围有206个（另含08-16两个范围外文件），与四个分周查询并集一致；burst用2015-10-16已知事件得到56文件作为阳性对照。
- CDAWeb存在49个原本未下载的fast文件，均为08-10至08-15，逐日7、12、5、2、12、11个。它们没有出现在本次SDC清单中。CDAWeb LastModified为09-29，但该字段不能确定首次公开时间，也不能据此反推昨天查询时已经可取得。
- 已从CDAWeb补齐49个原CDF至Z盘既有mms1/fpi/fast/l2/dis-moms/2026/08目录，总计95,966,211字节。下载中遇到不完整传输，由大小校验拦截并重试；最终49个大小与官网完全一致，CDF头及MATLAB完整读取均通过。未变更原155个CDF。
- MATLAB/IRFU基线：原155个CDF有219,482个有限速度样本，全部162个时段与旧绘图计数一致，速度数值和点数未因数据库读取或GSE→GSM而遗漏。
- 新49个CDF有71,349个有限速度样本，约新增89.187小时覆盖；新样本首末为08-10 08:51:03.758至08-15 23:09:19.183 UTC，中间保留空档。
- 实际受影响26个图件窗口/纯图片PPT页：W128–W140、W144、W150–W161。原图这些窗口的Vi都是0点；补齐后18段完整覆盖、8段局部空档。W134（08-11 06:21:28.920–10:21:28.920 UTC）可读到3200点，原图仍显示No available L2 data。
- 当前204个CDF合计290,831个有限速度样本。162个时段的可用数据状态为71段完整、42段局部空档、49段全空；原状态为53、34、75。全部CDF无读取错误，GSM转换丢点0。这里的覆盖核查未增加质量筛选、平滑或插值。
- 现有PNG和两份PPT仍是核查前版本。后续重画应只对上述26个WindowList设置SkipExisting=false，采用原绘图方法；之后更新对应PPT页并保留用户要求的75%图片尺寸。原science_complete及W*_overview记录仍描述旧图，不能将其当作本次204文件的完整清单。
- 用户核查表：Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815\MMS1_Vi_completeness_check_20260930.csv，列出162页原图点数、补齐后点数、剩余空档及需更新标记。
- 查询、下载校验、补齐前后MATLAB审计及来源证据：Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815\Vi_completeness_20260930；根summary.json汇总结果。当前162张PNG的哈希与制作PPT时的输入清单全部一致。
- 实际代码在MMS_fu：MMS1_Vi_sdc_audit_20260930.py、MMS1_Vi_cdaweb_audit_20260930.py、MMS1_Vi_CDAWeb_supplement_20260930.py、MMS1_Vi_local_audit_20260930.m、MMS1_Vi_audit_report_20260930.py。运行日志在TEMP。

### 2026-09-30：新增Vi数据的26张三panel图及指定75%PPT已更新

- 用户明确要求重画新增数据的图片，并只替换指定PPT内有变化的图片，其他内容不动。于北京时间19:01完成原文件覆盖交付。
- 已更新 B_Vi_AE 目录中的W128–W140、W144、W150–W161，共26张PNG；其余136张PNG逐文件哈希不变。原脚本的科学读取与处理、时间、坐标、配色及版式保持；26图Vi样本数均与补下载后的MATLAB审计一致，真实空档保留。26图新旧尺寸一致，磁场和AE区域逐像素相同。
- 已在原路径更新 MMS1_tail_20260720_0815\MMS1_B_Vi_AE_20260720_0815_75percent.pptx。以用户18:44保存的当前版本为源，包含其6处文字（第37、50、67、87、96、110页），全部保留。PPT仍162页、原图片尺寸和位置不变，仅26个PNG媒体部件改变，其余801个包内部件逐字节不变，包括页面XML、关系、文字、备注、母版及属性。
- PowerPoint原生渲染全部162页完成：136张未改页逐像素一致；26张改页的差异仅在原图片边界内，边界外差异0像素。26张更新页已逐页目视核对。最终PPT SHA256：adb51f0fe016ef77e4b39e33d86d13ca065bd13b3f2f289e4cb4357105f8643a。
- 本次未更新五panel图或100%尺寸PPT；它们仍为补下载前版本。原science_complete、原W*_overview及完整性核查CSV保留其历史状态，新三panel记录单独保存在下述redraw目录，后续勿将旧记录误作当前图件的数据清单。
- 实际代码：MMS_fu\Overview_MMS1_BViAE_refresh_20260930.m、MMS1_BViAE_PPT_replace_20260930.mjs、MMS1_BViAE_PPT_assets_20260930.py、MMS1_BViAE_PPT_publish_20260930.py。过程文件、完整原PPT和26张原PNG备份在TEMP\MMS1_BViAE_refresh_20260930。
- 验证与交付证据：Z:\SPART-WORK\Data\MMS\derived\MMS1_tail_20260720_0815\Vi_completeness_20260930\redraw，其中redraw_verification.json核对数据和图件，ppt_independent_verification.json核对PPT包及原生渲染，delivery_receipt.json记录已覆盖的最终路径及哈希。

### 2026-09-30：提取MMS1.pptx中标注1的页面

- 用户指定源文件为结果目录MMS1_tail_20260720_0815\MMS1.pptx，要求将标注1的页面单独做成新PPT。
- 已提取原第37、50、67、87、96、110、130、134、135页，共9页，按原顺序保存为同目录MMS1_标注1.pptx。原图、可编辑的标注1、尺寸、位置及所有页面内容保持，未增加文字。源文件未修改。
- 源SHA256：1eeeabceda87be48bf254ab9b0e990856ab6efe3b83c6c907317857b00ed67ed；输出SHA256：3d096dbad2814f109fce0ff9cb064edccc3de372538c7ba44466f87df6e323eb。
- 九页原slide XML及图片逐字节一致；PowerPoint原生渲染与对应原页逐页0像素差异。输出包只含所选9页和9张图片，无未选图页残留。
- 代码MMS_fu\MMS1_select_marked1_20260930.py、MMS1_marked1_20260930.mjs、MMS1_marked1_render_20260930.ps1；过程及验证记录TEMP\MMS1_marked1_20260930。

### 2026-10-01：对平均、插值及绘图抽稀的准确说明

- 用户询问此前这些图是否使用平均/插值。已核对原磁尾绘图脚本、26图重画脚本及实际IRFU调用链：官方CDF读取后的B、Vi、AE数值未额外做时间平均、平滑或插值补缺，也未先将三种物理量插值到统一时间轴。Vi有GSE→GSM坐标转换，IRFU按CDF积分时间信息校正FPI采样中心；缺测/断线按NaN处理。
- 必须同时说明显示处理：B、Vi的irf_plot调用含'reduce'。本地reduce_to_width按坐标轴像素宽度分箱，数据密集时选每箱最小值、最大值用于显示，并将显示横坐标放到像素时间网格。此处不取物理量均值，也不补造缺测物理量，但PNG并非将每个原始采样点及其原生时间戳逐点保留。AE调用未启用reduce，仍为原CDF1min采样。
- 轨道计算确有线性插值：MMS1_tail_20260720_0815_intervals.m用interp1估算GSM X=0的进入/离开时刻，用irf_resamp(...,'linear')计算每段中点的标题XYZ位置。因此以后不要笼统回答整个流程“没有任何插值或降采样”；需区分物理量时间序列、显示抽稀及轨道位置计算。
- 本次只核查并解释，未改动科学计算、图件或PPT。


### 2026-10-01：五个指定事件的四星磁场电流图

- 用户发来五张B/Vi/AE图，要求保留原三个panel，在Vi和AE之间插入电流强度、分量两个panel。附件对应W003、W025、W050、W051、W082，使用windows.json原精确时窗，边界保留到纳秒表示，不自行延长。
- 顺序为B、Vi、|J|/|J_parallel|/|J_perp|、GSM Jx/Jy/Jz、AE；单位nA/m^2。B/Vi仍为MMS1，标题仍为原图中点位置。AE继续用原CDF原生1min。08-02图原有Vi局部缺口保持。
- 参考用户Plot_TCS.m的8参数c_4_j调用，直接使用IRFU c_4_j、irf_dec_parperp、irf_abs；场向参考为c_4_j返回的四星平均磁场。四星FGM和MEC均为GSM，输入nT/km。
- 四星burst同时有效时优先burst电流，其他时刻用四星survey电流。新增电流对齐使用irf_resamp(...,'linear')，在有限观测连续段内对齐到本模式MMS1磁场时间，不外推。FGM连续性取自官方bdeltahalf采样宽度，MEC按本批原生30s；容差仅为epochUnix浮点精度。未平滑、未做时间平均，未添加divB/curlB或几何质量门槛，诊断值单独记录。
- 原B/Vi显示段沿用原图已有3*dt断线及irf_plot reduce，以保留附件显示流程。该旧阈值没有用于新增电流连续性判断。电流显示也采用reduce，AE未采用reduce。
- 正式结果目录：Recovery-Work_SMILE-MMS\MMS_current_5events_20261001，含5张PNG、5份独立矢量PDF及MMS_5events_B_Vi_J_AE_20261001.pdf（5页，附件顺序）。原图及PPT未修改。
- 实际主程序：MMS_fu\Overview_MMS_current_5events_20261001.m，直接脚本/%%分区，无新增科学算法function。辅助程序MMS_current_5events_download_20261001.py复用既有SDC查询下载接口，MMS_current_5events_pdf_20261001.py仅合并/检查PDF。下载候选104个原CDF，复用50个、补下载54个，按原产品层级归档Z盘，保留官网查询、大小及CDF头校验。
- 五段原B/Vi及AE点数均与原附件基线一致；总/场向/垂直分解恒等式通过。已知电流线性合成磁场验证单位和方向，最大误差1.31e-13 nA/m^2；缺失输入保持NaN。PNG和Poppler渲染的PDF五页逐页目视核对通过。
- 记录在Z盘derived\MMS_current_5events_20261001，包括windows、manifest、各事件_current.json、algorithm_verification.json及delivery_verification.json；过程日志/渲染页在TEMP\MMS_current_5events_20261001。


### 2026-10-01：五图第3栏总电流、第4栏带符号场向电流

- 用户随后指定每图第3个panel只画总电流强度，第4个panel画平行/垂直电流、去掉绝对值，以观察平行/反平行方向。已在同一主程序及同一五图结果路径更新完成。
- 当前顺序为B、Vi、|J|、J_parallel/J_perp、AE。第3栏黑色总电流；第4栏蓝色带符号J_parallel=J·Bhat，正值沿四星平均磁场、负值反向；红色J_perp为垂直矢量的模，天然非负。图例不加绝对值符号，并画零基线。带符号的垂直标量需要另指定垂直投影轴，本次未引入额外投影轴。
- 直接调用IRFU irf_dec_parperp的有符号场向输出，已删除原abs(J_parallel)。四星c_4_j计算、输入数据、对齐流程、单位、时间窗及burst/survey选择保持。
- 五张图原B/Vi区域及AE/UTC区域与本次修改前PNG逐像素相同；电流有效行数保持，五事件均有负J_parallel样本。已逐页检查五页PDF渲染，面板顺序、正负曲线、标签及缺口显示通过。
- 本次未补下载数据。原版本备份与过程文件在TEMP\MMS_current_5events_20261001\before_signed_parallel；最终记录更新于Z盘derived\MMS_current_5events_20261001\delivery_verification.json，包括面板对比、正负范围、输出及程序哈希。


### 2026-10-01：五事件磁场和电流的时间分辨率核对

- 用户询问磁场和电流精度。已用MATLAB mms.get_data重新核对原始时间戳：五事件survey主体16Hz（约62.5ms），W025、W050局部burst为128Hz（约7.8125ms），电流输出对齐到对应MMS1磁场时刻，无时间平均或平滑。
- W082的MMS4存在288个约125ms相邻间隔（局部8Hz）；该段四星输入的原生时间分辨率受MMS4限制，不能将16Hz电流输出当成新增独立四星16Hz观测。W082 MMS1–3还有少量非等间隔时间戳，按原数据保留。
- 原生MEC位置30s，通过已明确的线性插值对齐到磁场时刻；B和电流图均使用IRFU reduce显示抽稀。采样间隔与幅值测量误差须分开表述，当前五事件未计算逐点电流幅值不确定度。
- 核对记录：Z盘derived\MMS_current_5events_20261001\sampling_resolution_verification_20261001.json；只读核查程序和日志在TEMP，本次未改动图件、数据或正式计算程序。


### 2026-10-01：五图垂直电流曲线改为50%不透明度

- 用户要求降低垂直电流的不透明度，让平行电流在重叠处清楚可见。已将五图第4栏红色J_perp设为0.5不透明度，蓝色带符号J_parallel维持1.0。原图、各独立PDF及5页合并PDF在同一正式路径更新。
- 主程序继续用IRFU绘图及reduce产生显示顶点，全部坐标范围设定完成后，用MATLAB官方patch/EdgeAlpha接口将红色曲线半透明显示；使用rectangle裁剪，NaN保留断线，线宽不变。仅改变显示，科学计算记录逐字段与修改前一致，其他4个panel逐像素相同。
- PDF描边透明度确认为128/255（0.501960784，MATLAB 8位量化）；五页最终PDF渲染已逐页目视检查通过，红蓝重叠、裁剪、图例及UTC标签清楚。
- Overview_MMS_current_5events_20261001.m及PDF合并核对程序已同步更新。PNG和PDF先完整导出到TEMP，再替换正式文件。原版本备份位于TEMP\MMS_current_5events_20261001\before_perpendicular_alpha；最终核对与哈希在Z盘derived\MMS_current_5events_20261001\delivery_verification.json的opacityVerification。


### 2026-10-01：五图电流改为每1分钟一个平均值，并核对实时磁场方向

- 用户要求全部电流取1min平均，并核对平行/垂直分解是否相对于实时磁场方向。已更新同一MATLAB主脚本、5张PNG、5份单页PDF和5页合并PDF；正式结果路径保持MMS_current_5events_20261001。
- 先在每个原生电流时刻调用irf_dec_parperp(Bmean,J_B)进行分解，再分别平均总电流幅值、带符号J_parallel及垂直矢量幅值J_perp。Bmean是该时刻四星磁场的空间平均，方向逐时刻变化，未使用固定背景磁场；运行时检查其时间轴及数值逐行等于同一时刻四星空间平均，所有实际使用的数据模式均通过。
- 从每张原图起点开始分为不重叠60s时间窗，边界为(left,right]；直接调用IRFU irf_resamp(...,'mean','window',60,'thresh',0)，不做sigma排除。先按原burst优先规则选取有效样本，每个有限样本等权，分别平均三个已分解量。未先平均磁场再分解，也未取平均矢量的模替代平均总电流幅值。
- 前4段各240个分钟均值，W082为78个；W082末窗34.811311s保留实际数据均值，显示时间为实际窗的中点。空窗保留NaN，本批所有窗均有有效输入。当前电流仅有240/78点，绘图不再使用reduce；原B/Vi及AE的显示流程保持。
- 本次原生电流统计与修改前逐字段一致；五张PNG原B/Vi与AE/UTC区域逐像素相同。红色J_perp仍为0.5不透明度，蓝色J_parallel保留正负。W082两个电流panel的图例移到左上方，避开右侧峰值；最终PDF逐页目视核对通过。
- 每窗参数、数值、burst/survey样本数及实时参考检查保存在Z盘derived\MMS_current_5events_20261001各事件_current.json和delivery_verification.json（minuteAverageVerification）。未补下载数据。修改前版本备份及过程日志在TEMP\MMS_current_5events_20261001\before_1min_average及同级TEMP文件，正式结果目录仅保存交付图件。


### 2026-10-01：去掉五图电流图例中的1min标识

- 用户要求把图注中1min标识全部去掉。主脚本仅修改两处irf_legend文本，三种电流图例现在为|J|、J_parallel、J_perp；一分钟分箱平均、实时磁场分解及0.5垂直电流不透明度保留。
- 同一正式路径更新5张PNG、5份独立PDF及5页合并PDF。五事件科学记录与修改前完全一致，PNG在两栏电流图例区域以外逐像素相同；全部PDF页面提取文本确认无1min标识，最终五页渲染已逐页核对。
- 验证及更新哈希保存在Z盘derived\MMS_current_5events_20261001\delivery_verification.json（legendVerification）。修改前备份及运行/渲染过程文件在TEMP\MMS_current_5events_20261001\before_remove_1min_legend及同级TEMP目录；未补下载数据。


### 2026-10-01：框选EV01区间原始分辨率四星电流密度

- 用户发来07-21磁场/流速/AE拼接图，要求计算框选时段电流密度、原始精度。图中边界约00:04-03:48；用户随后明确选用整段2026-07-21 00:00-04:00 UTC。本次只对该新事件计算，已有五图一分钟平均结果保留。
- 新直接MATLAB脚本Overview_MMS_current_20260721_0000_0400_native.m复用Plot_TCS.m的8参数c_4_j及此前已核对的四星对齐、irf_dec_parperp、irf_abs。GSM磁场nT、MEC位置km，电流nA/m^2；逐时刻四星空间平均磁场作为分解参考，运行时对应时间轴与平均值检查通过。未取固定背景，未做时间平均/平滑或新增质量筛选。
- 四星burst原生128Hz（7.8125ms），survey约16Hz（62.5ms）；有效burst电流786822点，burst优先后使用survey132049点。原生MEC位置30s，沿用在观测连续段内的线性对齐，不外推或跨缺口补值。电流图完整原生时间戳绘制，不使用reduce；官方FGM采样宽度识别的时间间断仅在显示数组中插NaN断线，原生数值与时间轴核对未改变。
- 新结果目录Recovery-Work_SMILE-MMS\MMS_current_20260721_0000_0400_native_20261001，含EV01_20260721_0000_0400_MMS1_B_Vi_J_AE.png和同名矢量PDF（约46.1MB）。五栏顺序B、Vi、|J|、J_parallel/J_perp、AE；蓝色场向电流带符号、红色垂直幅值保持0.5不透明度，图例无平均标识；AE仍为原CDF原生240个分钟样本。
- 原生衍生电流、分量和逐时刻参考B保存到Z盘derived\MMS_current_20260721_0000_0400_native_20261001\EV01_20260721_0000_0400_native_current.mat。仅保存用户所求的衍生电流，不转换全仪器CDF。MATLAB重新读取该文件并核对原生时间范围、严格递增、有效行数和缺测对应关系通过；最终PNG及PDF渲染目视检查通过，PDF透明度128/255，验证/哈希在同目录delivery_verification.json。
- 已使用既有SDC查询/下载接口核对本批数据，最终计算清单192个原CDF均复用本地已有文件。官网日期结束边界返回了次日07-22日文件，初次清单多含8个FGM/MEC日文件，其中6个新下载、2个此前已有；已收紧候选筛选，本次计算清单排除这8个文件，保留官方原文件在Z盘归档供后续复用，记录于download_scope_verification.json。运行/QA文件在TEMP\MMS_current_20260721_native_20261001。
- 本次原生总电流幅值最高298.1966nA/m^2（burst，03:19:12.831 UTC）；该数值为未质量筛选的curlometer结果。用户图内文字仅为图件内容，未把其中KH/极光解释作为新的科学指令或验证结论。


### 2026-10-02：EV01原生电流图放大到03:07-03:15 UTC

- 用户要求把03:07-03:15这段时间单独绘图，按上文事件2026-07-21 UTC处理。最终直接截取既有四小时原生电流MAT中的J、|J|、J_parallel、J_perp、逐时刻参考B及采样宽度，保留原生时间与数值，不做平均/平滑或质量筛选。
- 初次尝试窄时段重算的逐位相等校验未通过，未交付该尝试。最终改为已验证四小时结果直接截取；MATLAB重新读取最终数据文件，以上六组数组与四小时结果在03:07-03:15范围逐元素isequaln检查全部通过。
- 新直接脚本为MMS_fu\Overview_MMS_current_20260721_0307_0315_native.m。B/Vi/AE仍用原CDF和IRFU绘图，AE为原生1min（8个样本，沿用[start,end)），中点位置标题为03:11 UTC MMS1位置。电流采用burst46856个有效点（128Hz）及其余survey1823点（约16Hz）；场向电流保留正负，垂直幅值0.5不透明度，无平均标识。
- 结果目录Recovery-Work_SMILE-MMS\MMS_current_20260721_0307_0315_native_20261002内保存EV01_20260721_0307_0315_MMS1_B_Vi_J_AE.png及同名矢量PDF，五栏B/Vi/|J|/J_parallel,J_perp/AE。最终PNG与PDF渲染均已检查通过。
- 衍生原生子集、参数、严格一致核对和最终哈希在Z盘derived\MMS_current_20260721_0307_0315_native_20261002；过程日志与QA在同名TEMP子目录。未下载新数据，四小时图及前五张一分钟平均图保持原结果。


### 2026-10-02：03:10-03:15原生/一分钟电流双版本，垂直电流alpha=0.7

- 用户指定2026-07-21 03:10-03:15 UTC，要求原始精度及1min两版，并把垂直电流透明度参数改为0.7。本次设置MATLAB EdgeAlpha=0.7（PDF量化为179/255=0.701960784）；蓝色带符号场向电流保持不透明。
- 原生电流直接截取已验证的四小时原生数据，本段38399个有效burst样本，128Hz；survey在burst优先后无有效选中点。最终保存子集与四小时原生值六组数组逐元素核对通过。
- 1min版沿用已确认流程，先逐时刻相对于实时四星空间平均磁场分解，再对|J|、带符号J_parallel及垂直幅值分别取IRFU irf_resamp(...,'mean','window',60,'thresh',0)算术平均。自03:10开始五个不重叠60s窗，显示在03:10:30、11:30、12:30、13:30、14:30；各窗有效点数7680、7680、7680、7680、7679。逐窗与保存原生样本直接算术平均核对通过。
- 两图五栏仍为B、Vi、|J|、J_parallel/J_perp、AE。B/Vi/AE不新增平均，两版原三栏像素完全相同；AE原生1min，沿用[start,end)的5个样本。图例不加1min标识。两份最终PDF渲染均检查通过。
- 实际直接脚本为MMS_fu\Overview_MMS_current_20260721_0310_0315_native_1min.m，依次导出两个版本，使用原脚本%%分区，无新增科学function。结果目录Recovery-Work_SMILE-MMS\MMS_current_20260721_0310_0315_native_1min_20261002，保存EV01_20260721_0310_0315_MMS1_B_Vi_J_AE_native及_1min的PNG/PDF各一份。
- 原生/分钟衍生数据、参数、一致性/均值核对及最终哈希在Z盘derived同名分类，日志和PDF渲染QA在TEMP同名目录。未下载新数据，既有03:07-03:15、四小时和五事件结果保持原文件。


### 2026-10-02：2019-08-05 16:24-16:25 UTC原生四星电流

- 用户指定2019.8.5/16:24-16:25，只算原始精度。已按UTC计算四星FGM curlometer电流；未做时间平均、平滑或新增质量筛选。此前事件图保持原结果。
- 本地四星FGM burst原CDF完整覆盖本段，均使用mms1–4_fgm_brst_l2_20190805162313_v5.202.0.cdf；原生时间分辨率128Hz（中位间隔7.8125ms）。电流对齐到MMS1磁场时间，7680个有效样本，首尾16:24:00.001919677至16:24:59.994823242 UTC。四星survey未参与该段电流。
- 本地MMS1 MEC已存在；补下载MMS2–4当日mec_srvy_l2_epht89d_v2.1.1原CDF，四星原生位置30s，经已确认的连续观测段内irf_resamp线性对齐，无外推或跨缺口补值。补下载2019-08 OMNI_HRO_1MIN原CDF；所有原数据按产品层级保存Z盘。
- 新直接MATLAB脚本MMS_fu\Overview_MMS_current_20190805_1624_1625_native.m，仿照原overview的%%分区，复用用户Plot_TCS.m的8参数c_4_j、IRFU irf_dec_parperp与irf_abs。GSM磁场nT/位置km，输出nA/m^2。场向参考为每个原生时刻的四星空间平均磁场，运行时核对其时间轴及逐行空间平均值通过；场向电流保留正负，垂直电流为垂直矢量模。原生电流绘图未启用reduce。
- 五栏顺序B、Vi、|J|、J_parallel/J_perp、AE；MMS1 B/Vi沿用原生数据和IRFU既有显示流程。垂直电流EdgeAlpha=0.7，PDF实际量化179/255=0.701960784。AE采用[start,end)原生1min，只有16:24一个有效值（1073nT），以点标记显示，未人为延长、平均或插值。
- 正式结果目录Recovery-Work_SMILE-MMS\MMS_current_20190805_1624_1625_native_20261002，含MMS_20190805_1624_1625_MMS1_B_Vi_J_AE.png及同名矢量PDF。总电流最大20.448976531nA/m^2，16:24:17.338064208 UTC；场向范围-11.882644330至14.053925856nA/m^2，垂直幅值范围0.006811127至20.437928849nA/m^2。上述为未加质量筛选的计算结果。
- 衍生J矢量、|J|、带符号J_parallel、垂直幅值和逐时刻参考B保存在Z盘derived\MMS_current_20190805_1624_1625_native_20261002\MMS_20190805_1624_1625_native_current.mat。MATLAB读回验证原始时间轴、有效样本和实时分解通过，最大分解恒等式残差2.2737e-13。最终PNG及PDF最新Poppler渲染目视核对通过，参数、下载清单、来源和哈希在同级derived记录；日志和渲染QA留TEMP同名目录。
### 2026-10-03：读取用户提供的两个 SMILE UVI L2 CDF

- 用户要求用本机 MATLAB 读取两个 UVI CDF，并明确告知任何假设。源文件仍在 Recovery-Work_SMILE-MMS 根目录：SMILE_UVI_L2_AURORA-GEO_20260720T230011(1).cdf、SMILE_UVI_L2_AURORA-GEO_20260720T230041(1).cdf，未改写或移动。它们属于 SMILE，未套用 MMS 数据目录约定。
- MATLAB R2026a 实跑成功。使用本地 IRFU 自带 spdfcdfinfo/spdfcdfread，KeepEpochAsIs=true；23变量/文件、每变量1记录全部读取。用 MATLAB cdflib.getVarRecordData 逐元素核对全部数值；本文件 ROW_MAJOR 导致其返回数组与 NASA 接口相差一次转置。正式 UVI 工作区保留 NASA 返回排列与数值类型。字符串 ITEN 的原生尾部 NUL 在 NASA 接口中为空格，文本一致。
- EPOCH按TT2000解码为2026-07-20T23:00:26.018023000Z与23:00:56.018377000Z；曝光分别23:00:11.018023000至23:00:41.018023000、23:00:41.018377000至23:01:11.018377000 UTC，均30s，EPOCH严格等于曝光中点，两帧间隔30.000354s。NASA与MATLAB时间分解一致；沿用与文件标记一致的2017-01-01最后更新闰秒表，未独立验证星上授时精度。
- 三类图像为SMILE_UVI_IMAGE 521x521（属性称corrected）、GRID_IMG_GEO 141x141、GRID_IMG_AACGM 501x501，单位属性均Rayleigh。投影高度变量110km，CADENCE=30s，CALIB_COEF=1；这些均直接来自文件，未重新投影、定标、归一化、平滑或插值。
- 本地 spdfcdfinfo.Variables 第10–12列对部分double变量返回异常FILLVAL/VALIDMIN/MAX，不能使用这些摘要值。VariableAttributes与MATLAB cdflib.getAttrEntry逐项一致，正式脚本使用直接读取的真实属性。常见FILLVAL实际为single(-1e31)，转为double是-9.9999998482432073e30；精确比较，仅在绘图副本转NaN。两帧GEO图分别686/692个FILLVAL，AACGM分别19846/19742个。IMAGE无FILLVAL，所有零值保留。
- PROJ_FLAG明示1-Earth/0-Space且FILLVAL也写0；保留全部0/1。两帧PROJ_GEO_LAT/LON的-999位置与FLAG=0逐像素相同。SC_POS描述写GSE、坐标属性写J2000；SC_ATT描述写GSE、坐标属性写Spacecraft Coordinate System，4分量却单位写deg；坐标与姿态含义未确认，未做转换。ROI区域定义、ITEN模式具体含义、keogram逐点坐标轴未提供，未自行补定。
- 实际直接脚本：SMILE_fu\Read_SMILE_UVI_CDF_20261003.m，%%分区、无自定义function。运行后工作区有UVI/UVIInfo/UVIAttr/TimeUTC，未另存全量MAT或CSV。结果目录Recovery-Work_SMILE-MMS\SMILE_UVI_readout_20261003含SMILE_UVI_two_frames_array_preview.png及SMILE_UVI_readout.txt。
- 预览为两行三列，按NASA读取后的数组行列显示，行号向下增加，未解释上方为北向或正午；各列两帧共享全范围线性色标，灰色为精确FILLVAL。图已目视核对，完整属性限制/假设写入TXT。诊断脚本、原始元数据、运行日志、交叉核对与源/交付SHA256留TEMP\SMILE_UVI_CDF_inspect_20261003。默认沙箱工具初始化失败，本次通过获准的沙箱外本机命令执行，未改科学计算工具。

### 2026-10-03：可配置路径和产品类型的 UVI 读取程序

- 用户明确要求可在程序中输入文件地址和文件类型，并确认“文件类型”为UVI产品类型，如AURORA-GEO；代码仿照原MMS overview的直接脚本/%%分区。
- 新主脚本SMILE_fu\Read_SMILE_UVI.m：顶部FilePath接受单个CDF完整路径或目录；FileType='AURORA-GEO'按标准产品文件名筛选目录，再核对GlobalAttributes.Logical_source，跳过空项及DOI。目录仅当前层、按文件名排序。ic选择当前提取常用变量的文件。
- MATLAB cdfinfo读取可靠元数据，避开本地spdfcdfinfo.Variables第10–12列问题；NASA spdfcdfread逐变量读取，KeepEpochAsIs=true，不统一不同变量记录数。UVI{n}保存各文件变量，CDFInfo{n}保存元数据；TT2000另生成TimeUTC{n}，原int64保留。其他时间类型仅保留原值，不擅自解码。
- Data=UVI{ic}。AURORA-GEO分区直接提取Epoch/EpochUTC、曝光起止、Cadence、Image/ImageGEO/ImageAACGM及各坐标、投影标志。当前两个文件结构已实测；其他产品通用读取保留原变量，未声称已验证其专用物理解释。
- 不更改原数组、类型、填充值、零值或-999，不平滑、不插值、不转置、不重新定标，不绘图或自动导出全量MAT/CSV。临时工作目录在TEMP；成功及产品不符退出均恢复原工作目录。字符串尾部NUL被NASA接口显示为空格的已知行为写入notes。
- MATLAB R2026a实跑：目录模式2文件共46变量与原生cdflib逐元素核对通过（考虑已验证的接口转置和字符串padding差异）；单文件模式读取第二文件通过，确认692/19742个图像FILLVAL和184411个空间标志0原样保留；错误产品类型被明确拒绝；checkcode无输出问题。测试副本和日志在TEMP\SMILE_UVI_CDF_inspect_20261003，正式旧读取脚本、图件和原CDF保持原文件。

### 2026-10-03：SMILE程序集中到SMILE_fu

- 用户要求在FWD_matlab下新建SMILE_fu，并把SMILE程序放入。已创建C:\Users\Administrator\Documents\FWD_matlab\SMILE_fu，将Read_SMILE_UVI.m和Read_SMILE_UVI_CDF_20261003.m从MMS_fu移动到该目录。逐文件SHA256核对一致，程序内容未变；今后SMILE卫星仪器新程序统一保存到SMILE_fu。
- 本次分类以实际处理对象为依据。名称含SMILE的历史MMS quicklook下载程序处理MMS图片，AE_Case_for_SMILE读取OMNI AE，因此归在MMS与辅助分析目录。联合项目承接记录继续位于MMS_fu\SMILE_MMS_project_context.md。
- AGENTS.md的代码目录约定及本记录中的UVI程序入口已同步更新。源CDF、已有图件和结果输出目录保持既定位置。程序没有依赖旧代码目录的内部路径，本次无需修改读取或计算逻辑。移动校验和说明文件备份位于TEMP\SMILE_UVI_program_move_20261003。
### 2026-10-03：清理读取脚本注释，新增逐文件三图主程序

- 用户要求Read_SMILE_UVI只保留大块%%注释，并新增调用该读取程序的主程序；对已提供两个CDF分别读取，每文件一张图、同图包含三类图像。
- Read_SMILE_UVI.m已删除全部普通整行与行尾注释，仅留%%分区。为主程序调用移除原clear;clc，FilePath/FileType/ic/IrfDir仅在不存在时赋默认值；原读取与时间解码保留。每次重置特定常用输出，避免切换产品后残留上次图像变量。
- 新入口SMILE_fu\Main_SMILE_UVI.m同为直接脚本/%%分区，无新增function。顶部设置InputDir、InputFiles、ProductType、OutputDir。循环分别设置FilePath并run同目录Read_SMILE_UVI.m，累计Data/Info/UTC保存到AllUVI/AllInfo/AllUTC，原数组保留。
- 主程序对每文件生成独立1x3图，依次Corrected image、GEO grid image、AACGM grid image，灰色为精确FILLVAL，数组行向下增加，不赋予地理方向；两图同类panel共用全范围线性色标，单位读取原属性。图题显示EPOCH与曝光区间UTC。仅绘图副本double化和FILLVAL转NaN，无新增平均、平滑、插值、背景扣除、重投影或定标。
- 已生成Recovery-Work_SMILE-MMS\SMILE_UVI_readout_20261003\SMILE_UVI_L2_AURORA-GEO_20260720T230011(1)_three_images.png及...230041(1)_three_images.png；原两帧合并预览保留。运行主程序会显示两张MATLAB figure并保存PNG。
- MATLAB实跑与独立核对通过：两文件46变量对原生cdflib逐元素相符（沿用已确认的数组排列及字符padding表示差异），6个图像CData与原数组的精确FILLVAL副本一致，2张图各3个image，AlphaData、YDir、共同CLim与两帧时间均核对通过。已目视检查两图标签、色标与布局。checkcode的主程序numel==1建议已改isscalar，未改变行为。
- 注释删除前读取脚本备份、运行和验证日志留TEMP\SMILE_UVI_main_20261003。主程序只支持本次每文件一帧三类图像的绘图结构，读取脚本仍按原变量读取；本次未改变原CDF或旧科研图。

### 2026-10-03：SMILE UVI 主程序改用 jet 色图

用户要求提高红蓝颜色对比度。Main_SMILE_UVI.m 中唯一修改为 colormap(Ax,jet(256))，已用本机 MATLAB 重新生成 SMILE_UVI_readout_20261003 内两个 *_three_images.png。已验证两张图的六个 panel 均使用 jet(256)，色标保持线性、原数值范围不变，坐标覆盖完整图像数组。数据读取、FILLVAL 显示规则和其他绘图设置保持原状，本次未增加科学假设。修改前的程序和图件备份、运行日志位于临时目录 SMILE_UVI_main_20261003。
### 2026-10-04：UVI 色标上限降为此前的 3/5

用户要求将当前色标上限调低至3/5。Main_SMILE_UVI.m在计算两帧共同有效值范围后增加 ColorLimits(:,2) = ColorLimits(:,2)*3/5;，各下限保持原值。三个新上限依次为35947.2、998.295703125、963.397106831 Rayleigh，两帧同类panel仍共用范围，使用jet(256)线性色图。已用MATLAB重新生成SMILE_UVI_readout_20261003内两张*_three_images.png，并核对六个panel上限、下限、色图和完整坐标范围通过。超过新上限的数据按最高端颜色显示，数组数值及FILLVAL显示规则沿用原设置，无新增科学假设。已目视检查两图；修改前程序、图件备份及日志保存在TEMP\SMILE_UVI_upper_limit_20261004。
### 2026-10-04：仅第一个 UVI panel 上限再降为当前的 3/5

用户要求只调整第一个panel。Main_SMILE_UVI.m增加 ColorLimits(1,2) = ColorLimits(1,2)*3/5;，该上限由35947.2降到21568.32 Rayleigh。GEO、AACGM上限仍为998.295703125、963.397106831 Rayleigh。已用MATLAB重新生成两张*_three_images.png，验证两图三个panel的上限并完成目视检查；读取、数据数值、完整坐标、jet及下限均沿用原设置。未新增科学假设。修改前程序、图件与运行日志留TEMP\SMILE_UVI_first_panel_limit_20261004。