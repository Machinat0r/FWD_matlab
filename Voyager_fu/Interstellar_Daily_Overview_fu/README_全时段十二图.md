# V1/V2 延长时段的十二张五面板图

入口：Run_Voyager_Extended_Overviews。默认参数 spacecraft=1:2，syncArchive=true；默认入口核对官方目录并下载缺失原始CDF，随后直接读取原始CDF计算。单颗可调用 Run_Voyager_Extended_Overviews(1,false)，false只跳过下载目录同步，科学数据仍从原始CDF读取。

程序目录：C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu。输出：C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/V1或V2。CDF及官方清单位于Z:/SPART-WORK/Data/Voyager，按航天器/仪器/级别/分辨率/年份分类。官方清单位于source_verification/V1_full_overview及V2_full_overview。图件和统计结果不写入数据归档。

## 科学定义

每颗航天器两种范围：2008-01-01至最新有效COHO观测；用户于2026-09-11指定的长时段图从1990-01-01开始，沿用all文件名。末端为最后有效日次日00:00的半开边界。原始CDF读数审计仍保留完整任务覆盖，绘图窗口单独记录。各仪器覆盖差异保留为空，不能把COHO最大日期当作所有仪器均有效的日期。

五面板依次为太阳黑子数、总磁场标量平均、COHO P1平均、P1中位数面板、LECP扇区求和。日图中位数为原始小时通量日中位数；三日图中位数对三日窗口内原始小时P1值直接求中位数；月图中位数面板对日中位数求月算术平均。其余面板在三日/月窗口内对有限日值等权平均。三日窗口从各图起点开始不重叠；月窗口按UTC自然月；首尾截取，不足完整窗口也保留，审计记录天数和有效样本数。不同变量在同窗口允许有效日期不同，不增加共同有效日筛选。

底部沿用此前各图：日图S1-S7包含S4；三日/月图S1,S2,S3,S5,S6,S7排除S4；始终排除S8。每一天要求所选扇区全部有限，避免缺测求和被当零。原始L1/L2完整S1-S7整格L1优先规则、负DeltaT舍弃和历史换算J=mean(R)/(0.44*(1.78-0.57))沿用原项目授权，未声称历史换算与官方L2等价；原始源能道元数据和源不确定度保留。扇区求和不涉及PA或姿态近似。

太阳黑子数直接读本地SILSO原始发布CSV：Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv，WDC-SILSO, Royal Observatory of Belgium，V2.0，第5列日总黑子数，-1留空，零保留；源元数据说明与日图入口一致。无传播时移、去趋势、平滑、插值或异常值剔除。

COHO采用原始ABS_B，文件无有效ABS_B时按旧入口回退F；P1使用protonFlux1_LECP。保留源CDF填充值/有效范围处理。若官方原始文件重复覆盖相同Epoch，要求所选磁场/P1值完全一致（含NaN位置）才保留一份，重复记录及来源写入审计；出现冲突立即停止，不隐式选优。每个源文件保留元数据和SHA256。每日计算采用accumarray等价UTC分组，避免全历史逐日遍历原始数组的额外开销。

## 图形与验证

太阳黑子数、磁场线性轴，三个粒子面板对数轴。点以0.4宽细线连接；NaN断线，对数显示隐藏非正值但审计保留。图内只保留正常标题、坐标轴、单位与panel编号，不加处理脚注。显示范围内的原始高场/高通量不新增剔除。

每颗航天器输出6组PNG/PDF/FIG及窗口CSV、extended_overviews_audit.mat、覆盖表、原始CDF派生与逐年扇区源审计。每个窗口保存实际/有效天数、三日中位数小时计数。有限样本计数守恒、窗口边界和代表性窗口的独立直接统计及图形YData检查在运行时断言。默认不从MAT/CSV派生值计算或重画。

开始实现前已检索IRFU的Voyager相关实现，沿用现有IRFU dataobj CDF读取器和项目原始扇区读数/优先级实现。代码风格沿用MMS_fu人工程序与当前日统计入口的路径、读数、计算、绘图分节。

## 2026-09-10 V2原六图磁场更新

用户要求将2008年起与全部时段的V2六张图全部重画补齐磁场。正式Run_Voyager_Extended_Overviews入口现对所有V2范围采用Voyager_Supplement_V2_MAG，同一MAG实测补充通路用于日球层顶外和长时段图。原CDF直接读数、48秒F1到UTC小时均值再到日均值，仅替换COHO完全缺失的日均值；缺少实测的日继续为空，不插值。整个范围内可用reviewed MAG原始文件均参与核查；原有效COHO日值不改变。三日/月磁场重新从补充后的日序列聚合。

六图PNG/PDF/FIG、窗口CSV、主审计和覆盖表原位更新；旧窗口CSV只作为验证基准保存在V2/mag_update_20260910/before_values，不用于科学输入。补充来源、逐原始CDF记录、小时/日样本数及替换掩码见V2/source_recompute/daily_with_MAG_supplement_audit.mat。验证报告在V2/mag_update_20260910/validation_audit.mat；校验其余面板和窗口不变、原有效磁场日值不变、日球层顶外重叠段与已核验图一致。

早期重复记录核查：reviewed MAG原始文件中发现2011-12-05 19:11:29（F1为0.13482和0.13407 nT）与2018-08-29 00:00:30（两条均0.37753 nT）的重复Epoch。两天原COHO日均均有效，按既定规则不采用MAG补充。现补充函数先按COHO缺测日选择MAG候选，候选中若存在重复Epoch仍停止；没有对冲突记录做平均或任意选优。非候选重复记录单独保存在UnusedDuplicateRecords，原始CDF不改写。

## 2026-09-11 长时段图改为1990年起

旅行者1号和2号原all组的日、三日、月图均从1990-01-01 00:00 UTC开始，终点沿用现有原始CDF覆盖。原位更新六组PNG/PDF/FIG、窗口CSV和审计。2008年起及日球层顶外图件保持原有版本。

仅重画这六张图的调用：`Run_Voyager_Extended_Overviews(1:2,false,'all')`。false跳过网站目录同步，科学计算仍直接读取原始CDF。默认extended入口中的all范围也已改为1990年起。三日窗口以1990-01-01重新锚定，三日中位数仍直接汇总原始小时P1；其余统计定义、S4选择、缺测规则和V2实测MAG补充保持原样。图内不新增处理说明。

旧统计仅作验证基准保存在结果目录start1990_update_20260911/V1与V2；不用于新图科学计算。主审计中的from2008条目原样保留为此前交付的结果记录。验证入口Validate_Voyager_From1990_Overviews检查日/月与旧数据在1990年后逐项一致、所有三日磁场/黑子数/P1平均及原始小时P1中位数、六图五面板数据/坐标/连线、V2已补回的340个磁场日值及2008组记录不变。另以SHA256核对2008组图件未改动。检索过IRFU irf_resamp，其自动插值行为不用于本任务；沿用原CDF读取器和明确UTC窗口分组。