# 日球层顶外六图与V2磁场缺测核查

入口：Run_Voyager_Extended_Overviews(1:2,false,'heliopause')。false仅跳过前一任务已完成的官方归档同步；科学计算和默认重画仍直接读原始CDF。正式程序沿用Interstellar_Daily_Overview_fu下的通用入口，新增Voyager_Supplement_V2_MAG原始磁场补充函数。

范围：V1 [2012-08-25,2025-07-01)；V2 [2018-11-05,2025-01-01)，均为UTC，起点采用日球层顶穿越日期，按日精度包含穿越当天。此处“出了太阳系之后”按项目语境解释为日球层顶外。日期来源NASA https://science.nasa.gov/mission/voyager/interstellar-mission/ 。三日窗口分别从各自穿越日锚定，月窗按自然月，首尾截取保留。日/三日/月定义与此前图一致；三日中位数池化原始小时P1，月图为日中位数算术平均。日图底部包含S4，三日/月图排除S4；始终排除S8。

## V2空缺原因与数据库核查

原COHO日图2020-12-31至2021-12-31连续缺少磁场日值。其中2021年COHO ABS_B和F为源CDF填充值；例如202106文件中两个字段各720个小时记录均为填充值，非有效范围阈值剔除。现有程序只使用COHO，导致独立MAG可用数据未显示。

官方网站存在独立reviewed VIM 48秒MAG CDF，2021年文件：
https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager2/magnetic_fields_cdaweb/vim_48secmag/voyager2_48s_mag-vim_20210101_v01.cdf

正式原始文件已在 Z:/SPART-WORK/Data/Voyager/voyager2/mag/48s/reviewed_vim/2021/ 内。本次重新下载官方文件，SHA256与本地完全一致：07F09BDF255B0790EDCDB71BA5B0ED9043E87D675907C1B29D70B589B85A26DB。无缺失文件或过期副本需要替换；核查元数据归档到数据目录source_verification/V2_MAG_gap_2021/official_MAG2021_verification.json。科学绘图补全读取通路，保留原始COHO文件内容不变。

2021年MAG CDF含141600条有效F1标量磁场记录，从2021-01-01 03:02:51至2021-11-26 20:56:42 UTC。采用既有CDF填充值/有效范围筛选，逐UTC小时对有效48秒F1求算术均值，再对有效小时等权求UTC日均；仅原COHO日均缺失时采用该值。保留全部有效COHO日值，且不对有效日里的缺失小时混合补值。三日和月图使用补充后的同一日序列。无插值、平滑、覆盖率阈值或推测值。source_recompute/daily_with_MAG_supplement_audit.mat保存每条原始48秒记录索引、源SHA256/元数据、小时和日均样本数及采用掩码。

V2日球层顶外共补回338天：2018年2天、2021年329天、2022年7天。2021年剩余36天为8月12日和11月27日至12月31日，已核对的COHO/MAG CDF无有效记录，继续留空。整个V2绘图窗口仍有159个磁场缺测日。没有将“所核对产品中无数据”等同于所有仪器从未观测。

图内遵循用户约定，不新增处理小字；方法与来源只在本说明、代码和审计中保存。保留五面板、薄线、线性磁场轴和对数粒子轴。图件、统计与日志统一保存在Recovery-Work-Voyager_betatron/Voyager_Heliopause_Overviews。
