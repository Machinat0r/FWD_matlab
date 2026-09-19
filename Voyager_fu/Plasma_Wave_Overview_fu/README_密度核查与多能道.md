# V1 密度核查与多能道日均图

正式 MATLAB 入口：Run_V1_Density_Multienergy_Overview。
验证入口：Validate_V1_Density_Multienergy_Overview。
密度审计函数：V1_Audit_Density_Completeness。
官方源文件核验：Check_V1_Density_Official_Archive.ps1。

时间固定为 [1990-01-01T00:00:00Z, 2025-07-01T00:00:00Z)。
输出到 C:\Users\Administrator\Documents\Recovery-Work-Voyager_betatron\V1_Density_Multienergy_Overview。
生成 LECP 三能道概览，以及 CRS 十五能道分三组概览，共四张 PNG/PDF/FIG。删除电场面板；此前五面板结果不覆盖。

所有磁场/粒子科学值直接读取 Z:\SPART-WORK\Data\Voyager 下原始 COHO CDF，复用 IRFU dataobj 和已有 Voyager_Read_CDF_Product。太阳黑子数沿用 Z:\SPART-WORK\Data\Solar_Indices 的 SILSO 官方日文件。密度使用用户于 2026-09-14 批准的官方 CSV 例外。

各质子通道分别对 UTC 日内有效 COHO 记录作算术平均，零值有效，NaN/填充值不参与。不插值、不补零、不做背景扣除、不合并能道。CRS 原生说明为 6 小时值，COHO Count 是合并文件有效条数。日志/审计记录每个 CDF 的原始能段和单位，读取时逐文件核验能段，所有源通量保持原单位。COHO 的 LECP P1 标签为 0.57–1.78 MeV；本任务不改写其他 LECP 扇区产品的能段元数据。

磁场沿用先 ABS_B、整文件没有有效 ABS_B 时使用 F 的已有规则。日通量用对数轴与细线连接，非正日均值仅在对数显示数组中隐藏；统计表原样保存。缺测日断线。

密度使用现行正式 PDS 版本全部原始 SCET/N_e。EPO/QTN 分别保留点序列，不做日均、不平滑、不插值，不重反演波形，不去重。版本审计保留重复次数，对旧 Annex 与现行 CSV 的 SCET + Source + N_e 多重集合进行一一比较。官方勘误未完整解释两版密度记录差异，正式图不合并旧版额外记录。

密度网站复核快照固定为 source_check_20260915，当前源文件 MD5 与配套标签相符。以后更新源产品需同步审核版本、源文件与标签；核验脚本不会静默替换当前源。其余图件方法、准确覆盖时间、空缺日期段、755 条原始点验证、新旧版差异、输出清单详见结果目录 README_密度核查与多能道图.md。

默认科学计算不读取派生 MAT/CSV 日统计。既有日均 CSV 仅用于原图前 3 条统计的保留性验证；验证入口读取本次审计 MAT 仅用于质量检查。
