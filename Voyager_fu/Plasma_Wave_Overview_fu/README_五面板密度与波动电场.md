# Voyager 1：日均概览、电子密度与波动电场

## 运行

在 MATLAB 中运行：

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Plasma_Wave_Overview_fu');
    result = Run_V1_Plasma_Wave_Overview(true,true);

第一个参数控制是否核对官方清单并补档；第二个参数控制 PWS 原始 CDF 的并行读取。默认补档并直接读取原始 CDF。已有完整数据库时可将第一个参数设为 false，仍重新读取原始 CDF 计算，不以输出 MAT/CSV 代替科学输入。并行工具箱不可用时串行运行。

补档程序 Sync_V1_PWS_Archive.ps1 仅负责 HTTPS 传输、官方文件清单和文件散列记录，科学计算及绘图全部由 MATLAB 完成。默认使用本机已验证可用的 PowerShell 7；在其他电脑运行时可在入口的 PowerShellExe 路径参数中指定 PowerShell 7。新下载文件只写入 Z 盘 Voyager 数据归档。

## 路径与时间

- 正式代码：C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Plasma_Wave_Overview_fu。
- Voyager 数据：Z:\SPART-WORK\Data\Voyager。
- 太阳黑子原文件：Z:\SPART-WORK\Data\Solar_Indices\raw\sunspot\SN_d_tot_V2.0.csv。
- 结果：C:\Users\Administrator\Documents\Recovery-Work-Voyager_betatron\V1_Plasma_Wave_Overview。
- 所有 panel 共用时间范围：1990-01-01T00:00:00Z 至 2025-07-01T00:00:00Z（右端不含）。
- 新图另存 V1_1990_20250630_daily_PWS_5panels.png/.pdf/.fig，保留已有概览及其穿越时间标记版本。
- PNG 与 PDF 均为 220 dpi；PDF 使用光栅化导出，避免约 20 万频谱色块产生超大矢量文件。FIG 保留可编辑图形。

## 五个 panel 的定义

(a) SILSO 原始日太阳黑子数，源值 -1 保留为缺测，0 为有效值。

(b) COHO 原始小时标量 ABS_B 的 UTC 日算术均值；沿用既有程序仅在该文件 ABS_B 无任何有限值时回退 F 的规则。单位 nT。

(c) COHO 原始小时 protonFlux1_LECP 的 UTC 日算术均值。P1 显示标签沿用 0.57–1.78 MeV；源通量及源能段元数据不重标定。绘图使用对数纵轴。

前三个 panel 的科学数值和缺测分布分别与上一版 all_daily_values.csv 核对。该旧表仅用于验证。日均点位于 UTC 12:00，细线连接相邻有效日；缺日以 NaN 断开。

(d) PWS 电子数密度。用户于 2026-09-14 明确同意直接采用官方原始 CSV，此授权限于该密度源格式。读取 SCET 和 N_e，保留每条原始记录。EPO 和 QTN 在同一 panel 中分别绘成点序列并设正常图例；不合并两类来源、不设置来源优先级、不插值、不自行从频谱反演密度，也不跨缺测时段连线。N_e_min 和 N_e_max 随原表保留在审计中，当前图没有误差棒。纵轴单位 cm^-3。

(e) PWS 低速谱 electric_field，单位 V/m。逐个原始 CDF 读取 epoch（CDF TT2000）、frequency 和 electric_field。使用源文件给出的 FILLVAL、VALIDMIN 和 VALIDMAX，有限有效原始样本按 UTC 日、按频道分别计算算术均值；不增加有效样本数或覆盖率阈值。空日、空频道留为 NaN。频率轴为真实 16 频道中心频率，显示色块边界取相邻中心的几何中点，仅用于显示，不代表仪器精确带宽。频率轴和电场色标采用对数显示，均值计算在线性 V/m 上完成。

电场色标范围 5e-7–1e-3 V/m 沿用已核验源 CDF 的 SCALEMIN/SCALEMAX，颜色仅在显示范围外饱和；源值与输出统计保持完整。输出电场为波动幅值，不提供直流电场矢量，不通过 -v×B 或模型补出未测参数。

V1 PLS 于 1980 年停止工作，本图时段缺少其直接测量的流速和温度；用户已选择省去这两个 panel。密度和电场覆盖不足时保持空白，不生成推测值。

## 官方源与归档

- NASA/JPL 仪器状态：
  https://www.jpl.nasa.gov/news/how-do-we-know-when-voyager-reaches-interstellar-space/
- PWS CDF 元数据（VG1_PWS_LR；DOI 10.17189/wp0z-1c51）：
  https://spase-metadata.org/NASA/NumericalData/Voyager1/PWS/CDF/PT4S
- SPDF 实际 CDF 目录：
  https://spdf.gsfc.nasa.gov/pub/data/voyager/voyager1/wave_spectra_pws/spectrum_analyzer_cdf/
- 官方密度（bundle 2.0；DOI 10.17189/zkkc-we77）：
  https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/
- 密度文件：
  https://pds-ppi.igpp.ucla.edu/data/voyager-pws-vlism-density/data/vg1-vlism-density-2012-2025.csv

实际归档的 Voyager 1 密度产品标签为版本 2.0，755 条记录；覆盖 2012-10-22T21:44:19.795Z 至 2025-04-14T17:27:35.917Z。官方 CSV MD5 为 9b2578494f7e1f07c99e1cdfe0297df4，读取时强制核对。该产品的真实终止时间取自 V1 标签和数据，不使用同时含 V2 数据的 bundle 总覆盖来代替。

PWS CDF 归档沿用：
voyager1/pws/calibrated/spectrum_analyzer/native/YYYY。

密度 CSV、标签及原始集合文档归档：
voyager1/pws/derived/electron_density/native/PDS_release_20260910。

官方逐年目录快照、选定文件清单、URL/大小/SHA256/下载状态及汇总归档：
voyager1/pws/source_verification/overview_1990_20250630。

补档以官方网站实际列出的文件为准。同日多版本选择编号最高版本；现有有效 CDF 不重复下载。程序校验 CDF 文件签名和 HTTP 长度，MATLAB 随后逐个解析真实 CDF；异常文件会停止运行。相同路径的更新文档或 CSV 若发生内容变化，旧内容保存在其 previous_versions/SHA256 子目录，不静默丢弃旧版本。

## 校验与结果记录

- 先检索并参考了 IRFU irf_cdf_read，使用既有 IRFU dataobj 与项目 CDF 读取器。
- MATLAB 路径、读数、统计、绘图分节沿用 MMS_fu 人工程序风格，未引用 codex 子目录。
- COHO 重复 Epoch 仅在所用 B/P1 数值完全一致时计一次；冲突停止并要求核查。
- PWS 检查每个文件的真实记录都位于文件标称 UTC 日内、Epoch 无重复、16 频道一致；异常情况停止，禁止静默丢弃或调整时间。
- first_three_preservation.csv：前三个面板的数值一致性与最大绝对差。
- coverage.csv：逐参数的首末有效时间及有效日数。
- V1_PWS_electric_field_daily_values.csv：16 频道每日均值、有效样本数、是否存在官方 CDF。
- V1_PWS_density_selected_records.csv：实际使用的原始密度记录、来源和原文件行号。
- V1_plasma_wave_overview_audit.mat：方法、完整来源、源属性、逐文件样本/排除数、散列、密度原表及标签、源清单、覆盖和输出清单。
- 原图 PNG/PDF/FIG 和旧日统计表在运行前后计算 SHA256，并断言完全未变。
- 本任务不新增图内方法说明、脚注或小字注释；仅保留正常标签、单位、图例与 panel 编号。

本次已完成 36 年度 SPDF/Iowa 日期及版本交叉核对，双方均列出 12,687 份本时段 CDF。原生补档程序清单解析也独立核验为相同的 12,687 份文件。源清单对比报告随数据归档，具体下载与全量绘图完成状态见结果目录最终运行记录。

## 本次完成记录

本次全部 12,687 份 CDF 已补齐并由 MATLAB 直接读取，共 28,098,881 条原始时间记录。日均电场有效 12,687 天，278 天按官网源产品缺测留空。密度 755 条（EPO 155、QTN 600），覆盖 742 个 UTC 日。最终五面板 FIG 与统计数据逐值一致，原图与旧日统计表 SHA256 完全不变，详见结果目录 README_结果与复现.md 和 final_validation.json。

Validate_V1_Plasma_Wave_Overview 可独立核对已生成的 FIG 与统计审计；此入口用于交付验证。V1_Format_Plasma_Wave_Figure 仅整理标题和图例位置。默认完整重画仍使用 Run_V1_Plasma_Wave_Overview，直接读取原始 CDF。
