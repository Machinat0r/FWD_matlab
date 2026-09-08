# 当前日均事件图：前后各 15 天、L1 优先

2026-09-04 用户要求将每个事件的 daily 平均图改为前后各 15 天。本入口覆盖当前全部 47 个 Voyager 1 事件：原 45 个事件，加 2020-07-30 与 2020-08-22 两个已建立的补充事件。

## 窗口与修改范围

单日事件 D 的绘图窗口为 [D-15日00:00, D+16日00:00)，共 31 个完整 UTC 日。所有 panel 共用这一窗口；原始事件日期和 EventID 不改。

仅原位更新 Epoch_daily 中的 47 张日平均总览图以及对应的 1d 派生数值审计。小时总览、五时刻 PAD、手工指定的六时刻 PAD、独立 Florinski Figure 4 图均保持不变。旧目录名中的 7d 是历史命名，窗口以图题、横轴和本次审计为准。

这里的 daily 图沿用现有六 panel 格式：磁场保留原小时测量及日均曲线，panel e 保留现有粒子通量曲线，底部为每日 PAD；本次不把所有 panel 改成仅显示日均值。

## 科学处理保持不变

- 完整有效的 S1–S7 L1 UTC 日格均值整格优先；L1 缺失或不完整时保留已有 L2 产品记录，不逐扇区混合产品。
- 原 Epoch 所属 UTC 日内，每个扇区有限非负计数率取算术均值；负 DeltaT 记录舍去，缺测不当零，真实零计数率参与平均。不增加样本数、时间覆盖、角跨度、磁场方向 RMS 或不确定度筛选阈值。
- L1 沿用已批准的历史换算 J=mean(R)/[0.44*(1.78-0.57)]，误差 sqrt(sum(sigma_i^2))/N 乘同一换算因子；参与误差缺失时误差保持缺测。这一历史换算未声称与官方 L2 标定等价。
- L2 保留原始 Epoch、通量、不确定度和元数据；仅舍去 DeltaT<0 的记录，不再重平均、拆分或扩展产品累计区间。
- L1 使用日格中点及同 UTC 日完整磁场矢量均值；L2 在原 Epoch 使用其所属日的磁场均值。磁场和通量不插值或人工补值。
- 投掷角沿用用户批准的官方预测姿态、日采样插值及官方标称安装/扇区几何，计算三维粒子方向与 B 的夹角。未施加 uN=0 或固定 RT 平面。预测姿态不等同于遥测姿态，实际安装偏差校准未知的限制保留。
- S8 不参与 S1–S7 PAD，不扣除背景；PAD 仅要求 S1–S7 正通量及有限投掷角。
- P1 显示能段继续为 0.57–1.78 MeV，源能段元数据保持原样。当前自动坐标和色限规则沿用原程序，在新窗口数据上求值。

## 数据覆盖与复现

所有正式卫星数据继续放在 Z:\SPART-WORK\Data\Voyager。扩展窗口整体涵盖 2012-12-25 至 2021-06-01；需 65 个既有 COHO 月 CDF，以及 2012–2021 年 LECP L1 和 L2 daily 原始 CDF。新补的 2012 年文件使用官方 SPDF 原始 bytes 和原文件名，分别放入 native/l1/sectored_rates/2012 和 1d/l2/sectored_flux/2012，不转存 CSV。

MATLAB 正式入口：
C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Case1_PPT_VerticalLine_Events_7d_fu\Run_Case1_Daily_Context15Days_L1_First.m

入口从原 45 事件函数及两个已有补充事件清单建立合并清单，通过既有 Run_Case1_PPT_VerticalLine_Events_7d / Voyager_Case1_Plot_Events 直接读取 CDF 出图。原 Case1_Event_Catalog 和小时入口不改。每个所需年度显式预检，避免漏掉整年源文件。

本轮审计保存在：
Z:\SPART-WORK\Data\Voyager\voyager1\lecp\validation\context15d_daily_l1_first

每次运行保留原日均数值基线、非日均输出文件哈希、输入源文件清单及哈希、全部 47 事件运行报告。验证逐个检查 31 日窗口、当前选项、旧七日窗口内的通量/误差/PA/磁场数值一致性、PNG 可读性，以及小时和独立 PAD 产物未变。

