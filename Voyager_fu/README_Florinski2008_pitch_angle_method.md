# Voyager 1 LECP 投掷角：Florinski et al. (2008) 处理方法

主程序：`Voyager_Interstellar_Monthly/Voyager1_Replot_Selected_Events_PaperStyle.m`

2020-E01 运行入口：`run_voyager1_2020e01_lecp_p1_pitch_angle_daily_from_hourly.m`

## 已实现的论文步骤

1. 使用 LECP Ch1 / P1 的 0.57–1.78 MeV 质子数据。
2. S1–S7 分别按 UTC 日形成 24 小时平均；S8 由遮阳板挡住，不作为有效观测方向。
3. 同一天的 COHOWeb 一小时 `BR、BT、BN` 分量分别取算术平均，得到日平均矢量 `<B>`。
4. LECP 扫描面采用论文给出的名义几何：扫描面近似平行 RT 平面，八个扇区中心间隔 45°，入射粒子方向与望远镜视轴相反。
5. 对每个有效扇区计算

   `mu_k = vhat_k dot <B> / |<B>|`

   `PA_k = acos(mu_k)`

   归一化使用完整三分量磁场，`BN` 不被丢弃。
6. 不插值，不做三日平均，不做二次多项式拟合，也不人工调整 MAG 方位角。
7. 当七个有效扇区不齐全、投掷角覆盖过窄、磁场方向在日内变化过大或相对统计误差过大时，该日会在 CSV 中标为不可用。

## 背景校正说明

论文说明“扣除了穿透宇宙线背景”，但没有公布逐日背景序列或自动校正公式。LECP 仪器资料说明，S8 可辅助估计背景，背景值还需要结合上下文确定。公开的自动化 LECP 产品也给出了未做背景校正的质量警告。

因此 2020-E01 默认设置为：

`LECPBackgroundMode = 'none'`

程序保留 NASA Level-2 的公开逐扇区日平均通量，S8 只写入诊断列。图中与 CSV 中都会标出 `public Level-2 background uncorrected`。这项限制影响通量的绝对值和各向异性幅度，不改变由扇区方向和磁场矢量确定的投掷角。

程序另提供：

`LECPBackgroundMode = 's8_daily_approx'`

该选项在日平均后用同日 S8 值减去 S1–S7，仅作为明确标注的近似方法。它不等同于论文作者使用的 LECP 团队背景产品。

## 依据

- Florinski, Decker, and le Roux (2008), JGR, doi:10.1029/2007JA012859
- NASA CDAWeb Voyager 1 LECP Level-1 / Level-2 数据说明
- Voyager LECP Data Analysis Handbook
