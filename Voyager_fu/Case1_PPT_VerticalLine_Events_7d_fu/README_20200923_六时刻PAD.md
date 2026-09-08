# 2020-09-23 事件：用户框选的六个小时 PAD

事件：Voyager 1，Case1-S01-L02。按照 2026-09-04 用户截图中框出的六根有效小时 PAD 色带，固定选择以下时刻，按时间从左到右排列：

1. 2020-09-22 22:30 UTC
2. 2020-09-22 23:30 UTC
3. 2020-09-23 00:30 UTC
4. 2020-09-23 01:30 UTC
5. 2020-09-23 02:30 UTC
6. 2020-09-23 03:30 UTC

这些时间是既有 Level-1 UTC 小时均值的格中点，对应 [22:00,23:00)、[23:00,00:00) 等半开小时区间。此图不进行峰值搜索，也不重新选择邻点。其后两条不完整 L2 记录不属于这六个有效 PAD。

## 绘图方式

沿用已有类 Florinski 点图的排版：六个并排面板，横轴投掷角 0–180°，黑色点和误差棒，每格保留 S1–S7 七个独立扇区点。纵轴为该小时的 J/Jmax，Jmax 为该小时 S1–S7 中最大通量。误差棒为 sigma/Jmax，归一化分母视为固定，不额外传播其误差。每格记录该格的未归一化 Jmax，单位 cm^-2 s^-1 sr^-1 MeV^-1。

不拟合，不连接散点，不合并投掷角，不进行新的平均、插值、磁场补值或背景扣除。S8 不进入七扇区点图。没有额外的科学筛选阈值。

## 数据处理与已有图的对应

直接选取原小时总览 MAT 的 pitchAngleTable 第 14–19 行，保证与截图使用同一数值。执行入口在绘图前直接读取原始官方 CDF，重新计算仅用于验证；绘图仍使用通过核对的原表数值。

LECP 原始文件：
Z:\SPART-WORK\Data\Voyager\voyager1\lecp\native\l1\sectored_rates\2020\voyager-1_lecp_lev-1-rates_20200101_v1.1.1-01.cdf

LECP 变量为 Epoch、DeltaT、FHDU_SectoredRates 与 FHDU_SectoredRateUncertainties，使用已有 P1 通道映射（FHDU 第 10 通道）。沿用获批处理：排除 DeltaT<0 的记录；按原 Epoch 所属 UTC 小时，对每个扇区的有限非负计数率求算术均值，缺测不作零，真实零计数率参与平均。各扇区样本数允许不同。L1 均值按历史换算 J=mean(R)/[0.44*(1.78-0.57)]，R 为计数率；这里沿用已有获批换算，未声称它与官方 L2 标定等价。

独立样本误差假设下，均值误差为 sqrt(sum(sigma_i^2))/N，再乘相同换算因子。任何参与平均的样本误差缺失时，该扇区均值误差保持缺失。

磁场原始文件：
Z:\SPART-WORK\Data\Voyager\voyager1\coho\1hr\l2\merged_mag_plasma\2020\09\voyager1_coho1hr_merged_mag_plasma_20200901_v01.cdf

使用 Epoch、BR、BT、BN，同一 UTC 小时内完整有限的 RTN 矢量均值。六个小时各有一条有效磁场记录。采用已有获批的官方预测姿态、日采样插值及官方标称 LECP 扇区/安装几何；预测姿态不等同于遥测重建姿态，实际安装偏差校准值缺失的限制沿用原方法说明。以粒子速度方向（仪器望向的负向）计算与 B 的三维夹角，未施加 uN=0 或固定 RT 平面。此轮重新核对了三分量方向和投掷角，未改变姿态通路。

显示能段沿用用户指定的 0.57–1.78 MeV；源 CDF 的能段元数据与文件保持不变。

## 复现与保存

MATLAB 入口：
C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Case1_PPT_VerticalLine_Events_7d_fu\Run_Case1_Selected6_20200923_PAD.m

通用固定时刻绘图函数：
C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Case1_PPT_VerticalLine_Events_7d_fu\Case1_Plot_Selected_PAD.m

图像：
C:\Users\Administrator\Documents\Recovery-Work-Voyager_betatron\Case1_PPT_VerticalLine_Events_7d\Selected_PAD_times_hourly\V1_Case1-S01-L02_20200923_PAD_selected6_20200922T2230_20200923T0330.png

数值与复现审计：
Z:\SPART-WORK\Data\Voyager\voyager1\lecp\1h\derived\selected_pitch_angle\2020\V1_Case1-S01-L02_20200923_PAD_selected6_20200922T2230_20200923T0330.mat

审计保留原始扇区通量、误差、投掷角、归一化分母及结果、原始记录索引、样本数、磁场和三维姿态核对结果、母图来源、代码与输入文件哈希、输出图哈希和像素尺寸。全部卫星数据留在 Z 盘分类目录，不导出 CSV。新图单独保存，既有日均、小时和五时刻图不改动。

