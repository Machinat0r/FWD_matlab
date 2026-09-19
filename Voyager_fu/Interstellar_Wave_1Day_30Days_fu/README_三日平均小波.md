# Voyager 1 三日平均磁场 Morlet 小波

## 运行

在 MATLAB 中运行：

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Wave_1Day_30Days_fu');
    result = Run_V1_ThreeDay_Morlet;

正式入口默认从 Z:/SPART-WORK/Data/Voyager/voyager1/coho/1hr/l2/merged_mag_plasma 的原始 CDF 读取。结果输出到 C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_ThreeDay_30Days。已有日均版入口和结果独立保留。代码依赖同目录的 V1_Wave_Read_Hourly、Plot_V1_ThreeDay_Morlet，现有 Case1_PPT_VerticalLine_Events_7d_fu 的 CDF 读取工具和本地 IRFU。

## 数据和用户指示

2026-09-11 用户要求用此前三天数据再做同样分析。本次沿用此前已要求的线性补缺、Morlet 小波、以天表示周期和色标范围 [-1,3]（用户最新修改）。执行前已说明三日平均定义、空三日窗口插值、六天 Nyquist 边界和末端不完整窗口处理。

时段为 2012-08-25 00:00 至 2021-12-17 00:00 UTC（后者不含）。读取 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA 的原始 Epoch 与 ABS_B，按既有源元数据处理填充值、有效范围。标量 |B| 单位为 nT。

先对每个 UTC 日内有效小时 |B| 求算术平均，再从 2012-08-25 开始每三天分一个不重叠窗口，对窗口内有限日均值等权求平均。这与 Run_V1_ThreeDay_Monthly_Overview 的磁场 panel a 完全相同；各天等权，不按小时样本数加权，不设置最低覆盖率。每个窗口记录实际有效天数及源小时数。

三日平均前保留原始日缺测；平均完成后仅对完全空的完整三日窗口做线性插值。通过 irf_resamp(...,'linear') 连接缺口两侧最近的原始有效三日值，禁止端点外推，不修改已有三日均值或观测计数。ObservedWindows 和 BMeanObserved_nT 保留原始 NaN，IsInterpolated 及插值审计表保留估计值、两端观测、权重、连续空窗口数。插值值不代表新增观测；三日平均及补缺会影响局部谱功率。

末尾 2021-12-15 至 2021-12-17 是两日的不完整窗口，保留原有实际中点和磁场曲线数值，不纳入等间隔三日小波输入。全部完整三日窗口形成一条连续序列。若序列长度为奇数，保留 irf_wavelet 原生舍去最后一个样本的行为，并记录返回样本标志，不改写实际时间戳。

## 小波和显示

直接调用用户 MMS_fu/Wave.m 所用的 irf_wavelet：Morlet 宽度 5.36，100 个对数频率，cutedge=1，采样率 fs=1/(3*86400) Hz。不额外做滑窗估计、去均值、去趋势、滤波或功率平滑。当前标量 |B| 功率保持 IRFU 归一化 P=2*pi*abs(W)^2/f，单位 nT²/Hz。

一个月按 30 天。三日采样的 Nyquist 频率为 1/(6*86400)=1.929012345679e-6 Hz，计算频率下限为 1/(30*86400)=3.858024691358e-7 Hz，显示周期 6–30 天，6 天在图上方。六天恰为采样边界，该处的振幅和相位识别受限。采用真实原生三日时间色块。纵轴标为“周期（天）”；色标为 log10 功率，范围 [-1,3]，功率单位仍为 nT²/Hz。两端保留 IRFU 随频率变化的原生边缘屏蔽，不新增图内方法注释。

沿用此前电场产品核查结论，所需低频段没有可用电场数据，本次输出磁场图。

## 核验和交付

正式运行直接读原始 CDF；旧 three_day_monthly_audit.mat 仅用于核验三日均值、时间戳、完整性及有效天数完全一致，不作为科学计算输入。程序独立核对分组均值、观测值不变、插值公式、采样率、Nyquist 边界及原生边缘缺测掩码，并用每三天采样的 12 天正弦信号检查周期换算。

输出 PNG/PDF/FIG、完整 MAT、原始日统计交付表、三日输入及标志表、插值审计、末端样本审计、CDF 路径和 SHA256、代码 SHA256、MATLAB 版本和运行日志。日统计 CSV 和 MAT 供交付及审计。可通过 DataRoot、OutputRoot、IRFURoot、ReferenceAuditFile、Visible 名值参数更改路径或可见性。

