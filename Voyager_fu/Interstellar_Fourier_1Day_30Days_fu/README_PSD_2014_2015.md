# Voyager 1：2014–2015 年日均磁场傅里叶 PSD 线图

用户确认使用 2014-01-01 至 2015-12-31；最新要求将横轴改为周期（天），纵轴仍为功率谱密度。

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Fourier_1Day_30Days_fu');
    result = Run_V1_Fourier_PSD_2014_2015;

结果目录：C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_PSD_2014_2015
图件：V1_Fourier_PSD_2014_2015_2day_90days.png / .pdf / .fig

## 时间与数据

计算区间为 [2014-01-01 00:00,2016-01-01 00:00) UTC，共 730 个完整日。仅读取该时段的 24 个原始 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA 月 CDF，使用 Epoch、标量 ABS_B（nT），按原有填充值和有效范围规则读数。

日均值继续按每个 UTC 日内有效小时样本的算术均值计算，不增加覆盖率或样本数阈值。沿用已经授权的内部日缺测线性插值，保留观测 NaN、有效小时计数、插值端点和权重。首尾禁止外推。本次日均值和插值记录与之前完整分析的 2014–2015 子集逐项核验；旧 MAT 只供核验，不作为科学输入。原始记录和结果表均仅包含这两年。

## 整段傅里叶 PSD

将两年日均数据作为一个完整的 730 点记录，减去该记录的算术均值一次，沿用周期 Hann 窗覆盖整段记录，并做 NFFT=730 的 FFT。已在执行前说明整段计算、Hann 窗及去均值。不使用 90 天滑动窗口，也不平均多个窗口的谱；不增加线性去趋势、滤波、频谱平滑或零填充。

复用 V1_Daily_Fourier 和 irf_psd，但输入长度与窗口长度完全相等，仅生成一个频谱。单边 PSD 为：

    abs(FFT((x-mean(x)).*w)).^2 / (Fs*sum(w.^2))

仅将正频率的非 Nyquist 点乘 2；DC 与 Nyquist 不翻倍。Fs=1/86400 Hz，单位 nT²/Hz。继续逐值核对独立 fft 公式和 Parseval 等式：

    sum(PSD)*(Fs/NFFT) = sum(((x-mean(x)).*w).^2)/sum(w.^2)

## 图形与频率范围

一条黑色 PSD 曲线，双对数坐标，横轴为周期（天），从左到右由 90 天递减至 2 天（XDir=reverse）；纵轴 PSD 单位仍为 nT²/Hz。横坐标按 T_days=1/f_Hz/86400 换算，并按周期从小到大排序，PSD 随对应频点一起重排。只改变横坐标，不转换为每单位周期的谱密度。对应频率范围保持：

    [1/(90*86400),1/(2*86400)] Hz
    = [1.286008230453e-7,5.787037037037e-6] Hz

FFT 原生频率间隔为 1/(730*86400) Hz。曲线仅连接原生频率点对应的周期，不人为制造恰好 90 天的点；实际最大周期约 81.11 天，横轴上限保留 90 天。完整频率格及 PSD（包括 DC 和显示范围外频点）均保存在 CSV 和 MAT 中。DC 的 PeriodDays 留为 NaN。

不添加图内方法脚注、拟合线、峰值标签或显著性结论。单次周期图相邻频率处可以有明显起伏；它表示两年数据的整体频谱，Hann 权重及已有线性补缺会影响估计。

## 输出与审计

保存 PNG、矢量 PDF、MATLAB FIG、完整分析 MAT、730 日输入、插值审计、全频率 PSD 表、一个整段窗口的均值和归一化审计、CDF 路径和 SHA256、代码 SHA256、MATLAB 版本及日志。图形中的周期按原生频率精确换算，PSD 按对应顺序与分析数组逐项核验，完整功率数组与修改前完全一致。


