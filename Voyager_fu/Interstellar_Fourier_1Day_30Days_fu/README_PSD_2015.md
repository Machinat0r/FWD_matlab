# Voyager 1：仅 2015 年的日均磁场傅里叶 PSD

用户要求只用 2015 年数据重画当前 PSD 线图。保留整段 Hann 窗、去均值和倒序周期横轴。

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Fourier_1Day_30Days_fu');
    result = Run_V1_Fourier_PSD_2015;

结果目录：C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_PSD_2015
图件：V1_Fourier_PSD_2015_2day_90days.png / .pdf / .fig

## 数据与计算

读取 [2015-01-01 00:00,2016-01-01 00:00) UTC 的 12 个原始月 CDF。保留既有 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA 标量 ABS_B 的日均定义：对每个完整 UTC 日的源有效小时样本求算术均值，时刻置于中午，不设置额外覆盖率阈值。

共 365 个日均点，其中 355 个为有效观测日均值，10 个缺测日沿用此前授权的线性插值。年初和年末均有真实观测，无端点外推。原始 NaN、有效小时计数、插值两端数据和权重保留，并与原 2014–2015 结果的 2015 年子集逐项核对。旧 MAT 仅供核验，不作为科学数据输入。

对整个 365 点序列减去一次算术均值，应用周期 Hann 窗并执行一次 NFFT=365 的 FFT。不滑窗、不平均多个频谱、不做额外线性去趋势、谱平滑或零填充。复用 V1_Daily_Fourier 和 irf_psd，单边 PSD 归一化为：

    abs(FFT((x-mean(x)).*w)).^2 / (Fs*sum(w.^2))

Fs=1/86400 Hz。DC 不翻倍；365 为奇数，所有保留的正频率都有独立负频率伙伴，因此均翻倍，没有恰好位于 Nyquist 的 FFT 点。逐值独立 FFT 检查与 Parseval 检查包含这一奇数长度端点处理。

## 图形和频率格

与上一张图相同：双对数坐标，横轴为周期（天），从左侧 90 天递减到右侧 2 天，纵轴为 PSD，单位 nT²/Hz。周期按 T_days=1/f_Hz/86400 换算，PSD 数值只随对应频点重排，不转换为每单位周期的谱密度。

显示频率范围保持 [1/(90*86400),1/(2*86400)] Hz。FFT 原生间隔为 1/(365*86400) Hz，范围内实际周期约为 73 天至 2.0055 天，共 178 个频点。绘图连接原生频点，不补造 90 天或 2 天的端点；完整 PSD 表保留 DC 和显示范围外的频率。旧两年图件及结果另行保留。

图内不添加方法脚注、拟合线、峰值标签或显著性判断。PSD 描述该年份整体的波动功率，Hann 权重和已有插值会影响估计。

## 输出

PNG、矢量 PDF、FIG、完整分析 MAT、2015 年日均输入、10 日插值审计、全频率 PSD 表、一个整年窗口的均值与功率归一化审计、原始 CDF 路径及 SHA256、代码哈希、MATLAB 版本和运行日志均保存在结果目录。
