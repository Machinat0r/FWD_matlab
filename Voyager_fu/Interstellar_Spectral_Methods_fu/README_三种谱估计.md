# Voyager 1：2014–2015 年三种方法的功率谱

## 运行
MATLAB 入口：

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Spectral_Methods_fu');
    result = Run_V1_PSD_ThreeMethods_2014_2015('Visible',false);

需要 MATLAB Signal Processing Toolbox 和本机 irfu-matlab-master。入口每次从原始月 CDF 重读数据；旧 MAT 仅用于核对是否与原图输入一致。

结果目录：
C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_PSD_ThreeMethods_2014_2015

## 数据与已授权处理
- 时间：2014-01-01 00:00 UTC 至 2016-01-01 00:00 UTC（右端不含），对应最新附件图的 2014–2015 两年。
- 输入：Z:/SPART-WORK/Data/Voyager/voyager1/coho/1hr/l2/merged_mag_plasma 下的 24 个官方月 CDF，数据集 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA，使用 Epoch、ABS_B，磁场单位 nT。
- 沿用同一套逐 UTC 日有限小时值算术平均：730 天，712 天有观测，18 天按此前明确授权作线性插值；不外推，不增加样本数阈值。
- 三种方法输入完全相同，并已与附件图所用日均序列逐值核对。
- 按本次提供的三个脚本直接使用 |B|，不额外去均值或去趋势。此前图采用去均值、周期 Hann 窗；因此本次 FFT 与旧图也会有变化。
- 横轴是周期（天），从左到右 90→2；换算 T_days=1/(f_Hz*86400)。纵轴仍是按频率定义的 PSD，单位 nT²/Hz，未转换成每单位周期的谱密度。

## 对所附三个脚本的适配
### FFT
参考 turbulence_power_law_fft.m 的 Hamming 加窗 FFT 核心。原脚本用于 Cluster 高频数据，给定一个事件时刻，从分窗频谱中抽取该时刻的曲线。本次请求是两年整体功率谱，因此用一段完整的 730 点 Hamming 窗。沿用 MATLAB hamming(N) 的对称窗，不补零。

本机 irf_wavefft 把 frame_length 和 frame_overlap 按毫秒换算成样本，原附件的 2048 数字不能作为日采样数据的固定点数直接传入；该旧函数循环还会跳过与输入等长的唯一窗口，并舍去偶数长度的 Nyquist 点。因此本适配直接采用相同的 fft(x.*hamming(N)) 核心，不修改 IRFU 原函数。

附件直接画 abs(FFT)^2，缺少每 Hz 的归一化。为保持所需 PSD 单位，使用：

    Y = fft(x.*w,N);
    PSD_two_sided = abs(Y).^2/(Fs*sum(w.^2));

保留非负频率，仅把正频率中非 Nyquist 点乘 2。DC 与 Nyquist 不倍增。FFT 全谱和 DC 都保留在数值输出中，图中显示 2–90 天范围。

### 小波时间平均
参考 turbulence_power_law_wavlet.m 的 irf_wavelet 和 nanmean(pB,1)。
调用本机 Morlet 小波，沿用默认宽度 5.36、200 个频点、cutedge=1，并把频率范围明确设为 2–90 天。每个频率对两年范围内有限的小波功率取算术平均，保持 IRFU 的 2*pi*abs(W)^2/f 归一化。
原函数边缘剔除照常保留，因此各频点实际参与平均的样本数为 551–725；逐频率计数及首尾时间见 wavelet_averaging_audit.csv。两年共有偶数个样本，不发生原函数的奇数末样本舍弃。此处只有磁场，无需重采样到电场时间轴。
小波曲线与 FFT、PMTM 的时间权重和频率分辨率不同，不强行缩放使曲线重合。

### PMTM
保留附件原调用：

    Fs = 1/86400;
    [P,f] = pmtm(x,3.5,[],Fs);

PMTM 是 Thomson 多窗功率谱估计。它给同一整段信号分别乘多个正交的 DPSS（Slepian）窗，分别求傅里叶功率，再用随频率变化的自适应权重合成估计。这里的 3.5 是无量纲的时间–半带宽积 NW；默认实际采用 2*NW-1=6 个窗。六个窗都覆盖同一段 730 天数据。

本数据长度下默认 NFFT=1024；补零只增加频率网格的密度。半带宽 W=NW/(730*86400)=5.5492136e-8 Hz，名义全带宽约为 1.1098427e-7 Hz。NW 增大会增加所用窗数并降低估计方差，同时使峰变宽，邻近频率成分更难区分。因此更平滑的谱线不自动意味着每个峰都更可信。

官方说明：
- https://www.mathworks.com/help/signal/ref/pmtm.html
- https://www.mathworks.com/help/signal/ug/nonparametric-methods.html

## 检查和输出
- FFT 与 MATLAB periodogram 计算一致，检查单边归一化和 Parseval 功率积分。
- PMTM 默认调用与显式传入六个 DPSS 窗、自适应权重的结果完全一致。
- 小波平均值及原函数剔除边缘后的逐频样本数已核对。
- 用已知 20 天正弦验证三种方法的频率/周期换算，峰的位置均在各自分辨率允许范围内。
- 三张图独立输出 PNG、矢量 PDF、可编辑 FIG；坐标范围一致，图内只保留常规标题和坐标标签。
- 每种方法各自输出原生频点的 CSV，不插值成共同频率网格。
- MAT 保存原始读数、缺测前后日均值、插值审计、完整小波功率、DPSS 窗、谱值、程序和附件哈希、CDF 来源及验证结果。
- 原始数据和以往图件保持原路径。没有进行功率律拟合或峰显著性推断。

上述日均、插值、保留均值及各方法的窗口/边缘权重都会影响谱值。不能仅凭曲线平滑程度判断哪个结果更接近真实谱。
