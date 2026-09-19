# Voyager 1 日均磁场傅里叶时频图

## 运行入口

正式程序目录：
C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Fourier_1Day_30Days_fu

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Fourier_1Day_30Days_fu');
    result = Run_V1_Daily_Fourier('WindowDays',90);

2026-09-11 用户在明确收到 90 天窗口、1 天步长、Hann 窗、逐窗去均值的方案后要求“画图”；本次按该方案执行。默认配置为 90 天窗口、每天移动 1 天、周期 Hann 窗、每窗减去自身算术均值。WindowDays 可配置。合成信号测试入口为 Validate_V1_Daily_Fourier。

默认输出目录：
C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Fourier_1Day_30Days

## 数据通路

时间范围为 2012-08-25 00:00 至 2021-12-17 00:00 UTC，后者不含。默认直接读取 Z:/SPART-WORK/Data/Voyager/voyager1/coho/1hr/l2/merged_mag_plasma 下的原始 CDF，使用 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA 的 Epoch 和标量 ABS_B（nT）。

复用日均小波程序中的 V1_Wave_Read_Hourly 和 V1_Interpolate_Daily_Gaps，以及现有 Voyager_Read_CDF_Product、irf_resamp。先按 UTC 日对源有效小时 |B| 求算术平均，再沿用用户此前授权的 79 个内部缺测日线性插值。日均值、插值值、原始观测计数及缺测审计均与日均小波版逐项核验。参考 MAT 仅用于核验，不作为科学输入。首尾禁止外推，不新增有效样本数阈值。

磁场曲线始终显示原有补齐后的日均值。分窗去均值仅作用于傅里叶计算副本，不修改磁场曲线或交付输入表。分析量为标量 |B|，无磁场矢量功率合并。

## 短时傅里叶方法

保留时间横轴的二维图需要分窗计算。每个完整窗口先减去该窗口的算术平均值，再乘周期 Hann 窗：

    w(n)=0.5-0.5*cos(2*pi*n/L), n=0,...,L-1

窗口每天移动 1 天；FFT 点数等于窗口实际日数，不做零填充。复用 irf_psd，每次输入恰好一个完整窗口，传入显式 Hann 向量、noverlap=0、dflag='mean'。外层循环移动窗口，各窗口谱保持独立输出，无额外跨窗平均。程序不做线性去趋势、频谱平滑或频谱插值。

单边 PSD 归一化：

    P2=abs(FFT((x-mean(x)).*w)).^2/(Fs*sum(w.^2))

严格正频率的非 Nyquist 点乘以 2。IRFU 当前 irf_psd 对全部点乘 2，因此本适配将 DC 和偶数 FFT 的 Nyquist 点除以 2，正确保留单边端点。单位为 nT²/Hz。频域积分逐窗核对时域 Hann 加权均方量：

    sum(P1)*(Fs/L) = sum(((x-mean(x)).*w).^2)/sum(w.^2)

irf_powerfft 和 irf_wavefft 已检索：前者含固定线性去趋势及额外数据整格行为；后者使用毫秒参数并遗漏 Nyquist 端点。当前选择显式分窗复用 irf_psd，以准确控制本次参数。

## 周期、窗口及边界

日采样 Fs=1/86400 Hz，Nyquist=1/(2*86400)=5.787037037037e-6 Hz。一个月按 30 天，下限为 1/(30*86400)=3.858024691358e-7 Hz；图中只显示 2–30 天。DC 和带外频点保留在数值审计中。

90 天窗口的 FFT 频率间隔为 1/(90*86400)=1.286008230453e-7 Hz，对应窗口内能覆盖三个 30 天周期。Hann 窗会使实际谱峰宽于一个 FFT 格点。60 天窗口提供较细的时间定位和较粗的频率区分，120 天反之。每天移动一次使输出更密，各谱点仍依赖整个窗口，重叠窗口也并非独立样本。

只计算完整窗口，不使用时段外数据，不以零值、镜像或外推补足首尾。时间戳为完整窗口实际时间区间的中点；图两端会保留约半个窗口的不可计算范围。

## 绘图和结果解释

沿用日均版两面板样式：上方日均 |B|，下方时间—周期傅里叶功率色图。周期轴单位为天，2 天在上、30 天在下；色标 log10 PSD 范围 [-2,3.5]。时间色块宽度是 1 天输出步长。FFT 原生频率格使用频率相邻中心的算术中点划分显示边界，不人为增加频率点。图内仅有正常坐标、单位、标题和 panel 标签。

傅里叶 PSD 和原 Morlet 图采用不同估计方法及归一化约定；色标单位相同，不代表所有数值可直接等同。线性补缺会影响局部谱功率，插值点不增加独立观测。窗口选择影响频率和时间分辨率，频谱本身不证明特定波模或显著性。

沿用此前 PWS 核验结论：所需低频段没有可用的电场观测，本次输出磁场图。

## 审计和验证

输出 PNG/PDF/FIG、完整 MAT、日均输入、79 日插值审计、逐窗起止时间/均值/原始与插值日数/归一化检查、实际频率格、CDF 清单及 SHA256、程序 SHA256、MATLAB 版本和日志。全部结果均保存到结果目录，原始数据归档不写入图或程序副本。

自动核验包括：与旧日均输入完全一致；所有窗口与独立 MATLAB fft 公式一致；单边谱 Parseval 等式；窗口时间中点和 1 天步长。合成信号检查 10 天正弦的周期和振幅功率、Nyquist 信号的端点权重、常量背景去除。该测试不使用 Voyager 数据。

