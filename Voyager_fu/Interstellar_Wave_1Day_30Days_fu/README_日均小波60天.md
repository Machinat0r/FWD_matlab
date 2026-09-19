# 日均磁场小波：周期范围扩展为 2–60 天

用户在日均小波色标 0–3 的图上要求将最低频率扩展至 60 天周期。本次重新计算 2–60 天范围；100 个对数频率覆盖整个新范围，色标仍为 log10 功率 0–3。

    addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Wave_1Day_30Days_fu');
    result = Run_V1_Wave_1Day_60Days;

结果目录：C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_60Days
图件：V1_daily_B_Morlet_2day_60days.png / .pdf / .fig

原有 30 天入口与结果保留。共享入口 Run_V1_Wave_1Day_30Days 新增 MaxPeriodDays 参数；默认仍为 30。60 天入口显式传入 MaxPeriodDays=60，并使用独立结果目录。

## 频率与方法

- 最低频率：1/(60*86400) = 1.929012345679e-7 Hz。
- 最高频率：1/(2*86400) = 5.787037037037e-6 Hz，为日采样 Nyquist。
- 纵轴“周期（天）”从顶部 2 天到底部 60 天；色标范围 0–3，功率单位 nT²/Hz。
- 时间范围仍为 [2012-08-25 00:00,2021-12-17 00:00) UTC。
- 默认直接读取同一批原始 CDF，按既有定义计算标量 |B| 日均值，沿用用户授权的 79 日线性插值。原始值、样本计数、插值记录保留。
- 继续直接调用 irf_wavelet Morlet，wavelet_width=5.36、nf=100、cutedge=1；保留原功率归一化，不新增去均值、滤波、背景扣除或功率平滑。
- 原生边界屏蔽随周期变长而扩展。60 天附近频点的前端屏蔽约 60 个日样本，末端再多 1 个；准确计数遵循 IRFU 的 floor(2*a)（保留原始浮点计算及取整），写入验证审计。3401 点输入仍沿用 IRFU 舍弃奇数末点的行为。
- 频率网格覆盖范围已改变，因此旧 2–30 天部分的离散频点也随之变化；全部数值均为新频点上的原生 Morlet 输出。

图内保留正常标签，不增加处理说明。原始 CDF 来源、SHA256、日均及补缺表、方法、代码哈希、MATLAB 版本、功率数组与输出清单保存在结果目录。旧 MAT 仅用于确认日均输入和插值完全一致，不参与科学计算。

## 验证

Validate_V1_Daily_Morlet_60Days 核对新频带与直接 irf_wavelet 输出逐值相同，60 天频率端点及原生边缘屏蔽正确，并用旧范围之外的 45 天正弦信号检查新增频带；同时验证原 30 天默认行为仍可复现。详细的原始数据、插值和电场核查见 README_波动分析.md。

