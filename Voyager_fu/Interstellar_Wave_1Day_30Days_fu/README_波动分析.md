# Voyager 1 日平均磁场 Morlet 小波分析

## 已完成与运行入口

2026-09-10 用户明确要求：仿照 MMS_fu/Wave.m 进行小波分析，直接基于所发图中的日平均磁场。
用户随后明确要求对缺测的 79 天插值后重画。正式入口现对这 79 个内部缺测日做线性插值，然后对完整日均标量 |B| 序列统一调用 irf_wavelet Morlet；图件原位更新。

```matlab
addpath('C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Wave_1Day_30Days_fu');
result = Run_V1_Wave_1Day_30Days('Visible',true);
```

正式代码：C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Wave_1Day_30Days_fu
正式结果：C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Wave_1Day_30Days
新图：V1_daily_B_Morlet_1day_30days.png / .pdf / .fig
分析与审计：V1_daily_morlet_analysis.mat、V1_daily_magnetic_input.csv、V1_daily_interpolation_audit.csv、V1_daily_morlet_segments.csv、source_manifest.csv、run_daily_morlet.log。
验证入口：Validate_V1_Daily_Morlet。

## 日均磁场与原始 CDF

- 时间为 [2012-08-25T00:00:00Z, 2021-12-17T00:00:00Z)，包含图中最后一天。
- 默认每次直接读取 113 个原始 VOYAGER1_COHO1HR_MERGED_MAG_PLASMA 月 CDF，使用 Epoch、标量 ABS_B（nT）。
- 原始目录：Z:/SPART-WORK/Data/Voyager/voyager1/coho/1hr/l2/merged_mag_plasma。
- 沿用所发图 panel a 的日统计：UTC [当日00:00,次日00:00) 内源有效小时 ABS_B 的算术均值，时刻置于12:00。每一天的有效样本数单独保留。
- 共有 3401 个日格，3322 个观测日均值，79 个原始缺测日。依用户最新授权，对全部 79 个缺测日做线性插值；3322 个观测值逐点不变。没有新增样本数或覆盖率阈值。
- CDF 填充值与有效范围处理复用现有 Voyager_Read_CDF_Product 和 IRFU dataobj。
- MAT/CSV 只用于交付、审计、交叉验证，默认科学入口不将它们当作原始输入。
- 本图的小波量为标量日均 |B|；不将它称为三个磁场分量的总波动功率。

## 频率与小波参数

- 1 个月定义为 30 天。
- 原请求频率范围：1/(30*86400) = 3.858024691358e-7 Hz 至 1/86400 = 1.157407407407e-5 Hz。用户随后要求纵轴改用天，顶部设为2天；当前实际显示周期为2–30天。
- 日均采样率 Fs = 1/86400 Hz，Nyquist 上限 Fs/2 = 5.787037037037e-6 Hz。
- 日均数据可分析的周期为 2–30 天；当前纵轴以“周期（天）”标注，顶部2天，底部30天；1–2天不在当前显示范围内。不能通过提高小波频点数或插值获得该部分真实信息。恰好 2 天位于采样极限附近。
- 小波函数直接使用本地 irf_wavelet。计算区间 2–30 天，100 个对数频点，Morlet wavelet_width=5.36，returnpower=1，cutedge=1，均沿用原函数约定。
- 没有 90 天滑窗，也没有 Lomb–Scargle 计算。
- 不做额外去均值、去趋势、预先滤波、背景扣除、沿场旋转、噪声扣除或时间压缩平均。
- 功率采用原 IRFU 定义 P=(2*pi)*|W|²/f，单位 nT²/Hz；图中显示 log10。
- 频率网格构造沿用原函数，完整连续输入的结果已验证与直接调用原函数完全一致。

## 已授权插值与边界

2026-09-10 用户明确要求：“请把缺的这79天插值后重新画这张图”。执行前已说明采用线性插值。

每个缺测日由前后最近的真实日均值连接成直线：B(t)=(1-alpha)*B_left+alpha*B_right，alpha=(t-t_left)/(t_right-t_left)。全部79日都被补齐，不按缺口长度额外筛选。
优先复用已核查的 irf_resamp(...,'linear')，显式方法参数关闭自动平均分支；仅查询缺测日。首尾均有真实观测，并通过程序断言禁止外推。保留已观测日均值不动，不改变它们的原始样本数。

V1_daily_magnetic_input.csv 与 result.Daily 保留：
- BMean_nT：小波与线图实际使用的完整序列。
- BMeanObserved_nT：含原始79个NaN的观测日均序列。
- IsInterpolated：79个插值日的标记。
- MAGSampleCount：原始有效小时样本数，插值日仍为0。

原始三列表另存 result.ObservedDaily。V1_daily_interpolation_audit.csv 与 result.InterpolationAudit 逐日记录左右真实观测日期、观测值、线性权重、插值值和连续缺口长度。
插值值属于估计，未当作新增独立观测。插值会影响缺口附近的波谱，长缺口内的小周期结构不能独立验证。

补齐后3401天构成一个连续输入，统一运行原 irf_wavelet，不再在原79个缺测日处拆段。
原生 cutedge 仅在整个时段两端应用：c=floor(Fs/f)，前c个样本和后c+1个样本置空。
保留原函数奇数输入舍弃末尾一个样本的行为：3401点传入，3400点参与变换，最后一天原始值仍保存在输入和上面板中。
内部缺口引起的大块边界空白消除；完整时段两端保留小波边界空白；纵轴当前裁至2–30天，1–2天频段不再显示。
原图及原分析结果保留在结果目录 before_daily_interpolation/ 中。

不使用姿态近似或额外仪器几何。波谱结构不自动构成传播波、特定模式或统计显著性的判据。

## 电场核查

Voyager 1 有 PWS 电场产品，其谱分析仪 VG1_PWS_LR 的 16 个频率通道覆盖 10–56200 Hz。
已下载并读取一个请求起始日的原始 CDF 样本：
Z:/SPART-WORK/Data/Voyager/voyager1/pws/calibrated/spectrum_analyzer/native/2012/vg1_pws_lr_20120825_v5.30.cdf
该文件有 2099 条记录，electric_field 单位 V/m，frequency 实测元数据范围为 10–56200 Hz。
这是频率核验样本，未声称下载整个 2012–2021 年 PWS 档案。
波形产品 VG1_PWS_WF 的接收机描述通带约 40 Hz–12 kHz，增益信息缺失，不能直接绝对标定。

目标超低频段没有对应电场观测，本任务仅绘制磁场。未把高频电场幅度的慢调制解释为超低频电场，未构造 -v×B。
下载元数据、SHA256 和产品文档随原始样本归档；核查结果见 electric_field_check.json、source_verification.mat。

官方来源（2026-09-10 核查）：
- [CDAWeb Voyager 产品说明](https://cdaweb.gsfc.nasa.gov/misc/NotesV.html)
- [NASA PDS PWS 产品说明](https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=VG1-J/S/SS-PWS-4-SUMM-SA1HOUR-V1.0)
- [NASA SPASE 电场波形说明](https://heliophysicsdata.gsfc.nasa.gov/WS/hdp/1/Spase?ResourceID=spase://NASA/NumericalData/Voyager1/PWS/Waveform-CDF/PT0.00003472S)

## 验证与显示

Validate_V1_Daily_Morlet 已通过：
1. 完整连续序列的小波数组与直接 irf_wavelet 输出完全一致。
2. 5 天合成信号的小波谱峰恢复为约 5.0692 天，符合原函数的离散频点与功率标度。
3. 未插值输入的底层小波函数仍正确保留缺口边界；授权插值后形成一个连续段，内部缺口处系数全部有限。
4. 线性插值与独立端点公式一致，真实观测值完全不变；整段小波与直接调用原函数完全一致。
5. 不计算超过 Nyquist 的频率，过短连续段保持空白。

图内只保留磁场、频率轴、单位、panel 编号与正常标题。频谱采用原始日格平色显示，不平滑或插值。
按用户最新要求，log10功率色标固定为[0,3]。色标之外的值仅在颜色上饱和，完整数值保留在MAT中。周期坐标使用T=1/f/86400换算；功率数值保持原样，色标单位仍为nT²/Hz。实际纵轴方向、周期范围与色标范围写入PlotAudit。


实测观测日均输入交叉验证：ObservedDaily 与原图的3401天观测值逐点完全一致（含NaN）。插值后仅79个原缺测值改变，观测值和样本数保持原样；最终审计见 daily_input_crosscheck.mat。旧滑动谱草稿仅留存在 superseded_lomb_draft 子目录，正式入口不调用。



