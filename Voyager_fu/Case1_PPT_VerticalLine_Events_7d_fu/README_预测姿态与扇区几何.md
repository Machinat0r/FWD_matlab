# Case 1：官方预测姿态与 LECP P1 扇区几何

## 授权与适用范围

2026-09-03，用户批准官方预测姿态、约 7 天姿态节点到每日的插值，以及官方标称 LECP 安装轴和扇区中心；明确要求不在图注中标注这些方法。本文件和逐图 MAT 审计承担溯源记录。所有推导仅用于 Voyager 1 当前 P1 低能孔径，不能直接推广到全部 FHDU 能道。用户最新要求新图 panel e/f 的对应能段按 Mosley 参考上限显示为 0.57–1.78 MeV。源 CDF `FHDU_EnergyRange` 中 P1 的 0.57–0.89 MeV 字段保持原样，与显示能段分开记录；当前核查的 V1/V2 官方相关元数据取值相同，差异原因尚未解释，不能据此判定真实仪器能段。标签变更不改变 P1 选取、原通量/不确定度、元数据或几何，也不涉及能宽缩放或重标定。

同日用户进一步明确：“如果出现负 deltaT，可以把该数据舍去。其余的直接按照 Epoch 时刻锚定即可。”当前据此舍去 `DeltaT < 0` 的记录，其他官方日/小时产品记录直接采用原始 Epoch 和原 Flux/Sigma；预测姿态在相同 Epoch 求值。磁场沿用该 Epoch 所属 UTC 日/小时内有效三分量记录的均值，不插值磁场或通量。此前累计区间待确认的统一暂停条件已解除，不再采用固定日中点/小时中点求值，也不因跨界或正时长过大而剔除或拆分记录。该规则不改变既有仪器几何和姿态来源限制。

当前原始 CDF 的 FHDU_Energy.CATDESC 按顺序列出 PL01–PL08、P32、P1、P10、P11、P16、P23、P27、P31；P1 是第 10 通道。程序保留最接近 0.73 MeV 的选取方式，并拒绝在新几何模式下把其他索引误当 P1。[Mosley et al. 2006 的 Table 2](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2006GL026732)用于核验 P1 的 D1/D2 符合、D3 反符合硬件逻辑，用户同时指定新图采用其参考上限作能段显示。源 CDF 数值和元数据另行保留，不由显示能段推算通量。孔径朝向另见 [Krimigis et al. 1977 的仪器说明](https://pds.nasa.gov/data/vg1-j-mag-4-rdr-hgcoords-9.60sec-v1.0/vg_1501/document/lecp/lecpinst.htm)。

## 姿态、时间与 RTN

- 原始文件与精确 URL 由 `Case1_Attitude_Files.m` 列出，下载和读取时保存 SHA256。CK 为 `vgr1_super.bc`，SCLK 为 `vg100051.tsc`，FK 为 `vg1_v02.tf`，LSK 为 `naif0012.tls`，SPK 为 `vgr1.x2100.bsp`，PCK 为 `pck00011.tpc`。
- [官方 CK 注释](https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/ck/vgr1_super_v2.cmt)明确 1990-04-16 起为预测段，节点平均间隔 604800 秒。此注释保留原始 CK 的说明；程序读取无附加角速度的原始 `vgr1_super.bc`。
- 在每条保留的官方 LECP 日/小时产品原始 Epoch 直接调用 `cspice_ckgp(-31000,clock,0,'J2000')`。官方 type-3 在相邻节点间按固定相对旋转轴和角速度插值；没有逐分量或欧拉角 `interp1`，没有把日中点向量再次插值到小时。`ckcov` 和零容差共同约束时间覆盖。
- `AttitudeDailyHourUTC=12` 仅为旧接口兼容字段，当前 Epoch 模式不使用；`PADCadence='day'/'hour'` 选择官方通量产品和磁场匹配时间格，不改变原始 Epoch。
- 返回的 C 是 J2000 到 SC_BUS，因此本体到 J2000 用 C 的转置。此方向按 [MICE ckgp 文档](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/mice/cspice_ckgp.html)实现，并与独立 FK/pxform 路径数值对照。
- 由同时刻几何星历（`NONE`，无光行时或像差修正）构造 R=太阳到 V1 的单位矢量、T=unit(太阳北自转轴×R)、N=R×T。太阳北极取官方 PCK 的 IAU_SUN，未以黄道北极替代。[COHOWeb RTN 定义](https://omniweb.gsfc.nasa.gov/coho/html/cw_data.html)

## 扇区中心与方向正负号

[SEDR Table G-1](https://voyager.ftecs.com/Handbook/DataFileDesc/SEDR/tabg1.html)给出 LECP 旋转轴本体 clock=200°、cone=90°。[官方 FK](https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/fk/vg1_v02.tf)将 VG1_CONE_CLOCK 转到 VG1_SC_BUS，得到：

```matlab
nSC = [-sind(15); cosd(15); 0];
hSC = [0; 0; -1];              % HGA 向外望向
qSC = cross(nSC, hSC);
beta = 22.5 + (0:7)*45;        % S1 ... S8，相对 HGA
lookSC = hSC*cosd(beta) + qSC*sind(beta);
particleSC = -lookSC;          % 粒子运动方向指向孔径内部
```

这里的本体轴零分量属于仪器安装坐标表达，没有把 RTN 的 uN 设为零。

| 扇区 | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 相对 HGA 的中心角 | 22.5° | 67.5° | 112.5° | 157.5° | 202.5° | 247.5° | 292.5° | 337.5° |

[Voyager 1 PDS 仪器目录](https://pds-ppi.igpp.ucla.edu/data/VG1-S-LECP-4-SUMM-SECTOR-15MIN-V1.0/CATALOG/VG1_LECP_INST.CAT)规定 S1/S8 分界与 HGA/roll 轴平行，S1、S8朝地球；Canopus 锁星时 S3 望向黄道北侧约 30°并逆轨道方向。用历史 CK 的 1980-01-01 和 1980-10-01 12:00 UTC 进行独立检验，正号方案符合两项条件，反号方案同时翻转南北和顺逆行。这是几何编号顺序，不代表电机始终按该方向扫描。NAVMAG 网页中“Y 沿 sector Z”的可疑转录没有被用作唯一符号依据。

计算次序为：原始 Epoch 的本体粒子方向 → J2000 → 与 COHO 一致的 RTN → 与该 Epoch 所属 UTC 日/小时的磁场均值方向做点积。PA=acosd(u_particle·Bhat)。S8 的几何可以计算，但其通量继续只作遮挡/校准相关诊断，不加入 S1–S7 正常 PAD。

## 仍需保留的解释限制

1. 姿态为预测值，在日产品或小时产品 Epoch 求值均不会恢复实际滚转、姿态机动或新的独立观测信息。
2. 使用官方标称安装中心线；未取得各扇区 as-built 安装偏差校准。有限视场通量被放在中心线对应的 PA 上，未进行响应积分或反卷积。
3. 使用原始 Epoch 的指向和其所属 UTC 日/小时的磁场均值，将官方平均通量放在该角度上；所得角度不等于逐曝光 PA 的时间平均，也未与逐扇区实际曝光时刻加权配对。较长或跨界的正 DeltaT 按用户要求保留，无法保证通量、姿态和磁场具有完全相同的时间支持区间。
4. 先按原始 Epoch 选择事件窗口并剔除负 DeltaT，随后只有 S1–S7 均有有效正通量和有限 PA 的记录进入 PAD。零或缺失 DeltaT 不新增排除条件；用户删去的样本覆盖、PA 跨度、方向 RMS、不确定度门限没有恢复。磁场矢量缺测或模为零时数学上无法求夹角，PA 保留 NaN。
5. 无背景扣除、无新磁场/通量插值；缺测保持缺测。概览时间色图保留已有 2° 中心合并显示规则，原扇区结果在 MAT 审计中不合并。新增五时刻黑色点图按用户要求完全保留七个扇区角度，不做此显示合并。
6. 新五时刻图按每条记录的七扇区最大原通量独立显示 `J/Jmax` 和 `sigma/Jmax`；分母作为固定值，不传播分母误差或协方差。窗口峰值按原通量搜索，峰值缺少 PA 时中格留空，不换用较弱峰值；两侧允许跳过缺测、各取最近两条可用 PAD，仍限制在原事件窗口内。该显示与选时授权不改变姿态或几何规则，详见数据处理说明第 5.5 节。

## 验证和复现

运行 `Case1_Test_Predicted_Attitude` 检查 2013-01-02 至 2021-05-25 的 3066 个日点、旋转矩阵正交/右手、FK 对照、RTN 径向、夹角不变、扇区间距与历史手性，以及缺口和调用者内核池保护。数值误差容限只是程序验证条件，不参与科学样本筛选。

逐图审计保存源 CDF 路径和原始记录行号、原始 Epoch/DeltaT、负 DeltaT 剔除明细、匹配磁场时间格及均值、各扇区原通量/不确定度/PA/完整 RTN 方向，以及原始 CK 等内核的 SHA256。新输出名包含 `predictedCK_1d` 或 `predictedCK_1h` 与 `nativeCDF_Epoch` 标识，现有历史图默认保留。图内不增加预测姿态或标称安装几何说明。

2026-09-03 早先日中点方案的历史验证：3066/3066 日、23/23 项数值回归通过；旋转正交误差最大 1.57e-15，独立 FK 对照误差为 0。测试审计文件为 `Z:\SPART-WORK\Data\Voyager\voyager1\attitude\spice\validation\predicted_attitude_regression_20260902_170107_834.mat`。当时 Case1-S01-L01 单图端到端验证在 15 个日格中有 14 个有效 PAD，完整 49 事件预检通过。后续 Epoch 日/小时概览图完成记录见项目上下文；本轮能段显示与五时刻点图仍需独立验证，实际输出以本轮运行清单为准。程序数值验证不构成预测姿态精度或安装误差的观测验证。Florinski 2008 原图视觉格式尚未完全核实，本轮五时刻图按用户明确的新 PA 横轴版式执行，不声称逐项复现该原图。
