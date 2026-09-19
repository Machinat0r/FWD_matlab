# 三日与月平均图新增太阳黑子数

运行 Run_V1_Averaged_Sunspots(false)。正式入口直接复用原始 Voyager CDF 统计通路，另读 Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv。太阳数据来源 WDC-SILSO, Royal Observatory of Belgium，Daily total sunspot number V2.0，https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv 。CSV 是该机构原始发布格式，非派生科学输入。

时间为 [2012-08-25,2021-12-17) UTC。三日窗口从起点起每三日分组且不重叠，末窗仅两天。自然月窗口首窗从8月25日起、末窗截至12月16日。直接使用原 Voyager 窗口起止，太阳黑子数在相同窗口对有效日值等权算术平均；-1 作为缺测，零值保留，空窗口为NaN。无平滑、插值、传播时间平移或最低覆盖率阈值。均值与独立 accumarray 分组核验，保存源行、源列、源文件校验值、分箱归属及有效天数。

新增太阳黑子数为a；原四面板顺延b–e。三日图d继续是原始小时P1样本的三日中位数，月图d继续是日中位数的月算术平均。e继续排除S4和S8，沿用原六扇区日均和再求窗口平均。太阳黑子数与磁场为线性轴，其余对数轴，细线0.4，缺测断线。图内不新增方法说明。

输出位于 Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/averaged_no_S4_sunspots，包含两套PNG/PDF/FIG、窗口值CSV和averaged_sunspots_audit.mat。源CDF审计与原磁场相关系数位于其source_recompute目录（配对编号仍对应原四面板）。原有四面板图未覆盖。日志位于上一级 averaged_sunspots_run_20260909.log。

复用已审查的显式UTC分组实现，IRFU irf_waverage 的权重/NaN规则与本任务不匹配。Run_V1_ThreeDay_Monthly_Overview 新增可选输出目录与makePlot参数，其默认调用行为保持原样。
