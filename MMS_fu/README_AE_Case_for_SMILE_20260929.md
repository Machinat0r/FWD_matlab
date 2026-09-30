# Case_for_SMILE 的 AE 指数图（2026-09-29）

程序：AE_Case_for_SMILE_20260929.m，直接MATLAB脚本，无新增function。

事件范围来自 Case_for_SMILE.pptx 的“事件时间”文字（第2、3、4、5、8页），全部UTC：
1. 2026-07-21 00:00–04:00
2. 2026-07-24 08:00–2026-07-25 05:00
3. 2026-07-28 05:10–06:00
4. 2026-08-04 01:50–03:00
5. 2026-08-11 06:50–08:30

数据：NASA CDAWeb OMNI_HRO_1MIN 的 AE_INDEX，nT，1分钟。该时段的AE由WDC Kyoto提供，为quicklook版本。来源说明：https://omniweb.gsfc.nasa.gov/html/omni_min_data.html 。

原CDF归档到 Z:\SPART-WORK\Data\MMS\ancillary\omni\hro_1min\2026，保持官方文件名和v01版本。使用IRFU dataobj、get_ts、irf.ts2mat、irf_tlim、irf_plot、irf_zoom及irf_subplot。CDF的AE为整数变量，因此根据原始FILLVAL=99999掩码在double绘图数组中明确保留NaN，并逐点核对原变量。没有自行构造AE指数，没有平滑、插值或重新平均。

时间选择沿用irf_tlim的[start,end)，保留每个原生分钟样本；图中坐标轴范围完全按PPT给定。5个事件分别有240、1260、50、70、100个有效分钟样本，所选时段没有缺测。8月CDF最后一天的缺测不在本次事件内。

运行：在MATLAB中将MMS_fu加入路径后直接执行 AE_Case_for_SMILE_20260929。程序复用已有CDF；文件不存在时才调用MATLAB websave下载。

输出目录：C:\Users\Administrator\Documents\KH\Case_SMILE_AE_20260929，含5张PNG和5页矢量PDF。原PPT未修改。每次运行覆盖这批AE图。运行记录位于Z盘derived\Case_SMILE_AE_20260929。按用户要求仅本地保存，不发送邮件。