# 20个事件的MMS四星位置图（2026-09-24）

入口：MMS_orbit_events_20260924.m

## 运行

在MATLAB中：
    addpath('C:\Users\Administrator\Documents\FWD_matlab\MMS_fu')
    MMS_orbit_events_20260924
    MMS_orbit_events_20260924([1 19],true)   % 指定事件并强制重画
    MMS_orbit_events_20260924(1:20,true)     % 全部重画

默认检测到已完成的PNG和记录后跳过。需要重画时第二参数设为true。
程序依赖本机 irfu-matlab-master 与已归档的数据；不需要转换整个CDF。
运行时原绘图函数会联网查询OMNI参数。无法取得参数时，沿用其默认边界模型。

## 复用的原程序

用户原稿：
C:\Users\Administrator\Documents\FWD_matlab\新建文件夹\MMS_orbit.m
原稿调用 mms.mms4_pl_conf(tint) 和 mms.mms4_pl_conf('gsm')。

本包中的 MMS_orbit_pl_conf_20260924.m 是本机IRFU原函数的兼容副本，
仅修改函数/回调名称，以及触发当前时间格式错误的毫秒字符串。
MMS_orbit_layout_20260924.m 将原8个坐标轴排为两行，放大字体，换行图例，
并将绝对位置X显示范围扩展到[-35,20] RE，便于PPT阅读。
原有科学算法、坐标轴物理量、卫星符号、边界模型保留。
用户原稿及IRFU原文件均未覆盖。

## 位置数据与时间

全部时间为UTC。单一时间直接使用，事件时间段取中点。
全部位置图使用GSM，RE=6372 km沿用原函数。
上排为地心位置投影，下排为相对四星几何中心的位置投影及原图例。
每个事件的图均包括MMS1–4。

EV01–EV15：MMS MEC L2 epht89d原始CDF。
EV16–EV20：NASA SSCWeb返回的原始GSE星历，原生60 s采样；
在有原生样本包围的区间内用irf_resamp线性对齐，然后由原IRFU绘图函数转换至GSM。
已下载的DEFEPH文件未用于最终图形，避免使用本机旧脚本中未完整处理的J2000转换。

数据全部保存在 Z:\SPART-WORK\Data\MMS。
MEC：mms1–mms4\mec\srvy\l2\epht89d\2026\07 或 08。
SSC：ancillary\sscweb\2026\EV15_gse.json 至 EV20_gse.json。
清单与位置审计：derived\event_orbits_20260924。
程序包manifests目录包含清单副本；本机运行读取上述Z盘derived目录中的清单。

MEC来源：https://lasp.colorado.edu/mms/sdc/public/about/how-to/
SSC来源：https://sscweb.gsfc.nasa.gov/WebServices/REST/
每张PPT的备注中记录原始数据路径、实际坐标、时间与来源URL。

## 边界曲线

原程序的磁层顶/弓激波曲线为模型示意。
本次EV06、EV07、EV08、EV09、EV11、EV12、EV19、EV20使用原函数的默认
P=2 nPa、Bz=0 nT，均已在图内标明“IMF using assumed model”。
其余事件使用函数成功取得的OMNI参数。
上述模型参数只影响边界曲线，卫星位置均来自星历。

## 输出

图形与可编辑MATLAB FIG：
C:\Users\Administrator\Documents\KH\MMS_event_overviews_PPT_20260924\orbits

更新后的PPT共127页：原107页与新增20页位置图，每页插在相应事件overview前。
原107页在PowerPoint中渲染后的像素与上一版完全一致。
图形以PNG嵌入PPT；修改科学图请编辑MATLAB程序或FIG后重绘。
