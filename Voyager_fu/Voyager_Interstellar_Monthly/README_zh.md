# 旅行者号星际空间逐月绘图程序

本程序把“太阳系以外”按**穿越太阳风层顶后进入局地星际介质（VLISM）**来处理：

- Voyager 1：从 2012-08-25 开始
- Voyager 2：从 2018-11-05 开始

程序读取 `Z:\SPART-WORK\Data\Voyager`，按航天器、年份和 UTC 自然月输出一张综合图。当前本地数据可生成 Voyager 1 的 2012-08 至 2025-12，以及 Voyager 2 的 2018-11 至 2022-12，共 211 个自然月。已有测量中的时间缺口会按数据等级补齐；整个变量完全没有测量基础时，面板会明确标为无测量。

## 快速运行

在 MATLAB 中打开并运行：

```matlab
run('C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Voyager_Interstellar_Monthly\run_voyager_interstellar_monthly.m')
```

默认输出到本程序目录下的 `monthly_plots`。首次运行需要把网络盘 CDF 文件复制到本地缓存并转换，后续运行会复用缓存。211 张图可能需要几十分钟，具体时间取决于 Z 盘速度和导出分辨率。

`examples` 文件夹中附有两张用 Z 盘真实数据生成的验收图：Voyager 1 的 2012-11 和 Voyager 2 的 2018-11。

## 每张月图的内容

1. 磁场强度 `|B|`
2. RTN 坐标系磁场分量 `B_R/B_T/B_N`
3. 等离子体体速度
4. 质子数密度
5. 质子温度
6. LECP 质子差分通量全部通道
7. CRS 质子差分通量全部能道，以时间—能量谱显示
8. 日心距离及月内代表性的 HGI 经纬度

太阳风层顶穿越所在月会画紫色虚线。程序同时输出逐月覆盖表，其中列出各产品文件、原始有效记录数、低精度补入量、插值量、状态和错误信息。

## 数据优先级与时间分辨率

磁场主图使用经过审校的 VIM 48 秒产品。它在星际阶段具有更可靠的定量标定，CDF 给出的典型不确定度约为 0.02–0.03 nT（`B_N` 通常为 0.1 nT）。本地还保存了 1.92 秒实验性、未审校磁场产品；其时间采样更密，产品说明提示它通常不适合作为科学质量数据，因此没有混入主图。

连续化按以下顺序执行：

1. 保留审校 VIM 48 秒磁场样本。
2. 用 COHO 1 小时磁场插值到 48 秒网格，填入审校数据缺口。
3. 两级产品仍有空点时，对月内已有值做线性插值，并用最近有效端点补齐月初/月末。
4. LECP/CRS 在对数通量空间按 1 小时网格插值，位置量按日网格插值。

覆盖表中的 `MAGLowPrecisionFilled`、`MAGInterpolated` 和 `COHOInterpolated` 会记录每个月各类补全的数值个数。

月图抽稀采用分时间箱保留极大值和极小值的方法。短时磁场尖峰、突变和可能的激波结构会被保留，常规等间隔抽样可能漏掉的峰值也能显示出来。

等离子体优先使用 PLS Level 3 高时间分辨率拟合矩，缺少时读取 COHO 1 小时合并产品。Voyager 1 穿越太阳风层顶后缺少直接 PLS 等离子体测量；Voyager 2 的高分辨率日球鞘 PLS 文件在 2018-11-05 附近结束。因此，许多星际月份的速度、密度和温度面板为空，这是仪器/产品覆盖造成的真实缺测。

## 关于激波判读

这些月图适合筛选候选事件：关注 `|B|` 和 RTN 分量的同步突变，并检查 LECP/CRS 粒子响应；有 PLS 覆盖时再核对速度、密度和温度。程序不自动给事件贴“激波”标签，因为可靠确认通常需要跨仪器证据、上下游平均和坐标分析。发现候选日期后，建议再做小时或分钟尺度的专题图。

## 缺数下载

将运行脚本中的：

```matlab
voyagerDownloadMissing = true;
```

程序会调用：

```text
C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Voyager_Download.m
```

下载范围会分别从两艘航天器的太阳风层顶穿越日期开始，并受 `StartTime/EndTime` 限制。下载器会检查已有文件大小并跳过完整文件。

## 常用配置

可以直接调用主函数，例如只画 Voyager 1 的 2012 年 11 月：

```matlab
report = Voyager_Plot_Interstellar_Monthly( ...
    'DataRoot', 'Z:\SPART-WORK\Data\Voyager', ...
    'OutputRoot', 'D:\Voyager_monthly', ...
    'Spacecraft', 1, ...
    'StartTime', '2012-11-01', ...
    'EndTime', '2012-11-30T23:59:59', ...
    'FillGaps', true, ...
    'Overwrite', true, ...
    'SaveFormats', {'png','pdf'});
```

主要参数：

- `Spacecraft`：`1`、`2` 或 `[1 2]`
- `StartTime`：空值表示从各自太阳风层顶穿越日期开始
- `EndTime`：`'auto'` 表示扫描本地最新数据
- `DownloadMissing`：是否调用现有下载器补齐缺失产品
- `FillGaps`：是否启用分级连续化，默认 `true`
- `MagCadenceSeconds`：连续磁场网格的时间步长，默认 48 秒
- `Overwrite`：是否覆盖已存在的月图
- `MakeEmptyFigures`：是否为完全缺数的月份保留空图
- `SaveFormats`：支持 PNG、PDF、JPG 和 TIFF
- `MaxPlotPoints`：每条曲线绘图点预算；抽稀会保留每箱极值

## 输出结构

```text
monthly_plots/
  Voyager1/YYYY/Voyager1_VLISM_YYYY-MM.png
  Voyager2/YYYY/Voyager2_VLISM_YYYY-MM.png
  Voyager_VLISM_monthly_coverage.csv
  Voyager_VLISM_monthly_report.mat
```

`Status` 可能为 `ok`、`empty`、`partial_error`、`error` 或 `plot_error`。单一产品读取失败时，其余可用产品仍会继续绘制，并在 `Notes` 中记录原因。

## CDF 兼容处理

这台机器上的 MATLAB R2026a 直接调用 `cdfinfo` 读取部分 NASA CDF 会触发访问冲突。程序附带 `voyager_cdf_bridge.py` 和 `cdflib 1.3.6`，先把所需变量转换成本地二进制缓存，再由 MATLAB 安全读取。无需修改原始 CDF，也不会向 Z 盘写缓存。

## 两组无插值输出

运行：

```matlab
run('C:\Users\Administrator\Documents\FWD_matlab\Voyager_fu\Voyager_Interstellar_Monthly\run_voyager_two_groups.m')
```

程序会在同一目录下生成两个扁平文件夹，文件夹内部不再设置年份或航天器子文件夹：

- `HighestPrecision_Segments`：仅使用审校 VIM MAG 48 秒记录。用 2013-01 图中的 59 个可见数据块进行标定，相邻记录间隔超过 7200 秒时开始新段。段内超过 240 秒的缺口会断线显示，且不会使用低精度数据连接。
- `UniformLowPrecision_Monthly_Raw`：统一使用 COHO 1 小时产品，逐月输出。程序直接连接实际存在的仪器记录点，不调用 `fillmissing`、`interp1`、`retime` 或重采样，也不会生成替代数值。

两个文件夹都包含 CSV 清单和 MAT 报告，可核对源 CDF、起止时间、记录数、通道数与状态。
