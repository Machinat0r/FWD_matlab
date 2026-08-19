# MMS KH FOTE / FOTE-V（MATLAB/IRFU模块化版本）

本程序在MATLAB中完成MMS四星数据读取、物理时间平均、FOTE/FOTE-V、质量筛选、绘图及结果输出。内部函数使用`+khfote` package管理，避免函数名冲突，并方便分别测试和修改。

## 程序结构

```text
kh_fote_matlab/
├─ kh_fote_event.m                 # 对外兼容入口
├─ kh_fote_event_irfu.m            # 单事件流程控制器
├─ run_kh005_matlab_irfu.m         # KH005示例
├─ run_kh005_eigenstability_compare.m # KH005有/无持续5 s对比
├─ kh_fote_batch_irfu.m            # 全事件批处理
└─ +khfote/
   ├─ parseOptions.m               # 参数和默认目录
   ├─ setupEnvironment.m           # IRFU/FOTE路径和MMS数据库
   ├─ loadMmsData.m                # 四星原生数据读取
   ├─ prepareData.m                # 物理时间平均、FOTE前irf_resamp、缺口保留
   ├─ runOriginalFote.m            # 调用原FOTE函数并整理零点类型
   ├─ computeFoteVError.m          # 计算alpha=|div(nV)|/|curl(nV)|
   ├─ applyQualityAndSummarize.m   # 40%、特征值稳定度和持续时间判据
   ├─ plotEvent.m                  # 7-panel图
   └─ writeResults.m               # PDF/PNG/CSV/MAT/JSON输出
```

每个模块文件的开头都说明职责、输入输出和单位；关键公式、缺口判定、平滑边界与类型编码位置也有行内注释。

## 调用的原程序包函数

- `FOTE\FOTE_function\FOTE_Taylor_Expansion.m`：分别输入四星磁场和四星流速，得到零点位置、类型、`eta`和`xi`。
- IRFU-MATLAB `c_4_grad.m`：计算FOTE-V连续性误差`alpha=|div(nV)|/|curl(nV)|`。

`FOTE_Taylor_Expansion`按8参数形式调用，其内部可选平滑分支不会执行。程序先在每种产品的原生时间戳上计算中心物理时间平均；进入FOTE和FOTE-V计算前，再调用IRFU `irf_resamp(...,'linear')`对齐到共同离子时间。`irf_resamp`按连续数据块分别调用，只保留源数据块起止时间以内的结果，因此不会跨缺口或使用区间外推值。原始场panel仍只选择真实观测样本。任一卫星缺数时，四星平均曲线和拓扑结果都会在相应时刻显示空缺。KH005的5 s平均窗对应约640个FGM点和33个DIS点。

当前逐时刻40%判据为：磁场要求`eta<=40%`且`xi<=40%`，流场要求`alpha<=40%`。As/Bs还要求特征值稳定度`M_lambda=|lambda_r|/max(|lambda_i|)>=0.5`，其中`lambda_r`是复共轭特征值组中的唯一实特征值。启用持续时间时，As和Bs分别要求连续合格至少5 s，类型翻转会切断连续区段；将`MinQualityDurationSeconds`设为0即可得到逐点版本。Poincaré index已从计算流程和图中移除。零点距离保存在CSV和MAT中供诊断，不参与质量筛选，图中也不保留距离panel。

## 运行方法

单事件：

```matlab
result = kh_fote_event('KH005', ...
  '2015-10-01T18:01:24Z','2015-10-01T18:09:00Z', ...
  'SmoothSeconds',5,'QualityPercent',40, ...
  'EigenvalueStabilityThreshold',0.5, ...
  'MinQualityDurationSeconds',5,'VelocityField','Vi');
```

复现KH005：

```matlab
result = run_kh005_matlab_irfu('on');
```

比较保留/去掉持续5 s判据：

```matlab
results = run_kh005_eigenstability_compare('off');
```

全事件批处理：

```matlab
report = kh_fote_batch_irfu('SmoothSeconds',5,'QualityPercent',40, ...
  'EigenvalueStabilityThreshold',0.5, ...
  'MinQualityDurationSeconds',5);
```

每个事件输出单页PDF、PNG、时间序列CSV、MAT和JSON摘要；批处理另行输出合并PDF与运行状态表。
