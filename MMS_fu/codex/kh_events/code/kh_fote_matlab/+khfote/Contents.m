% +KHFOTE  MMS Kelvin-Helmholtz FOTE/FOTE-V内部模块包
%
% 配置与环境
%   parseOptions               - 解析事件参数和默认目录
%   setupEnvironment           - 配置IRFU/FOTE路径与MMS数据库
%
% 数据处理
%   loadMmsData                - 读取四星B、V、N和R原生时间序列
%   prepareData                - 时间平均、FOTE前irf_resamp对齐和缺口保留
%
% 拓扑计算
%   runOriginalFote            - 调用原FOTE函数计算零点和类型
%   computeFoteVError          - 计算FOTE-V连续性误差alpha
%   applyQualityAndSummarize   - 误差、特征值稳定度与持续时间筛选
%
% 展示与输出
%   plotEvent                  - 绘制7-panel单事件图
%   writeResults               - 保存PDF、PNG、CSV、MAT和JSON
