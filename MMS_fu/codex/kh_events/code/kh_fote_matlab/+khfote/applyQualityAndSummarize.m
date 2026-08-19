function [mag,vel,agreement,summary] = ...
    applyQualityAndSummarize(mag,vel,t,opt)
%KH FOTE 质量控制与统计模块：应用误差、特征值稳定度及持续时间判据。
%
% 当前质量判据
%   磁场FOTE：eta<=阈值 且 xi<=阈值，并要求三个本征值有限；
%   流场FOTE-V：alpha<=阈值，并要求三个本征值有限；
%   As/Bs：M_lambda=|lambda_r|/max(|lambda_i|)>=指定阈值；
%   底部screw panel：每种类型分别连续成立至少指定时间。
%
% 零点距离不参与good筛选。near字段仅表示good且distance/L小于给定诊断阈值。
% good表示逐时刻40%误差判据；spiralQualified进一步包含M_lambda判据；
% sustainedGood表示最终通过持续时间判据的As/Bs点。
% agreement仅比较磁场和流场同时最终合格的点：+1同向，-1反向。

%% 1. FOTE和FOTE-V质量判据
mag.good = all(isfinite(mag.eigenvalues),2) & isfinite(mag.eta) & ...
    isfinite(mag.xi) & mag.eta<=opt.QualityPercent & ...
    mag.xi<=opt.QualityPercent;
vel.good = all(isfinite(vel.eigenvalues),2) & isfinite(vel.alpha) & ...
    vel.alpha<=opt.QualityPercent;

%% 2. As/Bs唯一实特征值的稳定度
mag.eigenvalueStability = spiralEigenvalueStability( ...
    mag.eigenvalues,mag.type);
vel.eigenvalueStability = spiralEigenvalueStability( ...
    vel.eigenvalues,vel.type);
mag.eigenvalueStable = isfinite(mag.eigenvalueStability) & ...
    mag.eigenvalueStability>=opt.EigenvalueStabilityThreshold;
vel.eigenvalueStable = isfinite(vel.eigenvalueStability) & ...
    vel.eigenvalueStability>=opt.EigenvalueStabilityThreshold;

mag.spiralQualified = mag.good & mag.eigenvalueStable & ...
    (mag.type=="As" | mag.type=="Bs");
vel.spiralQualified = vel.good & vel.eigenvalueStable & ...
    (vel.type=="As" | vel.type=="Bs");

% 对As和Bs分别施加持续时间判据，因此类型翻转会主动切断连续区段。
% MinQualityDurationSeconds=0时，函数原样返回逐点合格掩码。
mag.acceptedAs = sustainedQualityRuns( ...
    mag.spiralQualified & mag.type=="As",t, ...
    opt.MinQualityDurationSeconds);
mag.acceptedBs = sustainedQualityRuns( ...
    mag.spiralQualified & mag.type=="Bs",t, ...
    opt.MinQualityDurationSeconds);
vel.acceptedAs = sustainedQualityRuns( ...
    vel.spiralQualified & vel.type=="As",t, ...
    opt.MinQualityDurationSeconds);
vel.acceptedBs = sustainedQualityRuns( ...
    vel.spiralQualified & vel.type=="Bs",t, ...
    opt.MinQualityDurationSeconds);
mag.sustainedGood = mag.acceptedAs | mag.acceptedBs;
vel.sustainedGood = vel.acceptedAs | vel.acceptedBs;

% Panel 5/6中的As/Bs也应用M_lambda判据；其他类型仍只作为40%误差诊断。
mag.markerGood = mag.good & (~(mag.type=="As" | mag.type=="Bs") | ...
    mag.eigenvalueStable);
vel.markerGood = vel.good & (~(vel.type=="As" | vel.type=="Bs") | ...
    vel.eigenvalueStable);

% near保留在数据表中，方便未来检查距离对统计的影响。
mag.near = mag.good & mag.distanceL<=opt.ReliableDistanceL;
vel.near = vel.good & vel.distanceL<=opt.ReliableDistanceL;

%% 3. screw方向一致性
agreement = nan(size(mag.screwSense));
spiral = abs(mag.screwSense)==1 & abs(vel.screwSense)==1;
agreement(spiral & mag.screwSense==vel.screwSense) = 1;
agreement(spiral & mag.screwSense~=vel.screwSense) = -1;
agreement(~(mag.sustainedGood & vel.sustainedGood)) = NaN;

%% 4. 汇总计数
summary = struct;
summary.BQuality40 = nnz(mag.good);
summary.VQuality40 = nnz(vel.good);
summary.BAsQuality40 = nnz(mag.good & mag.type=="As");
summary.BBsQuality40 = nnz(mag.good & mag.type=="Bs");
summary.VAsQuality40 = nnz(vel.good & vel.type=="As");
summary.VBsQuality40 = nnz(vel.good & vel.type=="Bs");
summary.EigenvalueStabilityThreshold = opt.EigenvalueStabilityThreshold;
summary.BStableSpiral = nnz(mag.spiralQualified);
summary.VStableSpiral = nnz(vel.spiralQualified);
summary.BAsStable = nnz(mag.spiralQualified & mag.type=="As");
summary.BBsStable = nnz(mag.spiralQualified & mag.type=="Bs");
summary.VAsStable = nnz(vel.spiralQualified & vel.type=="As");
summary.VBsStable = nnz(vel.spiralQualified & vel.type=="Bs");
summary.MinQualityDurationSeconds = opt.MinQualityDurationSeconds;
summary.BQuality40Continuous = nnz(mag.sustainedGood);
summary.VQuality40Continuous = nnz(vel.sustainedGood);
summary.BAsQuality40Continuous = nnz(mag.acceptedAs);
summary.BBsQuality40Continuous = nnz(mag.acceptedBs);
summary.VAsQuality40Continuous = nnz(vel.acceptedAs);
summary.VBsQuality40Continuous = nnz(vel.acceptedBs);
summary.BAcceptedSpiral = nnz(mag.sustainedGood);
summary.VAcceptedSpiral = nnz(vel.sustainedGood);
summary.ScrewSame = nnz(agreement==1);
summary.ScrewOpposite = nnz(agreement==-1);
end


function stability = spiralEigenvalueStability(eigenvalues,type)
% M_lambda=|lambda_r|/max(|lambda_i|)，仅对As/Bs定义。
% lambda_r取虚部绝对值最小的特征值；实3x3矩阵的螺旋型谱含一个实根
% 和一对复共轭根。M_lambda接近0时，实根符号容易受扰动影响。
n = size(eigenvalues,1);
stability = nan(n,1);
spiralRows = find(type=="As" | type=="Bs");
for k = 1:numel(spiralRows)
    row = spiralRows(k);
    values = eigenvalues(row,:);
    if ~all(isfinite(values)), continue; end
    scale = max(abs(values));
    if ~isfinite(scale) || scale<=0, continue; end
    [~,realIndex] = min(abs(imag(values)));
    stability(row) = abs(real(values(realIndex)))/scale;
end
end


function keep = sustainedQualityRuns(pass,t,minDurationSeconds)
% 保留持续时间达到阈值的完整true区段；短区段全部置false。
pass = logical(pass(:));
t = double(t(:));
if numel(pass)~=numel(t)
    error('Quality mask and time vector must have the same length.');
end
keep = false(size(pass));
if isempty(pass) || minDurationSeconds<0
    return
end
if minDurationSeconds==0
    keep = pass;
    return
end

dt = diff(t);
normalDt = dt(isfinite(dt) & dt>0);
if isempty(normalDt)
    return
end
cadence = median(normalDt);
gapLimit = 3*cadence;

% link(k)=true表示第k和k+1个样本属于同一个连续合格区段。
link = pass(1:end-1) & pass(2:end) & isfinite(dt) & ...
    dt>0 & dt<=gapLimit;
starts = find(pass & [true; ~link]);
stops = find(pass & [~link; true]);
durationTolerance = max(1e-6,1e-6*cadence);
for k = 1:numel(starts)
    if t(stops(k))-t(starts(k))+durationTolerance >= minDurationSeconds
        keep(starts(k):stops(k)) = true;
    end
end
end
