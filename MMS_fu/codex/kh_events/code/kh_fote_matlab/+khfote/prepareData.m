function [rawAligned,aligned,smoothingInfo] = prepareData( ...
    raw,tStart,tStop,windowSeconds)
%KH FOTE 时间处理模块：原生时间平均、IRFU对齐和缺口保留。
%
% 模块职责
%   1. 生成未平滑的共同离子时间序列，用于图中Raw B和Raw V panel；
%   2. 在每颗卫星、每种产品的原生时间戳上计算中心时间平均；
%   3. 磁场和流场使用相同的物理时间窗，允许采样点数不同；
%   4. 原始场panel仅选择时间最近的真实观测样本；
%   5. FOTE/FOTE-V/PI输入在计算前使用irf_resamp对齐；
%   6. irf_resamp按连续数据块调用，不跨缺口，也不使用区间外推结果。
%
% windowSeconds=5时，KH005的FGM窗约640点，DIS窗约33点。
% 密度与位置不做时间平均；密度用于FOTE-V连续性误差，位置用于梯度。

%% 1. 未平滑共同时间序列
rawAligned = matchAtMomentCadence(raw,tStart,tStop);

%% 2. 原生时间平均
if windowSeconds>0
    [smoothedRaw,smoothingInfo] = smoothNativeVectorFields( ...
        raw,windowSeconds);
    analysisSource = smoothedRaw;
else
    smoothingInfo = struct([]);
    analysisSource = raw;
end

%% 3. FOTE前使用IRFU标准函数对齐到共同离子时间
aligned = resampleAtMomentCadence(analysisSource,tStart,tStop);

%% 4. 输出采样信息，便于检查物理时间窗
fprintf(['Raw grid samples=%d; analysis grid samples=%d; cadence=%.3f s; ' ...
    'irf_resamp B support=%d, V support=%d.\n'], ...
    numel(rawAligned.t),numel(aligned.t), ...
    median(diff(aligned.t),'omitnan'),nnz(aligned.support.B), ...
    nnz(aligned.support.V));
for k = 1:numel(smoothingInfo)
    fprintf(['MMS%d %s: native cadence %.9f s, %.3f Hz, ' ...
        '%g s window=%d samples.\n'], ...
        smoothingInfo(k).Spacecraft,smoothingInfo(k).Product, ...
        smoothingInfo(k).CadenceSeconds,smoothingInfo(k).SampleRateHz, ...
        windowSeconds,smoothingInfo(k).NominalWindowSamples);
end
end


function aligned = matchAtMomentCadence(raw,tStart,tStop)
% 以MMS1速度矩时刻为共同时间，只选取时间最近的真实观测样本。
% 最近样本与目标时刻之差若超过半个原生采样周期，输出保持NaN。
% 所有目标时刻都保留；程序不会删除含缺失值的行。
t = double(raw.V{1}(:,1));
t = t(isfinite(t) & t>=tStart & t<=tStop);
t = unique(t,'stable');
n = numel(t);
if n<20, error('Reference moment series has fewer than 20 samples.'); end

B = nan(n,3,4);
V = nan(n,3,4);
R = nan(n,3,4);
N = nan(n,4);
for sc = 1:4
    B(:,:,sc) = nearestObservedSample(raw.B{sc},t);
    V(:,:,sc) = nearestObservedSample(raw.V{sc},t);
    R(:,:,sc) = nearestObservedSample(raw.R{sc},t);
    N(:,sc) = nearestObservedSample(raw.N{sc},t);
end
aligned = assembleAligned(t,B,V,R,N);
end


function aligned = resampleAtMomentCadence(raw,tStart,tStop)
% FOTE输入统一到MMS1速度矩时间，所有数值对齐均调用irf_resamp。
% 显式指定linear可避免高频B在已完成5 s平均后再次被窗口平均。
% 每个连续数据块单独处理，避免irf_resamp默认的区间外推跨越缺口。
t = double(raw.V{1}(:,1));
t = t(isfinite(t) & t>=tStart & t<=tStop);
t = unique(t,'stable');
n = numel(t);
if n<20, error('Reference moment series has fewer than 20 samples.'); end

B = nan(n,3,4);
V = nan(n,3,4);
R = nan(n,3,4);
N = nan(n,4);
for sc = 1:4
    B(:,:,sc) = irfResampPreserveGaps(raw.B{sc},t);
    V(:,:,sc) = irfResampPreserveGaps(raw.V{sc},t);
    R(:,:,sc) = irfResampPreserveGaps(raw.R{sc},t);
    N(:,sc) = irfResampPreserveGaps(raw.N{sc},t);
end
aligned = assembleAligned(t,B,V,R,N);
end


function aligned = assembleAligned(t,B,V,R,N)
% 计算各FOTE输入的四星有效掩码，同时保留全部目标时刻。
supportB = isfinite(t);
supportV = isfinite(t);
supportNV = isfinite(t);
for sc = 1:4
    positionFinite = all(isfinite(R(:,:,sc)),2);
    supportB = supportB & all(isfinite(B(:,:,sc)),2) & positionFinite;
    supportV = supportV & all(isfinite(V(:,:,sc)),2) & positionFinite;
    supportNV = supportNV & all(isfinite(V(:,:,sc)),2) & ...
        isfinite(N(:,sc)) & positionFinite;
end
aligned = struct('t',t,'B',B,'V',V,'R',R,'N',N, ...
    'support',struct('B',supportB,'V',supportV,'NV',supportNV, ...
    'Joint',supportB & supportV & supportNV));
end


function values = irfResampPreserveGaps(source,targetTime)
% 在连续时间块内调用irf_resamp，并屏蔽其默认外推区域。
values = nan(numel(targetTime),size(source,2)-1);
source = source(isfinite(source(:,1)),:);
if isempty(source), return; end
[~,ia] = unique(source(:,1),'stable');
source = sortrows(source(ia,:),1);
if size(source,1)==1
    exact = abs(targetTime-source(1,1))<=1e-9;
    values(exact,:) = repmat(source(1,2:end),nnz(exact),1);
    return
end

dt = diff(source(:,1));
baseDt = median(dt(isfinite(dt) & dt>0),'omitnan');
if ~isfinite(baseDt) || baseDt<=0, return; end
gapLimit = max(10*baseDt,1);
starts = [1;find(dt>gapLimit)+1];
stops = [find(dt>gapLimit);size(source,1)];

for block = 1:numel(starts)
    rows = starts(block):stops(block);
    if numel(rows)<2, continue; end
    targetRows = targetTime>=source(rows(1),1) & ...
        targetTime<=source(rows(end),1);
    if ~any(targetRows), continue; end
    blockResult = irf_resamp( ...
        source(rows,:),targetTime(targetRows),'linear');
    if isempty(blockResult), continue; end
    values(targetRows,:) = blockResult(:,2:end);
end
end


function values = nearestObservedSample(source,targetTime)
% 无插值时间匹配：输出值始终来自source中的某一条真实记录。
%
% 采用双指针寻找每个目标时刻的最近源样本。允许的最大时间差为
% 0.51个原生采样周期；0.01余量仅用于浮点时间戳舍入。数据缺口内
% 最近样本会超过该阈值，因此对应输出自然保持NaN。
values = nan(numel(targetTime),size(source,2)-1);
source = source(isfinite(source(:,1)),:);
if isempty(source), return; end
[~,ia] = unique(source(:,1),'stable');
source = sortrows(source(ia,:),1);
if size(source,1)==1
    exact = abs(targetTime-source(1,1))<=1e-9;
    values(exact,:) = repmat(source(1,2:end),nnz(exact),1);
    return
end

dt = diff(source(:,1));
baseDt = median(dt(isfinite(dt) & dt>0),'omitnan');
if ~isfinite(baseDt) || baseDt<=0, return; end
maximumOffset = 0.51*baseDt;

sourceIndex = 1;
for targetIndex = 1:numel(targetTime)
    target = targetTime(targetIndex);
    while sourceIndex<size(source,1) && ...
            abs(source(sourceIndex+1,1)-target) <= ...
            abs(source(sourceIndex,1)-target)
        sourceIndex = sourceIndex+1;
    end
    if abs(source(sourceIndex,1)-target)<=maximumOffset
        % 整行复制真实样本；任何原始NaN会原样保留。
        values(targetIndex,:) = source(sourceIndex,2:end);
    end
end
end


function [out,info] = smoothNativeVectorFields(raw,windowSeconds)
% B与V分别在各自原生时间线上处理；N和R保持原始值。
out = raw;
info = repmat(struct('Spacecraft',0,'Product','', ...
    'CadenceSeconds',NaN,'SampleRateHz',NaN, ...
    'NominalWindowSamples',0),8,1);
counter = 0;
for sc = 1:4
    products = {'B','V'};
    for p = 1:numel(products)
        counter = counter+1;
        product = products{p};
        [out.(product){sc},diagnostics] = centeredTimeMean( ...
            raw.(product){sc},windowSeconds);
        info(counter).Spacecraft = sc;
        info(counter).Product = product;
        info(counter).CadenceSeconds = diagnostics.CadenceSeconds;
        info(counter).SampleRateHz = diagnostics.SampleRateHz;
        info(counter).NominalWindowSamples = ...
            diagnostics.NominalWindowSamples;
    end
end
end


function [out,diagnostics] = centeredTimeMean(source,windowSeconds)
% 按真实时间边界[t-window/2,t+window/2)计算中心算术平均。
% 要求窗口覆盖完整，且有限样本数至少达到名义样本数的80%。
out = source;
out(:,2:end) = NaN;
t = source(:,1);
dt = diff(t);
positive = dt(isfinite(dt) & dt>0);
if isempty(positive)
    diagnostics = struct('CadenceSeconds',NaN,'SampleRateHz',NaN, ...
        'NominalWindowSamples',0);
    return
end

baseDt = median(positive,'omitnan');
nominal = max(1,round(windowSeconds/baseDt));
diagnostics = struct('CadenceSeconds',baseDt, ...
    'SampleRateHz',1/baseDt,'NominalWindowSamples',nominal);
halfWindow = windowSeconds/2;
gapLimit = max(10*baseDt,1);
starts = [1;find(dt>gapLimit)+1];
stops = [find(dt>gapLimit);numel(t)];
minimumCount = max(1,floor(0.8*nominal));

for block = 1:numel(starts)
    rows = starts(block):stops(block);
    m = numel(rows);
    if m<2, continue; end
    ts = t(rows);
    left = ones(m,1);
    right = (m+1)*ones(m,1);
    il = 1;
    ir = 1;
    for i = 1:m
        while il<=m && ts(il)<ts(i)-halfWindow, il=il+1; end
        if ir<il, ir=il; end
        while ir<=m && ts(ir)<ts(i)+halfWindow, ir=ir+1; end
        left(i) = il;
        right(i) = ir;
    end
    complete = (ts-ts(1)>=halfWindow) & (ts(end)-ts>=halfWindow);
    for column = 2:size(source,2)
        y = source(rows,column);
        finite = isfinite(y);
        sums = [0;cumsum(fillNonfinite(y,0))];
        counts = [0;cumsum(double(finite))];
        localSums = sums(right)-sums(left);
        localCounts = counts(right)-counts(left);
        usable = complete & localCounts>=minimumCount;
        values = nan(m,1);
        values(usable) = localSums(usable)./localCounts(usable);
        out(rows,column) = values;
    end
end
end


function out = fillNonfinite(in,value)
out = in;
out(~isfinite(out)) = value;
end
