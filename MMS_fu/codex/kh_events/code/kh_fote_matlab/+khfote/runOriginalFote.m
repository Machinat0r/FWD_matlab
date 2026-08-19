function result = runOriginalFote(t,R,F)
%KH FOTE 核心计算模块：调用用户原有FOTE函数分析一个矢量场。
%
% 模块职责
%   1. 将Nx3x4数组转换为原FOTE函数要求的四个[t,x,y,z]矩阵；
%   2. 以8参数形式调用FOTE_Taylor_Expansion，关闭其内部smooth；
%   3. 调用IRFU c_4_grad获得每个时刻的3x3梯度矩阵和本征值；
%   4. 整理零点类型、screw方向、零点距离和FOTE误差。
%
% 输入单位
%   t : Unix秒
%   R : km，大小Nx3x4
%   F : 磁场nT或速度km/s，大小Nx3x4
%
% 输出关键字段
%   type        A/B/As/Bs/X/O/S+/S-
%   screwSense As=-1（screw-in），Bs=+1（screw-out）
%   distanceKm  零点到最近卫星的距离，单位km
%   distanceL   distanceKm除以四面体六条边的平均长度L
%   eta, xi     原FOTE函数给出的百分比误差
%
% 含NaN的时刻不会传入原FOTE函数。其结果字段保持NaN/type="missing"，
% 从而在图和CSV中保留真实数据缺口。

%% 1. 建立完整时间轴结果，并识别四星数据齐全的时刻
n = size(F,1);
valid = isfinite(t(:));
for sc = 1:4
    valid = valid & all(isfinite(R(:,:,sc)),2) & ...
        all(isfinite(F(:,:,sc)),2);
end

result = initializeResult(n,valid);
if ~any(valid), return; end

rows = find(valid);
tValid = t(valid);
RValid = R(valid,:,:);
FValid = F(valid,:,:);
[R1,R2,R3,R4] = tensorToFourSeries(tValid,RValid);
[F1,F2,F3,F4] = tensorToFourSeries(tValid,FValid);

%% 2. 仅对真实四星样本齐全的时刻调用原FOTE函数
% 某些远距离零点会产生病态梯度。保留结果，并暂时关闭重复警告。
warningState1 = warning('off','MATLAB:singularMatrix');
warningState2 = warning('off','MATLAB:nearlySingularMatrix');
cleanup = onCleanup(@()restoreWarnings(warningState1,warningState2));

% 仅传入8个必需参数，原函数内部的可选smooth分支不会执行。
[distanceInfo,typeInfo,errorInfo] = FOTE_Taylor_Expansion( ...
    R1,R2,R3,R4,F1,F2,F3,F4);
grad = c_4_grad(R1,R2,R3,R4,F1,F2,F3,F4,'grad');
clear cleanup

%% 3. 将有效结果放回完整时间轴
validEigenvalues = complex(nan(numel(rows),3));
for j = 1:numel(rows)
    G = reshape(grad(j,2:10),3,3);
    result.gradient(:,:,rows(j)) = G;
    validEigenvalues(j,:) = eig(G).';
end

result.eigenvalues(valid,:) = validEigenvalues;
result.type(valid) = decodeOriginalTypes(typeInfo,validEigenvalues);
result.screwSense(result.type=="As") = -1;
result.screwSense(result.type=="Bs") = 1;
result.distanceKm(valid) = distanceInfo.dRmin(:,2);
result.scaleL(valid) = tetrahedronMeanEdge(RValid);
result.distanceL = result.distanceKm./result.scaleL;
result.eta(valid) = errorInfo.err1(:,2);
result.xi(valid) = errorInfo.err2(:,2);
end


function result = initializeResult(n,valid)
% 预分配完整时间轴。缺数行保持NaN并明确标记为missing。
result = struct;
result.gradient = nan(3,3,n);
result.eigenvalues = complex(nan(n,3));
result.type = repmat("missing",n,1);
result.screwSense = nan(n,1);
result.distanceKm = nan(n,1);
result.scaleL = nan(n,1);
result.distanceL = nan(n,1);
result.eta = nan(n,1);
result.xi = nan(n,1);
result.validInput = valid;
result.engine = 'FOTE_Taylor_Expansion';
end


function restoreWarnings(state1,state2)
warning(state1);
warning(state2);
end


function type = decodeOriginalTypes(info,eigenvalues)
% 原程序以marker和face color编码类型，此处转换为可统计的字符串。
n = size(eigenvalues,1);
type = repmat("degenerate",n,1);
markers = char(info.type);
faces = char(info.f_color);
markers = markers(:);
faces = faces(:);
for i = 1:min([n,numel(markers),numel(faces)])
    switch markers(i)
        case '^'
            if faces(i)=='r', type(i)="As"; else, type(i)="A"; end
        case '>'
            if faces(i)=='b', type(i)="Bs"; else, type(i)="B"; end
        case 'X'
            type(i) = "X";
        case {'o','O'}
            type(i) = "O";
        case 's'
            % source/sink通过三个本征值实部的符号区分。
            ev = eigenvalues(i,:);
            scale = max(abs(ev));
            if isfinite(scale) && scale>0
                if all(real(ev)>1e-8*scale), type(i)="S+"; end
                if all(real(ev)<-1e-8*scale), type(i)="S-"; end
            end
    end
end
end


function [S1,S2,S3,S4] = tensorToFourSeries(t,tensor)
% 将Nx3x4张量转换成IRFU/FOTE常用的四个[t,x,y,z]矩阵。
t = t(:);
if numel(t)~=size(tensor,1)
    error('Time/tensor length mismatch.');
end
S1 = [t,tensor(:,:,1)];
S2 = [t,tensor(:,:,2)];
S3 = [t,tensor(:,:,3)];
S4 = [t,tensor(:,:,4)];
end


function scale = tetrahedronMeanEdge(R)
% L定义为四面体六条边长的算术平均，单位与R相同。
pairs = [1 2;1 3;1 4;2 3;2 4;3 4];
n = size(R,1);
separation = nan(n,6);
for k = 1:6
    delta = R(:,:,pairs(k,1))-R(:,:,pairs(k,2));
    separation(:,k) = sqrt(sum(delta.^2,2));
end
scale = mean(separation,2,'omitnan');
end
