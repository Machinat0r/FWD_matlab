function piValue = computePoincare(t,F)
%KH FOTE Poincaré模块：调用IRFU识别四面体内矢量场零点拓扑度。
%
% 输入
%   t  Unix秒；F为Nx3x4四星矢量场，可输入B或V。
%
% 输出
%   piValue通常接近-1、0或+1。|PI|>=0.5用于标记非零Poincaré index。
%   PI=0表示四面体内零点总拓扑度为0；它也可能对应偶数个符号相反零点。
% 含任意四星矢量缺失的时刻保持NaN，不会用周围时刻填补。

piValue = nan(numel(t),1);
validInput = isfinite(t(:));
for sc = 1:4
    validInput = validInput & all(isfinite(F(:,:,sc)),2);
end
if ~any(validInput), return; end

[F1,F2,F3,F4] = tensorToFourSeries(t(validInput),F(validInput,:,:));
raw = c_4_poincare_index(F1,F2,F3,F4);
if size(raw,2)>1
    validValue = raw(:,2);
else
    validValue = raw(:);
end

% 清理浮点舍入噪声，使CSV和图中显示稳定的整数拓扑度。
validValue(abs(validValue)<1e-8) = 0;
nearInteger = isfinite(validValue) & ...
    abs(validValue-round(validValue))<1e-6;
validValue(nearInteger) = round(validValue(nearInteger));
piValue(validInput) = validValue;
end


function [S1,S2,S3,S4] = tensorToFourSeries(t,tensor)
t = t(:);
if numel(t)~=size(tensor,1)
    error('Time/tensor length mismatch.');
end
S1 = [t,tensor(:,:,1)];
S2 = [t,tensor(:,:,2)];
S3 = [t,tensor(:,:,3)];
S4 = [t,tensor(:,:,4)];
end
