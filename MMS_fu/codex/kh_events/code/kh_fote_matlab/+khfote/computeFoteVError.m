function [alpha,fluxField] = computeFoteVError(t,R,V,N)
%KH FOTE-V误差模块：用IRFU梯度计算粒子通量连续性误差alpha。
%
% FOTE-V中速度场允许存在物理散度，因此程序使用粒子通量nV：
%   alpha = 100 * |div(nV)| / |curl(nV)|
%
% 输入
%   t  Unix秒；R为Nx3x4 km；V为Nx3x4 km/s；N为Nx4 cm^-3。
%
% 输出
%   alpha      Nx1百分比误差；alpha<=QualityPercent时通过当前判据。
%   fluxField  Nx3x4的nV，保存在MAT文件中供复核。
% 含NaN的时刻不会传入c_4_grad，alpha在对应位置保持NaN。

fluxField = V;
for sc = 1:4
    fluxField(:,:,sc) = V(:,:,sc).*N(:,sc);
end

alpha = nan(numel(t),1);
validInput = isfinite(t(:));
for sc = 1:4
    validInput = validInput & all(isfinite(R(:,:,sc)),2) & ...
        all(isfinite(fluxField(:,:,sc)),2);
end
if ~any(validInput), return; end

[R1,R2,R3,R4] = tensorToFourSeries(t(validInput),R(validInput,:,:));
[F1,F2,F3,F4] = tensorToFourSeries( ...
    t(validInput),fluxField(validInput,:,:));
divergence = c_4_grad(R1,R2,R3,R4,F1,F2,F3,F4,'div');
curlField = c_4_grad(R1,R2,R3,R4,F1,F2,F3,F4,'curl');
curlMagnitude = sqrt(sum(curlField(:,2:4).^2,2));

validRatio = isfinite(divergence(:,2)) & isfinite(curlMagnitude) & ...
    curlMagnitude>0;
alphaValid = nan(nnz(validInput),1);
alphaValid(validRatio) = 100*abs(divergence(validRatio,2))./ ...
    curlMagnitude(validRatio);
alpha(validInput) = alphaValid;
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
