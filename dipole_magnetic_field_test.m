clear;clc;close all
%------written by Wending Fu, Jan.2024 in Beijing------------

%% Parameter (two monopoles)
units = irf_units;
k0 = 1/(4*pi);
Q1 = 1e6;Q2 = -1e6;
r1 = [0,0,0.1];r2 = [0,0,-0.1];
x = -50:0.5:50;x=x';

a1 = 10*[0,0,3]; a2 = 10*[-2*sqrt(2)-1,2,-1]; a3 = 10*[sqrt(2)+2,sqrt(6)+1,-2]; a4 = 10*[sqrt(2)-1,-sqrt(6)+1,-3];
c_eval('a? = x.*[1,0,0] + ones(length(x),1).*a?;',1:4);
temp = 1:length(a1); temp = temp';
%% calculate B (two monopoles)
c_eval('d1? = sqrt((a?(:,1)-r1(1)).^2+(a?(:,2)-r1(2)).^2+(a?(:,3)-r1(3)).^2);',1:4);
c_eval('d2? = sqrt((a?(:,1)-r2(1)).^2+(a?(:,2)-r2(2)).^2+(a?(:,3)-r2(3)).^2);',1:4);
c_eval("Bt1? = transpose(k0*Q1./(d1?'.^2));",1:4);
c_eval("Bt2? = transpose(k0*Q2./(d2?'.^2));",1:4);
c_eval('Bt1?(d1?<=0,1)=0;',1:4);c_eval('Bt2?(d2?<=0,1)=0;',1:4);
c_eval('B?(:,1) = Bt1?.*((a?(:,1)-r1(1))./d1?) + Bt2?.*((a?(:,1)-r2(1))./d2?);',1:4);
c_eval('B?(:,2) = Bt1?.*((a?(:,2)-r1(2))./d1?) + Bt2?.*((a?(:,2)-r2(2))./d2?);',1:4);
c_eval('B?(:,3) = Bt1?.*((a?(:,3)-r1(3))./d1?) + Bt2?.*((a?(:,3)-r2(3))./d2?);',1:4);
c_eval('Bt? = irf_abs(B?);',1:4);
c_eval('Bt? = Bt?(:,4);',1:4);
c_eval("a? = [temp,a?];")
c_eval('B? = [temp,B?];',1:4);
c_eval('Bt? = [temp,Bt?];',1:4);
%%
% 定义磁场计算的函数
calculateMagneticField = @(r) k0 * (Q1 ./ sum((r - r1).^2, 2).^1.5) .* (r - r1) + ...
                                k0 * (Q2 ./ sum((r - r2).^2, 2).^1.5) .* (r - r2);

% 生成空间网格
[x, y, z] = meshgrid(-50:5:50, -50:5:50, -50:5:50);
gridPoints = [x(:), y(:), z(:)];

% 计算磁场在网格点上的值
B = calculateMagneticField(gridPoints);

% 画出磁力线
figure;
quiver3(gridPoints(:,1), gridPoints(:,2), gridPoints(:,3), B(:,1), B(:,2), B(:,3), 'LineWidth', 1, 'AutoScale', 'on');
streamline(x, y, z, reshape(B(:,1), size(x)), reshape(B(:,2), size(y)), reshape(B(:,3), size(z)),...
    gridPoints(1:5:end,1), gridPoints(1:5:end,2), gridPoints(1:5:end,3));

hold on;
% % % plot3(a1(:,1), a1(:,2), a1(:,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% % % plot3(a2(:,1), a2(:,2), a2(:,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% % % plot3(a3(:,1), a3(:,2), a3(:,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% % % plot3(a4(:,1), a4(:,2), a4(:,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Magnetic Field Lines near Magnetic Dipoles');
grid on;
axis equal;

% 可选：设置坐标轴范围，以适应您的数据
% xlim([-60, 60]);
% ylim([-60, 60]);
% zlim([-60, 60]);

hold off;
