close all; clear; clc
% 初始状态向量
x0 = [1; 1; 1];

% 模拟时间范围
tspan = [0 50];

% 模拟物理空间轨迹
[t, x] = ode45(@lorenz, tspan, x0);

% 绘制物理空间轨迹
figure;
plot3(x(:,1), x(:,2), x(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('物理空间轨迹');

% 绘制相空间轨迹
figure;
plot3(x(:,1), x(:,2), x(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
title('相空间轨迹');
hold on;

% 绘制庞加莱截面
poincare_section = @(x) x(3);
poincare_points = x(1:end-1, 3).* (x(2:end, 3) < x(1:end-1, 3));
scatter3(x(1:end-1, 1), x(1:end-1, 2), x(1:end-1, 3), 1, poincare_points, 'filled');
colormap jet;

% 绘制庞加莱截面平面
[xs, ys] = meshgrid(-30:30);
zs = 28*ones(size(xs));
surf(xs, ys, zs, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
view(60, 30);

% 绘制庞加莱截面的投影
figure;
scatter(x(1:end-1, 1), x(1:end-1, 2), 3, poincare_points,'filled','o','MarkerFaceColor','k');
xlabel('x');
ylabel('y');
title('庞加莱截面');


%%
function dxdt = lorenz(t, x)
    sigma = 10;
    r = 28;
    b = 8/3;
    
    dxdt = [sigma*(x(2)-x(1)); r*x(1)-x(2)-x(1)*x(3); -b*x(3)+x(1)*x(2)];
end
