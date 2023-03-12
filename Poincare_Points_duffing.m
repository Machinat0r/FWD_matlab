close all; clear; clc
% 初始状态向量
x0 = [0.5; 0];

% 参数
delta = 0.1;
alpha = -1;
beta = 1;
gamma = 0.5;
omega = 1.2;

% 模拟时间范围
tspan = [0 1000];

% 模拟物理空间轨迹
[t, x] = ode45(@(t, x) duffing(t, x, delta, alpha, beta, gamma, omega), tspan, x0);

% 绘制物理空间轨迹
figure;
plot(t, x(:,1));
xlabel('时间');
ylabel('位移');
title('物理空间轨迹');

% 绘制相空间轨迹和庞加莱截面
figure;
scatter(x(1,1),x(1,2),'r.');
hold on;
scatter(x(2:end,1),x(2:end,2),'.');
xlabel('位移');
ylabel('速度');
title('相空间轨迹和庞加莱截面');
for ii=1:10:length(x)-10
    if (x(ii,2) > 0 && x(ii+10,2) < 0)
        plot(x(ii+10,1),x(ii+10,2),'r.');
    end
end

% 绘制庞加莱截面图
figure;
poincare_points = x(1:end-1, 2).* (x(2:end, 2) < x(1:end-1, 2));
scatter(x(1:end-1, 1), x(1:end-1, 2), 1, poincare_points, 'filled');
xlabel('位移');
ylabel('速度');
title('庞加莱截面');
colormap jet;

%%
function dxdt = duffing(t, x, delta, alpha, beta, gamma, omega)
    dxdt = [x(2); -delta*x(2) - alpha*x(1) - beta*x(1)^3 + gamma*cos(omega*t)];
end