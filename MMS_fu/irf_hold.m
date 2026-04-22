function out = irf_hold(data, onlyHorizontal)
% data: [t, y1, y2, ...]，t 为 epochUnix 或 irf_plot 可识别的时间列
% onlyHorizontal = false: 画阶梯（含竖直跳变）
% onlyHorizontal = true : 只画水平线段，不画竖直跳变

if nargin < 2
  onlyHorizontal = false;
end

t = data(:,1);
y = data(:,2:end);
N = size(data,1);
M = size(y,2);

if N < 2
  out = data;
  return;
end

if ~onlyHorizontal
  % 阶梯：水平保持到下一时刻，时刻处竖直跳变
  t2 = zeros(2*N-1,1);
  y2 = zeros(2*N-1,M);
  t2(1)   = t(1);
  y2(1,:) = y(1,:);
  for k = 2:N
    t2(2*k-2)   = t(k);     y2(2*k-2,:) = y(k-1,:); % 水平段终点
    t2(2*k-1)   = t(k);     y2(2*k-1,:) = y(k,:);   % 竖直跳变
  end
else
  % 只水平：每段 [t(k), t(k+1)] 画 y(k)，段与段之间用 NaN 断开
  t2 = nan(3*(N-1),1);
  y2 = nan(3*(N-1),M);
  for k = 1:N-1
    idx = 3*(k-1) + (1:3);
    t2(idx)   = [t(k);   t(k+1); nan];
    y2(idx,:) = [y(k,:); y(k,:); nan(1,M)];
  end
end

out = [t2, y2];
end
