%% IGRF 地磁北极轨迹（1590–2025）地图绘制
clear;close all;
clc

% —— 读取数据 ——
load("NorthPole.mat");
lon = NorthPole(:,1);
lat = NorthPole(:,2);
yr  = NorthPole(:,3);

% 基本清洗
ok = isfinite(lon) & isfinite(lat) & isfinite(yr);
lon = lon(ok); lat = lat(ok); yr = yr(ok);

% 经度统一到 [-180, 180]
lon = mod(lon + 180, 360) - 180;

% 颜色映射
cmap = parula(256);
cmin = min(yr); cmax = max(yr);
cidx = round( (yr - cmin) / max(eps, (cmax - cmin)) * 255 ) + 1;
cidx = max(1, min(256, cidx));

% 绘图
figure('Color','w'); 
ax = axesm('stereo','Origin',[90 0 0],'MapLatLimit',[65 90]); % 北极立体投影
axis off; framem on; gridm on;


latlim = getm(ax,'MapLatLimit');      % [65 90]
mlabel_parallel = latlim(1) - 0.8;    % 比边界再小一点的纬度(在圆外)
setm(ax, ...
    'MLineLocation', 15, ...
    'PLineLocation', 5, ...
    'MLabelLocation', 90, ...
    'MLabelParallel', mlabel_parallel, ...
    'PLabelLocation', 20, ...
    'PLabelMeridian', 0, ...
    'LabelRotation','off', ...
    'FontSize', 9);
mlabel on;      % 打开经度标注
plabel on;      % 若不需要纬度可改为 off

% 为防止被裁剪，关闭经度文字的剪裁
hM = handlem('MLabel'); 
if ~isempty(hM), set(hM,'Clipping','off'); end

% 海岸线
load coastlines
coastlon = mod(coastlon + 180, 360) - 180;
geoshow(coastlat, coastlon, 'Color', [0.75 0.75 0.75]);

% 连接轨迹
plotm(lat, lon, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.5);

% 关键年份标记
ii = 1859 - yr(1) + 1;
scatterm(lat(1),  lon(1),  20, 'b', 'filled');
scatterm(lat(end),lon(end),20, 'b', 'filled');
if ii>=1 && ii<=numel(yr)
    scatterm(lat(ii), lon(ii), 20, 'r', 'filled');
    % === 关键：把“1859”放在红点左侧 ===
    textm(lat(ii), lon(ii), sprintf('%d', yr(ii)), ...
          'FontSize',10, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
end

% 起止年份标注（保持不变）
textm(lat(1),  lon(1),  sprintf('%d', yr(1)),  'FontSize',10, ...
      'VerticalAlignment','bottom');
textm(lat(end),lon(end),sprintf('%d', yr(end)),'FontSize',10, ...
      'VerticalAlignment','bottom');
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')