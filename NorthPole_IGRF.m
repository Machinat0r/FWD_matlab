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


%% figure 2
latLim = [66 74];
lonLim = [-105 -90];
inZoom = lat >= latLim(1) & lat <= latLim(2) & ...
         lon >= lonLim(1) & lon <= lonLim(2);

% 1859 年位置
idx1859 = find(yr == 1859, 1);
lat1859 = lat(idx1859);
lon1859 = lon(idx1859);

figure('Color','w'); hold on;

% 1. 画海岸线（转成 x=纬度, y=经度）
load coastlines
coastlon = mod(coastlon + 180, 360) - 180;
plot(coastlat, coastlon, 'Color',[0.85 0.85 0.85], 'LineWidth',1);  % x=lat, y=lon

% 2. 画轨迹散点（只画点）
scatter(lat(inZoom), lon(inZoom), 15, [0.2 0.2 0.2], 'filled'); % 其他年份，稍小

% 1859 年：蓝色五角星
scatter(lat1859, lon1859, 60, 'b', 'p', 'filled');   % 'p' 为五角星
text(lat1859, lon1859, ' 1859', ...
    'FontSize',10, 'VerticalAlignment','middle', ...
    'HorizontalAlignment','left', 'Color','b');

% 3. 坐标轴与样式
xlim(latLim);
ylim(lonLim);
set(gca,'YDir','reverse');          % 经度轴倒过来（上大下小）
box on; grid on;
set(gca,'GridLineStyle',':');

xt = 66:2:74;
yt = -105:5:-90;
xticks(xt);
yticks(yt);
xticklabels(strcat(string(xt),      char(176),'N'));
yticklabels(strcat(string(abs(yt)), char(176),'W'));

xlabel('Latitude');
ylabel('Longitude');

