clear; close all; clc

% ===== 读取 & 预处理数据 =====
load("NorthPole.mat");      % [lon, lat, yr]
lon = NorthPole(:,1);
lat = NorthPole(:,2);
yr  = NorthPole(:,3);

ok  = isfinite(lon) & isfinite(lat) & isfinite(yr);
lon = lon(ok); lat = lat(ok); yr = yr(ok);

% 经度统一到 [-180,180]
lon = mod(lon + 180, 360) - 180;

% 关键年份索引（如 1859）
idx1859 = find(yr == 1859, 1);


%% ================== 图 1：北极立体投影（m_map） ==================
figure('Color','w');

% 只画 65–90N 的北极立体投影，中央经线 0°，范围类似示意图
% rad = 90 - 最低纬度 = 90 - 65 = 25
m_proj('stereographic','lat',90,'long',0,'rad',25);

% 海洋背景色（浅蓝）
oceanColor = [214 232 242]/255;
set(gca,'Color',oceanColor);
hold on

% 陆地颜色（“土黄色”）
landColor  = [210 190 140]/255;

% 高精度海岸线 + 填色，白色描边
% 你也可以把 _h 换成 _i 或 _f，看你安装了哪一级精度
m_gshhs_h('patch',landColor,'edgecolor','w');

% 网格与标注（可根据喜好再调）
m_grid('linest','--', ...         % 虚线
       'tickstyle','dd', ...      % 度°格式
       'fontsize',9, ...
       'xaxislocation','bottom', ...
       'yaxislocation','left');

% 轨迹：粗红线
m_plot(lon, lat, '-', 'Color',[0.7 0 0], 'LineWidth',2);

% 起点、终点用橙色大点
m_plot(lon(1),  lat(1),  'o', ...
       'markerfacecolor',[1 0.4 0], 'markeredgecolor','none', 'markersize',7);
m_plot(lon(end),lat(end),'o', ...
       'markerfacecolor',[1 0.4 0], 'markeredgecolor','none', 'markersize',7);

% 1859 年位置（若存在）
if ~isempty(idx1859)
    m_plot(lon(idx1859), lat(idx1859), 'o', ...
           'markerfacecolor','b','markeredgecolor','none','markersize',7);
    % 把 “1859” 标在点的左侧
    m_text(lon(idx1859), lat(idx1859), sprintf(' %d', yr(idx1859)), ...
           'FontSize',10, ...
           'VerticalAlignment','middle', ...
           'HorizontalAlignment','right', ...
           'Color','b');
end

% 起止年份标注
m_text(lon(1),  lat(1),  sprintf(' %d', yr(1)),  ...
       'FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','center');
m_text(lon(end),lat(end),sprintf(' %d', yr(end)),...
       'FontSize',10,'VerticalAlignment','bottom','HorizontalAlignment','center');

set(gcf,'Renderer','painters');
set(gcf,'PaperPositionMode','auto');


%% ================== 图 2：放大局部（m_map） ==================
latLim = [69 69.5];
lonLim = [-97 -96.5];
inZoom = lat >= latLim(1) & lat <= latLim(2) & ...
         lon >= lonLim(1) & lon <= lonLim(2);

figure('Color','w', ...
       'Units','pixels', ...
       'Position',[200 200 600 600]);

% 这块区域用 Lambert 投影
m_proj('miller','lat',latLim,'long',lonLim);

% 海洋背景
set(gca,'Color',oceanColor);
hold on

% 高精度海岸线（同样：土黄色 + 白色描边）
m_gshhs_h('patch',landColor,'edgecolor','w');

% 网格
m_grid('linest',':', ...
       'tickstyle','dd', ...
       'fontsize',8, ...
       'box','fancy');

% 把窗口内的所有点连成线
idxZoom = find(inZoom);
m_plot(lon(idxZoom), lat(idxZoom), ...
       '-', 'Color',[0.6 0 0], 'LineWidth',1.5);
m_plot(lon(idxZoom), lat(idxZoom), ...
       'o', 'MarkerFaceColor',[0.6 0 0], ...
       'MarkerEdgeColor','w', 'MarkerSize',4);

% 1859 年位置（若存在）
if ~isempty(idx1859)
    m_plot(lon(idx1859), lat(idx1859), 'o', ...
           'markerfacecolor','b','markeredgecolor','none','markersize',4);
    % 把 “1859” 标在点的左侧
    m_text(lon(idx1859), lat(idx1859), sprintf(' %d', yr(idx1859)), ...
           'FontSize',10, ...
           'VerticalAlignment','middle', ...
           'HorizontalAlignment','right', ...
           'Color','b');
end

% === 给每一个点都标出年份 ===
for k = idxZoom(1:5:end)'
    % 文本稍微偏一点，避免遮住点
    m_text(lon(k), lat(k), sprintf(' %d', yr(k)), ...
           'FontSize',8, ...
           'VerticalAlignment','bottom', ...
           'HorizontalAlignment','left', ...
           'Color',[0.1 0.1 0.1]);
end

% 不需要经纬度坐标轴名称，m_grid 已经标出来；如果想加可以用 xlabel/ylabel
box on
