function h = cline_arrow(x, y, z, c, minc, maxc, cmap)
% cline_arrow3  用带箭头的小线段绘制三维彩色曲线
%
%   h = cline_arrow3(x,y,z,c,minc,maxc)
%   h = cline_arrow3(x,y,z,c,minc,maxc,cmap)
%
% 输入：
%   x,y,z    — 坐标向量（长度相同）
%   c        — 与各顶点对应的标量，用于着色
%   minc/maxc— 用于 colormap 归一化的最小/最大值
%   cmap     — （可选）colormap 名称或 M×3 数组，默认 'jet'
%
% 输出：
%   h        — 每段箭头的第一条线对象句柄数组

  if nargin<7
    cmap = 'jet';
  end

  % 1) 计算每一段的颜色
  if ischar(cmap) || isstring(cmap)
    cmap_mat = colormap(cmap);
  else
    cmap_mat = cmap;
  end
  % 把 c 插值到 [1,size(cmap_mat,1)] 的索引
  idx = interp1(linspace(minc, maxc, size(cmap_mat,1)), ...
                1:size(cmap_mat,1), ...
                c);
  % 确保索引在范围内并取整
  idx = round(max(min(idx, size(cmap_mat,1)), 1));
  colors = cmap_mat(idx, :);

  n = numel(x);
  % h = gobjects(n-1,1);

  % 2) 对每一段，用 arrow3 画箭头
  W = 0.8;     % 箭柄宽度，单位是 PlotBox 对角线的 1/72
  H = 1;   % 箭头高度（一般取 3 倍于宽度）
  % for i = 1:n-1
  dspan = 60;st = 40;
  P1 = [x(1:dspan:end-dspan)',   y(1:dspan:end-dspan)',   z(1:dspan:end-dspan)'];
    P0 = [x(st:dspan:end)', y(st:dspan:end)', z(st:dspan:end)'];
    if size(P1,1) > size(P0,1)
        P1 = P1(1:end-1,:);
    elseif size(P1,1) < size(P0,1)
        P0 = P0(1:end-1,:);
    end
    % close all;
    colors = colors(1:dspan:end,:);

    set(gca,'ColorOrder',colors);
    arrow3(P0, P1,'o-2',W,H);hold on;

    box_bd1 = 1000;box_bd2 = 250;box_bd3 = 300;
    set(gca, 'Ylim',[-box_bd1 box_bd1], 'Ylim',[-box_bd2 box_bd2], 'Zlim',[-box_bd3 box_bd3]);
    set(gca,'xtick',[-box_bd1:box_bd1:box_bd1], 'ytick',[-box_bd2:box_bd2:box_bd2], 'ztick',[-box_bd3:box_bd3:box_bd3],'fontsize',23);
    set(gca,'DataAspectRatio',[5 1.0 1]);
    set(gca,'view',[90 0])
    % set(gca,'YTickLabel',[]);set(gca,'ZTickLabel',[])
    xlabel('x');ylabel('y')
    light('position',[800,-100,-100],'style','infinite'), lighting gouraud
    % for hobj = hn(:)'
    %   if isprop(hobj, 'Color')
    %     hobj.Color     = colors(i,:);
    %   end
    %   if isprop(hobj, 'FaceColor')
    %     hobj.FaceColor = colors(i,:);
    %   end
    % end
    % h(i) = hn(1);
  % end

end
