%% 从图片中读取 Dst 曲线并复刻图像
% 1. 准备：把你给我的那张图保存为 dst_input.png
img = imread('/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/Carrington_fu/dst_input.png');      % 原图大小应为 476×1766×3 或 ×4

% 2. 裁剪出只包含坐标轴和曲线的区域
%    这些索引是针对你这张图（1766×476）的，如果将来图片尺寸变了需要手动微调
crop = img(27:368, 232:1741, :);    % [行, 列, RGB]
[H, W, ~] = size(crop);

% 3. 分离 R/G/B 通道，构造“蓝线”和“红线”的颜色掩膜
R = double(crop(:,:,1));
G = double(crop(:,:,2));
B = double(crop(:,:,3));

maxRG  = max(R, G);
maxRGB = max(maxRG, B);             % 每个像素的最大通道值

% 颜色判据（根据图像的实际颜色调出来的阈值）：
blueMask = (B > G + 10) & (B > R + 10) & (maxRGB < 240);      % 蓝色曲线像素
redMask  = (R > G + 20) & (R > B + 20) & (maxRGB < 240);      % 红色曲线像素

% 4. 对每一列寻找该列上蓝线/红线的平均行坐标（即 y 像素）
yb_pix = nan(1, W);   % 蓝线 y，0 基准坐标
yr_pix = nan(1, W);   % 红线 y，0 基准坐标

for j = 1:W
    rowsB = find(blueMask(:, j));
    if ~isempty(rowsB)
        yb_pix(j) = mean(rowsB - 1);    % Matlab 行号转为 0-based 像素坐标
    end

    rowsR = find(redMask(:, j));
    if ~isempty(rowsR)
        yr_pix(j) = mean(rowsR - 1);
    end
end

% 处理少量缺失点：线性插值补齐
x_col = 1:W;
yb_pix = fillmissing(yb_pix, 'linear', 'EndValues','nearest');
yr_pix = fillmissing(yr_pix, 'linear', 'EndValues','nearest');

% 5. 像素坐标 → 物理坐标的标定
%    水平方向：从 05:00 到 12:00，线性映射
t0 = 5;                 % 05:00
t1 = 12;                % 12:00
t  = t0 + ( (x_col - 1) / (W - 1) ) * (t1 - t0);   % 单位：小时 UT

%    垂直方向 Dst：根据网格线位置标定
%    在裁剪后的图中，网格线的行像素（0-based）大致为：
%    35, 109, 182, 254, 326  ↔  Dst = 0, -250, -500, -750, -1000 nT
y0 = 35;     % 对应 Dst =   0 nT 的网格线（0-based 行坐标）
yN = 326;    % 对应 Dst = -1000 nT 的网格线

dst_blue = -1000 * ( (yb_pix - y0) / (yN - y0) );
dst_red  = -1000 * ( (yr_pix - y0) / (yN - y0) );

% 可选：轻微裁剪到显示范围
dst_blue(dst_blue >  50)   =  50;
dst_blue(dst_blue < -1000) = -1000;
dst_red(dst_red   >  50)   =  50;
dst_red(dst_red   < -1000) = -1000;

% 现在 t, dst_blue, dst_red 就是从图中“读”出的数据
% 如需保存：
% save('dst_from_figure.mat','t','dst_blue','dst_red');

%% 6. 用提取的数据复刻原图
figure('Units','centimeters','Position',[2 2 22 6]);
hold on; box on;

plot(t, dst_blue, 'Color',[0 0.45 0.74], 'LineWidth',2);      % 蓝线
plot(t, dst_red , 'Color',[0.85 0.33 0.10], 'LineWidth',2);   % 红线

% 竖直虚线（按原图大致时间位置）
xline(5.0 , 'k--', 'LineWidth',2);    % 05:00
xline(6.5 , 'k--', 'LineWidth',2);    % 06:30
xline(8.0 , 'k--', 'LineWidth',2);    % 08:00

xlim([5 12]);
ylim([-1000 50]);
yticks(0:-250:-1000);

xticks(5:12);
xticklabels({'05:00','06:00','07:00','08:00', ...
             '09:00','10:00','11:00','12:00'});

ax = gca;
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridColor = [0.7 0.7 0.7];
ax.GridLineStyle = '-';
ax.LineWidth = 1;
ax.FontSize = 10;
ax.FontName = 'Times New Roman';

ylabel('Dst (nT)','FontSize',12);
xlabel('02 Sep.','FontSize',12);

set(gcf,'Color','w');
