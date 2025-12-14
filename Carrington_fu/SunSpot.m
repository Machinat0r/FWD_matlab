clear;clc
%% 0. 文件名（按你本地的名字修改）
monthFile = '/Users/fwd/Documents/Ti~mor~/M/Carrington_Pole/Data/SN_m_tot_V2.0.txt';   % 月平均数据（你第一次发的文件）
yearFile  = '/Users/fwd/Documents/Ti~mor~/M/Carrington_Pole/Data/SN_y_tot_V2.0.txt';   % 年平均数据（你第二次发的文件）

%% 1. 读取月平均数据
% 格式：year  month  decYear  monthlySSN  err  nObs
data_m = readmatrix(monthFile);

year_m  = data_m(:, 1);
month_m = data_m(:, 2);
ssn_m   = data_m(:, 4);   % 月平均太阳黑子数

% 选取 1840–1880 年
mask_m = (year_m >= 1840) & (year_m <= 1880);
year_m_sub  = year_m(mask_m);
month_m_sub = month_m(mask_m);
ssn_m_sub   = ssn_m(mask_m);

% 时间轴：每个月用 1 号表示
t_m = datetime(year_m_sub, month_m_sub, 1);

%% 2. 读取年平均数据（直接用你第二个文件）
% 格式：decYear  yearlySSN  err  nObs  （例如 1840.5  105.5  13.0  248）
data_y = readmatrix(yearFile);

decYear = data_y(:, 1);   % 十进制年份（如 1840.5）
ssn_y_all = data_y(:, 2); % 年平均太阳黑子数

% 选取 1840–1880 年
mask_y = (decYear >= 1840.5) & (decYear <= 1880.5);
decYear_sub = decYear(mask_y);
ssn_y       = ssn_y_all(mask_y);

% 把十进制年份转换成 datetime
% 这些年一般是 N+0.5，表示大约在该年的中间，我们就放在 7 月 1 日
year_int = floor(decYear_sub);
t_y = datetime(year_int, 7, 1);

%% 3. 作图：月平均用浅灰“误差条风格”，年平均用粗线
figure;
hold on;

% (1) 月平均：用 stem 模拟“误差bar”的浅灰色竖线，视觉上比较淡
% h_month = stem(t_m, ssn_m_sub, ...
%     'Marker', 'none', ...            % 不画点，只画竖线
%     'LineStyle', '-', ...
%     'LineWidth', 0.5, ...
%     'Color', [0.8 0.8 0.8]);         % 浅灰色

% 如果你更想要一条淡灰细折线而不是竖棍，把上面的 stem 换成：
h_month = plot(t_m, ssn_m_sub, ...
    'Color', [0.8 0.8 0.8], ...
    'LineWidth', 0.5);

% (2) 年平均：明显的粗线，这里用红色
h_year = plot(t_y, ssn_y, ...
    'r-o', ...
    'LineWidth', 2);

xlabel('年份');
ylabel('太阳黑子数');
title('1840–1880 年太阳黑子数：月平均（淡灰）与年平均（粗线）');

grid on;
legend([h_month, h_year], ...
       {'月平均（浅灰，淡化）', '年平均（粗线）'}, ...
       'Location', 'best');

xlim([datetime(1840,1,1) datetime(1880,12,31)]);
set(gcf, 'Color', 'w');
hold off;
