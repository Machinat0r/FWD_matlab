close all
% 读取数据
filename = '/Users/fwd/Documents/Ti~mor~/M/Carrington_Pole/Data/Digitised_Historical_Bfields/STPETERSBURG_DIGITISED.txt';
data = readtable(filename, 'Delimiter', '\t', 'Format', '%{yyyy-MM-dd HH:mm:ss}D%f');

% 提取日期和数值
dates = data{:, 1};
values = data{:, 2};

% 绘制数据
figure;
plot(dates, values, 'b-o', 'MarkerFaceColor', 'b', 'MarkerSize', 6); % 绿色圆点
xlabel('时间');
ylabel('Colaba B_H (nT)');
title('数据变化趋势');

xlim([min(dates) max(dates)]);
ylim([min(values)-100 max(values)+100]);

% 格式化日期
datetick('x', 'HH:MM', 'keepticks');
grid on;
