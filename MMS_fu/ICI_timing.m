clear; clc
% 假设文件名为 data.xlsx，第三列是你要处理的文本
filename = "C:\Users\pc\Documents\Timor\Recovery-Work-DF-MP\2-ICI_Magnetopause\Reply-1\velocity_timing.xlsx";

% 读取整个文件（也可以只读取第三列）
opts = detectImportOptions(filename);
T = readtable(filename, opts);

% 假设第三列是你想处理的列
textData = T{:,3};  % 第三列

% 初始化数组
N = length(textData);
Vmag = zeros(N,1);           % 速度大小
Vdir = zeros(N,3);           % 单位方向向量

for i = 1:N
    str = textData{i}; 
    
    % 使用正则表达式提取数值
    tokens = regexp(str, '([\d\.\-eE]+)\s*\*\s*\[([^\]]+)\]', 'tokens');
    
    if ~isempty(tokens)
        tokens = tokens{1};
        Vmag(i) = str2double(tokens{1});
        Vdir(i,:) = str2num(tokens{2}); %#ok<ST2NM> % 转换成向量
    else
        warning('第 %d 行格式不符合预期: %s', i, str);
        Vmag(i) = NaN;
        Vdir(i,:) = [NaN NaN NaN];
    end
end
