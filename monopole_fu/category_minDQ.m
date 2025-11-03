% scan_min_q.m
% 作用：遍历文件夹内所有 .mat 文件，读取变量 resQ（n×1 cell，部分为空），
%       对每个文件找出全局最小的 q，并把文件名、最小 q 和其所在 cell 的行号写入 txt。
% 使用方法：
%   1) 修改下方 folderPath 为你的 .mat 文件所在目录
%   2) 运行本脚本；输出将写到同目录下的 min_q_summary.txt

clear; clc;

%% ===== 需要修改：数据目录 =====
folderPath = 'D:\MMS\Monopole_Search\2020-01-21To2023-01-01\';

%% ===== 输出文件路径 =====
outTxt = fullfile(folderPath, 'min_q_summary.txt');

files = dir(fullfile(folderPath, '*.mat'));
[~, order] = sort({files.name});
files = files(order);


fid = fopen(outTxt, 'w');
if fid == -1
    error('无法创建输出文件：%s', outTxt);
end
fprintf(fid, 'filename\tmin_q\trow_index\n');

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    try
        S = load(fpath, 'resQ');

        if ~isfield(S, 'resQ')
            fprintf(fid, '%s\tNaN\t-1\n', files(k).name);
            warning('文件缺少变量 resQ：%s', files(k).name);
            continue;
        end
        resQ = S.resQ;

        if ~iscell(resQ)
            fprintf(fid, '%s\tNaN\t-1\n', files(k).name);
            warning('变量 resQ 不是 cell：%s', files(k).name);
            continue;
        end

        resQ = resQ(:);
        n = numel(resQ);

        cell_mins = inf(n, 1);

        for i = 1:n
            c = resQ{i};
            if isempty(c) || ~isnumeric(c)
                continue;
            end
            vals = c(:);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                cell_mins(i) = min(vals);
            end
        end

        [min_q, idxMin] = min(cell_mins);

        if isinf(min_q)
            fprintf(fid, '%s\tNaN\t-1\n', files(k).name);
            warning('文件中没有有效的 q 值：%s', files(k).name);
        else
            fprintf(fid, '%s\t%.9g\t%d\n', files(k).name, min_q, idxMin);
        end

    catch ME
        fprintf(fid, '%s\tNaN\t-1\n', files(k).name);
        warning('处理文件出错：%s\n原因：%s', files(k).name, ME.message);
    end
end

fclose(fid);