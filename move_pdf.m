% 主文件夹路径，包含88个子文件夹
source_folder = '/Users/fwd/Downloads/MIIT-papers-2023/input/files';  % 修改为你的主文件夹路径
% 目标文件夹路径
destination_folder = '/Users/fwd/Downloads/MIIT-papers-2023/output';  % 修改为你的目标文件夹路径

% 如果目标文件夹不存在，则创建该文件夹
if ~exist(destination_folder, 'dir')
    mkdir(destination_folder);
end

% 获取主文件夹下的所有子文件夹信息
subfolders = dir(source_folder);
subfolders = subfolders([subfolders.isdir]); % 只保留文件夹

% 遍历每个子文件夹
for i = 1:length(subfolders)
    % 跳过 "." 和 ".." 这两个系统目录
    if strcmp(subfolders(i).name, '.') || strcmp(subfolders(i).name, '..')
        continue;
    end
    
    % 当前子文件夹路径
    current_subfolder = fullfile(source_folder, subfolders(i).name);
    
    % 获取子文件夹中的所有PDF文件
    pdf_files = dir(fullfile(current_subfolder, '*.pdf'));
    
    % 将PDF文件移动到目标文件夹
    for j = 1:length(pdf_files)
        % 当前PDF文件的完整路径
        source_file = fullfile(current_subfolder, pdf_files(j).name);
        % 目标文件路径
        destination_file = fullfile(destination_folder, pdf_files(j).name);
        
        % 移动文件
        copyfile(source_file, destination_file);
        
        % 输出信息，显示已移动的文件
        fprintf('Moved: %s to %s\n', source_file, destination_file);
    end
end

disp('所有文件已成功移动。');
