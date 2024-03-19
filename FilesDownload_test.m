
file_url = 'https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?file=mms1_fpi_brst_l2_dis-dist_20170611175223_v3.3.0.cdf';
local_file_path = '/Users/fwd/Documents/Ti~mor~/M/dipolarization fronts/Search Case/201706111755/mms1_fpi_brst_l2_dis-dist_20170611175223_v3.3.0.cdf';
file_name = 'mms1_fpi_brst_l2_dis-dist_20170611175223_v3.3.0.cdf';

% 构建调用 Python 脚本的命令字符串
command = sprintf('python3 /Users/fwd/Documents/PyExercise/download_files.py "%s" "%s" "%s"', file_url, local_file_path, file_name);

try
    % 调用系统命令执行 Python 脚本
    status = system(command);
    
    % 检查命令是否成功执行
    if status == 0
        disp('Python 脚本成功执行。');
    else
        disp('执行 Python 脚本时出现错误。');
    end
    
catch
    % 处理错误
    disp('执行 Python 脚本时出现未捕获的错误。');
end
