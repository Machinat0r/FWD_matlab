close all
clear;clc

global ParentDir 
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);


TT = '2017-05-22T10:42:12.00Z/2017-05-22T10:42:14.00Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
filenames_fast = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms');
filenames = [filenames1,filenames2,filenames3,filenames4];

[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

% %SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 32, 'CheckSize', 0)
% % % SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 32, 'CheckSize', 0)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load P_a
c_eval('Pa_e?_ts = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_prestensor_gse_brst'',tint);',iic);
c_eval('Pa_i?_ts = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_prestensor_gse_brst'',tint);',iic);

% 将每个 3x3 矩阵展平成行向量,转换为double
c_eval('Pa_e? = double([Pa_e?_ts.time.epochUnix, reshape(Pa_e?_ts.data, size(Pa_e?_ts.data,1), [] )]);', iic);
c_eval('Pa_i? = double([Pa_i?_ts.time.epochUnix, reshape(Pa_i?_ts.data, size(Pa_i?_ts.data,1), [] )]);', iic);

%% 计算电子标量压力和偏应力张量
for idx = iic
    eval(sprintf('[Pa_e_sclar%d, Pa_e_deviatoric%d] = process_pressuretensor(Pa_e%d);', idx, idx, idx));
    eval(sprintf('res_e%d = Pa_e%d(:, 2:end) - Pa_e_deviatoric%d(:, 2:end);', idx, idx, idx));
    eval(sprintf('[Pa_i_sclar%d, Pa_i_deviatoric%d] = process_pressuretensor(Pa_i%d);', idx, idx, idx));
    eval(sprintf('res_i%d = Pa_i%d(:, 2:end) - Pa_i_deviatoric%d(:, 2:end);', idx, idx, idx));
end

function [Pa_sclar, Pa_deviatoric] = process_pressuretensor(Pa)
    % 提取时间戳和数据
    time = Pa(:, 1);
    data = Pa(:, 2:end);
    
    % 初始化处理后的数据矩阵
    processed_data1 = zeros(size(data));
    
    % 逐行处理每一行的 3x3 矩阵，计算标量压力
    for i = 1:size(data, 1)
        % 将展平的当前行数据恢复为 3x3 矩阵
        mat = reshape(data(i, :), [3, 3]);
        
        % 对对角线元素除以 3
        mat(1, 1) = mat(1, 1) / 3;
        mat(2, 2) = mat(2, 2) / 3;
        mat(3, 3) = mat(3, 3) / 3;
        
        % 将处理后的矩阵重新展平成行向量并存储
        processed_data1(i, :) = reshape(mat, [1, 9]);
    end
    
    % 将时间戳和处理后的数据部分重新拼接为原始大小
    Pa_sclar = [time, processed_data1];

    % 初始化存储进一步处理后的数据的矩阵
    processed_data2 = zeros(size(processed_data1));
    
    delta = eye(3);
    
    % 遍历每一行进行处理，计算偏应力张量
    for i = 1:size(processed_data1, 1)
        % 将展平的当前行数据恢复为 3x3 矩阵
        mat = reshape(processed_data1(i, :), [3, 3]);
        
        % 将非对角线元素置为零
        mat = mat .* delta; 
        
        % 将处理后的矩阵重新展平成行向量并存储
        processed_data2(i, :) = reshape(mat, [1, 9]);
    end
    
    % 将时间戳和处理后的数据部分重新拼接为原始大小
    Pa_deviatoric = [time, processed_data2];
end



