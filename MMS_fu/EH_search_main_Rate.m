%------written by Wending Fu, Nov.2023 in Beijing------------GSM坐标系
%------revised by Jiqiu Niu, Jan.2025 in Xi'an------------
%---【二轮-rate判据！！】---modified by Wending Fu, Mar.2025 in Beijing------------
%% GSE
clear; %清除工作空间中的所有变量
clc; %清空 MATLAB 命令窗口
global ParentDir OutputDir %全局变量
ParentDir = 'Z:/Data/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir); %形成临时目录路径

Date = '2017-05-28/2017-05-29';
splitDate = regexp(Date,'/','split');%基于斜杠 /，分割 Date 字符串，返回开始日期和结束日期
iic = 1:4;
ic = 1:4;%替换句
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');%返回日期内所有有的数据，命令行输出"文件目录爬取完成"
filenames2= SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames3 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');

filenames = [filenames1,filenames2,filenames3]; %文件名列表

% NameTags/FileGroups
expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,'EHsearch/',splitDate{1},'To',splitDate{2},'/'];%Z盘里新建文件夹记录xxxx-xx-xxToxxxx-xx-xx[起始终止日期]【DFsearch记得改名】
if ~isfolder([OutputDir,'OverviewFig/'])
    mkdir([OutputDir,'OverviewFig/']);
end
%分类和目录准备结束
units = irf_units; %%【常数统一】

for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
for i = 1:length(FileGroups{TDT})
    if ~isempty(strfind(FileGroups{TDT}{i},'fgm')) || ~isempty(strfind(FileGroups{TDT}{i},'edp')) ...
       || ~isempty(strfind(FileGroups{TDT}{i},'fpi')) 
        %检查当前文件名是否包含特定的字符串（如'fgm'、'des-moms'或'dis-moms'）。如果包含，则执行接下来的代码块。
    try
        SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 64, 'CheckSize', 0) %%% 替换句，不要直接把数据下到nas里，nas检索速度很慢
    catch 
        writematrix([NameTags{TDT}(2:end-2),'的Z盘数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
    end
    end
end
%下载文件到tempDir目录，并进行线程数设置、文件大小检查【后期可尝试删掉，影响运行速度】
SDCDataMove(TempDir,ParentDir)

%% get tint
mms.db_init('local_file_db',ParentDir);
formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
% 所得为字符串：两个NameTags元素的时间信息，并格式化为特定的日期时间格式
tempTint=irf.tint(tempDate);
% 字符串转换为时间间隔对象
try
    B1_ts=mms.get_data('B_gse_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));    
%     Pos = mms.get_data('R_gsm',tint);
%     Pos = Pos.gsmR1;
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    %将错误信息写入到OutputDir下的errorlog.txt文件中，使用追加模式。
    continue
end

%% Rate判据
flag = 0;
mu0 = units.mu0;ep0=units.eps0;
ic=1:4;
try
    c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
    c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
    c_eval('B_brst0?=irf.ts2mat(Bxyz?);',ic);
    c_eval('E_brst?=irf.ts2mat(Exyz?);',ic);
    Pos = mms.get_data('R_gse',tint);
    c_eval('Pos? = Pos.gseR?;',ic) %% xyz + total 四列
    c_eval('R?(1,2:4)=Pos?(1,1:3);',ic);
    c_eval('R?(1,1)=E_brst?(1,1);',ic); 
    c_eval('Ni?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',iic);
    c_eval('Ne?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',iic);
    c_eval('vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',iic);
    c_eval('ve?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',iic);
    c_eval('ni0?=irf.ts2mat(Ni?);',iic);
    c_eval('ne0?=irf.ts2mat(Ne?);',iic);
    c_eval('Vi0?=irf.ts2mat(vi?);',iic);
    c_eval('Ve0?=irf.ts2mat(ve?);',iic);
     %% 【dE/dT 】&【displacement current】
    dE_dt=[];
    c_eval('j_dis?=[];',ic)
    c_eval('j_v?=[];',iic)
    c_eval('err?=[];',iic)
    c_eval('err_rate?=[];',iic)
    c_eval('E_brst?=irf_resamp(E_brst?,E_brst1);',ic)
    for ic_idx = 1:length(ic)
        ic_val = ic(ic_idx);
        E_brst = eval(['E_brst' num2str(ic_val)]);  % 获取当前ic对应的E_brst数据
        dE_dt(:,1,ic_idx) = (E_brst(1:end-1,1) + E_brst(2:end,1)) / 2;% 计算时间中点 (t_i + t_{i+1})/2
        dt = diff(E_brst(:,1));  % 时间间隔
        dE_dt(:,2:4,ic_idx) = diff(E_brst(:,2:4)) ./ dt / 1e3;  % 计算电场分量的时间导数 
        eval(['dE_dt' num2str(ic_val) ' = dE_dt(:,:,ic_idx);']);
    end
    c_eval('j_dis?(:,1)=dE_dt?(:,1);',ic)
    c_eval('j_dis?(:,2:4)=dE_dt?(:,2:4).*mu0.*ep0.*1e12;',ic)%[T/m]*1e12=[nT/km]
    %% 【delta X B】
    c_eval('B_brst?=irf_resamp(B_brst0?,E_brst1);',ic)
    d_cros_B=c_4_grad('R?','B_brst?','curl'); % [nT/km]/1e9/1e3=[T/m]
    %% j=nqv mean
    c_eval('Ve?=irf_resamp(Ve0?,B_brst?);',iic);
    c_eval('Vi?=irf_resamp(Vi0?,B_brst?);',iic);
    c_eval('ni?=irf_resamp(ni0?,B_brst?);',iic);
    c_eval('ne?=irf_resamp(ne0?,B_brst?);',iic);
    c_eval('j_vi? = irf_multiply(1,[ni?(:,1) ni?(:,2)],1,[Vi?(:,1) Vi?(:,2:4)],1);',iic);
    c_eval('j_ve? = irf_multiply(1,[ne?(:,1) ne?(:,2)],1,[Ve?(:,1) Ve?(:,2:4)],1);',iic);
    c_eval('j_v?(:,1)=ni?(:,1);',iic);
    c_eval('j_v?(:,2:4) = (j_vi?(:,2:4)-j_ve?(:,2:4))*1.6/1e10;',iic);
    c_eval('j_v?(:,2:4)=j_v?(:,2:4)*mu0*1e12;',iic)%[T/m]*1e12=[nT/km]
    c_eval('j_v?(:,5)=vecnorm(j_v?(:,2:4),2,2);',iic)
    %% 【rate=位移电流/dxB and 夹角】【4颗卫星】
    c_eval('j_dis?=irf_resamp(j_dis?,d_cros_B);',iic)
    c_eval('j_dis?(:,5)=vecnorm(j_dis?(:,2:4),2,2);',iic)
    d_cros_B(:,5)=vecnorm(d_cros_B(:,2:4),2,2);
    c_eval('rate?=irf_multiply(1,[j_dis?(:,1) j_dis?(:,5)],1,[j_v?(:,1) j_v?(:,5)],-1);',iic)
    %% error
    c_eval('err_rate?(:,1)=j_v1(:,1);',iic)
    c_eval('err?(:,1)=j_v1(:,1);',iic)
    c_eval('err?(:,2:4)=j_v?(:,2:4)+j_dis?(:,2:4)-d_cros_B(:,2:4);',iic)
    c_eval('err?(:,5)=vecnorm(err?(:,2:4),2,2);',iic)
    c_eval('err_rate? = irf_multiply(1,[err?(:,1) err?(:,5)],1,[d_cros_B(:,1) d_cros_B(:,5)],-1);',iic);


    flagTime = []; flag = 0;validIdx=[];
%     rate_all = [rate1(:,2), rate2(:,2), rate3(:,2), rate4(:,2)];%替换句  % 合并所有rate变量（n×4矩阵）
    rate_all = [rate1(:,2), rate2(:,2), rate3(:,2)];%替换句 
%     err_rate_all = [err_rate1(:,2), err_rate2(:,2), err_rate3(:,2), err_rate4(:,2)];%替换句  % 合并所有error变量（n×4矩阵）
    err_rate_all = [err_rate1(:,2), err_rate2(:,2), err_rate3(:,2)];%替换句
    cond1 = any(rate_all > 0.25, 2);  % 条件1：至少一个rate > 0.25，% 按行判断（any(..., 2)）
    cond2 = any(err_rate_all < 0.3, 2);% 条件2：至少一个error < 0.3
    validIdx = find(cond1 & cond2);% 同时满足两个条件的索引

    if ~isempty(validIdx)
        flag=1;
        flagTime = [flagTime;rate1(validIdx, 1)];
        writematrix(['在',irf_time(rate1(validIdx(1), 1), 'epoch>utc'),'到',irf_time(rate1(validIdx(end), 1), 'epoch>utc'),...
            '期间找到了',num2str(length(validIdx)),'个jE'], [OutputDir, 'caselist.txt'], 'WriteMode', 'append', 'Encoding', 'UTF-8')
    end
    
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据检索出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
%% 符合判据的继续下载并出图 同1个时间段文件下载一次，逐个跑图
if flag == 1
try
    SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 32, 'CheckSize', 1)
    SDCDataMove(TempDir,ParentDir) %下载到tempDir的数据移动到ParentDir目录
    mms.db_init('local_file_db',ParentDir);
catch
     writematrix([NameTags{TDT}(2:end-2),'的画图-数据导入出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
    for ip = 1:length(flagTime)
        try
        PlotTint1 = irf_time([flagTime(ip)-0.02,flagTime(ip)+0.02],'epoch>epochTT'); %前后40ms时间
        PlotTint2 = irf_time([flagTime(ip)-0.1,flagTime(ip)+0.1],'epoch>epochTT');
        PlotTint3 = irf_time([flagTime(ip)-0.2,flagTime(ip)+0.2],'epoch>epochTT'); 
%         PlotTint2 = irf_time([flagTime(ip)-15,flagTime(ip)+15],'epoch>epochTT'); %前后30s时间
        % id_flagTime = SDCPlot(PlotTint,PlotTint2,ic,NameTags{TDT},flagTime(ip));需要返回值的话用这句
        SDCPlot(PlotTint1,PlotTint2,PlotTint3,ic,flagTime(ip),ip);
        catch
         writematrix([irf_time(flagTime(ip),'epoch>utc'),'的画图调用出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
        continue
        end    
    end
end
end