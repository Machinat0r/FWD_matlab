%------written by Wending Fu, Nov.2023 in Beijing------------GSM坐标系
%------revised by Jiqiu Niu, Jan.2025 in Xi'an------------
clear; %清除工作空间中的所有变量
clc; %清空 MATLAB 命令窗口
global ParentDir OutputDir %全局变量
ParentDir = 'Z:\Data\MMS\'; 
% ParentDir = 'Z:/Data/MMS/'; 原版
DownloadDir = 'E:/fu/data2/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir); %形成临时目录路径

Date = '2015-09-03/2015-09-04';
splitDate = regexp(Date,'/','split');%基于斜杠 /，分割 Date 字符串，返回开始日期和结束日期
ic = 1:4;iic = 1;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');%返回日期内所有有的数据，命令行输出"文件目录爬取完成"
filenames2= SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3,filenames4]; %文件名列表
% filenames = filenames1;

% NameTags/FileGroups
expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
desiredLength2 = 11;
desiredLength = 17;
NameTags2 = NameTags(cellfun(@(x) strlength(x) == desiredLength2, NameTags));
NameTags = NameTags(cellfun(@(x) strlength(x) == desiredLength, NameTags));
%%
FileGroups2 = cell(1,length(NameTags2)); 
for j = 1:length(NameTags2)
    FileGroups2{j} = filenames(contains(filenames,NameTags2{j}));
end
FileGroups2 = cellfun(@cellstr,FileGroups2,'UniformOutput',false);%按时间分类整理后的文件名组
%%
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

units = irf_units; %%【！！】【常数统一】

%%
% parpool;
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'/'];
clc;fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
%遍历当前NameTags元素对应的FileGroups组中的文件名。
for i = 1:length(FileGroups{TDT})
    if ~isempty(strfind(FileGroups{TDT}{i},'fgm')) || ~isempty(strfind(FileGroups{TDT}{i},'edp')) || ...
        ~isempty(strfind(FileGroups{TDT}{i},'scm'))
        %检查当前文件名是否包含特定的字符串（如'fgm'、'des-moms'或'dis-moms'）。如果包含，则执行接下来的代码块。
        try
            SDCFilesDownload_NAS(FileGroups{TDT}(i),tempDir, 'Threads', 32, 'CheckSize', 1)
        catch 
            writematrix([NameTags{TDT}(2:end-2),'的Z盘数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
            continue
        end
        %下载文件到tempDir目录，并进行线程数设置、文件大小检查【后期可尝试删掉，影响运行速度】
    end
end
SDCDataMove(tempDir,ParentDir)

mms.db_init('local_file_db',ParentDir);
    tempDate = [NameTags{TDT}(2:5),'-',NameTags{TDT}(6:7),'-',NameTags{TDT}(8:9),'T',...
        NameTags{TDT}(10:11),':',NameTags{TDT}(12:13),':',NameTags{TDT}(14:15),'.000Z/',...
        NameTags{TDT+1}(2:5),'-',NameTags{TDT+1}(6:7),'-',NameTags{TDT+1}(8:9),'T',...
        NameTags{TDT+1}(10:11),':',NameTags{TDT+1}(12:13),':',NameTags{TDT+1}(14:15),'.000Z'];
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


%% 强变化电场判据
flag = 0;
mu0 = units.mu0;ep0=units.eps0;
ic=1:4;
c_eval('E_fac?=[];',ic);
try
    c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
%     try
    c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
%     catch
%             for i = 1:length(FileGroups{TDT})Bscm
%                 ########先删掉再下载，再改
%                 if  ~isempty(strfind(FileGroups{TDT}{i},'edp')) 
%         %检查当前文件名是否包含特定的字符串（如'fgm'、'des-moms'或'dis-moms'）。如果包含，则执行接下来的代码块。
%                     SDCFilesDownload_NAS(FileGroups{TDT}(i),tempDir, 'Threads', 32, 'CheckSize', 1)
%         %下载文件到tempDir目录，并进行线程数设置、文件大小检查【后期可尝试删掉，影响运行速度】
%                 end
%             end
%         c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
%         continue
%     end

%     c_eval('Bscm?=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
    c_eval('B_brst0?=irf.ts2mat(Bxyz?);',ic);
    c_eval('E_brst?=irf.ts2mat(Exyz?);',ic);
    ic=2:4; c_eval('E_brst?=irf_resamp(E_brst?,E_brst1);',ic);
    ic=1:4;
    c_eval('B_brst?=irf_resamp(B_brst0?,E_brst?);',ic);
    c_eval('B_brst?(:,5)=vecnorm(B_brst?(:,2:4),2,2);',ic);
%     c_eval('E_fac?(:,1)=B_brst?(:,1);',ic);
    
    c_eval('EB_2? = irf_multiply(1,[B_brst?(:,1) B_brst?(:,2)],1,[E_brst?(:,1) E_brst?(:,2)],1);',ic);
    c_eval('EB_2? = irf_multiply(1,[EB_2?(:,1) EB_2?(:,2)],1,[B_brst?(:,1) B_brst?(:,5)],-1);',ic);
    c_eval('EB_3? = irf_multiply(1,[B_brst?(:,1) B_brst?(:,3)],1,[E_brst?(:,1) E_brst?(:,3)],1);',ic);
    c_eval('EB_3? = irf_multiply(1,[EB_3?(:,1) EB_3?(:,2)],1,[B_brst?(:,1) B_brst?(:,5)],-1);',ic);
    c_eval('EB_4? = irf_multiply(1,[B_brst?(:,1) B_brst?(:,4)],1,[E_brst?(:,1) E_brst?(:,4)],1);',ic);
    c_eval('EB_4? = irf_multiply(1,[EB_4?(:,1) EB_4?(:,2)],1,[B_brst?(:,1) B_brst?(:,5)],-1);',ic);
    
    c_eval('E_fac?(:,1)=EB_2?(:,1);',ic);
    c_eval('E_fac?(:,2)=EB_2?(:,2)+EB_3?(:,2)+EB_4?(:,2);',ic);

% %     for i= 1:length(B_brst1(:,1)) %生成四项平行电场，一会儿加验证【循环运行慢，已改为irf_multiply】
% %         c_eval('E_fac?(i,1)=B_brst?(i,1);',ic);
%         c_eval('E_fac?(i,2)=dot(E_brst?(i,2:4),B_brst?(i,2:4))/B_brst?(i,5);',ic);
% %     end
    
    flagTime = []; tempB = 1; flag = 0;
    deltaTime = 33; % frequency 8192Hz * duration 4ms=4/1e3s，4ms窗口不断滑动
    disp(['□□□□□□□□□□','开始检索✧(≖ ◡ ≖✿)'])
 while tempB <= size(E_fac1,1) - deltaTime %如果没有用到B1，该数据名需修改；作用为检测4ms间隔
    c_eval('[~,tempidx?_1] = max(E_fac?(tempB:tempB+deltaTime,2));',ic);
    c_eval('[~,tempidx?_2] = min(E_fac?(tempB:tempB+deltaTime,2));',ic);
    tempmid=mean(tempidx1_1,tempidx1_2);
    tempwindow_1=tempmid-17; tempwindow_2=tempmid+17;
    if max(E_fac1(tempB:tempB+deltaTime,2))-min(E_fac1(tempB:tempB+deltaTime,2)) >= 100  ...
            && abs(tempidx1_1-tempidx1_2)<= 16 ...
            if max(E_fac2(tempB:tempB+deltaTime,2))-min(E_fac2(tempB:tempB+deltaTime,2)) >= 100  ...
            && abs(tempidx2_1-tempidx2_2)<= 16  ...
            && tempwindow_2 >= tempidx2_1>= tempwindow_1 && tempwindow_2 >= tempidx2_2>= tempwindow_1 ...
            && max(E_fac3(tempB:tempB+deltaTime,2))-min(E_fac3(tempB:tempB+deltaTime,2)) >= 100  ...
            && abs(tempidx3_1-tempidx3_2)<= 16  ...
            && tempwindow_2 >= tempidx3_1>= tempwindow_1 && tempwindow_2 >= tempidx3_2>= tempwindow_1 ...
            && max(E_fac4(tempB:tempB+deltaTime,2))-min(E_fac4(tempB:tempB+deltaTime,2)) >= 100  ...
            && abs(tempidx4_1-tempidx4_2)<= 16  ...
            && tempwindow_2 >= tempidx4_1>= tempwindow_1 && tempwindow_2 >= tempidx4_2>= tempwindow_1 

            flagTime(end+1) = E_fac1(tempB,1); %E_fac1当前的时间存储
            writematrix([irf_time(flagTime(end),'epoch>utc'),'找到了EH'],...
                [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')

            tempB = tempB + deltaTime;
            flag = 1;
            end
    end
    tempB = tempB + 1;
    clc; disp(['૮₍ ˃ ⤙ ˂ ₎ა',repmat('■',1,round(10*tempB/size(E_fac1,1))),repmat('□',1,10-round(10*tempB/size(E_fac1,1))),'正在光速检索ing...'])
end
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据检索出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
%% 符合判据的继续下载并出图 同1个时间段文件下载一次，逐个跑图
if flag == 1
try
    goal=[NameTags{TDT}(1:9),NameTags{TDT}(16:17)];
    mask = cellfun(@(x) ~isempty(strfind(x, goal)), NameTags2); %#ok<STREMP> 
    TDT2 = find(mask);
    SDCFilesDownload_NAS(FileGroups2{TDT2},tempDir, 'Threads', 32, 'CheckSize', 1)
    SDCFilesDownload_NAS(FileGroups{TDT},tempDir, 'Threads', 32, 'CheckSize', 1)
    SDCDataMove(tempDir,ParentDir) %下载到tempDir的数据移动到ParentDir目录
    mms.db_init('local_file_db',ParentDir);
catch
     writematrix([NameTags{TDT}(2:end-2),'的画图-数据导入出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
    for ip = 1:length(flagTime)
        try
        PlotTint = irf_time([flagTime(ip)-0.15,flagTime(ip)+0.15],'epoch>epochTT'); %前后0.3s时间
        PlotTint2 = irf_time([flagTime(ip)-15,flagTime(ip)+15],'epoch>epochTT'); %前后30s时间
        % id_flagTime = SDCPlot(PlotTint,PlotTint2,ic,NameTags{TDT},flagTime(ip));需要返回值的话用这句
        SDCPlot(PlotTint,PlotTint2,ic,NameTags{TDT},flagTime(ip),ip);
        catch
         writematrix([irf_time(flagTime(ip),'epoch>utc'),'的画图调用出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
        continue
        end    
    end
end
%% 删除文件夹并生成记录文件
try
    cd(OutputDir)
    rmdir(tempDir,'s');    
catch
    writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
end
end