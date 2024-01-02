%------written by Wending Fu, Nov.2023 in Beijing------------
clear;clc;
Date = '2019-08-05/2019-08-06';
splitDate = regexp(Date,'/','split');
ic = 1;iic = 1;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3,filenames4];
% filenames = filenames1;

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

global OutputDir ParentDir
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,'DFSearch/',splitDate{1},'To',splitDate{2},'/'];
if ~isfolder([OutputDir,'OverviewFig/'])
    mkdir([OutputDir,'OverviewFig/']);
end

units = irf_units;

for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'/'];
clc;fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
for i = 1:length(FileGroups{TDT})
    if ~isempty(strfind(FileGroups{TDT}{i},'fgm')) || ~isempty(strfind(FileGroups{TDT}{i},'des-moms'))...
            || ~isempty(strfind(FileGroups{TDT}{i},'dis-moms'))
        SDCFilesDownload_NAS(FileGroups{TDT}(i),tempDir)
    end
end
SDCDataMove(tempDir,ParentDir)

mms.db_init('local_file_db',ParentDir);
tempDate = [NameTags{TDT}(2:5),'-',NameTags{TDT}(6:7),'-',NameTags{TDT}(8:9),'T',...
    NameTags{TDT}(10:11),':',NameTags{TDT}(12:13),':',NameTags{TDT}(14:15),'.000Z/',...
    NameTags{TDT+1}(2:5),'-',NameTags{TDT+1}(6:7),'-',NameTags{TDT+1}(8:9),'T',...
    NameTags{TDT+1}(10:11),':',NameTags{TDT+1}(12:13),':',NameTags{TDT+1}(14:15),'.000Z'];
tempTint=irf.tint(tempDate);

try
    B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));    
    Pos = mms.get_data('R_gsm',tint);
    Pos = Pos.gsmR1;
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end


%% 偶极化锋面判据
flag = 0;
if Pos(1,1) <= -10*units.RE/1e3 && abs(Pos(1,2)) <= 12*units.RE/1e3 %GSM坐标系X向小于10Re
try
    B1_ts=mms.get_data('B_gsm_brst',tint,1);
    Bt1_ts=B1_ts.abs; Bt1=irf.ts2mat(Bt1_ts);
    B1 = irf.ts2mat(B1_ts);

    Ne1_ts=mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint);
    Ne1=irf.ts2mat(Ne1_ts);
    Ni1_ts=mms.db_get_ts('mms1_fpi_brst_l2_dis-moms','mms1_dis_numberdensity_brst',tint);
    Ni1=irf.ts2mat(Ni1_ts);
    c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',ic);
    c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],ic);
    c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',ic);
    c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],ic);
    c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],ic);

    c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',tint);',ic);
    c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
    c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',tint);',ic);
    c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
    c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);
    
    mu0 = units.mu0; kB = units.kB;
    c_eval('Ni?_res=irf_resamp(Ni?,Bt?);',ic)
    c_eval('Ne?_res=irf_resamp(Ne?,Bt?);',ic)
    c_eval('Te?_res=irf_resamp(Te?,Bt?);',ic)
    c_eval('Ti?_res=irf_resamp(Ti?,Bt?);',ic)
    c_eval('Pb?=[Bt?(:,1) ((Bt?(:,2).^2))/(2*mu0)*1e-9];',ic);%nPa
    c_eval('Pti? = irf_multiply(11604.505*kB*1e6*1e9,[Ni?_res(:,1) Ni?_res(:,2)],1,[Ti?_res(:,1) Ti?_res(:,2)],1);',ic);
    c_eval('Pte? = irf_multiply(11604.505*kB*1e6*1e9,[Ne?_res(:,1) Ne?_res(:,2)],1,[Te?_res(:,1) Te?_res(:,2)],1);',ic);
    c_eval('Pt? = [Pti?(:,1) Pti?(:,2)+Pte?(:,2)];',ic);
    c_eval('beta? = [Bt?(:,1) Pt?(:,2)./Pb?(:,2)];',ic);

    c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); 
    c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
    c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
    c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);


    flagTime = []; tempB = 1; flag = 0;
    deltaTime = 128*5; % frequency 128Hz * duration 5s
    disp(['□□□□□□□□□□','开始检索✧(≖ ◡ ≖✿)'])
while tempB <= size(B1,1) - deltaTime %如果没有用到B1，该数据名需修改；检测5s间隔
    NeTime = irf_time(Ne1_ts.time,'epochTT>epoch');
    ViTime = irf_time(Vi1_ts.time,'epochTT>epoch');
    [~,tempNeTime1] = min(abs(NeTime-B1(tempB,1))); [~,tempNeTime2] = min(abs(NeTime-B1(tempB+deltaTime,1)));
    [~,tempViTime1] = min(abs(ViTime-B1(tempB,1))); [~,tempViTime2] = min(abs(ViTime-B1(tempB+deltaTime,1)));

    if abs(B1(tempB+deltaTime,4)) - abs(B1(tempB,4)) >= 5 && max(abs(B1(tempB:tempB+deltaTime,4))./Bt1(tempB:tempB+deltaTime,2)) >= 0.5 ...
            && Ne1(tempNeTime1,2) >= Ne1(tempNeTime2,2) ...
            && mean(Vit1(tempViTime1:tempViTime2,2)) >= 100 ...
            && min(beta1(tempB:tempB+deltaTime,2)) >= 0.5

            flagTime(end+1) = B1(tempB,1);
            writematrix([irf_time(flagTime(end),'epoch>utc'),'找到了DF'],...
                [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')

            tempB = tempB + deltaTime;
            flag = 1;
    end
    tempB = tempB + 1;
    clc; disp(['૮₍ ˃ ⤙ ˂ ₎ა',repmat('■',1,round(10*tempB/size(B1,1))),repmat('□',1,10-round(10*tempB/size(B1,1))),'正在光速检索ing...'])
end
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end
end
%% 获得des-moms的文件名
for tempnum = 1:length(FileGroups{TDT})
    if strfind(FileGroups{TDT}{tempnum},'des-moms') > 0
        desmoms = [ParentDir,'mms1/fpi/brst/l2/des-moms/',...
            NameTags{TDT}(2:5),'/',NameTags{TDT}(6:7),'/',...
            NameTags{TDT}(8:9),'/',FileGroups{TDT}{tempnum}];
    end
end
%% 符合判据的继续下载并出图
if flag == 1
try
    SDCFilesDownload(FileGroups{TDT},tempDir) 
    SDCDataMove(tempDir,ParentDir)
    mms.db_init('local_file_db',ParentDir);
    PlotTint = irf_time([flagTime(end)-20,flagTime(end)+20],'epoch>epochTT');
    id_flagTime = SDCPlot(PlotTint,desmoms,ic,NameTags{TDT},flagTime(end));
catch
    writematrix([irf_time(flagTime(end),'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
    continue
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