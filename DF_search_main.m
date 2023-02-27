%main
%gsm
clear;
clc;
Date = '2016-07-01/2016-08-30';
splitDate = regexp(Date,'/','split');
ic = 1;
filenames1 = SDCFilenames(Date,ic,'inst','fgm,edp','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3];
% filenames = filenames1;

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags));
for j = 1:length(NameTags)
    temp = 0;
    for  i = 1:length(filenames)
        if strfind(filenames{i},NameTags{j}) > 0
            temp = temp + 1;
            FileGroups{j}{1,temp} = filenames{i};
        end
    end
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);

global OutputDir ParentDir
ParentDir = 'D:\MMS\'; 
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,splitDate{1},'To',splitDate{2},'\'];
if ~isfolder([OutputDir,'OverviewFig\'])
    mkdir([OutputDir,'OverviewFig\']);
end

for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'\'];
fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
for i = 1:length(FileGroups{TDT})
    if ~isempty(strfind(FileGroups{TDT}{i},'fgm')) || ~isempty(strfind(FileGroups{TDT}{i},'des-moms'))
        SDCFilesDownload(FileGroups{TDT}(i),tempDir)
    end
end
for j = 1:length(filenames_srvy)
    if ~isempty(strfind(filenames_srvy{j},NameTags{TDT}(2:9)))
        SDCFilesDownload(filenames_srvy(j),tempDir)
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
if Pos(1,1) < -63720 && abs(Pos(1,2)) < 95580 %GSM坐标系X向小于10Re
    try
        B1_ts=mms.get_data('B_gsm_brst',tint,1);
        Bt1_ts=B1_ts.abs;
        Bt1=irf.ts2mat(Bt1_ts);
        B1=irf.ts2mat(B1_ts);
        Ne1_ts=mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint);
        Ne1=irf.ts2mat(Ne1_ts);
        energy_e1 = mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_energyspectr_omni_brst',tint);
        e1_time = struct('t',irf_time(energy_e1.DEPEND_0.data,'ttns>epoch'));
    catch
        writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],...
            'WriteMode','append','Encoding','UTF-8')
        continue
    end

    deltadot = round((e1_time.t(1,1)-Ne1(1,1))/0.03); %一般会前后少检验最多15个点左右，最多共约1s
    flagTime = [];e_enhance = [];
for dotB = 1:length(B1) - 200 %如果没有用到B1，该数据名需修改；200个磁场点,共1.56s
for dotNe = abs(deltadot)+1:length(Ne1) - 52 - abs(deltadot)  %52个电子点，共1.56s 
    if B1(dotB,1) - Ne1(dotNe,1) <= 0.03 && ...
            abs(B1(dotB+200,4) - B1(dotB,4)) >= 5 && abs(max(B1(dotB:dotB+200,4))) > abs(0.5*max(Bt1(dotB:dotB+200,2))) ...
            && abs(mean(B1(dotB:dotB+200,4))) > 0.5*abs(mean(Bt1(dotB:dotB+200,2))) && Ne1(dotNe,2) - Ne1(dotNe+52,2) >= 0.2
        %这组判据寻找的是每个时间中通量最大的能道
% % %         tempidx1 = find(energy_e1.data(dotNe+deltadot,:) == max(energy_e1.data(dotNe+deltadot,4:end)));
% % %         tempidx2 = find(energy_e1.data(dotNe+deltadot+52,:) == max(energy_e1.data(dotNe+deltadot+52,4:end)));
        %这组判据是寻找每个时间中通量变化最大的能道
% % % % %         tempidx1 = find(diff(energy_e1.data(dotNe+deltadot,4:end)) == max(diff(energy_e1.data(dotNe+deltadot,4:end)))) +4;
% % % % %         tempidx2 = find(diff(energy_e1.data(dotNe+deltadot+52,4:end)) == max(diff(energy_e1.data(dotNe+deltadot+52,4:end)))) +4;
% % % % % 
% % % % %         tempchid1 = energy_e1.DEPEND_1.data(dotNe+deltadot,tempidx1);
% % % % %         tempchid2 = energy_e1.DEPEND_1.data(dotNe+deltadot+52,tempidx2);
        %这组判据是寻找每个时间中的通量加权平均能量
% % %         tempidx1 = energy_e1.data(dotNe+deltadot,4:end).*energy_e1.DEPEND_1.data(dotNe+deltadot,4:end);
% % %         tempidx2 = energy_e1.data(dotNe+deltadot+52,4:end).*energy_e1.DEPEND_1.data(dotNe+deltadot+52,4:end);
% % %             
% % %         tempchid1 = sum(tempidx1)/sum(energy_e1.data(dotNe+deltadot,4:end));
% % %         tempchid2 = sum(tempidx2)/sum(energy_e1.data(dotNe+deltadot+52,4:end));
        
% % %         if tempchid2/tempchid1 >= 2.5
% % %             e_enhance(end+1) = tempchid2/tempchid1;
            flagTime(end+1) = B1(dotB,1);
            flag = 1;
            break %只跳出一层循环，继续寻找该时间段是否有其他符合判据的时间点
% % %         end   
    end
end
end
end
%% 获得des-moms的文件名
for tempnum = 1:length(FileGroups{TDT})
    if strfind(FileGroups{TDT}{tempnum},'des-moms') > 0
        desmoms = ['D:\MMS\mms1\fpi\brst\l2\des-moms\',...
            NameTags{TDT}(2:5),'\',NameTags{TDT}(6:7),'\',...
            NameTags{TDT}(8:9),'\',FileGroups{TDT}{tempnum}];
    end
end
%% 符合判据的继续下载并出图
    if flag == 1
        SDCFilesDownload(FileGroups{TDT},tempDir) 
        SDCDataMove(tempDir,ParentDir)
        mms.db_init('local_file_db',ParentDir);
        Units=irf_units;
        me=Units.me;
%         irf.tint(B1(dotB,1))
        id_flagTime = SDCPlot(tint,desmoms,ic,NameTags{TDT},flagTime);
    end
%% 删除文件夹并生成记录文件
try
    cd(OutputDir)
    rmdir(tempDir,'s');    
    if flag == 0
        fprintf([NameTags{TDT}(2:end-2),'中无事件\n'])
    else 
        fprintf('找到事件啦φ(≧ω≦*)♪\n')
        writematrix([NameTags{TDT}(2:end-2),'找到了符合条件的事件，电子通量增强了',num2str(e_enhance(id_flagTime)),'倍'],...
            [OutputDir,'caselist.txt'],'WriteMode','append','Encoding','UTF-8')
    end
catch
    writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
end
end