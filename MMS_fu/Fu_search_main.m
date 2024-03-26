%main
%gsm
clear;
clc;
Date = '2019-05-01/2019-07-01';
% Date = '2018-09-16/2018-09-17';
splitDate = regexp(Date,'/','split');
ic = 1;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames3 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3];

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

global OutputDir ParentDir
ParentDir = 'E:\MMS\'; 
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,splitDate{1},'To',splitDate{2},'\'];
if ~isfolder([OutputDir,'OverviewFig\'])
    mkdir([OutputDir,'OverviewFig\']);
end

for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'\'];
fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
for i = 1:length(FileGroups{TDT}) %这里是为了在判断的时候不全下载，只下需要的数据
    if ~isempty(strfind(FileGroups{TDT}{i},'fgm')) || ~isempty(strfind(FileGroups{TDT}{i},'edp'))
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


%% By bipolr && Epara
flag = 0;
if Pos(1,1) > -6372*14 && abs(Pos(1,2)) < 6372*10 %GSM坐标系X向小于10Re
    try
        c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
        c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
        c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
        c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
        c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
        c_eval('E_res = irf_resamp(E?,B?);',ic);
        c_eval(['Efac?=irf_convert_fac(E_res,B?,[1,0,0]);'],ic);
    catch
        writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],...
            'WriteMode','append','Encoding','UTF-8')
        continue
    end
    
    flagTime = [];
    for dotB = 1:length(B1)-64
        if mean(abs(B1(dotB:dotB+64,4))) >  mean(abs(B1(dotB:dotB+64,3)))  && ...
                mean(abs(B1(dotB:dotB+64,4))) >  mean(abs(B1(dotB:dotB+64,2))) && ...
                mean(abs(Efac1(dotB:dotB+64,4))) > 5 
                
            fprintf(num2str(abs(mean(Efac1(dotB:dotB+64,2:4)))))
            flagTime(end+1) = B1(dotB,1);
            flag = 1;
        end
    end
    
end
%% 获得des-moms的文件名
for tempnum = 1:length(FileGroups{TDT})
    if strfind(FileGroups{TDT}{tempnum},'des-moms') > 0
        desmoms = ['E:\MMS\mms1\fpi\brst\l2\des-moms\',...
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
        try
            id_flagTime = SDCPlot(tint,desmoms,ic,NameTags{TDT},flagTime);
        catch
            writematrix([NameTags{TDT}(2:end-2),'画图出现问题，但该时间段内有事件'],[OutputDir,'errorlog.txt'],...
            'WriteMode','append','Encoding','UTF-8')
        end
    end
%% 删除文件夹并生成记录文件
try
    cd(OutputDir)
    rmdir(tempDir,'s');    
    if flag == 0
        fprintf([NameTags{TDT}(2:end-2),'中无事件\n'])
    else 
        fprintf('找到事件啦φ(≧ω≦*)♪\n')
    end
catch
    writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
end
end