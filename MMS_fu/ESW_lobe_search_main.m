%------written by Wending Fu, July.2025 in Singapore------------
clear;clc;
global ParentDir OutputDir
ParentDir = '/Volumes/SPART-NAS/Data/MMS/'; 
DownloadDir = '/Volumes/SPART-NAS/Data/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

Date = '2015-09-01/2016-01-01';
splitDate = regexp(Date,'/','split');
ic = 1;iic = 1;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2, filenames3, filenames4];
% filenames = filenames1;

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
OutputDir = [ParentDir,'ESWSearch/',splitDate{1},'To',splitDate{2},'/'];
if ~isfolder([OutputDir,'OverviewFig/'])
    mkdir([OutputDir,'OverviewFig/']);
end

units = irf_units;
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'/'];
clc;fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])
SDCFilesDownload_NAS(FileGroups{TDT},tempDir)
SDCDataMove(tempDir,ParentDir)

formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
% tempDate = [char(NameTags(TDT)) '/' char(NameTags(TDT+1))];
tempTint=irf.tint(tempDate);

try
    B1_ts=mms.get_data('B_gsm_brst',tempTint,ic);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));    
    Pos = mms.get_data('R_gsm',tint);
    Pos = Pos.gsmR1;
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据导入出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end


%% lobe location
flag = 0;
if Pos(1,1) <= -10*units.RE/1e3 && abs(Pos(1,3)) >= 5*units.RE/1e3 % X<-10, |Z|>5
try
    % FLAG 1, discrete location
    Ne1_ts=mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint);
    Ne1=irf.ts2mat(Ne1_ts);
    
    if min(Ne1(:,2)) > 0.1, continue; else, flag = 1; end

    % FLAG 2, discrete E zero
    B1_ts=mms.get_data('B_gsm_brst',tint,ic);
    B1 = irf.ts2mat(B1_ts);

    c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic); 
    c_eval(['E?_ts = irf_gse2gsm(E?_ts);'],ic);
    c_eval('E? = irf.ts2mat(E?_ts);',ic);
    c_eval(['Efac=irf_convert_fac(E?,B?,[1,0,0]);'],ic);

    Ecri = Efac(1:end-1,4) .* Efac(2:end,4); 
    index_Ezero = find(Ecri <= 0);  % E过零点

    Efre = 1/median(diff(Efac(:,1)));
    index_Elen = index_Ezero + [-1 * Efre, 1 * Efre];
    index_Elen(index_Elen(:,1) < 1 , 1) = 1;
    index_Elen(index_Elen(:,2) > size(index_Elen,1) , 2) = size(index_Elen,1);

    for i_indexE = 1:size(index_Elen, 1)
        tempE = Efac(index_Elen(i_indexE , 1):index_Elen(i_indexE,2),4);
        if max(tempE) >= 10
            flag = 2;
            break
        else
            continue
        end
    end

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
else
    continue
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