clear;clc
clear;clc

global ParentDir 
ParentDir = 'E:\MMS\'; 
TempDir = 'E:\MMS\temp\';mkdir(TempDir);
TT = '2017-08-07T17:02:43.300Z/2017-08-07T17:02:43.800Z';
% TT = '2017-08-07T16:37:17.690Z/2017-08-07T16:37:17.820Z';
% TT = '2017-08-04T09:01:08.070Z/2017-08-04T09:01:08.140Z';
% TT = '2016-09-27T01:19:37.489Z/2016-09-27T01:19:37.529Z';
% TT = '2016-10-19T18:36:11.100Z/2016-10-19T18:36:11.700Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:4;
iic = 1:4;
filenames = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
% filenames = [filenames1,filenames2,filenames3,filenames4];

expr1 = '_\d+\_v';
NameTags = regexp(filenames,expr1,'match');
NameTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NameTags)),'UniformOutput',false);

TTlist = strjoin(regexp(TT,'\d+','match'),'');
i = 1;flag = 0;%若flag=0，说明整段时间都在第i-1个Tag里
while str2double(TTlist(18:31)) > NameTags{i}  % 如果时间段刚好仅在某天的最后一个文件里会出bug，可以把时间往前调1ms
    if str2double(TTlist(1:14)) < NameTags{i}
        flag=1; break  %若flag=1，说明时间段的开始在第i-1个Tag里，结束在第i个里
    else
        i=i+1;
        if i > length(NameTags) break; end
    end
end

if flag == 0
    tempTag = num2str(NameTags{i-1});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'\fpi\brst\l2\des-moms\',tempTag(1:4),'\',tempTag(5:6),'\',...
            tempTag(7:8),'\',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
    desmoms1 = desmoms; desmoms2 = desmoms1;
else
    if i == 1
        errordlg('时间起始处无brst数据，请检查时间范围,或使用Overview_srvydownload程序')
    end
    tempTag1 = num2str(NameTags{i-1});
    tempTag2 = num2str(NameTags{i});
    filenames1 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag1)));
    filenames2 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag2)));
    filenames = [filenames1,filenames2];
    desmoms1 = [ParentDir,'mms',num2str(ic),'\fpi\brst\l2\des-moms\',tempTag1(1:4),'\',tempTag1(5:6),'\',...
            tempTag1(7:8),'\',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'\fpi\brst\l2\des-moms\',tempTag2(1:4),'\',tempTag2(5:6),'\',...
            tempTag2(7:8),'\',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
end
SDCFilesDownload(filenames,TempDir)
SDCFilesDownload(filenames_srvy(1:2:7),TempDir)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);


% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_fast_l2'',tint,?);'],ic);
c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
c_eval('dfE? =1/median(diff(E?_gsm.time.epochUnix));',ic);
c_eval('Ebf? = E?_gsm.filt(0,1024,dfE?,3);',ic);
c_eval(['E?_gsm=irf.ts2mat(Ebf?);'],ic);

% c_eval(['E?_gsm=irf.ts2mat(E?_gsm);'],ic);
c_eval('E?=irf_abs(E?_gsm);',ic);

Pos = mms.get_data('R_gsm',tint);
c_eval('R?_gsm = Pos.gsmR?;',ic);
%% 卷积
load monopole_example
% conv_E = E2(475:545,2);
c= conv(E2(:,4),conv_E,'same');

findchangepts(c)

corr_E = ones(length(E3)-71,1);
for tt = 1:length(corr_E)
    temp = corrcoef(E2(tt:tt+70,4),conv_E);
    corr_E(tt) = temp(1,2);
end
%% plot
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 200; ySize = 60; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

subplot(2,2,1)
plot(E3(:,1),E3(:,2))

subplot(2,2,2)
plot(conv_E)

subplot(2,2,3)
plot(E3(:,1),c)

subplot(2,2,4)
plot(corr_E)