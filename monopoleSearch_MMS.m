%main
%gsm
clear;
clc;
%%
<<<<<<< HEAD
% Date = '2015-09-01/2017-04-30';
% Date = '2015-09-01/2023-04-30';
Date = '2018-06-16/2018-06-17';
=======
Date = '2015-09-01/2017-04-30';
% Date = '2017-04-30/2017-05-01';
>>>>>>> parent of 6233c67 (Ver23.9.11)
% Date = '2017-01-01/2021-01-01';
% Date = '2022-06-09/2022-07-02';

splitDate = regexp(Date,'/','split');
ic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames3 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %To get loaction
filenames = filenames1;

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
%The dir of "SDCFilesDownload" to "datamove" must be the ParentDir!
OutputDir = [ParentDir,splitDate{1},'To',splitDate{2},'/'];
if ~isfolder([OutputDir,'meanFig/'])
    mkdir([OutputDir,'meanFig/']);
end
%%
units = irf_units;
NameTags{end+1} = ['_' strrep(splitDate{2},'-','') '235959_v'];
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'/'];
clc
fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])

tempDate = [NameTags{TDT}(2:5),'-',NameTags{TDT}(6:7),'-',NameTags{TDT}(8:9),'T',...
    NameTags{TDT}(10:11),':',NameTags{TDT}(12:13),':',NameTags{TDT}(14:15),'.000Z/',...
    NameTags{TDT+1}(2:5),'-',NameTags{TDT+1}(6:7),'-',NameTags{TDT+1}(8:9),'T',...
    NameTags{TDT+1}(10:11),':',NameTags{TDT+1}(12:13),':',NameTags{TDT+1}(14:15),'.000Z'];
tempTint=irf.tint(tempDate);

%% Poincare Index  
% srvyIdx = find(contains(filenames_srvy,NameTags{TDT}(2:9))==1);
flag = 0;flag2 = 0;
if length(FileGroups{TDT}) == 4 
try
    SDCFilesDownload_NAS(FileGroups{TDT},tempDir);
    SDCDataMove(tempDir,ParentDir); mms.db_init('local_file_db',ParentDir);
    B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));  
    c_eval("B?_ts=mms.get_data('B_gsm_brst',tint,?);");
    c_eval('B?_gsm = irf.ts2mat(B?_ts);'); 
    c_eval('B? = irf_abs(B?_gsm);');
    c_eval('B? = irf_resamp(B?,B1);',2:4);

    PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
    PI(PI>=0.5) = 1;
    PI(PI<=-0.5) = -1;
    PI(abs(PI)<0.5) = 0;
%% solve monopole
    if ~isempty(find(PI ~= 0,1))
%         SDCFilesDownload(filenames_srvy(srvyIdx),tempDir);
%         SDCDataMove(tempDir,ParentDir); mms.db_init('local_file_db',ParentDir);
        Pos = mms.get_data('R_gsm',tint);
        c_eval('R? = Pos.gsmR?;')
        c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
        c_eval('R? = irf_resamp(R?,B1);')
        CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
        c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');
        
    LocPoint = zeros(length(PI),3)*nan;
    LocRes = cell(length(PI),1);
    Q = zeros(length(PI),1)*nan;
    resQ = cell(length(PI),1);
    
    Qerror = ones(length(PI),1)*1000;
    Locerror = ones(length(PI),1)*200;
    dLoc = ones(length(PI),15)*5;
    
    % div
    gradB=c_4_grad('R?','B?_gsm','grad');
    divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度

    PI_id = find(PI~=0)';
    for i = PI_id
    flag_m = 0;
    time_flagm = 0;
    clc;
    disp(['当前日期:' NameTags{TDT}(2:end-2)])
    disp(['current calculate:',num2str(i),'/',num2str(length(PI))]);
    
    MultiPower = ceil(max([log10(abs(R1(i,2:4))),log10(abs(R2(i,2:4))),log10(abs(R3(i,2:4))),log10(abs(R4(i,2:4)))]));
    
    if MultiPower > 3
        continue
    end
    
    RR_mean = zeros(1,4);
    for ii = 1:3 
    c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
        'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
    c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
    end
    RR_mean = RR_mean(4)/6;
    
    % solve
    [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?_gsm',i,i*sign(divB(i,2)),RR_mean,1);
    
    id = nchoosek(1:6,2);
    c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
    tempd = [];
    c_eval('tempd = [tempd,tempd?(4)/RR_mean];',1:15);
    dLoc(i,:) = tempd;
    
    
    if mean(tempd)<0.3
        flag_m = 1;
        time_flagm = B1(i,1);
        tempidx_B1 = i;
    elseif mean(tempd)<0.6
        flag_m = 2;
        time_flagm = B1(i,1);
        tempidx_B1 = i;
    elseif mean(tempd) < 1
        flag_m = 3;
        time_flagm = B1(i,1);
        tempidx_B1 = i;
    end
    
    switch flag_m
        case 1
        writematrix(['Flag 03 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
            [OutputDir,'case_mean_03.txt'],'WriteMode','append','Encoding','UTF-8')
        writematrix(['mean RR = ', num2str(mean(tempd),5),],...
            [OutputDir,'case_mean_03.txt'],'WriteMode','append','Encoding','UTF-8')
        case 2
        writematrix(['Flag 06 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
            [OutputDir,'case_mean_06.txt'],'WriteMode','append','Encoding','UTF-8')
        writematrix(['mean RR = ', num2str(mean(tempd),5),],...
            [OutputDir,'case_mean_06.txt'],'WriteMode','append','Encoding','UTF-8')
        case 3
        writematrix(['Flag 06 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
            [OutputDir,'case_mean_1.txt'],'WriteMode','append','Encoding','UTF-8')
        writematrix(['mean RR = ', num2str(mean(tempd),5),],...
            [OutputDir,'case_mean_1.txt'],'WriteMode','append','Encoding','UTF-8')
    end
<<<<<<< HEAD
%% Plot
if flag_m ~=0
try
    monopoleSearch_MMS_plot('tint',[OutputDir,'meanFig/'])
=======
    
    tof05 = tempd<0.5;tof1 = tempd<1;
    if sum(tof05) == 15
        flag_ang = 3;
        time_flag = B1(i,1);
%         tempidx_B1 = i;
    elseif sum(tof1) == 15
        flag_ang = 2;
        time_flag = B1(i,1);
%         tempidx_B1 = i;
    end
%%
    %error
% % %     Qerror(i) = 1000*std(resQ{i})/Q(i);
% % %     tri_a = delaunayTriangulation([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)]);
% % %     [~,volume_a] = convexHull(tri_a);
% % %     tri = delaunayTriangulation(LocRes{i});%%delaunay三角剖分
% % %         if size(tri.Points,1)==1
% % %             volume = 0;
% % %         else
% % %         [~,volume] = convexHull(tri);%%计算多面体体积
% % %         end
% % %     Locerror(i) = 100*volume/volume_a;
% % %     
% % %     if Locerror(i)<=50 && Qerror(i)<1000
% % %         flag_ang = 2;
% % %         time_flg = B1(i,1);
% % %     end
    end
    end
%% Calculate Error
% % %     Qerror = zeros(length(PI),1);
% % %     Locerror = zeros(length(PI),1);
% % % 
% % %     for i = 1:length(PI)
% % %     if ~isempty(resQ{i})
% % %         Qerror(i) = 100*std(resQ{i})/Q(i);
% % %     else
% % %         Qerror(i) = 200;
% % %     end
% % % 
% % %     if ~isempty(LocRes{i})
% % %     tri_a = delaunayTriangulation([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)]);
% % %     [~,volume_a] = convexHull(tri_a);
% % %     tri = delaunayTriangulation(LocRes{i});%%delaunay三角剖分
% % %         if size(tri.Points,1)==1
% % %             volume = 0;
% % %         else
% % %         [~,volume] = convexHull(tri);%%计算多面体体积
% % %         end
% % %     Locerror(i) = 100*volume/volume_a;
% % %     else
% % %     Locerror(i) = 200;
% % %     end
% % % 
% % %     if Locerror(i)<=50 && Qerror(i)<=100
% % %         flag_ang = 2;
% % %     end
% % %     end
>>>>>>> parent of 6233c67 (Ver23.9.11)
catch
    writematrix([NameTags{TDT}(2:end-2),'的画图出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
end
break
end

    end
    end
catch
    writematrix([NameTags{TDT}(2:end-2),'的数据下载或读取出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
end
%% Delete folder
try
    cd(OutputDir)
    rmdir(tempDir,'s');
    fclose all;
catch
    fprintf(['删除文件夹',NameTags{TDT}(2:end-2),'失败\n'])
end
%% continue
continue 
end

%% Plot
switch flag
case 1
writematrix(['mean 1 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flag/86400,'yyyymmdd HH:MM:SS.FFF')],...
                [OutputDir,'meancase.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
case 2
writematrix(['mean 06 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flag/86400,'yyyymmdd HH:MM:SS.FFF')],...
                [OutputDir,'meancase.txt'],'WriteMode','append','Encoding','UTF-8')
fprintf('找到事件啦φ(≧ω≦*)♪\n')
try                
%% FOTE err & div
gradB=c_4_grad('R?','B?_gsm','grad');
divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度
%% Init figure 1
n=9;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Btotal
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2(:,1) B2(:,5)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3(:,1) B3(:,5)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,5)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,5)) min(B2(:,5)) min(B3(:,5)) min(B4(:,5))])-5 ...
    max([max(B1(:,5)) max(B2(:,5)) max(B3(:,5)) max(B4(:,5))])+5]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('|B| [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bx
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gsm(:,1) B1_gsm(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gsm(:,1) B2_gsm(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gsm(:,1) B3_gsm(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gsm(:,1) B4_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
    max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bx [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% By
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gsm(:,1) B1_gsm(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gsm(:,1) B2_gsm(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gsm(:,1) B3_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gsm(:,1) B4_gsm(:,3)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,3)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,3))])-100 max([max(B?_gsm(:,3))])+100]);",ic);
set(gca,'Ylim',[min([min(B1_gsm(:,3)) min(B2_gsm(:,3)) min(B3_gsm(:,3)) min(B4_gsm(:,3))])-5 ...
    max([max(B1_gsm(:,3)) max(B2_gsm(:,3)) max(B3_gsm(:,3)) max(B4_gsm(:,3))])+5]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('By [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bz
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gsm(:,1) B1_gsm(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gsm(:,1) B2_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gsm(:,1) B3_gsm(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gsm(:,1) B4_gsm(:,4)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,4))])-100 max([max(B?_gsm(:,4))])+100]);",ic);
set(gca,'Ylim',[min([min(B1_gsm(:,4)) min(B2_gsm(:,4)) min(B3_gsm(:,4)) min(B4_gsm(:,4))])-5 ...
    max([max(B1_gsm(:,4)) max(B2_gsm(:,4)) max(B3_gsm(:,4)) max(B4_gsm(:,4))])+5]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bz [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% η
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([err_4C(:,1) err_4C(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,4))])-3 max([max(B?_gsm(:,4))])+3]);",ic);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % c_eval("set(gca,'Ylim',[min(err_4C(:,2))-0.1 max(err_4C(:,2))+0.1]);",ic);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('η','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% ξ
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([eigVal_err_v2(:,1) eigVal_err_v2(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(eigVal_err_v2(:,2))-0.1 max(eigVal_err_v2(:,2))+0.1]);",ic);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('ξ','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% divB
h(i)=irf_subplot(n,1,-i);
irf_plot([divB(:,1) divB(:,2)], 'color','k', 'Linewidth',0.75); hold on;

% c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
% set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('divB [nT/km^2]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% PI
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) PI], 'color','k', 'Linewidth',0.75); hold on;

c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('PI ','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% monopole_index
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) monopole_index], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
% % % set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('monopole index ','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Location Error
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) Locerror], 'color','k', 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca,'Ylim',[0 max(Locerror)]);
% % % set(gca,'Ylim',[0,200], 'ytick',[0 100 200],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Loc Err [%]','fontsize',12);
% % % i=i+1;
%% dLoc Error
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) max(dLoc,[],2)], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Locerror)]);
set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('Loc Err','fontsize',12);
i=i+1;
%% meandLoc Error
h(i)=irf_subplot(n,1,-i);

irf_plot([B1(:,1) mean(dLoc,2)], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Locerror)]);
set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('mean','fontsize',12);
i=i+1;
%% Q Error
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) Qerror], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Qerror)]);
set(gca,'Ylim',[0,1000], 'ytick',[0 500 1000],'fontsize',9);
pos1=get(gca,'pos');
ylabel('Q Err [%]','fontsize',12);
i=i+1;
%% Adjust the position
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
figname = [OutputDir,'meanFig\',NameTags{TDT}(2:end-2)];    
print(gcf, '-dpng', [figname '.png']);
clf

%% id
[~,tempidx_R] = sort(abs(R1(:,1)-B1(tempidx_B1,1)));
tempidx_R = tempidx_R(1);
%% Init Figure 2
figure(2)
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Coordinate
c_eval('R?(:,2:4) = R?(:,2:4)-LocPoint(tempidx_B1,:);');
LocRes{tempidx_B1} = LocRes{tempidx_B1}-LocPoint(tempidx_B1,:);
LocPoint = LocPoint - LocPoint(tempidx_B1,:);
%% Location 
cor = 'krgb';
c_eval("plot3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),'s' ,'color',cor(?),'linewidth',5);hold on;")
%% Tetrahedron configuration
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(tempidx_R,2),R',num2str(ii),'(tempidx_R,3),R',num2str(ii),'(tempidx_R,4);',...
    'R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
plot3(RR12(:,1),RR12(:,2),RR12(:,3),'--k');hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'--k');hold on;  
plot3(RR14(:,1),RR14(:,2),RR14(:,3),'--k');hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'--k');hold on;  
plot3(RR34(:,1),RR34(:,2),RR34(:,3),'--k');hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'--k');hold on;  
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
xlabel('e_1 [km]','fontsize',12);
ylabel('e_2 [km]','fontsize',12);
zlabel('e_3 [km]','fontsize',12);

%% Quiver
c_eval('B?_gsm = irf_abs(B?_gsm);');
maxB = max([B1_gsm(tempidx_B1,5),B2_gsm(tempidx_B1,5),B3_gsm(tempidx_B1,5),B4_gsm(tempidx_B1,5)]);
c_eval("quiver3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),RR_mean*B?_gsm(tempidx_B1,2)/maxB,RR_mean*B?_gsm(tempidx_B1,3)/maxB,RR_mean*B?_gsm(tempidx_B1,4)/maxB,'color',cor(?));hold on;")

%% Loc res
plotPolyhedron(LocRes{tempidx_B1}(:,1),LocRes{tempidx_B1}(:,2),LocRes{tempidx_B1}(:,3),'#FFB8CE',0.3);
% idx = [166:198];
% plot3(LocPoint(idx,1),LocPoint(idx,2),LocPoint(idx,3),'*','color','#FFBAF1');
% line(LocPoint(173:195,1),LocPoint(173:195,2),LocPoint(173:195,3),'color','#FFBAF1');
plot3(LocPoint(tempidx_B1,1),LocPoint(tempidx_B1,2),LocPoint(tempidx_B1,3),'*','color','m');
plot3(LocRes{tempidx_B1}(:,1),LocRes{tempidx_B1}(:,2),LocRes{tempidx_B1}(:,3),'*','color','#FFB8CE');
%%
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname = [OutputDir,'meanFig\',NameTags{TDT}(2:end-2),'-Configuration'];
print(gcf, '-dpng', [figname '.png']);
clf
catch
writematrix([NameTags{TDT}(2:end-2),'画图出现问题，但该时间段内有事件'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
end
case 0
    fprintf([NameTags{TDT}(2:end-2),'中无事件\n'])
end
end