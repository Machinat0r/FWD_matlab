clear;clc
close all

%% load data
cd G:\Cluster
ParentDir = 'G:\Cluster\';
ic=1:4;

% TT = '2002-03-17\2003-01-01';
TT = '2001-01-01\2004-06-01';
Datelist = regexp(TT,'\d+-\d+-\d+','match');
TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'\']; mkdir(TaskDir)
Datelist = datenum(Datelist,'yyyy-mm-dd');
Datelist = datestr(Datelist(1):Datelist(2),'yyyy-mm-dd');

mkdir([TaskDir 'meanPlot\'])
mkdir([TaskDir 'maxPlot\'])
% load monopole_example

for tempDate = 1:length(Datelist)-1
clc
Tsta = [Datelist(tempDate,:) 'T00:00:000Z'];
Tend = [Datelist(tempDate+1,:) 'T00:00:000Z'];
tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH
disp(['当前日期:' Datelist(tempDate,:)])

try 
    caa_load_changed_by_fwd('CL_SP_AUX',Tsta,Tend);
catch
%     caa_download(tint,'CL_SP_AUX');
    caa_load_changed_by_fwd('CL_SP_AUX',Tsta,Tend);
end

clear RR
R = c_caa_var_get('sc_r_xyz_gse__CL_SP_AUX','mat');
tint = [R(1,1) R(end,1)];
Tsta = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(1))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];
Tend = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(2))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];

% if min(sqrt(R(:,2).^2+R(:,3).^2+R(:,4).^2)) < 20000+6371
%     continue
% end

try
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
catch
    try
    %    Magnetic fields
%     c_eval("caa_download(tint,'C*_CP_FGM_FULL')",ic);
%     c_eval("caa_download(tint,'C*_CP_AUX_POSGSE_1M')",ic);  % position & velocity for each sc
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
    catch
        continue
    end
continue
end
% caa_download(tint,'CL_SP_AUX')% position,attitude.. for all sc
% caa_download(tint,'C2_CP_FGM_FULL');
% caa_download(tint,'C4_CP_FGM_FULL');

% caa_load C  %load data from datebase form C

%background magnetic field
try
% % % dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
% % % varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
c_eval('B?_gsm = irf_gse2gsm(B?_gse);',ic);
c_eval('B?=irf_abs(B?_gsm);',ic);
c_eval('B? = irf_resamp(B?,B1);',2:4);
c_eval('B?_gsm = irf_resamp(B?_gsm,B1);',2:4);
catch
        writematrix(['Data incompleted at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
end
%% PI
flag_ang = 0;
% PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
PI(abs(PI)<0.5) = 0;
%% solve monopole
if ~isempty(find(PI ~= 0,1))
    flag_ang = 1;
%         SDCFilesDownload(filenames_srvy(srvyIdx),tempDir);
%         SDCDataMove(tempDir,ParentDir); mms.db_init('local_file_db',ParentDir);
try
    c_eval('R?_gse = c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');',ic);
    c_eval('R? = irf_gse2gsm(R?_gse);',ic);
    c_eval('R? = irf_resamp(R?,B1);')
    CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
    c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

    LocPoint = zeros(length(PI),3)*nan;
    LocRes = cell(length(PI),1);
    Q = zeros(length(PI),1)*nan;
    resQ = cell(length(PI),1);
    
    Qerror = ones(length(PI),1)*1000;
    Locerror = ones(length(PI),1)*200;
    dLoc = ones(length(PI),15)*10;
    
    
    % div
    gradB=c_4_grad('R?','B?_gsm','grad');
    divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度

    PI_id = find(PI~=0)';
    if length(PI_id)>1e4
        writematrix([Datelist(tempDate,:),'需计算数据点超过1w'],[TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
    end
    for i = PI_id
    flag_m = 0;
    clc;
    disp(['当前日期:' Datelist(tempDate,:)])
    disp(['current calculate:',num2str(i),'/',num2str(length(PI))]);
    
    if isnan(divB(i)), continue;end
    
    MultiPower = ceil(max([log10(abs(R1(i,2:4))),log10(abs(R2(i,2:4))),log10(abs(R3(i,2:4))),log10(abs(R4(i,2:4)))]));
    if MultiPower > 3, continue; end
    
    RR_mean = zeros(1,4);
    for ii = 1:3 
    c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
        'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
    c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
    end
    RR_mean = RR_mean(4)/6;
    
    % solve
    [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?_gsm',i,i*sign(divB(i,2)),RR_mean,10);
    
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
                    [TaskDir,'case_mean_03.txt'],'WriteMode','append','Encoding','UTF-8')
        writematrix(['mean RR = ', num2str(mean(tempd),5),],...
            [TaskDir,'case_mean_03.txt'],'WriteMode','append','Encoding','UTF-8')
        case 2
        writematrix(['Flag 06 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
            [TaskDir,'case_mean_06.txt'],'WriteMode','append','Encoding','UTF-8')
        writematrix(['mean RR = ', num2str(mean(tempd),5),],...
            [TaskDir,'case_mean_06.txt'],'WriteMode','append','Encoding','UTF-8')
        case 3
        writematrix(['Flag 1 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
            [TaskDir,'case_mean_1.txt'],'WriteMode','append','Encoding','UTF-8')
        writematrix(['mean RR = ', num2str(mean(tempd),5),],...
            [TaskDir,'case_mean_1.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    
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
catch
    writematrix([Datelist(tempDate,:),'的数据下载或读取出现问题'],[TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue
end
%% continue
continue 
%% flag_ang
switch flag_ang
    case 0
%         continue
    case 1
        writematrix(['Flag1 find at: ',datestr(datenum(1970,1,1,0,0,0)+B1(find(PI~=0,1))/86400,'yyyymmdd HH:MM:SS.FFF')],...
                        [TaskDir,'caselog.txt'],'WriteMode','append','Encoding','UTF-8')
%         continue
    case 2
        writematrix(['Flag2 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flag2/86400,'yyyymmdd HH:MM:SS.FFF')],...
                        [TaskDir,'case1.txt'],'WriteMode','append','Encoding','UTF-8')
%         continue
    case 3
        writematrix(['Flag3 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flag3/86400,'yyyymmdd HH:MM:SS.FFF')],...
                        [TaskDir,'case05.txt'],'WriteMode','append','Encoding','UTF-8')
%         continue
end

%% flag_m
if flag_m ==1
    writematrix(['Flag 06 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
                        [TaskDir,'case_mean.txt'],'WriteMode','append','Encoding','UTF-8')
elseif flag_m == 2
    writematrix(['Flag 1 find at: ',datestr(datenum(1970,1,1,0,0,0)+time_flagm/86400,'yyyymmdd HH:MM:SS.FFF')],...
                        [TaskDir,'case_mean.txt'],'WriteMode','append','Encoding','UTF-8')
else
    continue
end

%% Plot
try
%% Init figure
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
c_eval("irf_plot([B?(:,1) 0*B?(:,5)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,5)) min(B2(:,5)) min(B3(:,5)) min(B4(:,5))])-10 ...
    max([max(B1(:,5)) max(B2(:,5)) max(B3(:,5)) max(B4(:,5))])+10]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
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
set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-10 ...
    max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+10]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
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
set(gca,'Ylim',[min([min(B1_gsm(:,3)) min(B2_gsm(:,3)) min(B3_gsm(:,3)) min(B4_gsm(:,3))])-10 ...
    max([max(B1_gsm(:,3)) max(B2_gsm(:,3)) max(B3_gsm(:,3)) max(B4_gsm(:,3))])+10]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
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
set(gca,'Ylim',[min([min(B1_gsm(:,4)) min(B2_gsm(:,4)) min(B3_gsm(:,4)) min(B4_gsm(:,4))])-10 ...
    max([max(B1_gsm(:,4)) max(B2_gsm(:,4)) max(B3_gsm(:,4)) max(B4_gsm(:,4))])+10]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
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
ylabel('divB [nT/km^2] ','fontsize',12);
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
%% Annotation
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
figname = [TaskDir 'meanPlot\' datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd')];
print(gcf, '-dpng', [figname '.png']);
clf
%%
if flag_m ~= 0
%% Index id
% % % tempidx_B1 = find(Locerror<=50 & Qerror<=1000, 1);
% % % %     [~,tempidx_B] = max(abs(divB(:,2)));
% % % c_eval('[~,tempidx_B?] = sort(abs(B?_gsm(:,1)-B1_gsm(tempidx_B1,1)));',2:4);
% % % c_eval('tempidx_B? = tempidx_B?(1);')
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
axis equal 
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname = [TaskDir 'meanPlot\' datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd')];
print(gcf, '-dpng', [figname '.png']);
clf
end

catch
    writematrix(['Plot Failed at: ',datestr(datenum(1970,1,1,0,0,0)+B1(tempidx_B1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
end
end
end