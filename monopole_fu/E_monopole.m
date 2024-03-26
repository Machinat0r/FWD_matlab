clear;clc

global ParentDir 
ParentDir = 'E:\MMS\'; 
TempDir = 'E:\MMS\temp\';mkdir(TempDir);
% TT = '2017-08-07T17:02:43.300Z/2017-08-07T17:02:43.800Z';
% TT = '2017-08-07T16:37:17.690Z/2017-08-07T16:37:17.820Z';
% TT = '2017-08-04T09:01:08.070Z/2017-08-04T09:01:08.140Z';
% TT = '2016-09-27T01:19:37.489Z/2016-09-27T01:19:37.529Z';
TT = '2016-10-19T18:36:11.100Z/2016-10-19T18:36:11.700Z';

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
        if i > length(NameTags), break; end
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
%% FOTE 误差
gradE=c_4_grad('R?_gsm','E?_gsm','grad');
eigVal_err_v2=E1(:,1);

for ii=1:length(E1(:,1))  
deltE_null=reshape(gradE(ii,2:end),3,3);
[V,D] = eig(deltE_null);
% Figure 1o    以最大特征值归一化 
% 百分比
% eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;  
% 小数形式
eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]);  
end
    

[j,divE,~,jxB,divTshear,divPb] = c_4_j('R?_gsm','E?_gsm');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divE,1,jmag,-1);          %% Figure 1n    受背景电流影响大
% err_4C(:,2)=abs(err_4C(:,2))*100;   
err_4C(:,2)=abs(err_4C(:,2)); 

%% 以上两个实际上都不是很合理，建议按照如下计算
c_eval('E? = irf_resamp(E?,E1);',2:4);
E_mean=mean([E1(:,5) E2(:,5) E3(:,5) E4(:,5)],2);

divE=[gradE(:,1) sum([gradE(:,2) gradE(:,6) gradE(:,10)],2)];      %% 未归一化散度
units = irf_units;
divE(:,2) = divE(:,2)*units.eps0;
divE(:,3)=divE(:,2)./E_mean;  

%% 卷积
load monopole_example
% conv_E = E3(2168:2470,2);
% load tripolar_example
c= E1(:,2)/max(E1(:,2));
[~,~,dist] = findsignal(c,conv_E,'Timealignment','dtw','Normalization','zscore')
findsignal(c,conv_E,'Timealignment','dtw','Normalization','zscore')

%% Poincare Index
indices=c_fgm_poincare_index(E1(:,2:4),E2(:,2:4),E3(:,2:4),E4(:,2:4));
% indices=c_4_poincare_index(E1(:,2:4),E2(:,2:4),E3(:,2:4),E4(:,2:4));
indices(abs(indices)<0.5) = 0;
indices(indices>=0.5) = 1;
indices(indices<=-0.5) = -1;
%% Init figure
n=7;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Ex
h(i)=irf_subplot(n,1,-i);
irf_plot([E1_gsm(:,1) E1_gsm(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([E2_gsm(:,1) E2_gsm(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([E3_gsm(:,1) E3_gsm(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([E4_gsm(:,1) E4_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([E?_gsm(:,1) 0*E?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[min([min(E1_gsm(:,2)) min(E2_gsm(:,2)) min(E3_gsm(:,2)) min(E4_gsm(:,2))])-10 ...
    max([max(E1_gsm(:,2)) max(E2_gsm(:,2)) max(E3_gsm(:,2)) max(E4_gsm(:,2))])+10]);
% c_eval("set(gca,'Ylim',[min([min(E?_gsm(:,2))])-100 max([max(E?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ex [mV/m]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Ey
h(i)=irf_subplot(n,1,-i);
irf_plot([E1_gsm(:,1) E1_gsm(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([E2_gsm(:,1) E2_gsm(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([E3_gsm(:,1) E3_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([E4_gsm(:,1) E4_gsm(:,3)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([E?_gsm(:,1) 0*E?_gsm(:,3)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(E?_gsm(:,3))])-100 max([max(E?_gsm(:,3))])+100]);",ic);
set(gca,'Ylim',[min([min(E1_gsm(:,3)) min(E2_gsm(:,3)) min(E3_gsm(:,3)) min(E4_gsm(:,3))])-10 ...
    max([max(E1_gsm(:,3)) max(E2_gsm(:,3)) max(E3_gsm(:,3)) max(E4_gsm(:,3))])+10]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ey [mV/m]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Ez
h(i)=irf_subplot(n,1,-i);
irf_plot([E1_gsm(:,1) E1_gsm(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([E2_gsm(:,1) E2_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([E3_gsm(:,1) E3_gsm(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([E4_gsm(:,1) E4_gsm(:,4)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([E?_gsm(:,1) 0*E?_gsm(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(E?_gsm(:,4))])-100 max([max(E?_gsm(:,4))])+100]);",ic);
set(gca,'Ylim',[min([min(E1_gsm(:,4)) min(E2_gsm(:,4)) min(E3_gsm(:,4)) min(E4_gsm(:,4))])-10 ...
    max([max(E1_gsm(:,4)) max(E2_gsm(:,4)) max(E3_gsm(:,4)) max(E4_gsm(:,4))])+10]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ez [mV/m]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% η
h(i)=irf_subplot(n,1,-i);
irf_plot([err_4C(:,1) err_4C(:,2)], 'color','k', 'Linewidth',0.75); hold on;

% c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(err_4C(:,2))]) max([max(err_4C(:,2))])]);",ic);
% set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('η','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% ξ
h(i)=irf_subplot(n,1,-i);
irf_plot([eigVal_err_v2(:,1) eigVal_err_v2(:,2)], 'color','k', 'Linewidth',0.75); hold on;

% c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(eigVal_err_v2(:,2))]) max([max(eigVal_err_v2(:,2))])]);",ic);
% set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('ξ','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;

%% divE
h(i)=irf_subplot(n,1,-i);
irf_plot([divE(:,1) divE(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([divE(:,1) 0*divE(:,2)], 'k--', 'Linewidth',0.75); hold on;
% c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(divE(:,2))]) max([max(divE(:,2))])]);",ic);
% set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('\rho ','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% normalized divE
% % h(i)=irf_subplot(n,1,-i);
% % irf_plot([divE(:,1) divE(:,3)], 'color','k', 'Linewidth',0.75); hold on;
% % irf_plot([divE(:,1) 0*divE(:,3)], 'k--', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % grid off;
% % c_eval("set(gca,'Ylim',[min([min(divE(:,3))]) max([max(divE(:,3))])]);",ic);
% % % set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
% % pos1=get(gca,'pos');
% % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % ylabel('\rho/|E| ','fontsize',12);
% % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % i=i+1;
%% PI
h(i)=irf_subplot(n,1,-i);
irf_plot([E1(:,1) indices], 'Linewidth',0.75, 'color','k', 'Linewidth',0.75); hold on;
grid off;
pos1=get(gca,'pos');
ylabel('PI','fontsize',10);
set(gca, 'Ylim',[-1 1], 'ytick',[-1 0 1],'fontsize',9);
% set(gca,'xtick',[])
i = i+1;
%%%
%% conv
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([E2(:,1) c], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([E2(:,1) 0*c], 'k--', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(c) max(c)]);",ic);
% % % % set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('c ','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Adjust the position
% irf_adjust_panel_position;
irf_zoom(tint,'x',h(1:n));
irf_plot_axis_align(h)

%% Annotation
% irf_zoom(tint,'x',h)
% set(h(1:end-1),'XTickLabe','');
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;

% leg = 'abcdefghi';
% for ii=1:i_subplot-1
%     set(h(ii),'ColorOrder',[0 0 0]);
%     irf_legend(h(ii),{['(' leg(ii) ')']},[0.02, 0.05]);
% end

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
% figname = ['E:\Cluster\' 'xmp'];
% print(gcf, '-dpng', [figname '.png']);    