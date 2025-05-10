% gsm
clear;close all
clc
global ParentDir 
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
%% load data
ic=1:4;
hh = 20;
hh1= 20;
thresold = 0.4;  
% Tsta = '2017-07-06T01:44:50.000Z';
% Tend = '2017-07-06T01:45:10.000Z';

Tsta = '2017-07-12T11:54:33.600Z';
Tend = '2017-07-12T11:54:35.400Z';
% Tsta = '2018-07-06T12:37:04.500Z';
% Tend = '2018-07-06T12:37:05.250Z';
% Tsta = '2017-05-22T08:23:52.000Z';
% Tend = '2017-05-22T08:23:55.000Z';
% Tsta = '2017-05-25T04:40:29.500Z';
% Tend = '2017-05-25T04:40:33.000Z';
% Tsta = '2017-07-06T17:31:55.000Z';
% Tend = '2017-07-06T17:32:05.000Z';
tint = irf.tint(Tsta,Tend);
% tint=irf.tint('2017-01-27T12:05:42.50Z/2017-01-27T12:05:43.80Z');
TT = '2017-07-12T11:54:00.000Z/2017-07-12T11:55:00.000Z';
% % % TT = '2017-07-06T17:46:30.000Z/2017-07-06T17:48:00.000Z';
% TT = '2018-07-06T12:36:30.000Z/2018-07-06T12:38:00.000Z';
% TT = '2017-05-22T08:23:30.000Z/2017-05-22T08:24:30.000Z';
% TT = '2017-05-25T04:40:00.000Z/2017-05-25T04:41:00.000Z';

tint2=irf.tint(TT);

Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames = [filenames1, filenames2,filenames4];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 48, 'CheckSize', 0)
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db','D:/MMS/')

%% load B Fields data
c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
c_eval('B? = irf_resamp(B?, B1);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
c_eval('Bt? = irf.ts2mat(Bt?);')
%% smooth
c_eval('BS?(:,2)=smooth(B?(:,2),hh);',ic);
c_eval('BS?(:,3)=smooth(B?(:,3),hh);',ic);
c_eval('BS?(:,4)=smooth(B?(:,4),hh);',ic);
c_eval('BS?(:,1)=smooth(B?(:,1),hh);',ic);
c_eval('BS?=irf_abs(BS?);',ic);

c_eval('Ne?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval('ne?=irf.ts2mat(Ne?);',ic);

% c_eval('Ni?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
% c_eval('ni?=irf.ts2mat(Ni?);',ic);
c_eval('ni? = ne?;',ic)
c_eval('Vi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); 
% c_eval('Vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval('vi?=irf.ts2mat(Vi?);',ic);
c_eval('vi? = irf_gse2gsm(vi?);',ic);
%% Te
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);%GSE坐标下的B
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);%GSM坐标下的B
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)%GSE坐标下的Te
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)%FAC坐标下的Te
%GSE坐标系下的等离子体温度张量数据旋转到了磁力线坐标系（FAC，Field Aligned Coordinates）
% c_eval('matParTe? = facTe?.xx.resample(gsmB?.time).data;',ic)
% %FAC坐标系下等离子体温度张量的平行（Parallel）分量，.xx是张量的平行分量。
% c_eval('matPerTe?= (facTe?.yy.resample(gsmB?.time).data + facTe?.zz.resample(gsmB?.time).data)/2;',ic)
% %FAC坐标系下等离子体温度张量的垂直（Perpendicular）分量，取yy和zz分量的平均值计算垂直分量。
%温度先不重采样
c_eval('matParTe?(:,2) = facTe?.xx.data;',ic)
%FAC坐标系下等离子体温度张量的平行（Parallel）分量，.xx是张量的平行分量。
c_eval('matPerTe?(:,2)= (facTe?.yy.data + facTe?.zz.data)/2;',ic)
c_eval('matParTe?(:,1)= ne?(:,1);',ic)
c_eval('matPerTe?(:,1)= ne?(:,1);',ic)

c_eval('matParTe? = irf_resamp(matParTe?, matParTe1);',ic);
c_eval('matPerTe? = irf_resamp(matPerTe?, matPerTe1);',ic);
%% load Ve Fields data
c_eval('Vexyz?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval('Vexyz? = irf_gse2gsm(Vexyz?);',ic);
% L=[0.05 -0.21 0.98]; % Vi
% M=[-0.23 0.95 0.22];
% N=[-0.97 -0.24 -0.01];

c_eval('Ve?=irf.ts2mat(Vexyz?);',ic);
c_eval('Ve?=irf_resamp(Ve?,Ve1);',ic);
c_eval('ne? = irf_resamp(ne?, Ve1);',ic);

c_eval('Rgsm?=mms.get_data(''R_gsm'',tint2,?);',ic);
c_eval('R?=irf.ts2mat(Rgsm?);',ic);
c_eval('R? = irf_resamp(R?, B?);',ic);

c_eval('e_r? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energy_brst'',tint);',ic);
%% smooth
% smooth_step1
kk = length (e_r1.data);
if mod(kk,2) == 1
   le = kk-2; lo = kk-3;%even偶数（这里当它是奇数），odd奇数
else
   le = kk-3; lo = kk-2;%even偶数，odd奇数
end
for ic =1:4
    c_eval('e_r = e_r?;',ic)
    if e_r.data(1)>e_r.data(2)
        for ii = 1:2:le
            c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('Ve?(ii+1,2:4)=(Ve?(ii+2,2:4)+Ve?(ii,2:4))/2;',ic)
        end
    else
        for ii = 2:2:lo
            c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('Ve?(ii+1,2:4)=(Ve?(ii+2,2:4)+Ve?(ii,2:4))/2;',ic)
        end
    end
end
ic = 1:4;
c_eval('Ve?=irf_resamp(Ve?,B?);',ic);
c_eval('ne?=irf_resamp(ne?,B?);',ic);

c_eval('ne?(:,2)=smooth(ne?(:,2),hh1);')
c_eval('Ve?(:,2)=smooth(Ve?(:,2),hh1);')
c_eval('Ve?(:,3)=smooth(Ve?(:,3),hh1);')
c_eval('Ve?(:,4)=smooth(Ve?(:,4),hh1);')

c_eval('ne?(:,2)=smooth(ne?(:,2),hh1);')
c_eval('Ve?(:,2)=smooth(Ve?(:,2),hh1);')
c_eval('Ve?(:,3)=smooth(Ve?(:,3),hh1);')
c_eval('Ve?(:,4)=smooth(Ve?(:,4),hh1);')

c_eval('ne?(:,2)=smooth(ne?(:,2),hh1);')
c_eval('Ve?(:,2)=smooth(Ve?(:,2),hh1);')
c_eval('Ve?(:,3)=smooth(Ve?(:,3),hh1);')
c_eval('Ve?(:,4)=smooth(Ve?(:,4),hh1);')

c_eval('ne?(:,2)=smooth(ne?(:,2),hh1);')
c_eval('Ve?(:,2)=smooth(Ve?(:,2),hh1);')
c_eval('Ve?(:,3)=smooth(Ve?(:,3),hh1);')
c_eval('Ve?(:,4)=smooth(Ve?(:,4),hh1);')
%% calculate
% c_eval('j?(:,1)=Ve?(:,1);',ic);
% c_eval('j?(:,2:4)=(ni?(:,2).*vi?(:,2:4)-ne?(:,2).*Ve?(:,2:4))*1.6/1e10;',ic);

gradVe=c_4_grad('R?','Ve?','grad');
%find null position
for ii=1:length(Ve1(:,1))
    dBeach=reshape(gradVe(ii,2:end),3,3);
    dR1(ii,2:4)=Ve1(ii,2:4)*inv(dBeach');
    dR2(ii,2:4)=Ve2(ii,2:4)*inv(dBeach');
    dR3(ii,2:4)=Ve3(ii,2:4)*inv(dBeach');
    dR4(ii,2:4)=Ve4(ii,2:4)*inv(dBeach');
end
dR1(:,1)=Ve1(:,1);
dR2(:,1)=Ve1(:,1);
dR3(:,1)=Ve1(:,1);
dR4(:,1)=Ve1(:,1);

dRmag1=irf_abs(dR1);
dRmag2=irf_abs(dR2);
dRmag3=irf_abs(dR3);
dRmag4=irf_abs(dR4);

% dRmin(:,2)=min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,2) = dRmag1(:,5);
dRmin(:,1)=dRmag1(:,1);
dRmin(dRmin(:,2)>=1000,2) = nan;

dxVey = gradVe(:,3);
dyVex = gradVe(:,5);
c_eval('dR?_2D = sqrt((Ve?(:,3)./dxVey).^2 + (Ve?(:,2)./dyVex).^2);')
% c_eval('dR?_2D(dR?_2D >= 100) = nan;')

d_cros_Ve=c_4_grad('R?','Ve?','curl');
units = irf_units;
c_eval('delta_B?(:,1)=Ve?(:,1);',ic);
% c_eval('delta_B?(:,2) = 1e9* -units.mu0 * units.e .* 1e6 .* ne?(:,2) .* d_cros_Ve(:,4) .* 1e6 .* dR?_2D.^2;')% dBz
c_eval('delta_B?(:,2) = 1e9* -units.mu0 * units.e .* 1e6 .* ne?(:,2) .* d_cros_Ve(:,4) .* 1e6 .* dRmin(:,2).^2;')% dBz


[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','Ve?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

delta_B1(err_4C(:,2)>40,2) = nan;
%% Init figure
% ic = 1;
n=4;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% delta x Ve x plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([d_cros_Ve(:,1) 0*d_cros_Ve(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;

set(gca,'Ylim',[0 10])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
ylabel('B [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% delta x Ve x plot
h(i)=irf_subplot(n,1,-i);
d_cros_Ve = irf_abs(d_cros_Ve);
c_eval("irf_plot([d_cros_Ve(:,1) d_cros_Ve(:,5)], 'color','k', 'Linewidth',0.75);",ic); hold on;

%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([d_cros_Ve(:,1) 0*d_cros_Ve(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;

% set(gca,'Ylim',[-5 80])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
irf_legend(gca,{'J_x','J_y','J_z'},[0.97 0.92]);
ylabel('J [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% delta x Ve z plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([delta_B1(:,1) delta_B1(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;

%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([d_cros_Ve(:,1) 0*d_cros_Ve(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;

set(gca,'Ylim',[-10 10])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
set(gca,'xtick',[])
irf_legend(gca,{'J_x','J_y','J_z'},[0.97 0.92]);
ylabel('/delta B [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%%
h(i)=irf_subplot(n,1,-i);
hold(h(i),'on');
a1=irf_plot(h(i),[err_4C(:,1) err_4C(:,2)],'k');%散度/旋度
fill(h(i), [a1.XData,fliplr(a1.XData)],[zeros(1,length(a1.YData)),fliplr(a1.YData)],[0.8706 0.9216 0.9804]);%浅灰色
irf_plot(h(i),[err_4C(:,1) err_4C(:,2)*0+40], 'k--'); hold off;
set(gca,'Ylim',[0 100])
% a1=irf_plot(h(9),[err_4C(:,1) err_4C(:,2)],'k','Linewidth',lnwid);
% a2=irf_plot(h(9), eigVal_err_v2, 'c','Linewidth',lnwid); 
% 
% % fill(h1(i_plot),[a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
% fill(h(9),[a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
% fill(h(9), [a1.XData,fliplr(a1.XData)],[zeros(1,length(a1.YData)),fliplr(a1.YData)],[0.7 0.7 0.7]);


% irf_plot(h1(i_plot), [err_4C(:,1) err_4C(:,2)*0+40], 'k--','Linewidth',lnwid); 
% irf_plot(h(9),[eigVal_err(:,1) eigVal_err(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;

grid(h(i),'off'); 
% set(h1(i_plot),'Ylim',[0 160]);
%ylabel('\fontcolor{k}|\nabla\cdot\bf{B}|/|\nabla\times\bf{B}|,\fontcolor{b}|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}| [%]');
ylabel(h(i),{'\alpha',' (%)'},'Interpreter','tex');
set(h(i),'Ylim',[0 130],'Ytick',[0 40 100]);
irf_legend(h(i),'(j)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h(i),{'|\nabla\cdot\bf{nVe}|/|\nabla\times\bf{nVe}|'},[0.85 0.9],'color','k','fontsize',10)
