clc
clear
%% B_data
Tsta='2017-03-19T19:32:00.000Z';
Tend='2017-03-19T20:25:00.000Z';
% Tmia='/';
% Tints=[Tsta,Tmia,Tend];%投掷角
% 
tint1=irf.tint(Tsta,Tend);
tint = irf_time(tint1,'epochtt>epoch');
tintS=tint1;

% mvu_mag_ss highres 
dataBR=importdata("/Users/fwd/Documents/MATLAB/Code/wcq/20170319_Rolling-pin/overview/fgm_jno_l3_2017078ss_r1s_v01 - 1.txt");
Tsmag='2017-03-19T00:00:01.715';
Temag='2017-03-20T00:00:02.820';
T1=iso2epoch(Tsmag);
T2=iso2epoch(Temag);
TT = linspace(T1,T2,length(dataBR(:,1)))';

%% 磁场 insitu(4s)\mst(1s)\MAG(30ms) MSO
B_JSO=[TT,dataBR(:,8),dataBR(:,9),dataBR(:,10)];% nT，磁场的三个分量BX BY BZ
B=irf_tlim(B_JSO,tint);
%% |B|
Bt=[];
Bt(:,2)=sqrt(B(:,2).^2+B(:,3).^2+B(:,4).^2);%用三个磁场分量的平方和开平方计算磁场的绝对值大小
Bt(:,1)=B(:,1);
%% PA_data
load( "/Users/fwd/Documents/MATLAB/Code/wcq/20170319_Rolling-pin/overview/JED_data_180.mat");
load ("/Users/fwd/Documents/MATLAB/Code/wcq/20170319_Rolling-pin/overview/JED_data_270.mat");
%% JED180
time = JED_data_180.t;
energy180 = JED_data_180.energy;
pitch_angle = JED_data_180.pa;
flux = JED_data_180.flux;
%% Eedata
c_eval('flux?=flux(:,:,?)',1:24);
c_eval('Ee_flux_?=(flux?(:,1)+flux?(:,2)+flux?(:,3)+flux?(:,4)+flux?(:,5)+flux?(:,6))/6',1:24);
Time1=iso2epoch('2017-03-19T19:30:00Z'); 
Time2=iso2epoch('2017-03-19T20:35:00Z'); 
[~, idx] = min(abs(time - Time1)); % 找到距离a最近的元素位置
[~, idx1] =min(abs(time - Time2)); % 找到距离b最近的元素位置
Ee_t=time(idx:idx1,1);
Ee_flux=zeros(length(Ee_t),24);
% % Ee_flux(:,1)=time(:,1);
c_eval('Ee_flux(:,?)=Ee_flux_?(idx:idx1,1)',1:24);
Ee_f=energy180(1,5:24);
%% smooth
n=floor(length(Ee_t)/30);
for i = 1:n
    idx = (i-1)*30+1 : i*30;
    p180_123(i,:) = nanmean(Ee_flux(idx,:));
    t_123(i,1) = mean(Ee_t(idx,:));         
end
p180need=zeros(length(p180_123),20); 
for i=1:20
   p180need(:,i)=p180_123(:,i+4);
end

Ee = struct('t',t_123,'f',Ee_f,'p',p180need);
specEe=struct('t',Ee.t);
specEe.p=Ee.p;
specEe.f=Ee.f;
specEe.f_label='';
specEe.p_label={' ','1/s/cm^2/sr/kev'};
% cal pa spectra
flux_e0 = (flux(:,:,6)+flux(:,:,7)+flux(:,:,8)+flux(:,:,9)+flux(:,:,10)+flux(:,:,11)+flux(:,:,12)+flux(:,:,13)+flux(:,:,14)+flux(:,:,15)+flux(:,:,16)+flux(:,:,17)+flux(:,:,18)+flux(:,:,19)+flux(:,:,20)+flux(:,:,21)+flux(:,:,22)+flux(:,:,23)+flux(:,:,24))/19;
% flux_e0 = (flux(:,:,16)+flux(:,:,17))/2;
flux_e5 = flux_e0;
flux_et = flux_e0;
[pitch_angle2, ind] = sort(pitch_angle, 2, 'ascend');
ag_int = 10;
for ii = 1: length(pitch_angle2)
    flux_et(ii,:) = flux_e0(ii,ind(ii,:));
    pa_tmp = pitch_angle2(ii,:);
    for ia = 1:18
        an1 = (ia-1)*ag_int; an2 = ia*ag_int;
        id_aa = find((an1 < pa_tmp) & (pa_tmp < an2));
        flux_e5(ii,ia) = nanmean(flux_et(ii,id_aa));
    end
end

specepad2=struct('t',time);
specepad2.p = flux_e5;
specepad2.p_label={' ',' cm-2 s-1 sr-1 keV-1 '};
specepad2.f_label={''};
specepad2.f = [5:10:175];

%% JED270
time = JED_data_270.t;
energy270 = JED_data_270.energy;
pitch_angle = JED_data_270.pa;
flux = JED_data_270.flux;

% cal pa spectra
flux_e0 =(flux(:,:,6)+flux(:,:,7)+flux(:,:,8)+flux(:,:,9)+flux(:,:,10)+flux(:,:,11)+flux(:,:,12)+flux(:,:,13)+flux(:,:,14)+flux(:,:,15)+flux(:,:,16)+flux(:,:,17)+flux(:,:,18)+flux(:,:,19)+flux(:,:,20)+flux(:,:,21)+flux(:,:,22)+flux(:,:,23)+flux(:,:,24))/19;
% flux_e0 = (flux(:,:,16)+flux(:,:,17))/2;
flux_e5 = flux_e0;
flux_et = flux_e0;
[pitch_angle2, ind] = sort(pitch_angle, 2, 'ascend');
ag_int = 10;
for ii = 1: length(pitch_angle2)
    flux_et(ii,:) = flux_e0(ii,ind(ii,:));
    pa_tmp = pitch_angle2(ii,:);
    for ia = 1:18
        an1 = (ia-1)*ag_int; an2 = ia*ag_int;
        id_aa = find((an1 < pa_tmp) & (pa_tmp < an2));
        flux_e5(ii,ia) = nanmean(flux_et(ii,id_aa));
    end
end

specepad3=struct('t',time);
specepad3.p = flux_e5;
specepad3.p_label={' ',' cm-2 s-1 sr-1 keV-1 '};
specepad3.f_label={''};
specepad3.f = [5:10:175];

%% combine
% Time1=iso2epoch('2022-01-11T19:40:00Z'); 
% Time2=iso2epoch('2022-01-07T23:59:00Z'); 
Tsta='2017-03-19T19:30:00Z'; 
Tend='2017-03-19T20:35:00Z'; 
tint_cut = [iso2epoch(Tsta) iso2epoch(Tend)];

p180 = specepad2.p;p270 = specepad3.p; 
t180 = specepad2.t;t270 = specepad3.t;

f180 = irf_tlim([t180,p180],tint_cut);
f270 = irf_tlim([t270,p270],tint_cut);
f_comb = zeros(length(f180),18);

for i = 1:length(f180)
    for j=1:18
        f_comb(i,j) = nanmean([f180(i,j+1),f270(i,j+1)]);
%         f_comb(i,j) = nansum([f180(i,j+1),f270(i,j+1)]);
    end
end

specepad4=struct('t',f180(:,1));
specepad4.p  = f_comb;
specepad4.p_label={' ',' cm-2 s-1 sr-1 keV-1 '};
specepad4.f_label={''};
specepad4.f =[5:10:175];
%specepad是98.4keV

%% smooth
ista = 15;
f_comb2 = zeros(floor(length(f180)/30)-1,18);
t_comb2 = zeros(length(f_comb2),1);
for i = 1:floor(length(f180)/30)-1
    for j=1:18
        f_comb2(i,j) = nanmean(f_comb(ista+(i-1)*30:ista+i*30,j));
    end
    t_comb2(i,1) = f180(ista+(i-1)*30+15, 1);
end

specepadnum=struct('t',t_comb2);
specepadnum.p  = f_comb2;
specepadnum.p_label={' ',' cm-2 s-1 sr-1 keV-1 '};
specepadnum.f_label={''};
specepadnum.f = [5:10:175];



%% plot
n_subplots=4;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters');
xSize = 19; ySize = 29; coef=floor(min(600/xSize,600/ySize));
xLeft = (20-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[20 20 400 800]);
colormap(jet);
%% B plot
%% plot B (1-5000ev)
%第一幅图 总磁场
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bt(:,1) Bt(:,2)], 'color','k', 'Linewidth',1.5); hold on;
irf_plot([B(:,1) B(:,2)*0],'k--', 'Linewidth',1);hold off;
grid off;
y=ylabel({'|B|','(nT)'},'fontsize',8);
% irf_legend(h(1),'|B|',[0.05 0.95],'color','k','fontsize',12)
% y=ylabel({'Bx，By','(nT)'},'fontsize',8);
% set(y,'position',get(y,'position')-[0.1,0,0]); %将y标签向左移动0.1
% set(h(1),'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% set(gca,'Ylim',[-600 1500]);
[maxB,I1]=max(Bt(:,2));
BB1=B(:,2);BB2=B(:,3);BB3=B(:,4);
BBB=[BB1 BB2 BB3];
minB=min(min(BBB));
set(gca,'Ylim',[-5 1.2*maxB]);
% set(gca,'Ylim',[1.2*minB 1.2*maxB]);
% pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[0 1 0];[1 0 0];[0 0 1]]);%这一行不同的代号代表着不同的颜色
irf_legend(h(1),'(c)',[0.99 0.98],'color','k','fontsize',12)

%第二幅图 磁场三分量
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B(:,1) B(:,2)], 'color','b', 'Linewidth',1.5); hold on;
irf_plot([B(:,1) -B(:,3)], 'color','g', 'Linewidth',1.5); hold on;
irf_plot([B(:,1) -B(:,4)], 'color','r', 'Linewidth',1.5); hold on;
%irf_plot([Bt(:,1) Bt(:,2)], 'color','b', 'Linewidth',1); hold on;
irf_plot([B(:,1) B(:,2)*0],'k--', 'Linewidth',1);hold off;
grid off;
y=ylabel({'B','(nT)'},'fontsize',9);
% y=ylabel({'Bz','(nT)'},'fontsize',9);
% set(y,'position',get(y,'position')-[0.1,0,0]); %将y标签向左移动0.1
% set(h(1),'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% set(gca,'Ylim',[-1200 -600]);
[maxB,I1]=max(Bt(:,2));
BB1=B(:,2);BB2=B(:,3);BB3=B(:,4);
BBB=[BB1 BB2 BB3];
minB=min(min(BBB));
% set(gca,'Ylim',[1.2*minB 1.2*maxB]);
set(gca,'Ylim',[1.2*minB 10]);
% pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[0 1 0];[1 0 0];[0 0 1]]);%这一行不同的代号代表着不同的颜色
%irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.1 0.12],'fontsize',8);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
irf_legend(h(2),'(d)',[0.99 0.98],'color','k','fontsize',12)
irf_legend(h(2),'Bx',[0.05 0.95],'color','b','fontsize',12)
irf_legend(h(2),'By',[0.15 0.95],'color','g','fontsize',12)
irf_legend(h(2),'Bz',[0.25 0.95],'color','r','fontsize',12)
% irf_legend(h(2),'(d)',[0.99 0.98],'color','k','fontsize',12)
% irf_legend(h(2),'B\theta',[0.25 0.95],'color','r','fontsize',12)
% irf_legend(h(2),'B\phi',[0.05 0.95],'color','b','fontsize',12)
% irf_legend(h(2),'Br',[0.15 0.95],'color','g','fontsize',12)


%% plot PA
%第三幅图 电子投掷角
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
[h(3) hcb3]=irf_spectrogram(h(3),specepadnum,'log');
set(h(3),'ytick',[45 90 135]);
set(h(3),'ylim',[0 180]);
%   caxis(h(4),[1.7 3]);
ylabel({'PA [\circ]'},'fontsize',8)
poscbarx=get(hcb3,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb3,'pos',poscbarx);
set(hcb3,'fontsize',8);
irf_legend(h(3),'(e)',[0.99 0.98],'color','k','fontsize',12)
irf_legend(h(3),'30keV-834keV',[0.02 0.75],'color','k','fontsize',10)


%% plot Ee (1-5000ev)
%第四幅图 电子能谱
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
[h(4) hcb4]=irf_spectrogram(h(4),specEe);
ylabel({'Ee','(keV)'},'fontsize',8)
% caxis(h(4),[3 6.3])
set(gca,'yscale','log');
set(gca,'ytick',[1e2 1e3 1e4 1e5]);
set(gca,'ylim',[0 1000]);
% set(gca,'ylim',[10 400]);
poscbarx=get(hcb4,'pos');
poscbarx(3)=poscbarx(3)*0.5;
set(hcb4,'pos',poscbarx);
set(hcb4,'fontsize',8);
irf_legend(h(4),'(f)',[0.99 0.98],'color','k','fontsize',12)
Tsta='2017-03-19T19:32:00Z'; 
Tend='2017-03-19T20:25:00Z'; 
tint1=irf.tint(Tsta,Tend);%irf.tint语法：irf.tint(起始时间,终止时间）

irf_zoom(tintS,'x',h(1:end));
irf_plot_axis_align;

%% save figure
 % 保存图形为矢量图
% print(gcf,'-depsc','figname') % 使用print函数，gcf表示当前图形窗口，'-depsc'表示eps格式，figname表示文件名
% print(gcf,'-dpdf','overview_figure4') % 使用print函数，gcf表示当前图形窗口，'-depsc'表示eps格式，figname表示文件名
%   print(gcf, '-dpdf', [figname '.pdf'])
%   print(gcf, '-dpng', [figname '.png'])
% irf_zoom(h(1),'x',tint);
% irf_zoom(h(2),'x',tint);
% irf_zoom(h(3),'x',tint);
% irf_zoom(h(4),'x',tint);