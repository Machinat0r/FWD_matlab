clear; close all
clc

ic=2;

% tint=irf.tint('2017-08-23T15:38:20Z/2017-08-23T15:39:15Z');
% tint=irf.tint('2019-07-19T13:46:50.000Z/2019-07-19T13:47:20.000Z');
% tint=irf.tint('2019-08-16T09:31:57.500Z/2019-08-16T09:31:59.000Z');
tint=irf.tint('2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.900Z');

c_eval('fpiFile? = dataobj(''D:\MMS\mms2\fpi\brst\l2\des-dist\2020\08\03\mms2_fpi_brst_l2_des-dist_20200803014513_v3.4.0.cdf'');',ic);
c_eval('diste_struct = get_variable(fpiFile?,''mms?_des_dist_brst'');',ic);
% c_eval('diste_struct = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
thetae=diste_struct.DEPEND_2.data;
energye0=diste_struct.DEPEND_3.data(1,1:end);
energye1=diste_struct.DEPEND_3.data(2,1:end);

c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
% c_eval('energy0=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy0_brst'',Tintr);',ic);
% c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
c_eval('phie=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_phi_brst'',tint);',ic);
% c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);
c_eval('stepTablee=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_steptable_parity_brst'',tint);',ic);


diste.data = diste.data*1e30; 
% data = diste.data(:,:,:,[1:3,14:16]);
% thetae = thetae([1:3,14:16]);
data = diste.data;

energyspec = ones(length(diste.time),1)*energye0;
for ii = 1:length(diste.time)
    if stepTablee.data(ii)
        energyspec(ii,:) = energye1;
    end
end


% define angles
dangle = pi/16;
lengthphi = 32;

z2 = ones(lengthphi,1)*sind(thetae);
solida = dangle*dangle*z2;
allsolide = zeros(size(data));

for ii = 1:length(diste.time)
    for jj=1:length(energye0)
        allsolide(ii,jj,:,:) = solida;
    end
end

distes = data.*allsolide;

% Electron analysis - OMNI
PSDomni = zeros(length(diste.time),length(energye0));
for ii = 1:length(diste.time)
    disttemp = squeeze(distes(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
end

timeomni=irf_time(phie.time,'epochTT>epoch');
% time1='2019-08-16T09:31:57.500Z';
time1='2020-08-03T01:45:27.500Z';
time1=iso2epoch(time1);
% % % 


idt1=find(timeomni>time1);
t=idt1(1);
dt = 40;
chan=energyspec(t,:);
psd=mean(PSDomni(t:t+dt,:));
psd=transpose(psd);
% loglog(chan,psd,'b+');
hold on;



% % % for i=2:4
% % % c_eval('idt?=find(timeomni>time?);',i);
% % % c_eval('t?=idt?(1);',i);
% % % c_eval('chan?=mean(energyspec(t?:t?+dt,:));',i);
% % % 
% % %     c_eval('psd?=mean(PSDomni(t?:t?+dt,:));',i);
% % %     c_eval('psd?=transpose(psd?);',i);
% % % end



%% 设置图形属性并画图
n_subplots=1;

i_subplot=1;
set(0,'DefaultAxesFontSize',9);
set(0,'DefaultLineLineWidth', 1);

set(gcf,'PaperUnits','centimeters')
xSize = 6; ySize = 16; coef=floor(min(1000/xSize,1000/ySize));
xLeft = 5; yTop = -1;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
lnwid=0.3;
m=9.1*10^-31;

k=1.38e-23; 
e=1.6e-19;

%% 拟合前九点
% Ne=0.1;
% Te=50;
% X1=[14.63,19.15,25.07,32.810001,42.950001,56.23,73.599998,96.339996,126.12,165.09,216.11,282.89001,370.31,484.73999,634.53998,830.63,1087.3101]%此处X1取能道在14.63eV（包含）之后
% F1=Ne*1e24*power(m/(2*pi*k*Te*11605),3/2)*exp(-X1*e/(k*Te*11605));%此处X1取能道在14.63eV（包含）之后
Ne=0.10;
Te=15;
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);
% X1=[32.810001,42.950001,56.23,73.599998,96.339996,126.12,165.09,216.11,282.89001,370.31,484.73999];%此处X1取能道在165eV（包含）之后
X1 = energye0(6:12);
F1=Ne*1e24*power(m/(2*pi*k*Te*11605),3/2)*exp(-X1*e/(k*Te*11605));%此处X2取能道在165eV（包含）之后

%% 拟合后九点
Ne=0.23;
Te=1500;
% Ne=0.13;
% Te=500;
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);
% X2=[830.63,1087.3101,1423.3199,1863.16,2438.9199,3192.6101,4179.2002,5470.6802,7161.25,9374.25,12271.12,16063.2,21027.109,27525];%此处X1取能道在165eV（包含）之后
X2 = energye0(10:28);
F2=Ne*1e24*power(m/(2*pi*k*Te*11605),3/2)*exp(-X2*e/(k*Te*11605));%此处X2取能道在165eV（包含）之后

%% 画图
xSize = 60; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%%
plot(X1,F1,'-.');% 画前九个点拟合图像
hold on

plot(X2,F2,'-.');% 画后九个点拟合图像
hold on

plot(chan,psd,'color','k');hold on;
% loglog(chan,psd,'^','color','#007CB3','MarkerSize',4,'MarkerFaceColor','#007CB3');hold on;
% plot(chan2,psd2,'b');hold on;
% loglog(chan2,psd2,'bv','MarkerSize',4,'MarkerFaceColor','b');hold on;
% loglog(chan3,psd3,'ko');
% loglog(chan4,psd4,'r^');

set(h(1),'Xlim',[30 1e5]);
set(h(1),'Ylim',[1e-4 1e3]);
set(h(1),'yscale','log');set(h(1),'xscale','log');
% end

ylabel('PSD (s^3 km^{-6})','fontsize',12);
xlabel('E (ev)','fontsize',12);
% set(gca,'Xlim',[1700 6000])
% irf_legend({'^ 16:24:10.00'},[0.9 0.9],'color','r')
% irf_legend({'v 16:24:30.00'},[0.9 0.8],'color','b')
% irf_legend({'○ 01:59:45.00'},[0.9 0.7],'color','k')
% irf_legend({'Δ 01:59:50.00'},[0.9 0.6],'color','r')


