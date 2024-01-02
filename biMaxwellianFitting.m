%------written by Wending Fu, Sept.18.2023 in Beijing------------
%%
clear; close all; clc
global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
TempDir = [ParentDir,'temp/'];mkdir(TempDir);

ic=3;
% tint = irf.tint('2019-07-19T13:47:05.200Z/2019-07-19T13:47:06.500Z');
TT = '2016-11-21T07:22:21.800Z/2016-11-21T07:22:21.900Z';
tint = irf.tint(TT);

inpath = '/Users/fwd/Documents/MATLAB/Code/XIE_HS/pdrk/pdrk_master_v181027/pdrk-master/BEW/input/pdrk.in';

Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
filenames = [filenames1,filenames2];

[filenames,~,~] = findFilenames(TT,filenames,'brst',ic);
SDCFilesDownload_NAS(filenames,TempDir)
SDCDataMove(TempDir,ParentDir)
%% Load in data
units = irf_units;
indata = importdata(inpath).data;
n = indata(:,3); % m-3
Tpara = indata(:,4) * units.e; % J
Tperp = indata(:,5) * units.e; % J
Vd = indata(:,8)*units.c; % m/s
%% Load data
mms.db_init('local_file_db',ParentDir);
Bxyz1=mms.get_data('B_gse_brst',tint,ic);

c_eval('Ne_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval(['Ne=irf.ts2mat(Ne_ts);'],ic);
Ne = Ne(1,2);

c_eval('diste1 = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
theta1=diste1.DEPEND_2.data;
c_eval('diste1 = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
length_diste1=length(diste1.time.epoch);
timeUTC1=irf_time(diste1.time,'epochtt>utc');
%% Produce PAD at all selected time
% tint1 = irf_time(timeUTC1(1,:),'utc>epochTT');
tint1 = tint;
[paddist10,thetapad1,energypad1,tintpad1] = mms_get_pitchangledist_my_change(diste1,Bxyz1,tint1);
energypad1 = energypad1(1,:);
paddist11 = paddist10*1e30;
paddist11 = mean(paddist11.data,1);
paddist11 = reshape(paddist11, 32, 30);

psd1=nanmean(paddist11(:,1:30),2);
psd1_par=nanmean(paddist11(:,1:10),2);
psd1_perp=nanmean(paddist11(:,11:20),2);
psd1_antipar=nanmean(paddist11(:,21:30),2);
psd1 = psd1_par;
psd2 = psd1_perp;
psd3 = psd1_antipar;
%% parameters
% % % m=units.me;
% % % k=units.kB; 
% % % e=units.e;
%% curve fitting function
% % % energy_section = {8:10;11:14;15:25};
% % % % low energy electron maxwell
% % % X1=energypad1(energy_section{1});
% % % Y1=double(psd1(energy_section{1}));
% % % p01 = [0.1,5];
% % % 
% % % % mid energy electron maxwell
% % % X2=energypad1(energy_section{2});
% % % Y2=double(psd1(energy_section{2}));
% % % p02 = [0.1,60];
% % % 
% % % % high energy electron maxwell
% % % X3=energypad1(energy_section{3});
% % % Y3=double(psd1(energy_section{3}));
% % % p03 = [0.01,6000];
% % % 
% % % % constraint condition
% % % fun4 = @(p,Ne)p(1)+p(3)+p(5)-Ne;
% % % p0 = [p01,p02,p03];
% % % p = lsqnonlin(@(p)MaxwellFunc(p, X1', X2', X3', Y1, Y2, Y3, Ne), p0);
% % % 
% % % % maxwell function ouput
% % % X1=energypad1(min(energy_section{1}):32);
% % % F1=p(1)*1e24*power(m/(2*pi*k*p(2)*11605),3/2)*exp(-X1*e/(k*p(2)*11605));
% % % 
% % % % X2=energypad1(min(energy_section{2}):32);
% % % F2=p(3)*1e24*power(m/(2*pi*k*p(4)*11605),3/2)*exp(-X2*e/(k*p(4)*11605));
% % % 
% % % X3=energypad1(min(energy_section{3}):32);
% % % F3=p(5)*1e24*power(m/(2*pi*k*p(6)*11605),3/2)*exp(-X3*e/(k*p(6)*11605));
%% Flux from biMaxwell function
idx = 1; % 1 for parallel flux, and 2 for perpendicular flux
FbiMaxwell1 = zeros(size(energypad1));
FbiMaxwell2 = FbiMaxwell1; FbiMaxwell3 = FbiMaxwell1;
for Eid = 1:length(energypad1)
    E = energypad1(Eid);
    for Sid = 1:length(n)
        tempflux1 = biMaxwellFunc(E, n(Sid), Tpara(Sid), Tperp(Sid), Vd(Sid), 1);
        FbiMaxwell1(Eid) = FbiMaxwell1(Eid) + tempflux1 * 1e18;
        tempflux2 = biMaxwellFunc(E, n(Sid), Tpara(Sid), Tperp(Sid), Vd(Sid), 2);
        FbiMaxwell2(Eid) = FbiMaxwell2(Eid) + tempflux2 * 1e18;
        tempflux3 = biMaxwellFunc(E, n(Sid), Tpara(Sid), Tperp(Sid), Vd(Sid), 3);
        FbiMaxwell3(Eid) = FbiMaxwell3(Eid) + tempflux3 * 1e18;
    end
end
%%  plot figure
fn1=figure(1);
xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(fn1,'PaperPosition',[xLeft yTop xSize ySize]);
set(fn1,'Position',[10 10 xSize*coef ySize*coef]);

h(1)=axes;
ymin = 10^-6;
yrange=[10^-4 10^6];
% ymax = ceil(max(max(log10(paddist_aver1(1:32,:)))));
% yrange = [ymin 10^ymax];
% yrange=[ymin 10^4];


plot(h(1),energypad1,psd1,'ko');hold on;
plot(h(1),energypad1,FbiMaxwell1,'k',LineWidth=2); hold on
% plot(h(1),X1,F1,'-.m','LineWidth',1);hold on;
% plot(h(1),X2,F2,'-.c','LineWidth',1);hold on;
% plot(h(1),X3,F3,'-.b','LineWidth',1);hold on;

ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)');
set(h(1),'yscale','log');
set(h(1),'xscale','log');

irf_zoom(h(1),'y',yrange);
irf_zoom(h(1),'x',[10 1e4]);

% irf_legend(h(1),{'23:22:22.14Z/23:22:22.58Z'},[0.91 0.96],'color','k','Fontsize',10);
% irf_legend(h(1),{'19:40:30.00-19:40:55.00'},[0.91 0.88],'color','r','Fontsize',10);
% irf_legend(h(1),{'19:40:55.00-19:41:10.00'},[0.91 0.80],'color','g','Fontsize',10);
% irf_legend(h(1),'Ne=0.4cm^{-3},Te=5eV',[0.1 0.1]);
% irf_legend(h(1),'Ne=0.31cm^{-3},Te=75eV',[0.3 0.1]);
% irf_legend(h(1),'PSD \propto \epsilon^{-3.48}',[0.8 0.1]);

% title(h(1),strcat('F_{||}'));
% set(h(1),'Fontsize',12);

% set(gcf,'render','painters');
% figname=['R2fitomini0720'];
%% figure 2
% % % fn2=figure(2);
% % % xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
% % % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % % set(fn2,'PaperPosition',[xLeft yTop xSize ySize]);
% % % set(fn2,'Position',[10 10 xSize*coef ySize*coef]);
% % % 
% % % h(1)=axes;
% % % ymin = 10^-6;
% % % yrange=[10^-4 10^6];
% ymax = ceil(max(max(log10(paddist_aver1(1:32,:)))));
% yrange = [ymin 10^ymax];
% yrange=[ymin 10^4];


plot(h(1),energypad1,psd2,'ro');hold on;
plot(h(1),energypad1,FbiMaxwell2,'r',LineWidth=2);hold on;
ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)');
set(h(1),'yscale','log');
set(h(1),'xscale','log');
% title(h(1),strcat('F_{anti-para}'));

irf_zoom(h(1),'y',yrange);
irf_zoom(h(1),'x',[10 1e4]);
%% figure 3
% % % fn3=figure(3);
% % % xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
% % % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % % set(fn3,'PaperPosition',[xLeft yTop xSize ySize]);
% % % set(fn3,'Position',[10 10 xSize*coef ySize*coef]);
% % % 
% % % h(1)=axes;
% % % ymin = 10^-6;
% % % yrange=[10^-4 10^6];
% ymax = ceil(max(max(log10(paddist_aver1(1:32,:)))));
% yrange = [ymin 10^ymax];
% yrange=[ymin 10^4];


plot(h(1),energypad1,psd3,'bo');hold on;
plot(h(1),energypad1,FbiMaxwell3,'b',LineWidth=2);hold on;
ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)');
set(h(1),'yscale','log');
set(h(1),'xscale','log');
% title(h(1),strcat('F_{\perp}'));
legend({'','0','','90','','180'})

irf_zoom(h(1),'y',yrange);
irf_zoom(h(1),'x',[10 1e4]);
%% Maxwell function
function MaxwellFun = MaxwellFunc(p, X1, X2, X3, Y1, Y2, Y3, Ne)
MaxwellFun = [p(1)+p(5)-Ne;...
    p(1)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(2)*11605),3/2)*exp(-X1*(1.6e-19)/((1.38e-23)*p(2)*11605))-Y1;...
    p(3)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(4)*11605),3/2)*exp(-X2*(1.6e-19)/((1.38e-23)*p(4)*11605))-Y2;...
    p(5)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(6)*11605),3/2)*exp(-X3*(1.6e-19)/((1.38e-23)*p(6)*11605))-Y3];
end
%% convert E to V
% eV to m/s
function V = convertEtoV(E)
units = irf_units;
V = sqrt(2*E*units.e/units.me);
end
%% convert V to E
% c to eV
function E = convertVtoE(V)
units = irf_units;
V = V/units.c;
E = 0.5*units.me*V^2/units.e;
end
%% biMaxwell function
function flux = biMaxwellFunc(E, n, Tpara, Tperp, Vd, idx)
units = irf_units;
m=units.me;
k=units.kB; 
e=units.e;
switch idx
case 1
Vpara  = convertEtoV(E);
Vperp = 0;
case 2
Vperp  = convertEtoV(E);
Vpara = 0;
case 3
Vpara  = convertEtoV(E);
Vperp=0;
Vd = -Vd;
end
Vtpara = sqrt(2*Tpara/units.me); % m/s
Vtperp = sqrt(2*Tperp/units.me); % m/s
flux = n * (sqrt(pi)*Vtpara)^-3 * Tpara/Tperp * exp(-(Vpara-Vd)^2/(Vtpara)^2) * exp(-Vperp^2/(Vtperp^2));
end