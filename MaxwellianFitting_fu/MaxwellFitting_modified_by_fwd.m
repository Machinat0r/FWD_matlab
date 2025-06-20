%------written by Yue Yu in Beijing------------
% modified the curve fit sesction to ensure the density can be identical to
% the observation
%------modified by Wending Fu, Sept.18.2023 in Beijing------------
%%
close all
clear;clc

global ParentDir 
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

% TT = '2019-08-05T16:24:00.000Z/2019-08-05T16:25:00.000Z';
TT = '2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.900Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 2;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,des-dist');
filenames = [filenames1, filenames2];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 64, 'CheckSize', 0)
%% Load data
% % % Bxyz1=mms.get_data('B_gse_brst',tint,ic);
% % % 
% % % c_eval('Ne_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
% % % c_eval(['Ne=irf.ts2mat(Ne_ts);'],ic);
% % % Ne = mean(Ne(:,2));
% % % 
% % % % c_eval('diste1 = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
% % % c_eval('diste1 = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
% % % % theta1=diste1.DEPEND_2.data;
% % % % length_diste1=length(diste1.time.epoch);
% % % % length_diste1 = size(diste1.data,1);
% % % % timeUTC1=irf_time(diste1.time,'epochtt>utc');
%% Produce PAD at all selected time
% % % % tint1 = irf_time(timeUTC1(1,:),'utc>epochTT');
% % % tint1 = tint;
% % % [paddist10,thetapad1,energypad1,tintpad1] = mms_get_pitchangledist_my_change(diste1,Bxyz1,tint1);
% % % paddist10 = paddist10.data;
% % % paddist11 = paddist10*1e30;
% % % 
% % % 
% % % psd1=nanmean(paddist11(:,1:30),2);
% % % psd1_par=nanmean(paddist11(:,1:6),2);
% % % psd1_perp=nanmean(paddist11(:,15:16),2);
% % % psd1_antipar=nanmean(paddist11(:,25:30),2);
% % % psd1 = psd1;
% % % psd2 = psd1_par;
% % % psd3 = psd1_antipar;
%%
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
psd1=transpose(psd);
% loglog(chan,psd,'b+');
hold on;

%% parameters
m=9.1*10^-31;
k=1.38e-23; 
e=1.6e-19;
%% curve fitting function
energy_section = {6:10;10:25;25:32};
% low energy electron maxwell
energypad1 = double(energye1);
X1=energypad1(energy_section{1});
Y1=double(psd1(energy_section{1}));
p01 = [1,5];

% mid energy electron maxwell
X2=energypad1(energy_section{2});
Y2=double(psd1(energy_section{2}));
p02 = [0.1,4000];

% high energy electron maxwell
X3=energypad1(energy_section{3});
Y3=double(psd1(energy_section{3}));
p03 = [1e-6,-4];

% constraint condition
% fun4 = @(p,Ne)p(1)+p(3)+p(5)-Ne;
p0 = [p01,p02,p03];
% p = lsqnonlin(@(p)MaxwellFunc(p, X1', X2', X3', Y1, Y2, Y3), p0, 0, inf);
p1 = lsqnonlin(@(p)MaxwellFunc1(p, X1', X2', X3', Y1, Y2, Y3), p0, 0, inf);
p2 = lsqnonlin(@(p)MaxwellFunc2(p, X1', X2', X3', Y1, Y2, Y3), p0, 0, inf);
p3 = lsqnonlin(@(p)MaxwellFunc3(p, X1', X2', X3', Y1, Y2, Y3), p0, 0, inf);
kk=polyfit(log(X3'),log(Y3),1);


% maxwell function ouput
X1=energypad1(min(energy_section{1}):32);
p1(1)=0.10;
p1(2)=15;
F1=p1(1)*1e24*power(m/(2*pi*k*p1(2)*11605),3/2)*exp(-X1*e/(k*p1(2)*11605));

X2=energypad1(10:32);
p2(3) = 0.2503;
p2(4) = 1200;
F2=p2(3)*1e24*power(m/(2*pi*k*p2(4)*11605),3/2)*exp(-X2*e/(k*p2(4)*11605));

X3=energypad1(min(energy_section{3}):32);
% F3=p3(5)*1e24*power(m/(2*pi*k*p3(6)*11605),3/2)*exp(-X3*e/(k*p3(6)*11605));
% F3 = p3(5)*X3.^(p3(6));
b1=exp(1)^kk(2);
F3=b1*X3.^kk(1);
%%  plot figure
fn1=figure(1);
xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(fn1,'PaperPosition',[xLeft yTop xSize ySize]);
set(fn1,'Position',[10 10 xSize*coef ySize*coef]);

h(1)=axes;
yrange=[10^-4 10^4];
% ymax = ceil(max(max(log10(paddist_aver1(1:32,:)))));
% yrange = [ymin 10^ymax];
% yrange=[ymin 10^4];


plot(h(1),energypad1,psd1,'k.-','LineWidth',1);hold on;
% plot(h(1),energypad1,psd2,'b.-','LineWidth',1);hold on;
% plot(h(1),energypad1,psd3,'r.-','LineWidth',1);hold on;

plot(h(1),X1,F1,'-.m','LineWidth',1);hold on;
plot(h(1),X2,F2,'-.c','LineWidth',1);hold on;
plot(h(1),X3,F3,'-.b','LineWidth',1);hold on;
hold(h(1),'off');

ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)');
set(h(1),'yscale','log');
set(h(1),'xscale','log');

set(h(1),'ylim',yrange);
set(h(1),'xlim',[3e1 3e4]);

% irf_legend(h(1),{'23:22:22.14Z/23:22:22.58Z'},[0.91 0.96],'color','k','Fontsize',10);
% irf_legend(h(1),{'19:40:30.00-19:40:55.00'},[0.91 0.88],'color','r','Fontsize',10);
% irf_legend(h(1),{'19:40:55.00-19:41:10.00'},[0.91 0.80],'color','g','Fontsize',10);
% irf_legend(h(1),'Ne=0.4cm^{-3},Te=5eV',[0.1 0.1]);
% irf_legend(h(1),'Ne=0.31cm^{-3},Te=75eV',[0.3 0.1]);
% irf_legend(h(1),'PSD \propto \epsilon^{-3.48}',[0.8 0.1]);

% title(h(1),strcat('omni'));
% set(h(1),'Fontsize',12);

colormap(jet)
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
%% Maxwell function
function Maxwellfun = MaxwellFunc(p, X1, X2, X3, Y1, Y2, Y3)
Maxwellfun = [...
    % p(1)+p(3)+p(5)-Ne;...
    p(1)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(2)*11605),3/2)*exp(-X1*(1.6e-19)/((1.38e-23)*p(2)*11605))-Y1;...
    p(3)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(4)*11605),3/2)*exp(-X2*(1.6e-19)/((1.38e-23)*p(4)*11605))-Y2;...
    p(5)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(6)*11605),3/2)*exp(-X3*(1.6e-19)/((1.38e-23)*p(6)*11605))-Y3...
    ];
end

function Maxwellfun = MaxwellFunc1(p, X1, X2, X3, Y1, Y2, Y3)
Maxwellfun = ...
    log10(p(1)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(2)*11605),3/2)*exp(-X1*(1.6e-19)/((1.38e-23)*p(2)*11605)))-log10(Y1);
end

function Maxwellfun = MaxwellFunc2(p, X1, X2, X3, Y1, Y2, Y3)
Maxwellfun = ...
log10(p(3)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(4)*11605),3/2)*exp(-X2*(1.6e-19)/((1.38e-23)*p(4)*11605)))-log10(Y2);
end

function Maxwellfun = MaxwellFunc3(p, X1, X2, X3, Y1, Y2, Y3)
Maxwellfun = ...
log10(p(5)*X3.^(p(6)))-log10(Y3);
end