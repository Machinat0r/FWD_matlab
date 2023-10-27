%------written by Yue Yu in Beijing------------
% modified the curve fit sesction to ensure the density can be identical to
% the observation
%------modified by Wending Fu, Sept.18.2023 in Beijing------------
%%
clear; 

ic=1;
tint = irf.tint('2015-12-15T23:43:20.000Z/2015-12-15T23:44:21.000Z');
mms.db_init('local_file_db','/Users/fwd/Documents/MATLAB/MMS/');
%% Load data
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
tint1 = irf_time(timeUTC1(1,:),'utc>epochTT');
[paddist10,thetapad1,energypad1,tintpad1] = mms_get_pitchangledist_my_change(diste1,Bxyz1,tint1);
paddist11 = paddist10*1e30;

psd1=nanmean(paddist11(:,1:30),2);
psd1_par=nanmean(paddist11(:,1:6),2);
psd1_perp=nanmean(paddist11(:,12:21),2);
psd1_antipar=nanmean(paddist11(:,25:30),2);
psd1 = psd1_par;
%% parameters
m=9.1*10^-31;
k=1.38e-23; 
e=1.6e-19;
%% curve fitting function
energy_section = {5:20;11:14;15:25};
% low energy electron maxwell
X1=energypad1(energy_section{1});
Y1=double(psd1(energy_section{1}));
p01 = [10,100];

% mid energy electron maxwell
X2=energypad1(energy_section{2});
Y2=double(psd1(energy_section{2}));
p02 = [0.1,60];

% high energy electron maxwell
X3=energypad1(energy_section{3});
Y3=double(psd1(energy_section{3}));
p03 = [0.01,6000];

% constraint condition
fun4 = @(p,Ne)p(1)+p(3)+p(5)-Ne;
p0 = [p01,p02,p03];
p = lsqnonlin(@(p)MaxwellFunc(p, X1', X2', X3', Y1, Y2, Y3, Ne), p0);

% maxwell function ouput
X1=energypad1(min(energy_section{1}):32);
F1=p(1)*1e24*power(m/(2*pi*k*p(2)*11605),3/2)*exp(-X1*e/(k*p(2)*11605));

% X2=energypad1(min(energy_section{2}):32);
F2=p(3)*1e24*power(m/(2*pi*k*p(4)*11605),3/2)*exp(-X2*e/(k*p(4)*11605));

X3=energypad1(min(energy_section{3}):32);
F3=p(5)*1e24*power(m/(2*pi*k*p(6)*11605),3/2)*exp(-X3*e/(k*p(6)*11605));
%%  plot figure
fn1=figure(1);
xSize = 80; ySize = 50; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(fn1,'PaperPosition',[xLeft yTop xSize ySize]);
set(fn1,'Position',[10 10 xSize*coef ySize*coef]);

h(1)=axes;
ymin = 10^-6;
yrange=[10^-3 10^5];
% ymax = ceil(max(max(log10(paddist_aver1(1:32,:)))));
% yrange = [ymin 10^ymax];
% yrange=[ymin 10^4];


plot(h(1),energypad1,psd1,'k.-','LineWidth',1);
hold on;
plot(h(1),X1,F1,'-.m','LineWidth',1);hold on;
% plot(h(1),X2,F2,'-.c','LineWidth',1);hold on;
% plot(h(1),X3,F3,'-.b','LineWidth',1);hold on;
hold(h(1),'off');

ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)');
set(h(1),'yscale','log');
set(h(1),'xscale','log');

irf_zoom(h(1),'y',yrange);
irf_zoom(h(1),'x',[1 20000]);

% irf_legend(h(1),{'23:22:22.14Z/23:22:22.58Z'},[0.91 0.96],'color','k','Fontsize',10);
% irf_legend(h(1),{'19:40:30.00-19:40:55.00'},[0.91 0.88],'color','r','Fontsize',10);
% irf_legend(h(1),{'19:40:55.00-19:41:10.00'},[0.91 0.80],'color','g','Fontsize',10);
% irf_legend(h(1),'Ne=0.4cm^{-3},Te=5eV',[0.1 0.1]);
% irf_legend(h(1),'Ne=0.31cm^{-3},Te=75eV',[0.3 0.1]);
% irf_legend(h(1),'PSD \propto \epsilon^{-3.48}',[0.8 0.1]);

% title(h(1),strcat('omni'));
% set(h(1),'Fontsize',12);

set(gcf,'render','painters');
figname=['R2fitomini0720'];

%% Maxwell function
function Maxwellfun = MaxwellFunc(p, X1, X2, X3, Y1, Y2, Y3, Ne)
Maxwellfun = [p(1)-Ne;...
    p(1)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(2)*11605),3/2)*exp(-X1*(1.6e-19)/((1.38e-23)*p(2)*11605))-Y1;...
    % p(3)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(4)*11605),3/2)*exp(-X2*(1.6e-19)/((1.38e-23)*p(4)*11605))-Y2;...
    % p(5)*1e24*power((9.1*10^-31)/(2*pi*(1.38e-23)*p(6)*11605),3/2)*exp(-X3*(1.6e-19)/((1.38e-23)*p(6)*11605))-Y3...
    ];
end