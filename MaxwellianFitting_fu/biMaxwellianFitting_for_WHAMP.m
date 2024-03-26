clear;
clc;
close all

global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
mms.db_init('local_file_db', ParentDir);
inpath = '/Users/fwd/Documents/MATLAB/Code/XIE_HS/pdrk/pdrk_master_v181027/pdrk-master/ECW/input/yy.in';
ic=1;

% Tintr = irf.tint('2021-07-22T05:20:55.00Z/2021-07-22T05:20:58.00Z');
Tintr = irf.tint('2019-07-19T13:47:05.200Z/2019-07-19T13:47:06.500Z');
%% Load plasma parameters
units = irf_units;
indata = importdata(inpath).data;
n = indata(:,3); % m-3
T = indata(:,4); % eV
Tps = indata(:,5); % eV
% Q = indata(:,5); % eV
Q= T./Tps;
Vd = indata(:,8) * units.c; % m/s
%% Load data
%
tic;
c_eval('diste = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',Tintr);',ic);
%theta=diste.DEPEND_2.data;
c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',Tintr);',ic);
c_eval('energy=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy_brst'',Tintr);',ic);

% c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
c_eval('phi=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_phi_brst'',Tintr);',ic);
% c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);

c_eval('stepTable=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_steptable_parity_brst'',Tintr);',ic);

c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',Tintr,?);',ic);
c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',Tintr);',ic);
c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',Tintr);',ic);


toc;

Bxyz=mms.get_data('B_gse_brst',Tintr,ic);
Bxyz = irf_gse2gsm(Bxyz);
% Bxyz = Bxyz.resample(diste);
%
% c_eval('diste = mms.db_get_variable(''mms?_fpi_fast_l2_des-dist'',''mms?_des_dist_fast'',Tintr);',ic);
% theta=diste.DEPEND_2.data;
% c_eval('diste = mms.db_get_ts(''mms?_fpi_fast_l2_des-dist'',''mms?_des_dist_fast'',Tintr);',ic);
% c_eval('energy=mms.db_get_variable(''mms?_fpi_fast_l2_des-dist'',''mms?_des_energy_fast'',Tintr);',ic);
% 
% % c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
% % c_eval('phi=mms.db_get_ts(''mms?_fpi_fast_l2_des-dist'',''mms?_des_phi_fast'',Tintr);',ic);
% % c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);
% 
% c_eval('stepTable=mms.db_get_ts(''mms?_fpi_fast_l2_des-dist'',''mms?_des_steptable_parity_fast'',Tintr);',ic);
% toc;
% 
% Bxyz=mms.get_data('B_gse_brst',Tintr,ic);

%% Produce a single PAD at a selected time
% tint1 = irf_time('2017-06-11T01:59:25.00Z','utc>epochTT');
tint1=irf.tint('2019-07-19T13:47:05.200Z/2019-07-19T13:47:06.500Z');
[paddist01,thetapad,energypad,tintpad1] = mms_get_pitchangledist_my_change(diste,Bxyz,tint1); 
paddist01 = paddist01*1e30; 
%
paddist1=irf.nanmean(paddist01.data);
paddist1=squeeze(paddist1);
%

for ii=1:32,
    sum_paddist1(ii,:)=nansum(paddist1(ii,:));
end
paddist_1=paddist1(17:32,:);

paddist11 = paddist01;
paddist11 = mean(paddist11.data,1);
paddist11 = reshape(paddist11, 32, 30);

psd1=nanmean(paddist11(:,1:30),2);
psd1_par=nanmean(paddist11(:,1:10),2);
psd1_perp=nanmean(paddist11(:,11:20),2);
psd1_antipar=nanmean(paddist11(:,21:30),2);
psd1 = psd1_par;
psd2 = psd1_perp;
psd3 = psd1_antipar;
%% calculate flux
% % % Vd = 5e5;
% % % c_eval('n = Ne?_ts.tlim(tint1);',ic);
% % % n = irf.ts2mat(n);
% % % n = mean(n(:,2));
% % % c_eval('Tpara = Te_para?_ts.tlim(tint1);',ic);
% % % c_eval('Tperp = Te_perp?_ts.tlim(tint1);',ic);
% % % Tpara = irf.ts2mat(Tpara);
% % % Tperp = irf.ts2mat(Tperp);
% % % T = [Tpara(:,1),(Tpara(:,2)+2*Tperp(:,2))/3.0];
% % % T = mean(T(:,2));
FbiMaxwell1 = zeros(3, size(energypad,2));
for Eid = 1:size(energypad,2)
    E = energypad(1,Eid);
    theta = [0,90,180];
    for Sid = 1:length(n)
        tempn = n(Sid);
        tempT = T(Sid);
        tempQ = Q(Sid);
        tempVd = Vd(Sid);
        tempflux1 = biMaxwellFunc(tempT, tempQ, tempn, E, tempVd, theta(1));
        tempflux2 = biMaxwellFunc(tempT, tempQ, tempn, E, tempVd, theta(2));
        tempflux3 = biMaxwellFunc(tempT, tempQ, tempn, E, tempVd, theta(3));
        FbiMaxwell1(1, Eid) = FbiMaxwell1(1, Eid) + tempflux1 * 1e18;
        FbiMaxwell1(2, Eid) = FbiMaxwell1(2, Eid) + tempflux2 * 1e18;
        FbiMaxwell1(3, Eid) = FbiMaxwell1(3, Eid) + tempflux3 * 1e18;
    end
end


Flux_calc = zeros(size(energypad, 2), size(thetapad, 2));
for k=1:32
E=energypad(1,k);
theta=thetapad;
for j = 1:length(n)
tempn = n(j);
tempT = T(j);
tempQ = Q(j);
tempVd = Vd(j);

Flux_calc(k,:) = Flux_calc(k,:) + biMaxwellFunc(tempT, tempQ, tempn, E, tempVd, theta)*1e18;
end
end

% Flux_calc = Flux_calc * 1e30;
%% Plot PAD

fn=figure;
set(fn,'Position',[100 20 1000 600])
h(1)=axes('position',[0.06 0.19 0.35 0.65]); 
h(2)=axes('position',[0.60 0.19 0.35 0.65]); 

ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',1.5); 
%%
ymin = 10^-6;
yrange=[10^-4 10^6];
% ymax = ceil(max(max(log10(paddist_aver1(1:32,:)))));
% yrange = [ymin 10^ymax];
% yrange=[ymin 10^4];

hold(h(1),'on');
plot(h(1),energypad(1,:),psd1,'ko');
plot(h(1),energypad(1,:),FbiMaxwell1(1,:),'k',LineWidth=1);

plot(h(1),energypad(2,:),psd2,'ro');
plot(h(1),energypad(2,:),FbiMaxwell1(2,:),'r',LineWidth=1);

plot(h(1),energypad(3,:),psd3,'bo');
plot(h(1),energypad(3,:),FbiMaxwell1(3,:),'b',LineWidth=1);
hold(h(1),'off');
% plot(h(1),X1,F1,'-.m','LineWidth',1);hold on;
% plot(h(1),X2,F2,'-.c','LineWidth',1);hold on;
% plot(h(1),X3,F3,'-.b','LineWidth',1);hold on;

ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)');
set(h(1),'yscale','log');
set(h(1),'xscale','log');
set(h(1),'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(h(1),{'0 deg','90 deg','180 deg'},[0.97 0.92]);

irf_zoom(h(1),'y',yrange);
irf_zoom(h(1),'x',[10 1e4]);
%%
jetcolor = colormap('jet');

% c_eval('plot(h(2),thetapad,paddist(?,:),''color'',jetcolor(?+9,:));',1);
hold(h(2),'on');
% c_eval('plot(h(2),j,f?,''color'',jetcolor(?+2,:));',[2:18]);
% paddist_1=paddist1(13:32,:);
c_eval('scatter(h(2),thetapad,paddist_1(?,:),''color'',jetcolor(15*?+12,:));',1:16);
c_eval('plot(h(2),thetapad,Flux_calc(16+?,:),''color'',jetcolor(15*?+12,:));',1:16);
% hold(h(1),'on');
% c_eval('plot(h(1),thetapad,paddist_2(?,:),''--'',''color'',jetcolor(3*?+12,:));',1:16);
% for i=1;32,
% plot(h(2),thetapad,paddist(i,:),'color',jetcolor(2*i,:) );
% end
% c_eval('plot(h(2),thetapad,paddist(?,:));',[1:32]);
hold(h(2),'off');

% 
% for i=1:length(energypad(1:32))
%     a=roundn(energypad(1,:)/1000,-3);
% c_eval('irf_legend(h(2),''? keV'',[1.02 1.03-i*0.033],''color'',jetcolor(2*i,:));',a(i));
% end

ylabel(h(2),'f_e (s^3 km^{-6})');
xlabel(h(2),'\theta (deg.)')
set(h(2),'yscale','log');
set(h(2),'Ylim',[1e-4 5*1e1]);
irf_zoom(h(2),'x',[0 180]);
box on;
% irf_zoom(h(2),'y',yrange);
tintutc1 = tintpad1.utc;
 
% title(h(1),strcat('MMS',num2str(ic),'-',tintutc1(12:23),'UT'));
title(h(2),strcat('MMS',num2str(ic),'-',tintutc1(23:2:38),'UT'));
%%

set(gcf,'color','w');
set(gcf,'render','painters');
set(h,'box','on')

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf, '-dpdf', [figname '.pdf']);
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
%% bi-Maxwell function
function flux = biMaxwellFunc(T, Q, n, E, Vd, theta)
Vpara = convertEtoV(E)*cosd(theta);
Vperp = convertEtoV(E)*sind(theta);
Vth = convertEtoV(T); % m/s

flux = n .* (sqrt(pi).*Vth).^-3 .* Q .* exp(-(Vpara-Vd).^2./(Vth).^2) .* exp(-Q.*Vperp.^2./(Vth.^2));
end