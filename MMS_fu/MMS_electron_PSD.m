% % Script to plot electron PSD around pitch angles 0, 90, and 180 deg 
% and PSD versus pitch angle L1b brst data 
%
% Written by D. B. Graham
clear;
clc;
mms.db_init('local_file_db','D:/MMS/');
ic = 2; % Spacecraft number

% Tintr = irf.tint('2019-08-16T09:31:45.0Z/2019-08-16T09:32:15.0Z');
% 

Tintr = irf.tint('2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.900Z');
% Tintr = irf.tint('2020-08-03T01:45:24.000Z/2020-08-03T01:45:25.000Z');
% Tintr = irf.tint('2020-07-15T23:22:21.100Z/2020-07-15T23:22:21.510Z');
% Tintr = irf.tint('2020-07-15T23:22:21.800Z/2020-07-15T23:22:22.800Z');
% Tintr = irf.tint('2019-08-16T09:31:04.700Z/2019-08-16T09:32:04.960Z');
%% Load data

tic;
c_eval('ePDist = mms.get_data(''PDe_fpi_brst_l2'',Tintr,?);',ic)
c_eval('Bxyz=mms.get_data(''B_gsm_brst_l2'',Tintr,?);',ic);
c_eval('SCpot=mms.get_data(''V_edp_brst_l2'',Tintr,?);',ic);
toc;
ePDist = ePDist.convertto('s^3/km^6');
SCpot = SCpot.resample(ePDist);



%% 连续时间平均分布

Paddist(32,24)=0;
%% Produce a single PAD at a selected time
% for ii=1:length(ePDist.time.epoch)
time = irf_time(ePDist.time,'epochTT>utc');
for ii=1:length(time(:,1))
tint = irf_time(time(ii,:),'utc>epochTT');
% tint = irf_time('2018-08-22T15:34:33.28Z','utc>epochTT');
clear paddist;
clear thetapad;
clear energypad;
clear tintpad;
clear energypad;
[paddist,thetapad,energypad,tintpad] = mms.get_pitchangledist(ePDist,Bxyz,tint,'angles',24); 
[~,idx] = min(abs(SCpot.time-tint));
energypad = energypad-SCpot.data(idx);

Paddist=Paddist+paddist

end

Paddist_aver=Paddist/9

paddistsmooth =[];
for i=1:32
%   paddisti = (smooth(paddist(i,:),180))'
  paddisti = (smoothts(Paddist_aver(i,:),'b',1))'
  paddistsmooth(i,:)= paddisti
end
% paddistsmooth(9,:)=smoothts(paddistsmooth(9,:),'b',15);



%% Plot PAD

fn=figure;
set(fn,'Position',[200 200 600 400])
%     h(1)=axes('position',[0.08 0.12 0.3 0.80]); 
%     h(1)=axes('position',[0.30 0.10 0.30 0.80]); 
 h(1)=axes('position',[0.10 0.10 0.30 0.80]); 
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 
% 
ymin = 3*10^-2;
% ymax = ceil(max(max(log10(paddist))));
ymax = 70;
yrange = [ymin ymax];
% plot(h(1),energypad,paddist(:,1),'k',energypad,mean(paddist(:,[8:9]),2),'r',energypad,paddist(:,16),'b');
% ylabel(h(1),'f_e (s^3 km^{-6})');
% xlabel(h(1),'E (eV)')
% set(h(1),'yscale','log');
% set(h(1),'xscale','log');
% axis(h(1),[5 3e4 yrange])
% set(h(1),'xtick',[1e0 1e1 1e2 1e3 1e4 1e5])
% irf_legend(h(1),{'0 deg'},[0.91 0.92],'color','k')
% irf_legend(h(1),{'90 deg'},[0.91 0.84],'color','r')
% irf_legend(h(1),{'180 deg'},[0.91 0.76],'color','b')    

jetcolor = colormap('jet');
lll = length(jetcolor(:,1));
vcolors = floor(lll/length(energypad));
vcolors = (1:length(energypad))*vcolors;

% paddist=smooth(paddist(:,:),4);
c_eval('plot(h(1),thetapad,squeeze(paddist(?,:)),''color'',jetcolor(vcolors(?),:));',5);
hold(h(1),'on');

% c_eval('plot(h(1),thetapad,mean(paddist(?,:)),''color'',jetcolor(vcolors(?),:));', 21:32);

c_eval('plot(h(1),thetapad,squeeze(paddistsmooth(?,:)),''color'',jetcolor(vcolors(?),:));', 14:21);
hold(h(1),'off')
ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'\theta (deg.)')
set(h(1),'yscale','log');
axis(h(1),[0 180 yrange])
set(h(1),'xtick',[0 45 90 135 180])
irf_zoom(h(1),'y',yrange);
tintutc = tintpad.utc;

for i=5:length(energypad)
    a=roundn(energypad(1,:)/1000,-3);
c_eval('irf_legend(h(1),''? keV'',[1.02 1-(i-4)*0.035],''color'',jetcolor(vcolors(i),:));',a(i));
end
% a=roundn(energypad(1,:)/1000,-3);
% c_eval('irf_legend(h(1),''? keV'',[1.02 1-1*0.03],''color'',jetcolor(vcolors(1),:));',a(1));     


% title(h(1),strcat(tintutc(12:23),'UT'));
title(h(1),strcat(tintutc(12:23),'UT'));
set(gcf,'color','w');

%  set(gcf,'render','painters');
% % set(gcf,'render','opengl');
% % set(gcf,'paperpositionmode','auto')
% % figname='Figure1';
% 
% figname=['Electron_PSD_ver2' '_' num2str(ic) '_' time(ii,1:4) time(ii,6:7) time(ii,9:10) '_' time(ii,12:13) time(ii,15:16) time(ii,18:19) time(ii,21:23)];
% print(gcf, '-dpng', [figname '.png']);
% close;
%%
paddistsmooth2=paddistsmooth
save energypad.mat energypad
save thetapad.mat thetapad
save paddist1.mat paddistsmooth