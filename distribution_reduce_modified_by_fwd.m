% A routine to compute and plot reduced electron distributions from FPI
% 
% Compare with: Wilder, F. D., et al. (2016), GRL, 43, 5909?5917,
% doi:10.1002/2016GL069473.
%
% The Example is fairly slow. Approx 4 min.
%
% Written by A. Johlander
%------modified by Wending Fu, Sept.21.2023 in Beijing------------
clear;
clc;
mms.db_init('local_file_db','/Volumes/172.17.190.41/Data/MMS/');

%% Set parameters and get data
% time interval
% tint = irf.tint('2019-08-02T10:19:05/2019-08-02T10:19:30');
tint=irf.tint('2019-07-19T13:46:50.000Z/2019-07-19T13:47:10.000Z');
% tint=irf.tint('2015-09-23T09:25:40.00Z/2015-09-23T09:26:00.00Z');
% times to make lines
% t1 = irf.time_array('2015-09-23T09:25:48.00Z');
% t2 = irf.time_array('2015-09-23T09:25:53.00Z');
t1 = irf.time_array('2019-07-19T13:46:51.000Z');
t2 = irf.time_array('2019-07-19T13:47:06.000Z');

ic = 1;
%%
% get distribution function
c_eval('ePDist = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint));',ic)
% c_eval('ePDist = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint));',ic)
ePDist = ePDist.tlim(tint);

% get magnetic field in DMPA (since ePDist is in DMPA)
c_eval('B = mms.get_data(''B_dmpa_fgm_brst_l2'',tint,?);',ic)

% Get the spacecraft potential (in this example, it is not so important)
c_eval('scPot = mms.get_data(''V_edp_brst_l2'',tint,?);',ic)

% remove flux from bottom two energy levels to make it more like in the
% paper
ePDist.data(:,1:2,:,:) = 0;

c_eval('diste_struct = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
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

energyspec = ones(length(diste.time),1)*energye0;
for ii = 1:length(diste.time);
    if stepTablee.data(ii),
        energyspec(ii,:) = energye1;
    end
end


% define angles
dangle = pi/16;
lengthphi = 32;

z2 = ones(lengthphi,1)*sind(thetae);
solida = dangle*dangle*z2;
allsolide = zeros(size(diste.data));

for ii = 1:length(diste.time);
    for jj=1:length(energye0);
        allsolide(ii,jj,:,:) = solida;
    end
end

distes = diste.data.*allsolide;

% Electron analysis - OMNI
PSDomni = zeros(length(diste.time),length(energye0));
for ii = 1:length(diste.time);
    disttemp = squeeze(distes(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
end
%%
% color/y-limit
clim = 10.^[-3.5,-1]; % s m^-4

% define velocity grid
% vg = linspace(-45e3,45e3,100); % km/s  
vg = ePDist.omni.depend{1}(1,:);

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 1e5;

% velocity limit in plot
vlim = 50e3; % km/s
%% get indicies for the times
it1 = interp1(ePDist.time.epochUnix,1:length(ePDist),(t1.epochUnix),'nearest');
it2 = interp1(ePDist.time.epochUnix,1:length(ePDist),(t2.epochUnix),'nearest');
if isnan(it1) && isnan(it2)
    it1=1;it2=length(ePDist.time);
end
% if just want to plot the line, ues this⬇️
ePDist.data([1:it1-1,it1+1:it2-1,it2+1:end],:,:,:) = 0;
%% Reduce distribution
tic
% reduced distribution along B
f1D = ePDist.reduce('1D',B,'vg',vg,'nMC',nMC,'scpot',scPot); 
% fOMNI = ePDist.omni.data;
toc
%% plot distribution as lines for the two selected lines
% matlab colors
col = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980];

% initiate figure
hca = irf_plot(1,'newfigure');
hold(hca,'on')
% plot the lines
plot(hca,f1D(it1).depend{1},f1D(it1).data,'Color',col(2,:),'linewidth',2)
plot(hca,f1D(it2).depend{1},f1D(it2).data,'Color',col(1,:),'linewidth',2)
hca.YScale = 'log';
hca.YLim = [10^(-2.5),10^(-1)];
% legends show time centers
hca.ColorOrder = flipud(col(1:2,:)); % set colr order for legends
irf_legend(hca,{ePDist(it1).time.toUtc;ePDist(it2).time.toUtc;},[0.98,0.98])
hca.XLim = [min(vg),2e4];
%labels
xlabel(hca,'V_{||} [km/s]')
ylabel(hca,'F [s m^-^4]')

