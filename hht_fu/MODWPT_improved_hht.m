%------written by Wending Fu, Jan.2024 in Beijing------------
% Reference: http://dx.doi.org/10.1063/1.4887415
close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

% TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
% % % TT = '2017-06-25T05:06:50.00Z/2017-06-25T05:07:10.00Z'; % weak magnetic field whistler from xzy [0.1-10]
% TT = '2015-09-19T07:43:29.00Z/2015-09-19T07:43:32.00Z';
TT = '2016-01-07T09:34:27.00Z/2016-01-07T09:34:37.00Z'; % ECH from 10.1029/2023JA031865 [10-4096]
% TT = '2017-07-17T07:53:05.00Z/2017-07-17T07:53:06.00Z';

mode = 'E';
comp = 'x';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:4;
iic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% % % filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
switch mode
case 'B'
    filenames = [filenames1,filenames3];
case 'E'
    filenames = [filenames1,filenames4];
end
[filenames,~,~] = findFilenames(TT,filenames,'brst',ic);
SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)

% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);
% SDCFilesDownload_NAS(filenames_srvy,TempDir)
% % % id_flagTime = 3OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval('Bbg?_ts=mms.get_data(''B_gsm_brst'',tint,?);',ic);
c_eval('Bbg?=irf.ts2mat(Bbg?_ts);',ic);

switch mode
case 'B'

c_eval('B?_ts=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
% c_eval('B?_ts=mms.get_data(''B_gsm_brst'',tint,?);',ic);
% % % c_eval(['B?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
% % % c_eval(['B?_ts=irf_gse2gsm(B?_ts);'],ic);
c_eval('Bt?_ts=B?_ts.abs;',ic);
c_eval('B?=irf.ts2mat(B?_ts);',ic);
c_eval('B?=irf_gse2gsm(B?);',ic);
c_eval('Bt?=irf.ts2mat(Bt?_ts);',ic);

c_eval(['Bfac?=irf_convert_fac(B?,Bbg?,[1,0,0]);'],ic);

c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',ic);
c_eval('Bfac? = irf_resamp(Bfac?,B1);',ic);
c_eval('Bt? = irf_resamp(Bt?,B1);',ic);

case 'E'
% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?_ts=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?_ts);'],ic);
c_eval(['Bt?=irf.ts2mat(Et?_ts);'],ic);

c_eval(['Bfac?=irf_convert_fac(E?,Bbg?,[1,0,0]);'],ic); % call 'E' as 'B' ^_^

c_eval('dfB? =1/median(diff(E?_ts.time.epochUnix));',ic);
end
%% MODWPT
% lev = floor(log2(numel(Bfac1(:,1)))); % default, about 16 for 10 sec
lev = 4; % 2^lev frequency bands
c_eval(['[wpt_x?, packetlevels_x?, Falign_x?] = modwpt(Bfac?(:,2), lev, TimeAlign = true);'],ic);
c_eval(['[wpt_y?, packetlevels_y?, Falign_y?] = modwpt(Bfac?(:,3), lev, TimeAlign = true);'],ic);
c_eval(['[wpt_z?, packetlevels_z?, Falign_z?] = modwpt(Bfac?(:,4), lev, TimeAlign = true);'],ic);
c_eval(['[wpt_t?, packetlevels_t?, Falign_t?] = modwpt(Bt?(:,2), lev, TimeAlign = true);'],ic);
%% emd and choose IMFs
for wpt_i = 1:2^lev
c_eval('[imf_x?{wpt_i}, res_x?{wpt_i}, info_x?{wpt_i}] = emd(wpt_x?(wpt_i,:)'', ''MaxNumIMF'', 20);',ic);
c_eval('[imf_y?{wpt_i}, res_y?{wpt_i}, info_y?{wpt_i}] = emd(wpt_y?(wpt_i,:)'', ''MaxNumIMF'', 20);',ic);
c_eval('[imf_z?{wpt_i}, res_z?{wpt_i}, info_z?{wpt_i}] = emd(wpt_z?(wpt_i,:)'', ''MaxNumIMF'', 20);',ic);
c_eval('[imf_t?{wpt_i}, res_t?{wpt_i}, info_t?{wpt_i}] = emd(wpt_t?(wpt_i,:)'', ''MaxNumIMF'', 20);',ic);

% c_eval('[imf_x?, res_x?, info_x?] = vmd(Bfac?(:,2), ''NumIMF'', 10);',ic);
% c_eval('[imf_y?, res_y?, info_y?] = vmd(Bfac?(:,3), ''NumIMF'', 10);',ic);
% c_eval('[imf_z?, res_z?, info_z?] = vmd(Bfac?(:,4), ''NumIMF'', 10);',ic);
% c_eval('[imf_t?, res_t?, info_t?] = vmd(Bt?(:,2), ''NumIMF'', 10);',ic);

c_eval('corr_x?{wpt_i} = corr(imf_x?{wpt_i},Bfac?(:,2));',ic);
c_eval('corr_y?{wpt_i} = corr(imf_y?{wpt_i},Bfac?(:,3));',ic);
c_eval('corr_z?{wpt_i} = corr(imf_z?{wpt_i},Bfac?(:,4));',ic);
c_eval('corr_t?{wpt_i} = corr(imf_t?{wpt_i},Bt?(:,2));',ic);

c_eval('u_threshold_x?(wpt_i) = threshold_cal(corr_x?{wpt_i});',ic);
c_eval('u_threshold_y?(wpt_i) = threshold_cal(corr_y?{wpt_i});',ic);
c_eval('u_threshold_z?(wpt_i) = threshold_cal(corr_z?{wpt_i});',ic);
c_eval('u_threshold_t?(wpt_i) = threshold_cal(corr_t?{wpt_i});',ic);

c_eval('imf_x_phy?{wpt_i} = imf_x?{wpt_i}(:,corr_x?{wpt_i} >= u_threshold_x?(wpt_i));',ic);
c_eval('imf_y_phy?{wpt_i} = imf_y?{wpt_i}(:,corr_y?{wpt_i} >= u_threshold_y?(wpt_i));',ic);
c_eval('imf_z_phy?{wpt_i} = imf_z?{wpt_i}(:,corr_z?{wpt_i} >= u_threshold_z?(wpt_i));',ic);
c_eval('imf_t_phy?{wpt_i} = imf_t?{wpt_i}(:,corr_t?{wpt_i} >= u_threshold_t?(wpt_i));',ic);
end
%% choose component
% c_eval(['imf? = imf_',comp,'?;'],ic);
c_eval(['corr? = corr_',comp,'?;'],ic);
c_eval(['u_threshold? = u_threshold_',comp,'?;'],ic);
c_eval(['imf_phy? = imf_',comp,'_phy?;'],ic);
c_eval('imf? = [];');
for wpt_i = 1:2^lev
    c_eval('imf? = [imf?, imf_phy?{wpt_i}];')
end

%% inst frequencies
% % % c_eval(['instf? = instfreq(imf?(:,1),8192,''method'',''hilbert'');'],ic);
% % % c_eval('instf? = abs(instf?);',ic);
% % % % c_eval(['instf?(instf?<=0)=0;'],ic)
% % % dspan = 128;
% % % c_eval('instf? = smooth(instf?, dspan);', ic);
% % % plot(instf1);hold on;plot(instf2);hold on;plot(instf3);hold on;plot(instf4);hold on;


lf = 2e2; hf = 4096;
% hht(imf, dfB1, 'FrequencyLimits', [lf, hf]);
% set(gca,'yscale','log');
%% Init figure
n=length(ic);
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Plot sub IMF
for ii = ic
h(ii)=irf_subplot(n,1,-ii);
c_eval('[p1, f1, t1, imfinsf, imfinse, PhaseAngle] = hht_modified_by_fwd(imf?, dfB1, ''FrequencyLimits'', [lf, hf]);',ii);
% % c_eval(['[p,f,t,imfinsf,imfinse] = hht(imf?(:,',panel,'),dfB?,''FrequencyLimits'',[lf, hf]);'], i);
% c_eval(['[p?,f?,t?,imfinsf,imfinse, PA?] = hht_modified_by_fwd(imf?(:,',panel,'),dfB?,''FrequencyLimits'',[lf, hf]);'], i);
% [p,f,t,imfinsf,imfinse] = hht(imf(:,i),dfB1,'FrequencyLimits',[0,4e3],'FrequencyResolution',4e3/1e3);
imfinse(imfinsf <= 0) = 0; imfinse(imfinse <= 0) = 0;
imfinsf(imfinsf <= 0) = 0;

PRange = log10([min(imfinse,[],'all'), max(imfinse,[],'all')]);
hht_plot(t1, imfinsf, imfinse, 'TimeAxis', 0, 'FRange', [lf, hf], 'PRange', PRange);hold on;

% c_eval('PA? = PA?*180/pi;',i);
% c_eval('irf_plot([Bfac?(:,1),PA?]);',i);

%%% HILBERT SPECTRUM plot
% % % for j = size(imfinse,2):-1:1
% % % % smooth
% % % dspan = 3;
% % % meanf = mean(imfinsf(:,j));
% % % tspan = dspan*round(dfB1/meanf);
% % % t = smooth(t, tspan);
% % % tempimfinsf = smooth(imfinsf(:,j), tspan);
% % % tempimfinse = smooth(imfinse(:,j), tspan);
% % % 
% % % 
% % % hht_plot(t, tempimfinsf, tempimfinse, 'TimeAxis', 0, 'FRange', [lf, hf], 'PRange', PRange);hold on;
% % % end

i=1;
%%% IRF-SPECTROGRAM plot
% % % c_eval('B?_hht=struct(''t'',B?(:,1));',i);
% % % c_eval('B?_hht.f=f?;',i);%energy levels
% % % c_eval('p? = full(p?'');',i);
% % % c_eval('p?(p? == 0) = nan;',i);
% % % c_eval('B?_hht.p=p?;',i);%data matrix
% % % c_eval('B?_hht.f_label=''f [Hz]'';',i);
% % % c_eval('B?_hht.p_label={'' '',''''};',i);
% % % c_eval('[h(i), hcb8]=irf_spectrogram(h(i),B?_hht);',i);
% % % set(gca,'yscale','log');
% % % colormap(gca,jet)

set(gca,"XTickLabelRotation",0)
end
irf_timeaxis(gca,'date');
