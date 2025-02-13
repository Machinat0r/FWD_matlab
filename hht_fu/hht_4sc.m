%------written by Wending Fu, Nov.2023 in Beijing------------
close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
TempDir = [ParentDir,'temp/'];mkdir(TempDir);

% TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
% % % TT = '2017-06-25T05:06:50.00Z/2017-06-25T05:07:10.00Z'; % weak magnetic field whistler from xzy [0.1-10]
% TT = '2015-09-19T07:43:29.00Z/2015-09-19T07:43:32.00Z';
TT = '2016-01-07T09:34:27.00Z/2016-01-07T09:34:37.00Z'; % ECH from 10.1029/2023JA031865 [10-4096]
% TT = '2017-07-17T07:53:05.00Z/2017-07-17T07:53:06.00Z';

mode = 'B';

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
SDCFilesDownload_NAS(filenames,TempDir)

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
%% emd 
c_eval('[imf_x?, res_x?, info_x?] = emd(Bfac?(:,2), ''MaxNumIMF'', 20);',ic);
c_eval('[imf_y?, res_y?, info_y?] = emd(Bfac?(:,3), ''MaxNumIMF'', 20);',ic);
c_eval('[imf_z?, res_z?, info_z?] = emd(Bfac?(:,4), ''MaxNumIMF'', 20);',ic);
c_eval('[imf_t?, res_t?, info_t?] = emd(Bt?(:,2), ''MaxNumIMF'', 20);',ic);
% c_eval('[imf_x?, res_x?, info_x?] = vmd(Bfac?(:,2), ''NumIMF'', 10);',ic);
% c_eval('[imf_y?, res_y?, info_y?] = vmd(Bfac?(:,3), ''NumIMF'', 10);',ic);
% c_eval('[imf_z?, res_z?, info_z?] = vmd(Bfac?(:,4), ''NumIMF'', 10);',ic);
% c_eval('[imf_t?, res_t?, info_t?] = vmd(Bt?(:,2), ''NumIMF'', 10);',ic);
c_eval('imf? = imf_t?;',ic);
%% inst frequencies
% % % c_eval(['instf? = instfreq(imf?(:,1),8192,''method'',''hilbert'');'],ic);
% % % c_eval('instf? = abs(instf?);',ic);
% % % % c_eval(['instf?(instf?<=0)=0;'],ic)
% % % dspan = 128;
% % % c_eval('instf? = smooth(instf?, dspan);', ic);
% % % plot(instf1);hold on;plot(instf2);hold on;plot(instf3);hold on;plot(instf4);hold on;


lf = 1; hf = 4096;
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
panel = '1';
panelnum = 4;
for i = ic
h(i)=irf_subplot(n,1,-i);

% hht(imf_z(:,i),dfB1,'FrequencyLimits',[0,10])
% % c_eval(['[p,f,t,imfinsf,imfinse] = hht(imf?(:,',panel,'),dfB?,''FrequencyLimits'',[lf, hf]);'], i);
c_eval(['[p?,f?,t?,imfinsf,imfinse, PA?] = hht_modified_by_fwd(imf?(:,',panel,'),dfB?,''FrequencyLimits'',[lf, hf]);'], i);
% [p,f,t,imfinsf,imfinse] = hht(imf(:,i),dfB1,'FrequencyLimits',[0,4e3],'FrequencyResolution',4e3/1e3);
imfinse(imfinsf <= 0) = 0; imfinse(imfinse <= 0) = 0;
imfinsf(imfinsf <= 0) = 0;

PRange = log10([min(imfinse,[],'all'), max(imfinse,[],'all')]);

c_eval('PA? = PA?*180/pi;',i);
c_eval('irf_plot([Bfac?(:,1),PA?]);',i);

%%% HILBERT SPECTRUM plot
% % % for j = panelnum:-1:1
% % % % smooth
% % % dspan = 3;
% % % meanf = mean(imfinsf(:,j));
% % % c_eval('tspan = dspan*round(dfB?/meanf);', i);
% % % c_eval('t = smooth(t, tspan);', i);
% % % tempimfinsf = smooth(imfinsf(:,j), tspan);
% % % tempimfinse = smooth(imfinse(:,j), tspan);
% % % 
% % % 
% % % hht_plot(t, tempimfinsf, tempimfinse, 'TimeAxis', 0, 'FRange', [lf, hf], 'PRange', [-4, -2]);hold on;
% % % end

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
% % % 
set(gca,"XTickLabelRotation",0)
end
irf_timeaxis(gca,'date');
%% shift by frequency & phase angle
c_eval('[tempfbox?,tempcol?,~] = find(p?);',ic);
c_eval('fbox? = zeros(size(B1,1),1);',ic);
c_eval('fbox?(tempcol?,1) = tempfbox?;',ic);
PAbox = -180:5:180;
shift_length = zeros(size(B1,1),3);

% shift to SC1
for i = 1:size(B1,1)
    c_eval('tempid = find(PAbox>PA?(i),1);',ic);
    c_eval('PA?(i) = mean(PAbox(tempid-1:tempid));',ic);
end

for tempt = 1:size(B1,1)
    c_eval('df? = sqrt((fbox1(tempt) - fbox?).^2 + (PA1(tempt) - PA?).^2);',2:4);
    c_eval('id? = find(df? == 0);',2:4);
    c_eval('[~,id_id?] = min(abs(id? - tempt));',2:4)
    if ~isempty(id2(id_id2)), shift_length(tempt,1) = id2(id_id2); else, shift_length(tempt,1) = nan; end
    if ~isempty(id3(id_id3)), shift_length(tempt,2) = id3(id_id3); else, shift_length(tempt,2) = nan; end
    if ~isempty(id4(id_id4)), shift_length(tempt,3) = id4(id_id4); else, shift_length(tempt,3) = nan; end
end