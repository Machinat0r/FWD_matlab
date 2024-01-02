%------written by Wending Fu, Nov.2023 in Beijing------------
close all
clear;clc

global ParentDir 
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
TempDir = [ParentDir,'temp/'];mkdir(TempDir);

% TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
% TT = '2017-06-25T05:06:50.00Z/2017-06-25T05:07:10.00Z';
% TT = '2015-09-19T07:43:29.00Z/2015-09-19T07:43:32.00Z';
TT = '2016-01-07T09:34:27.00Z/2016-01-07T09:34:37.00Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% % % filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames = [filenames1,filenames3, filenames4];
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
SDCFilesDownload_NAS(filenames,TempDir)

% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);
% SDCFilesDownload_NAS(filenames_srvy,TempDir)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval('B1_ts=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
% c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],iic);
% % % c_eval(['B?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
% % % c_eval(['B?_ts=irf_gse2gsm(B?_ts);'],ic);
Bt1_ts=B1_ts.abs;
B1=irf.ts2mat(B1_ts);
B1=irf_gse2gsm(B1);
Bt1=irf.ts2mat(Bt1_ts);

dfB1 = 1/median(diff(B1_ts.time.epochUnix));
%% emd 
[imf_x, res_x, info_x] = emd(B1(:,2), MaxNumIMF=20);
[imf_y, res_y, info_y] = emd(B1(:,3), MaxNumIMF=20);
[imf_z, res_z, info_z] = emd(B1(:,4), MaxNumIMF=20);
[imf_t, res_t, info_t] = emd(Bt1(:,2), MaxNumIMF=20);
imf = imf_t;

lf = 500; hf = 3000;
%% Init figure
n=size(imf,2);
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Plot sub IMF
for i = 1:n
f1(i)=irf_subplot(n,1,-i);
% hht(imf_z(:,i),dfB1,'FrequencyLimits',[0,10])
[p,f,t,imfinsf,imfinse] = hht(imf(:,i),dfB1,'FrequencyLimits',[lf, hf]);
% [p,f,t,imfinsf,imfinse] = hht(imf(:,i),dfB1,'FrequencyLimits',[0,4e3],'FrequencyResolution',4e3/1e3);

imfinse(imfinsf <= 0) = 0; imfinse(imfinse <= 0) = 0;
imfinsf(imfinsf <= 0) = 0;

PRange = log10([min(imfinse,[],'all'), max(imfinse,[],'all')]);

%%% HILBERT SPECTRUM plot
% smooth
dspan = 3;
meanf = mean(imfinsf);
c_eval('tspan = dspan*round(dfB?/meanf);', ic);
c_eval('t = smooth(t, tspan);', ic);
tempimfinsf = smooth(imfinsf, tspan);
tempimfinse = smooth(imfinse, tspan);
hht_plot(t, imfinsf, imfinse, 'TimeAxis', 0, 'FRange', [lf, hf]);


%%% IRF-SPECTROGRAM plot
% B_hht=struct('t',B1(:,1));
% B_hht.f=f;%energy levels
% p = full(p'); p(p == 0) = nan;
% B_hht.p=p;%data matrix
% B_hht.f_label='f [Hz]';
% B_hht.p_label={' ',''};
% [h(i), hcb8]=irf_spectrogram(h(i),B_hht);
% set(gca,'yscale','log');
% colormap(gca,jet)
end
set(gca,"XTickLabelRotation",0)
%% Plot total IMF
% % % figure(2)
% % % for i = 1:n
% % % % hht(imf_z(:,i),dfB1,'FrequencyLimits',[0,10])
% % % [p,f,t,imfinsf,imfinse] = hht(imf(:,n+1-i),dfB1,'FrequencyLimits',[lf, hf]);
% % % % [p,f,t,imfinsf,imfinse] = hht(imf(:,i),dfB1,'FrequencyLimits',[0,4e3],'FrequencyResolution',4e3/1e3);
% % % 
% % % %%% HILBERT SPECTRUM plot
% % % hht_plot(t, imfinsf, imfinse, 'EdgeAlpha', 0.2, 'FRange', [lf, hf]);
% % % hold on;
% % % end
%% Wavelet transforms
nf = 400;
Bwavelet1 = irf_wavelet(Bt1,'returnpower',1,'cutedge',1,'w',5.36,'nf',400,'f',[lf hf]);
% Bwavelet1 = irf_wavelet(B1,'returnpower',1,'cutedge',1,'w',5.36*2,'linear', [500:2.5:3000]);
             
%compress wavelet transform data 10 point average
nc = 20;
idx1 = [nc/2:nc:length(Bwavelet1.t)-nc/2];
Bwavelettimes1 = Bwavelet1.t(idx1);
Bwaveletx1 = zeros(length(idx1),nf);
Bwavelety1 = zeros(length(idx1),nf);
Bwaveletz1 = zeros(length(idx1),nf);
for ii = [1:length(idx1)]
    Bwaveletx1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,1}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
    % Bwavelety1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,2}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
    % Bwaveletz1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,3}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
end

specB1=struct('t',Bwavelettimes1);
specB1.f=Bwavelet1.f;
% specB1.p=Bwaveletx1+Bwavelety1+Bwaveletz1;
specB1.p=Bwaveletx1;
% specB1.p=Bwaveletz1;
specB1.f_label='';
specB1.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};
figure(3);
f3(1)=irf_subplot(1,1,-1);
[f3(1), ~] = irf_spectrogram(f3(1),specB1);
set(gca,'yscale','log');
colormap jet
%% wavelet for each imf
f4 = figure(4);
f4(1)=irf_subplot(1,1,-1);
% [mra,cfs,wfb,info] = ewt(Bt1(:,2));

Bwavelet1 = irf_wavelet([B1(:,1), imf],'returnpower',1,'cutedge',1,'w',5.36*2,'nf',nf,'f',[lf hf]);
% Bwavelet1 = irf_wavelet([B1(:,1), mra],'returnpower',1,'cutedge',1,'w',5.36*2,'nf',nf,'f',[lf hf]);
             
%compress wavelet transform data 10 point average
nc = 20;
idx1 = [nc/2:nc:length(Bwavelet1.t)-nc/2];
Bwavelettimes1 = Bwavelet1.t(idx1);
Bwaveletx1 = zeros(length(idx1),nf,n);
for ii = [1:length(idx1)]
    c_eval('Bwaveletx1(ii,:,?) = squeeze(irf.nanmean(Bwavelet1.p{1,?}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));',1:n);
end

specB1=struct('t',Bwavelettimes1);
specB1.f=Bwavelet1.f;
specB1.p= sum(Bwaveletx1,3);
specB1.f_label='';
specB1.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};

[f4(1), ~] = irf_spectrogram(f4(1),specB1);
set(gca,'yscale','log');
colormap jet
caxis([-9.75,-4.15]);