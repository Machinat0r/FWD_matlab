%% 数据库初始化
clear;
clc;
mms.db_init('local_file_db','D:\MMS\');
% mms.db_init('local_file_db','D:\MATLAB\matlab1\work_zhang\Dispersion Relation Solver\20180416_data')
% Tint = irf.tint('2018-04-16T10:20:00.00Z/2018-04-16T11:00:00.00Z');
Tint = irf.tint('2021-07-21T13:19:57.000Z/2021-07-21T13:20:01.000Z');
ic=1;

%R
Tintlong = Tint+[-60 60];
R  = mms.get_data('R_gse',Tintlong);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?(:,1:3));',ic);

flag = 'B';
switch flag
    case 'B'
%         c_eval('Bxyz?=mms.get_data(''B_gse_fgm_srvy_l2'',Tint,?);',ic);
        c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',Tint);',ic);
        c_eval('Bscm_ts=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);
        Bscm1=irf_gse2gsm(Bscm_ts);
        Bscmfac1 = irf_convert_fac(Bscm1,Bxyz1,[1 0 0]);

        % c_eval('Bxyz?=mms.get_data(''B_gsm_fgm_srvy_l2'',Tint,?);',ic);
        Bav = irf.ts_vec_xyz(Bscmfac1.time,double(Bscmfac1.data)); % 直流 磁场平均
        % Bav = irf.ts_vec_xyz(Bxyz1.time,double(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4); % 直流 磁场平均
%         Bav=Bav.filt(0,0.05,16,5);
        c_eval('dfB? =1/median(diff(Bscmfac1.time.epochUnix));',ic);
        Bav=Bav.filt(2,0,dfB1,5);        
        
%         Bav = irf.ts_vec_xyz(Bxyz1.time,Bav);
        
        Bav=[Bav.time.epochUnix double(Bav.data)];
        Bxyz=[Bscmfac1.time.epochUnix double(Bscmfac1.data)];
%         [Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(Bxyz,Bav,1.0e-7,61440);
        [Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(Bxyz,Bav,1.0e-7,256);
%         [Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_means(Bxyz,Bav,1.0e-7,256);
    case 'E'
        c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
        Eav = irf.ts_vec_xyz(Exyz1.time,double(Exyz1.data)); % 直流 磁场平均
        c_eval('dfE? =1/median(diff(Exyz?.time.epochUnix));',ic);
        % Bav = irf.ts_vec_xyz(Bxyz1.time,double(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4); % 直流 磁场平均
        Eav=Eav.filt(2,0,dfE1,5);
        Eav=[Eav.time.epochUnix double(Eav.data)];
        Exyz=[Exyz1.time.epochUnix double(Exyz1.data)];
        [Epsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(Exyz,Eav,1.0e-7,256);
%         [Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_means(Exyz,Eav,1.0e-7,256);
end

%%
n_subplots=1;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;           
set(gcf,'PaperUnits','centimeters');
xSize = 40; ySize = 25; coef=floor(min(600/xSize,600/ySize));
xLeft = (20-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[20 20 400 800]);
colormap(jet);
%%

h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(gca,jet);
specWaveangle=struct('t',waveangle.t);
specWaveangle.f=waveangle.f;
specWaveangle.p=waveangle.p;
specWaveangle.f_label='';
specWaveangle.p_label={' ','k-B angle'};
[h(1) hcb1]=irf_spectrogram(h(1),specWaveangle,'lin');
ylabel({'f','(Hz)'},'fontsize',12);
set(gca,'ylim',[2 1e4]);
set(gca,'yscale','log');





