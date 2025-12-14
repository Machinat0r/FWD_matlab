clear;clc;close all
% cd 'Z:\SPART-WORK\Data\MMS\'
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
% ParentDir = 'Z:\SPART-WORK\Data\MMS\';
OutputDir = '/Users/fwd/Documents/MATLAB/MMS/FEEPS_srvy_plot/';
TT = '2015-10-07T01:50:26.00Z/2015-10-07T01:56:26.00Z';
PlotTint1=irf.tint(TT);
ic=1;

  units = irf_units; %单位
  markr = mms.get_data('R_gsm',PlotTint1,1);
  RE = units.RE;
  markpos1=markr.data(1,1)/RE*1000
  markpos2=markr.data(1,2)/RE*1000
  markpos3=markr.data(1,3)/RE*1000
%% load B
c_eval(['B?_ts=mms.get_data(''B_gsm_srvy'',PlotTint?,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
 % c_eval(['B?=irf_gse2gsm(B?);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);

%% R
units = irf_units; % 标准物理常数
% Pos = mms.get_data('R_gsm',tint);
% c_eval('R? = Pos.gsmR?;',iic)
% c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];',iic)
% c_eval('R? = irf_resamp(R?,Bt?);',iic)

eEcorr = [14.0, -1.0, -3.0, -3.0];
iEcorr = [0.0, 0.0, 0.0, 0.0];
eGfact = [1.0, 1.0, 1.0, 1.0];
iGfact = [0.84, 1.0, 1.0, 1.0];
specie = 'electron';
% csv_file = '20170531'; csv_dir = [cd, '/sun/'];
%% load FEEPS data
sensor_id = [3:5,11:12];
c_eval('e?_Tit! = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_top_intensity_sensorid_!'',PlotTint?);', ic, sensor_id);
c_eval('e?_Bit! = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_bottom_intensity_sensorid_!'',PlotTint?);', ic, sensor_id);

c_eval('ePA? = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_pitch_angle'',PlotTint?);', ic);
% % % c_eval('eLow = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_lower_bound'',tint);',ic);
% % % c_eval('eUp = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_upper_bound'',tint);',ic);
% % % c_eval('energies? = (eLow.data + eUp.data)/2. + eEcorr(?);',ic)
c_eval('Flux_e_feeps = mms.get_data(''Omnifluxelectron_epd_feeps_srvy_l2'',PlotTint?,?);',ic);
c_eval('energies? = Flux_e_feeps.depend{1, 1}(1,:);',ic);

%% calculate Omni flux
% sensors = [1:5, 9:12];
sensors = [3:5, 11:12];
nSensors = length(sensors);
c_eval('dTmp= e?_Bit3;', ic)
% omniD = NaN([size(dTmp.data) nSensors*2]);
for iSen = 1:nSensors
    c_eval(['omniD?(:,:,iSen) = ' specie(1) '?_Tit!.data;'...
        'omniD?(:,:,nSensors+iSen) = ' specie(1) '?_Bit!.data;'],ic,sensors(iSen))
end
c_eval([specie(1) 'Omni? = dTmp; ' specie(1) 'Omni?.data =' ...
    'mean(double(omniD?),3,''omitnan'')*' specie(1) 'Gfact(?);'],ic)

%% calculate PAD 
sensors = [3:5,11:12];
nSensors = length(sensors);
bin_size = 16.3636; % default 
n_pabins = round(180./bin_size);
pa_bins = 180*(0:1:n_pabins)./n_pabins; 
pa_label = 180*(0:1:n_pabins-1)/n_pabins+bin_size/2;
dAngResp = 21.4; %default
delta_pa = (pa_bins(2)-pa_bins(1))/2.0;

c_eval('dTmp= e?_Bit3;', ic)
c_eval('[sensor?_PA,sensor?_PA_id] = sort(ePA?.data,2,''ascend'');',ic)

PAD_raw_Data = NaN([size(dTmp.data) nSensors*2]);
Size_sp = size(PAD_raw_Data);
for it = 1:Size_sp(1)
    c_eval('PAD_Data?(it,:,:) = omniD?(it,:,sensor?_PA_id(it,:));',ic)
end

c_eval('PAD_spec_data?=NaN([size(dTmp.data) n_pabins]);',ic)
 
for it = 1:Size_sp(1)
    c_eval('PAD_mat_tmp = squeeze(PAD_Data?(it,:,:));',ic)
    for ipa = 1:n_pabins
        %     ind = where((dpa[it,*] + dAngResp ge pa_label[ipa]-delta_pa) and (dpa[it,*] - dAngResp lt pa_label[ipa]+delta_pa))
        %     c_eval('ind = find((((sensor?_PA(it,:)+dAngResp) > (pa_label(ipa)-delta_pa)) .* ((sensor?_PA(it,:)-dAngResp)<(pa_label(ipa)+delta_pa)))==1);',ic)
        c_eval('ind = find(((sensor?_PA(it,:)+dAngResp) > (pa_label(ipa)-delta_pa)) & ((sensor?_PA(it,:)-dAngResp)<(pa_label(ipa)+delta_pa)));',ic)
        if ~isempty(ind)
        c_eval('PAD_spec_data?(it,:,ipa)= mean(PAD_Data?(it,:,ind),3,''omitnan'') ;',ic)   
        end
               
    end
end

% find plot range
c_eval('idx_findout = find(eOmni?.time.epochUnix > irf_time(PlotTint1(1),''epochtt>epoch''));',ic)
idx_start = idx_findout(1);
c_eval('idx_findout = find(eOmni?.time.epochUnix < irf_time(PlotTint1(2),''epochtt>epoch''));',ic)
idx_end = idx_findout(end);


c_eval('PAD_size = size(PAD_Data?);', ic)
c_eval('PAD_resamp? = PAD_Data?;', ic)
anglea = linspace(5,175,18)';

for ix = idx_start: idx_end
    for i_cha = 1:PAD_size(2)
        c_eval('flux_tmp = squeeze(PAD_Data?(ix,i_cha,:));',ic)
        c_eval('ang_tmp = sensor?_PA(ix,:);',ic)
%         c_eval('anglea = linspace(sensor?_PA(ix,1),sensor?_PA(ix,end),18)'';',ic)
        
        tsin_tmp = timeseries(flux_tmp,ang_tmp);        
        tsout_tmp = resample(tsin_tmp,anglea);
        c_eval('PAD_resamp?(ix,i_cha,:) = tsout_tmp.data'';',ic)
    end  
end
c_eval('PAD_Data? = PAD_resamp?;', ic)

%2 PAD
channel = [14:16];
% anglea = linspace(5,175,18)';
% c_eval('spe_PAD?_c! = struct(''t'', eOmni?.time.epochUnix);',ic,channel)
% c_eval('spe_PAD?_c!.p = squeeze(PAD_Data?(:,!,:)); spe_PAD?_c!.p_label = {[''log('' eOmni?.units '')'']};',ic,channel)
% c_eval('spe_PAD?_c!.f = anglea; spe_PAD?_c!.f_label = {''Energy''};',ic,channel)

anglea = pa_label';
c_eval('spe_PAD?_c! = struct(''t'', eOmni?.time.epochUnix);',ic,channel)
c_eval('spe_PAD?_c!.p = squeeze(PAD_spec_data?(:,!,:)); spe_PAD?_c!.p_label = {['' '']};',ic,channel)
c_eval('spe_PAD?_c!.f = anglea; spe_PAD?_c!.f_label = {''Energy''};',ic,channel)
%% Init figure
% ic = 1;
n=5;
fonts=8;
h1=irf_plot(n,'newfigure');
i=1;
set(0,'DefaultAxesFontSize',8); %字号
set(0,'DefaultLineLineWidth', 0.5); %线宽
% fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80;
coef=floor(min(800/xSize,800/ySize)); %计算适应屏幕的最大缩放因子，并取其中较小的一个，确保图形在显示时不会超出屏幕范围
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef]) %位置和大小
%% B plot
h(i)=irf_subplot(n,1,-i);%创建子图
c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off; % 零参考线
grid off;
c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);%自动设置Y轴范围 
% set(gca,'Ylim',[-10 15], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% set(gca,'Ylim',[-16 30])
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]); % 设置颜色顺序
set(gca,'xtick',[])
set(gca,'fontsize',9)
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);%在指定位置创建图例
ylabel('B [nT]','fontsize',8); % Y轴标签
xlabel(' ');
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;

% % % %% N plot
% % % h(i)=irf_subplot(n,1,-i);
% % % 
% % % %滤波
% % % %     irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
% % % %     irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;
% % % 
% % % %非滤波
% % % % c_eval('irf_plot(h(i),Ne?,''Linewidth'',0.75,''color'',''b'');',ic);
% % % % % c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
% % % % c_eval('irf_plot(h(i),Ni?,''Linewidth'',0.75,''color'',''g'');',ic);
% % % c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;%绘制电子密度数据
% % % c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % if ~isnan(Ne1(:,2));
% % % c_eval("set(gca,'Ylim',[max([0 min([min(Ne?(:,2)) min(Ni?(:,2))])]) max([max(Ne?(:,2)) max(Ni?(:,2))])]);",ic);
% % % else
% % % end
% % % % set(gca,'Ylim',[0.1 0.4])
% % % % set(gca,'Ylim',[0.15 0.45], 'ytick',[0.1 0.2 0.3 0.4],'fontsize',9);
% % % % pos1=get(h(1),'pos');
% % % %  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
% % % %  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
% % %   set(gca,'ColorOrder',[[0 0 1];[0 1 0]]); % 图例颜色顺序
% % %   set(gca,'xtick',[])
% % %   set(gca,'fontsize',9)
% % %  irf_legend(gca,{'Ne','Ni'},[0.97 0.92]); % 图例位置
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
% % % ylabel('N [cm^{-3}]','fontsize',8);
% % % i=i+1;
% % % 
% % % 
% % % %% Vi plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % if ~isnan(Ne1(:,2));
% % % c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+30],'fontsize',9);",ic); % 自动调Y轴显示范围else
% % % else
% % % end
% % % % set(gca,'Ylim',[-100 200], 'ytick',[0 200 400]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
% % % set(gca,'xtick',[])
% % % ylabel('Vi [km/s]','fontsize',8);
% % % i=i+1;

%% plot FEEPS electron flux spectrom
h(i)=irf_subplot(n,1,-i);
% colormap(h(i),jet)
c_eval('speOmni? = struct(''t'', eOmni?.time.epochUnix);',ic)
c_eval('speOmni?.p = eOmni?.data; speOmni?.p_label = {[''log('' eOmni?.units '')'']};',ic)
c_eval('speOmni?.f = double(energies?); speOmni?.f_label = {''Energy''};',ic)
c_eval('spec_pl1 = speOmni?;',ic)
[a5,b5]=irf_spectrogram(h(i),spec_pl1,'donotfitcolorbarlabel');
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

% irf_legend(h1(i_plot),'e',[0.99 0.98],'color','w','fontsize',fonts);
grid off;
set(h(i),'yscale','log');
set(h(i),'Ylim',[2e5 6e5]);
% set(h(i),'ytick',[1e4 2e5 4e5 6e5 1e6]);
set(h(i),'ytick',[2e5 4e5]);
% ytickformat('%g');  % 通用数值格式
h(i).YRuler.Exponent = 0;
ytickformat('%.0e'); 
% set(h(i),'ytick',[1e4 1e5]);
set(gca,'fontsize',9)
ylabel(h(i),{'\color{black}E(eV)','\color{black}'},'Interpreter','tex','fontsize',8);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a5,'position'); Pos_s = get(h(i),'position');
set(a5,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a5,'fontsize',8); set(b5,'fontsize',8);
% colormap(h(i),jet)
i=i+1;
caxis(gca,[-0.5, 4])   

% % % %% plot e energy spectrom 绘制电子全向能谱
% % % c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',PlotTint?,?);',ic);
% % % if    isempty(Ne1_ts)
% % %    h(i)=irf_subplot(n,1,-i);
% % %    c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;%绘制电子密度数据
% % % c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % i=i+1;
% % % else
% % % c_eval('Omni_flux_e?(:,1)=irf_time(energy_e?.DEPEND_0.data,''ttns>epoch'');',ic)
% % % c_eval('Omni_flux_e?(:,2:33)=energy_e?.data;;',ic)
% % % c_eval('channel?=transpose(energy_e?.DEPEND_1.data(1,1:32));',ic)
% % % 
% % % c_eval('spec_fpi_e?=struct(''t'',Omni_flux_e?(:,1));',ic)
% % % c_eval('spec_fpi_e?.f=transpose(energy_e?.DEPEND_1.data(1,1:32));',ic) %energy levels
% % % c_eval('spec_fpi_e?.p=Omni_flux_e?(:,2:33);',ic) %data matrix
% % % c_eval('spec_fpi_e?.f_label=''Energy'';',ic)
% % % c_eval('spec_fpi_e?.p_label={'' '',''(cm^2 s sr keV)^{-1}''};',ic)
% % % 
% % % c_eval('spec_pl2 = spec_fpi_e?;',ic)
% % % 
% % % h(i)=irf_subplot(n,1,-i);
% % % % colormap(h(i),jet)
% % % [a5,b5]=irf_spectrogram(h(i),spec_pl2,'donotfitcolorbarlabel');
% % % 
% % % % hold on;
% % % % irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
% % % grid off;
% % % set(h(i),'yscale','log');  % 对数刻度
% % % set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
% % % ylabel('Ee(ev)','fontsize',8)
% % % set(gca,'Ylim',[1.5e1 3e4]);
% % % % caxis(gca,[7.1 7.4])
% % % % set(gca,'xtick',[])
% % % % irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
% % % xwidth=0.02; ywidth=0.154;
% % % Pos_c = get(a5,'position'); Pos_s = get(h(i),'position');
% % % set(a5,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
% % % set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
% % % set(a5,'fontsize',8); set(b5,'fontsize',8);
% % % % colormap(h(i),jet)
% % % 
% % % i=i+1;
% % % end

%% PAD channel14
h(i)=irf_subplot(n,1,-i);
c_eval('spec_pl2 = spe_PAD?_c!;',ic, 14)
[a6,b6]=irf_spectrogram(h(i),spec_pl2,'donotfitcolorbarlabel');
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

% irf_legend(h1(i_plot),'e',[0.99 0.98],'color','w','fontsize',fonts);
irf_legend(h(i),'450keV',[0.02 0.98],'color','k','fontsize',13);
grid off;
set(h(i),'yscale','lin');
set(h(i),'ytick',[0:90:180]);
ylabel(h(i),{'PA(deg)'},'Interpreter','tex','fontsize',fonts);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a6,'position'); Pos_s = get(h(i-2),'position');
set(a6,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b6,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a6,'fontsize',fonts); set(b6,'fontsize',fonts);
colormap(h(i),jet)
caxis(h(i),[-1, 2])
i=i+1;
%% PAD channel15
h(i)=irf_subplot(n,1,-i);
c_eval('spe_PAD?_c!.p_label = {[''log('' eOmni?.units '')'']};',1,15);
c_eval('spec_pl2 = spe_PAD?_c!;',ic, 15)
[a6,b6]=irf_spectrogram(h(i),spec_pl2,'donotfitcolorbarlabel');
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

% irf_legend(h1(i_plot),'e',[0.99 0.98],'color','w','fontsize',fonts);
irf_legend(h(i),'523keV',[0.02 0.98],'color','k','fontsize',13);
grid off;
set(h(i),'yscale','lin');
set(h(i),'ytick',[0 90 180]);
ylabel(h(i),{'PA(deg)'},'Interpreter','tex','fontsize',fonts);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a6,'position'); Pos_s = get(h(i-1),'position');
set(a6,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b6,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a6,'fontsize',fonts); set(b6,'fontsize',fonts);
colormap(h(i),jet)
caxis(h(i),[-1, 2])
i=i+1;
%% PAD channel16
h(i)=irf_subplot(n,1,-i);
c_eval('spec_pl2 = spe_PAD?_c!;',ic, 16)
[a6,b6]=irf_spectrogram(h(i),spec_pl2,'donotfitcolorbarlabel');
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

% irf_legend(h1(i_plot),'e',[0.99 0.98],'color','w','fontsize',fonts);
irf_legend(h(i),'575keV',[0.02 0.98],'color','k','fontsize',13)
grid off;
set(h(i),'yscale','lin');
set(h(i),'ytick',[0 90 180]);
ylabel(h(i),{'PA(deg)'},'Interpreter','tex','fontsize',fonts);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a6,'position'); Pos_s = get(h(i-1),'position');
set(a6,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b6,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a6,'fontsize',fonts); set(b6,'fontsize',fonts);
colormap(h(i),jet)
caxis(h(i),[-1, 2])
i=i+1;

%%  出图保存部分
title(h(1),sprintf('X = %.2f  Y = %.2f  Z = %.2f', markpos1,markpos2, markpos3),'fontsize',10);
set(h(1:n),'fontsize',12);
irf_zoom(PlotTint1,'x',h(1:n));
    % irf_adjust_panel_position;
    % %   irf_plot_axis_align(h)
    irf_plot_axis_align(h)
colormap(jet) % 全局颜色映射设置
set(gca,"XTickLabelRotation",0) % 确保X轴刻度标签方向水平
set(gcf,'render','painters'); % 设置图形矢量渲染
set(gcf,'paperpositionmode','auto') % 打印或导出文件的尺寸与屏幕上显示的图形尺寸完全一致

% str_time = irf_time(flagTime,'epoch>utc');
fig_dir0 = strcat( OutputDir,'OverviewFig4\');
% fig_dir = strcat(fig_dir0, '\',str_time(1:4),'\',str_time(6:7),'\',str_time(9:10));
fig_dir = fig_dir0;
if ~exist( fig_dir, 'dir')
    mkdir(fig_dir);
end
fig_name = strcat(t_center1(1:4),t_center1(6:7),t_center1(9:10),t_center1(12:13),t_center1(15:16),t_center1(18:19));
fig_path_name = strcat(fig_dir,'\',fig_name,'.png');
% fig_path_name = strcat(fig_dir,'\',fig_name,'.pdf');

% set(gcf, 'InvertHardCopy', 'off');
% set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% if(eps)
% print('-depsc','-painters','-r150',fig_path_name);
% else
print('-dpng','-painters','-r300',fig_path_name);
% end
close all
% catch
%     writematrix([t_center1(1:end-1),'的数据导入2出现问题需要手动下载'],[OutputDir,'errorlog4.txt'],...
%         'WriteMode','append','Encoding','UTF-8')
%     continue
% end