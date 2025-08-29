clear
clc
thresold=0.5; 
mms.db_init('local_file_db','/Volumes/SPART-NAS/Data/MMS/')
% tint=irf.tint('2017-07-26T07:10:00Z/2017-07-26T07:50:00Z');
% % % tint=irf.tint('2017-07-26T07:25:00Z/2017-07-26T07:31:00Z');
% % % tint3=irf.tint('2017-07-26T07:27:29Z/2017-07-26T07:27:39Z');

tint=irf.tint('2019-08-05T16:24:00.000Z/2019-08-05T16:25:00.000Z');
tint3 = tint;
% tint3=irf.tint('2017-07-26T07:37:15Z/2017-07-26T07:37:35Z');
% from SPEDAS: mms_feeps_omni.pro - correction factors per SC
eEcorr = [14.0, -1.0, -3.0, -3.0];
iEcorr = [0.0, 0.0, 0.0, 0.0];
eGfact = [1.0, 1.0, 1.0, 1.0];
iGfact = [0.84, 1.0, 1.0, 1.0];
specie = 'electron';
% csv_file = '20170531'; csv_dir = [cd, '/sun/'];
%% load data
ic = 1;
% c_eval('spinsect_num? = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_spinsectnum'',tint);', ic);
% c_eval('mask_csv? = [csv_dir,''MMS'',num2str(ic),''_FEEPS_ContaminatedSectors_'',csv_file,''.csv'']; ',ic)
% c_eval('mask_sec? = csvread(mask_csv?);',ic)

sensor_id = [1:5,9:12];
c_eval('e?_Tit! = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_!'',tint);', ic, sensor_id);
c_eval('e?_Bit! = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_!'',tint);', ic, sensor_id);

c_eval('mask_sun_con_Tit! = mms.db_get_variable(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_top_sun_contamination_sensorid_!'',tint);', ic, sensor_id);
c_eval('mask_sun_con_Bit! = mms.db_get_variable(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_bottom_sun_contamination_sensorid_!'',tint);', ic, sensor_id);
c_eval('mask_bad_eye_Tit! = mms.db_get_variable(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_top_sector_mask_sensorid_!'',tint);', ic, sensor_id);
c_eval('mask_bad_eye_Bit! = mms.db_get_variable(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_bottom_sector_mask_sensorid_!'',tint);', ic, sensor_id);

% c_eval('e?_Tit!.data(e?_Tit!.data == 0) = NaN;',ic, sensor_id)
% c_eval('e?_Bit!.data(e?_Bit!.data == 0) = NaN;',ic, sensor_id)

% c_eval('e?_Tit!.data(e?_Tit!.data > 1e6) = NaN;',ic, sensor_id)
% c_eval('e?_Bit!.data(e?_Bit!.data > 1e6) = NaN;',ic, sensor_id)
% c_eval('daT_check?(!) = mean(e?_Tit!.data(:,1),''omitnan'');',ic, sensor_id)
% c_eval('daB_check?(!) = mean(e?_Bit!.data(:,1),''omitnan'');',ic, sensor_id)

c_eval('ePA? = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_pitch_angle'',tint);', ic);
% % % c_eval('eLow = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_lower_bound'',tint);',ic);
% % % c_eval('eUp = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_upper_bound'',tint);',ic);
% % % c_eval('energies? = (eLow.data + eUp.data)/2. + eEcorr(?);',ic)
c_eval('Flux_e_feeps = mms.get_data(''Omnifluxelectron_epd_feeps_brst_l2'',tint,?);',ic);
c_eval('energies? = Flux_e_feeps.depend{1, 1}(1,:);',ic);
%% clean up 
% First, remove bad eyes
% MMS1 
c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN; e?_Tit1.data(:,:) = NaN;',ic);

% MMS2
% % % c_eval('e?_Tit5.data(:,:) = NaN; e?_Tit7.data(:,:) = NaN; e?_Tit12.data(:,:) = NaN;',ic);
% % % c_eval('e?_Bit7.data(:,:) = NaN;',ic);
% MMS3
% c_eval('e?_Tit2.data(:,:) = NaN; e?_Tit12.data(:,:) = NaN;',ic);
% c_eval('e?_Bit2.data(:,:) = NaN; e?_Bit5.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN;',ic);

% MMS4
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit2.data(:,:) = NaN; e?_Tit7.data(:,:) = NaN;',ic);
% c_eval('e?_Bit2.data(:,:) = NaN; e?_Bit4.data(:,:) = NaN; e?_Bit5.data(:,:) = NaN;',ic);
% c_eval('e?_Bit10.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN;',ic);

% Second, remove bad channels
% channel 0
% MMS1
sensor_id = [2, 5];
c_eval('e?_Tit!.data(:,1) = NaN;',ic,sensor_id)
sensor_id = [2, 3, 4, 5, 6, 8, 9, 11, 12];
c_eval('e?_Bit!.data(:,1) = NaN;',ic, sensor_id)
% MMS2, 3 & 4
% % % sensor_id = [1, 2, 3, 4, 5, 9, 10, 11, 12];
% % % c_eval('e?_Tit!.data(:,1) = NaN;',ic,sensor_id)
% % % sensor_id = [1, 2, 3, 4, 5, 9, 10, 11, 12];
% % % c_eval('e?_Bit!.data(:,1) = NaN;',ic, sensor_id)
% channel 1
% % MMS1
c_eval('e?_Tit6.data(:,2) = NaN;',ic);
c_eval('e?_Bit6.data(:,2) = NaN; e?_Bit11.data(:,2) = NaN;',ic);
% MMS2
% c_eval('e?_Tit8.data(:,2) = NaN;',ic);
% c_eval('e?_Bit7.data(:,2) = NaN; e?_Bit12.data(:,2) = NaN;',ic);
% % MMS3
% c_eval('e?_Tit1.data(:,2) = NaN;e?_Tit6.data(:,2) = NaN;e?_Tit7.data(:,2) = NaN;',ic);
% c_eval('e?_Bit6.data(:,2) = NaN; e?_Bit7.data(:,2) = NaN;',ic);
% % MMS4
% c_eval('e?_Tit1.data(:,2) = NaN;e?_Tit6.data(:,2) = NaN;',ic);
% c_eval('e?_Bit6.data(:,2) = NaN; e?_Bit7.data(:,2) = NaN;e?_Bit8.data(:,2) = NaN; e?_Bit0.data(:,2) = NaN;',ic);
%% sun contamination
% % % for isen = [1:5, 9:12]
% % %     c_eval('mask_tmp_T = mask_sec?(:,isen);',ic,isen)
% % %     idx_tmp1 = find (mask_tmp_T == 1);
% % %     idx_tt = [];
% % %     for aaa = 1:length(idx_tmp1)
% % %         c_eval('idx_ttt = find(double(spinsect_num?.data) == idx_tmp1(aaa)-1);',ic)
% % %         idx_tt = [idx_tt;idx_ttt];
% % %     end
% % %     c_eval('A = zeros(size(spinsect_num2.data));',ic)
% % %     A(idx_tt) = 1;
% % %     idx_tmp2 = logical(A);
% % %     c_eval('e?_Tit!.data(idx_tmp2,:) = NaN;',ic,isen);
% % % 
% % %     c_eval('mask_tmp_B = mask_sec?(:,isen+12);',ic,isen)
% % %     idx_tmp1 = find (mask_tmp_B == 1);
% % %     idx_tt = [];
% % %     for aaa = 1:length(idx_tmp1)
% % %         c_eval('idx_ttt = find(double(spinsect_num?.data) == idx_tmp1(aaa)-1);',ic)
% % %         idx_tt = [idx_tt;idx_ttt];
% % %     end
% % %     A = zeros(size(spinsect_num2.data));
% % %     c_eval('A(idx_tt) = 1;',ic)
% % %     idx_tmp2 = logical(A);
% % %     c_eval('e?_Bit!.data(idx_tmp2,:) = NaN;',ic,isen);
% % % end

%% calculate Omni flux
sensors = [1:5, 9:12];
nSensors = length(sensors);
c_eval('dTmp= e?_Bit2;', ic)
% omniD = NaN([size(dTmp.data) nSensors*2]);
for iSen = 1:nSensors
    c_eval(['omniD?(:,:,iSen) = ' specie(1) '?_Tit!.data;'...
        'omniD?(:,:,nSensors+iSen) = ' specie(1) '?_Bit!.data;'],ic,sensors(iSen))
end
c_eval([specie(1) 'Omni? = dTmp; ' specie(1) 'Omni?.data =' ...
    'mean(double(omniD?),3,''omitnan'')*' specie(1) 'Gfact(ic);'],ic)

%% calculate PAD 
sensors = [1:5, 9:12];
nSensors = length(sensors);
bin_size = 16.3636; % default 
n_pabins = round(180./bin_size);
pa_bins = 180*(0:1:n_pabins)./n_pabins; 
pa_label = 180*(0:1:n_pabins-1)/n_pabins+bin_size/2;
dAngResp = 21.4; %default
delta_pa = (pa_bins(2)-pa_bins(1))/2.0;

c_eval('dTmp= e?_Bit2;', ic)
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
c_eval('idx_findout = find(eOmni?.time.epochUnix > irf_time(tint3(1),''epochtt>epoch''));',ic)
idx_start = idx_findout(1);
c_eval('idx_findout = find(eOmni?.time.epochUnix < irf_time(tint3(2),''epochtt>epoch''));',ic)
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

%%
%% make spectrum 
%1 omni
% eval(['speOmni = struct(''t'', ' specie(1) 'Omni.time.epochUnix);'...
%     'speOmni.p = ' specie(1) 'Omni.data;'...
%     'speOmni.p_label = {[''log('' ' specie(1) 'Omni.units '')'']};']);
% 
% speOmni.f_label = {'Energy'};
% speOmni.f = double(energies1);

c_eval('speOmni? = struct(''t'', eOmni?.time.epochUnix);',ic)
c_eval('speOmni?.p = eOmni?.data; speOmni?.p_label = {[''log('' eOmni?.units '')'']};',ic)
c_eval('speOmni?.f = double(energies?); speOmni?.f_label = {''Energy''};',ic)

%2 PAD
channel = [1:16];
% anglea = linspace(5,175,18)';
% c_eval('spe_PAD?_c! = struct(''t'', eOmni?.time.epochUnix);',ic,channel)
% c_eval('spe_PAD?_c!.p = squeeze(PAD_Data?(:,!,:)); spe_PAD?_c!.p_label = {[''log('' eOmni?.units '')'']};',ic,channel)
% c_eval('spe_PAD?_c!.f = anglea; spe_PAD?_c!.f_label = {''Energy''};',ic,channel)

anglea = pa_label';
c_eval('spe_PAD?_c! = struct(''t'', eOmni?.time.epochUnix);',ic,channel)
c_eval('spe_PAD?_c!.p = squeeze(PAD_spec_data?(:,!,:)); spe_PAD?_c!.p_label = {[''log('' eOmni?.units '')'']};',ic,channel)
c_eval('spe_PAD?_c!.f = anglea; spe_PAD?_c!.f_label = {''Energy''};',ic,channel)



%% plot
% set figure
h1=irf_plot(17,'newfigure');
lnwid = 1; fonts = 9;
set(0,'DefaultAxesFontSize',9);
set(0,'DefaultLineLineWidth',lnwid);
xSize=400; ySize=1200;
set(gcf,'Position',[100 100 xSize ySize]);

% Omniflux_1
i_plot = 1;
% [h(i_plot), hcb1]=irf_spectrogram(h(i_plot),specrec_p_e,);
c_eval('spec_pl1 = speOmni?;',ic)
[a5,b5]=irf_spectrogram(h1(i_plot),spec_pl1,'donotfitcolorbarlabel');
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

irf_legend(h1(i_plot),'e',[0.99 0.98],'color','w','fontsize',fonts);
grid off;
set(h1(i_plot),'yscale','log');
set(h1(i_plot),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h1(i_plot),{'\color{black}E(eV)','\color{black}(MMS1)'},'Interpreter','tex','fontsize',fonts);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a5,'position'); Pos_s = get(h1(1),'position');
set(a5,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a5,'fontsize',fonts); set(b5,'fontsize',fonts);
colormap(h1(i_plot),jet)
caxis(gca,[-0.5, 4.5])   


for i_plot = 2:17
% [h(i_plot), hcb1]=irf_spectrogram(h(i_plot),specrec_p_e,);
c_eval('spec_pl2 = spe_PAD?_c!;',ic, i_plot-1)
[a6,b6]=irf_spectrogram(h1(i_plot),spec_pl2,'donotfitcolorbarlabel');
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;

irf_legend(h1(i_plot),'e',[0.99 0.98],'color','w','fontsize',fonts);
grid off;
set(h1(i_plot),'yscale','lin');
% set(h1(i_plot),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h1(i_plot),{'\color{black}E(eV)','\color{black}(MMS1)'},'Interpreter','tex','fontsize',fonts);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a6,'position'); Pos_s = get(h1(i_plot-1),'position');
set(a6,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b6,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a6,'fontsize',fonts); set(b6,'fontsize',fonts);
colormap(h1(i_plot),jet)
caxis(gca,[1, 4])
end
irf_zoom(h1,'x',tint3);