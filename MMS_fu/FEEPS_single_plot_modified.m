% -----------modified by Wending Fu, Dec.2025 in Beijing------------
% ----------------------HAPPY NEW YEAR (๑´ڡ`๑)----------------------
%%
clear;clc;close all
global ParentDir
% cd 'Z:\SPART-WORK\Data\MMS\'
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Volumes/SPART-WORK/Data/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
% ParentDir = 'Z:\SPART-WORK\Data\MMS\';
OutputDir = '/Volumes/SPART-WORK/Data/MMS/';
% % % ic=1;
% % % 
% % % % TT = '2002-03-17\2003-01-01';
% % % % TT = '2001-05-01\2008-01-01';
% % % TT = '2015-10-01\2016-10-02';
% % % % % % TT = '2002-01-01\2003-01-01';
% % % Datelist = regexp(TT,'\d+-\d+-\d+','match');
% % % TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'/']; mkdir(TaskDir)
% % % 
% % % load([TaskDir,'FEEPSdata.mat'],  'all_RecordTime','all_channel15','all_channel16','all_Rx','all_Ry','all_Rz');%加载之前保存的数据
% % % % load(['/Volumes/SPART-NAS/Data/Cluster/Geomagnetic_Field_Statistic/2001-01-01T2004-06-02/','BRdata.mat'], 'all_MLT','all_MLat','all_Bphi','all_L');
% % combined
% file_list = {
%     % ['E:\fu\data2\Juno\sta2017-01-01T2018-01-01\BRdata.mat']
%     % ['E:\fu\data2\Juno\sta2018-01-01T2019-01-01\BRdata.mat'] 
%     % ['E:\fu\data2\Juno\sta2019-01-01T2020-01-01\BRdata.mat']
% %     ['D:\MMS_FEEPS_SURY\Feeps_statistic\2015-10-01To2016-10-02\FEEPSdata.mat']
% %     ['D:\MMS_FEEPS_SURY\Feeps_statistic\2016-10-02To2017-12-31\FEEPSdata.mat']
% %     ['D:\MMS_FEEPS_SURY\Feeps_statistic\2017-08-18To2018-01-01\FEEPSdata.mat']
% %     ['D:\MMS_FEEPS_SURY\Feeps_statistic\2018-01-01To2019-01-01\FEEPSdata.mat']
% %     ['D:\MMS_FEEPS_SURY\Feeps_statistic\2019-01-01To2020-01-01\FEEPSdata.mat']
% %     ['D:\MMS_FEEPS_SURY\Feeps_statistic\2020-01-01To2021-01-01\FEEPSdata.mat']
%       ['D:\MMS_FEEPS_SURY\Feeps_statistic\2021-01-01To2022-01-01\FEEPSdata.mat']
%       ['D:\MMS_FEEPS_SURY\Feeps_statistic\2022-01-01To2023-01-01\FEEPSdata.mat']
% };
% 
% all_RecordTime_combined = [];
% all_channel14_combined = [];
% all_channel15_combined = [];
% all_channel16_combined = [];
% all_Rx_combined = [];
% all_Ry_combined = [];
% all_Rz_combined = [];
% 
% for i = 1:length(file_list)
%     % 加载当前文件
%     data = load(file_list{i});
% 
%     % 拼接数据
%     all_RecordTime_combined  = [all_RecordTime_combined; data.all_RecordTime];
%     all_channel14_combined  = [all_channel14_combined; data.all_channel14];
%     all_channel15_combined  = [all_channel15_combined; data.all_channel15];
%     all_channel16_combined  = [all_channel16_combined; data.all_channel16];
%     all_Rx_combined = [all_Rx_combined; data.all_Rx];
%     all_Ry_combined = [all_Ry_combined; data.all_Ry];
%     all_Rz_combined = [all_Rz_combined; data.all_Rz];
%     
% end
% 
% all_RecordTime  = []; all_channel14 = []; all_channel15 = []; all_channel16 = [];all_Rx = []; all_Ry = []; all_Rz = [];
% 
% all_RecordTime = all_RecordTime_combined ;
% all_channel14 = all_channel14_combined;
% all_channel15 = all_channel15_combined;
% all_channel16 = all_channel16_combined;
% all_Rx =  all_Rx_combined;
% all_Ry =  all_Ry_combined;
% all_Rz =  all_Rz_combined;
% 
% %% 挑选特定时间范围数据
% % all_RecordTime_str = string(all_RecordTime);
% % % start
% % search_str1 = '2015-12-30T';
% % rows_with_search_str1 = contains(all_RecordTime_str, search_str1);
% % matching_rows1 = find(rows_with_search_str1);
% % marktime1 = matching_rows1(1,:)
% % % end
% % search_str2 = '2020-12-31T';
% % rows_with_search_str2 = contains(all_RecordTime_str, search_str2);
% % matching_rows2 = find(rows_with_search_str2);
% % marktime2 = matching_rows2(1,:)
% % % X = all_Rx(marktime1:marktime2);
% % % Y = all_Ry(marktime1:marktime2);
% % all_RecordTime = all_RecordTime(marktime1:marktime2,:);
% % all_channel15 = all_channel15(marktime1:marktime2);
% % all_channel16 = all_channel16(marktime1:marktime2);
% % all_Rx =  all_Rx(marktime1:marktime2);
% % all_Ry =  all_Ry(marktime1:marktime2);
% % all_Rz =  all_Rz(marktime1:marktime2);
% 
% %% 降采样 每13个数据取一个样点（约30s间隔）
% sel_idx =1:3:length(all_RecordTime);
% 
% all_RecordTime = [all_RecordTime(sel_idx,:)];
% all_channel14 = [all_channel14(sel_idx)];
% all_channel15 = [all_channel15(sel_idx)];
% all_channel16 = [all_channel16(sel_idx)];
% all_Rx =  [all_Rx(sel_idx)];
% all_Ry =  [all_Ry(sel_idx)];
% all_Rz =  [all_Rz(sel_idx)];
% 
% %% 筛选错误数据
% sel_idx =(all_channel16 >= 0);
% 
% all_RecordTime = [all_RecordTime(sel_idx,:)];
% all_channel14 = [all_channel14(sel_idx)];
% all_channel15 = [all_channel15(sel_idx)];
% all_channel16 = [all_channel16(sel_idx)];
% all_Rx =  [all_Rx(sel_idx)];
% all_Ry =  [all_Ry(sel_idx)];
% all_Rz =  [all_Rz(sel_idx)];
% 
% %% 筛选条件：位置的筛选
% idx = (all_Rx > 10) | (all_Rx < -12);
% % idx = (all_Rx >=-8 & all_Rx <= -7) & (all_Ry >= -9 & all_Ry <= -8);
% % idx = (all_Rx >=-4 & all_Rx <= -3) & (all_Ry >= -10 & all_Ry <= -9);
% 
% all_RecordTime = all_RecordTime(idx,:) ;
% all_channel14 = all_channel14(idx);
% all_channel15 = all_channel15(idx);
% all_channel16 = all_channel16(idx);
% all_Rx =  all_Rx(idx);
% all_Ry =  all_Ry(idx);
% all_Rz =  all_Rz(idx);
% 
% %% 通量的筛选
% idx_ch14 = (all_channel14 > 180)
% 
% all_RecordTime = all_RecordTime(idx_ch14,:);
% all_channel14 = all_channel14(idx_ch14);
% all_channel15 = all_channel15(idx_ch14);
% all_channel16 = all_channel16(idx_ch14);
% all_Rx =  all_Rx(idx_ch14);
% all_Ry =  all_Ry(idx_ch14);
% all_Rz =  all_Rz(idx_ch14);
% % save([OutputDir, 'FEEPSdata2.mat'], 'all_RecordTime')%
% % load('D:\FEEPS_plot_4sc\FEEPSdata1.mat')
% Nrec = size(all_RecordTime,1);
% for k = 1:Nrec
%  clc;clear B1 B2 B3 B4 R1 R2 R3 R4 Pos omniD1 energy_e1 Omni_flux_e1 omniD2 energy_e2 Omni_flux_e2 filenames_fast1 omniD3 energy_e3 Omni_flux_e3 omniD4 energy_e4 Omni_flux_e4;

%  t_center1 = all_RecordTime(k,:);
%  t_center =irf_time(t_center1,'utc>epoch');
%  t_start = t_center-150;
%  t_end = t_center+150;
% %  PlotTint1 = irf.tint(t_start,t_end);
%  t_start2 = irf_time(t_start,'epoch>utc');
%  t_end2 = irf_time(t_end,'epoch>utc');
% TT = [t_start2(1:21) 'Z/' t_end2(1:21) 'Z'];
ic = 1;
% TT='2016-03-25T17:33:07Z/2016-03-25T17:38:00Z';
% TT='2019-08-31T12:10:07Z/2019-08-31T12:25:00Z';
% TT='2017-09-07T00:20:07Z/2017-09-07T00:25:00Z';
% TT='2017-07-26T07:10:07Z/2017-07-27T23:00:00Z';
% TT='2024-05-11T16:00:00.000Z/2024-05-11T16:10:00.000Z';
TT='2024-05-10T16:00:00.000Z/2024-05-11T19:00:00.000Z';
% TT='2018-06-24T17:00:00.000Z/2018-06-24T17:10:00.000Z';
% TT='2019-03-03T14:20:07Z/2019-03-03T14:26:00Z';
% TT='2019-09-29T18:22:07Z/2019-09-29T18:27:00Z';
c_eval('PlotTint?=irf.tint(TT);',ic);
 
 %%
 Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
% iic = 1:4;
iic=1;
try
% filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames5 = SDCFilenames(Date,ic,'inst','feeps','drm','brst','dpt','electron');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); 
% % filenames_fast1 = SDCFilenames(Date,ic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms');
% % filenames_fast2 = SDCFilenames(Date,ic,'inst','edp','drm','fast');
% filenames = [filenames1, filenames2, filenames3, filenames4,filenames5];
% filenames_fast = [filenames_fast1, filenames_fast2];
% % % 
% [filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% [fileames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',ic);
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',ic);

% SDCFilesDownload_NAS(filenames,TempDir,'CheckSize', 0)
% SDCFilesDownload(filenames,TempDir)
% % % 
% SDCFilesDownload_NAS(filenames_fast,TempDir, 'Threads', 64, 'CheckSize', 0)
% SDCFilesDownload_NAS(filenames_srvy,TempDir, 'Threads', 64, 'CheckSize', 0)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
% % % filenames_srvy3 = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); 
% filenames_srvy2 = SDCFilenames(Date,ic,'inst','scm','drm','srvy','dpt','scsrvy');
% filenames_srvy1 = SDCFilenames(Date,iic,'inst','mec','drm','srvy','dpt','epht89d');
% % % filenames_srvy2 = SDCFilenames(Date,iic,'inst','feeps','drm','srvy','dpt','electron');
% % % filenames_srvy = [filenames_srvy2,filenames_srvy3];
% filenames_fast = SDCFilenames(Date,iic,'inst','fpi','drm','fast','dpt','des-moms,dis-moms');
% filenames_fast = [filenames_fast1];
% [filenames_fast,~,~] = findFilenames(TT,filenames_fast,'fast',iic);
% [filenames_srvy,~,~] = findFilenames(TT,filenames_srvy,'srvy',iic);

% % % SDCFilesDownload_NAS(filenames_srvy ,TempDir,  'CheckSize', 0)
% filenames_srvy = [filenames_srvy1];
% SDCFilesDownload_NAS(filenames_fast,TempDir, 'CheckSize', 0)
catch
    warning('no files have been downloaded')
end
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

%%
% try
  units = irf_units; %单位
  c_eval('markr = mms.get_data(''R_gsm'',PlotTint?,?);',ic);
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

% % load E
% c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',PlotTint?,?);'],ic);
% %%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
% c_eval(['Et?_ts=E?_ts.abs;'],ic); 
% c_eval(['E?=irf_gse2gsm(E?_ts);'],ic); % 坐标系转换
% c_eval(['E?=irf.ts2mat(E?);'],ic);
% c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
% c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);%电场数据重采样
% 
% c_eval(['Bt?_res=irf_resamp(Bt?,Et?);'],ic);
% 
% c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);%电场数据转换到以当地磁场为基准的场向FAC坐标系
% 
% c_eval('E?_err_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',PlotTint?);',ic);
% c_eval('E?_err=irf.ts2mat(E?_err_ts);',ic);%电场误差数据


%% load FPI
% c_eval('Ne?_ts = mms.get_data(''Ne_fpi_fast_l2'',PlotTint?,?);',ic);
% c_eval('Ni?_ts = mms.get_data(''Ni_fpi_fast_l2'',PlotTint?,?);',ic);
% % c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
% c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
% % c_eval('if isempty(Ne?), Ne? = NaN; end', ic);
% % c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
% c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);
% % c_eval('if isempty(Ni?), Ni? = NaN; end', ic);
% 
% 
% c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_fast_l2'',PlotTint?,?);',ic); 
% % c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
% c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
% c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
% c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],ic);
% c_eval(['gsmVi?=irf.ts2mat(gsmVi?_ts);'],ic);
% c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);
% % merge data/time from 2 cdf files
% % c_eval('energy_low?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_lowen_brst'',tint);',ic)
% % c_eval('energy_mid?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_miden_brst'',tint);',ic)
% % c_eval('energy_high?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_highen_brst'',tint);',ic) % 获取电子投掷角分布分析的三个不同能量区间的能谱信息
% c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_fast_l2_des-moms'',''mms?_des_energyspectr_omni_fast'',PlotTint?);',ic)
% c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_energyspectr_omni_fast'',PlotTint?);',ic) % 获取电子和离子的全向能量谱信息
% 


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
sensor_id = [1:5,9:12];
c_eval('e?_Tit! = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_top_intensity_sensorid_!'',PlotTint?);', ic, sensor_id);
c_eval('e?_Bit! = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_bottom_intensity_sensorid_!'',PlotTint?);', ic, sensor_id);
c_eval('ePA? = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_pitch_angle'',PlotTint?);', ic);
% c_eval('eLow = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_lower_bound'',tint);',ic);
% c_eval('eUp = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_upper_bound'',tint);',ic);
% c_eval('energies? = (eLow.data + eUp.data)/2. + eEcorr(?);',ic)
c_eval('Flux_e_feeps = mms.get_data(''Omnifluxelectron_epd_feeps_srvy_l2'',PlotTint?,?);',ic);
c_eval('energies? = Flux_e_feeps.depend{1, 1}(1,:);',ic);
%MMS1
c_eval('energies? =[48305.5000000000,67250,86075,104781.000000000,122185.500000000,140179,161732.500000000,186598.500000000,215486.500000000,249350.500000000,288661.500000000,334833.000000000,388701.500000000,451215,524269,575411.500000000]/1000;',ic);
%MMS2
% c_eval('energies? =[32488,50952.0000000000,69661.0000000000,88742.5000000000,106710.500000000,124924.500000000,145865.500000000,170279.000000000,199153.500000000,233105.500000000,272627.500000000,318592,372373,435080.000000000,508309.500000000,579540.000000000]/1000;',ic);
%MMS3
% c_eval('energies? =[31776.5000000000,50426.5000000000,69080,87986,105621.000000000,123762.000000000,145076.000000000,169944.000000000,199002.500000000,232623.500000000,271954,318007.000000000,371672.500000000,434475.500000000,507688.000000000,586024.500000000]/1000;',ic);
%MMS4
% c_eval('energies? =[28490,47126.5000000000,65763.0000000000,84522,102166,120309.500000000,141682.500000000,166535.500000000,195362.500000000,229152,268663.500000000,314637.000000000,368550.500000000,431543.500000000,504604,577058]/1000;',ic);
c_eval('Quality?sensor!top = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_top_quality_indicator_sensorid_!'',PlotTint?);', ic, sensor_id);
c_eval('Quality?sensor!bottom = mms.db_get_ts(''mms?_feeps_srvy_l2_electron'',''mms?_epd_feeps_srvy_l2_electron_bottom_quality_indicator_sensorid_!'',PlotTint?);', ic, sensor_id);
%%
% 假设每组数据存储在类似 Quality1sensorX 和 e1_BitX 这样的变量中
% 例如 Quality1sensor10bottom 和 e1_Bit10 表示第 10 个传感器的质量指示和数据
 % 需要处理的传感器 ID

% 循环遍历所有传感器
for sensor_idx = [1:5, 9:12]
    % 生成当前传感器的 Quality 和 Bit 数据的变量名
    quality_var_name = sprintf('Quality1sensor%dbottom', sensor_idx);  % 例如 'Quality1sensor10bottom'
    bit_var_name = sprintf('e1_Bit%d', sensor_idx);  % 例如 'e1_Bit10'
    
    % 获取 Quality 和 Bit 数据
    quality_data = eval([quality_var_name, '.data']);  % 例如 'Quality1sensor10bottom.data'
    bit_data = eval([bit_var_name, '.data']);  % 例如 'e1_Bit10.data'
    
    % 获取数据的大小
    [num_rows, num_cols] = size(bit_data);
    
    % 处理数据：检查 Quality 数据，根据条件修改 Bit 数据
    for i = 1:num_rows
        if quality_data(i) ~= 0 && quality_data(i) ~= 1
            % 如果 Quality 数据不是 0 或 1, 将对应行的 Bit 数据设为空值
            bit_data(i, :) = NaN;
        end
    end
    
    % 将处理后的结果保存回原数据结构中
    eval([bit_var_name, '.data = bit_data;']);  % 将修改后的数据保存回 'e1_Bit10.data' 等变量中
end

for sensor_idx = [1:5, 9:12]
    % 生成当前传感器的 Quality 和 Bit 数据的变量名
    quality_var_name = sprintf('Quality1sensor%dtop', sensor_idx);  % 例如 'Quality1sensor10bottom'
    tit_var_name = sprintf('e1_Tit%d', sensor_idx);  % 例如 'e1_Bit10'
    
    % 获取 Quality 和 Bit 数据
    quality_data = eval([quality_var_name, '.data']);  % 例如 'Quality1sensor10bottom.data'
    tit_data = eval([tit_var_name, '.data']);  % 例如 'e1_Bit10.data'
    
    % 获取数据的大小
    [num_rows, num_cols] = size(tit_data);
    
    % 处理数据：检查 Quality 数据，根据条件修改 Bit 数据
    for i = 1:num_rows
        if quality_data(i) ~= 0 && quality_data(i) ~= 1
            % 如果 Quality 数据不是 0 或 1, 将对应行的 Bit 数据设为空值
            tit_data(i, :) = NaN;
        end
    end
    
    % 将处理后的结果保存回原数据结构中
    eval([tit_var_name, '.data = tit_data;']);  % 将修改后的数据保存回 'e1_Bit10.data' 等变量中
end

%%
% MMS1:20170818前
% c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit10.data(:,:) = NaN;e?_Bit2.data(:,:) = NaN; e?_Bit9.data(:,:) = NaN;',ic);
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit10.data(:,:) = NaN;e?_Tit2.data(:,:) = NaN; e?_Tit9.data(:,:) = NaN;',ic);

% MMS1:20170818后
c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN;e?_Bit12.data(:,:) = NaN; e?_Bit3.data(:,:) = NaN;',ic);
c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit11.data(:,:) = NaN;e?_Tit2.data(:,:) = NaN; e?_Tit4.data(:,:) = NaN;',ic);

% MMS2:20170818前
% c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit10.data(:,:) = NaN; e?_Bit2.data(:,:) = NaN;e?_Bit9.data(:,:) = NaN;',ic);
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit10.data(:,:) = NaN; e?_Tit2.data(:,:) = NaN;e?_Tit9.data(:,:) = NaN;',ic);

% MMS2:20170818后
% c_eval('e?_Bit10.data(:,:) = NaN; e?_Bit12.data(:,:) = NaN; e?_Bit2.data(:,:) = NaN;e?_Bit3.data(:,:) = NaN;',ic);
% c_eval('e?_Tit12.data(:,:) = NaN;e?_Tit4.data(:,:) = NaN;e?_Tit2.data(:,:) = NaN;e?_Tit9.data(:,:) = NaN;',ic);

% MMS3前
% c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit10.data(:,:) = NaN;e?_Bit2.data(:,:) = NaN;e?_Bit9.data(:,:) = NaN;',ic);
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit10.data(:,:) = NaN;e?_Tit2.data(:,:) = NaN;e?_Tit9.data(:,:) = NaN;',ic);

% MMS3后
% c_eval('e?_Bit11.data(:,:) = NaN; e?_Bit12.data(:,:) = NaN;e?_Bit4.data(:,:) = NaN;e?_Bit5.data(:,:) = NaN;',ic);
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit11.data(:,:) = NaN;e?_Tit2.data(:,:) = NaN;e?_Tit4.data(:,:) = NaN;',ic);

% MMS4前
% c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit10.data(:,:) = NaN; e?_Bit2.data(:,:) = NaN;e?_Bit9.data(:,:) = NaN;',ic);
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit10.data(:,:) = NaN; e?_Tit2.data(:,:) = NaN;e?_Tit9.data(:,:) = NaN;',ic);

% MMS4后
% c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN; e?_Bit2.data(:,:) = NaN;e?_Bit4.data(:,:) = NaN;',ic);
% c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit12.data(:,:) = NaN; e?_Tit2.data(:,:) = NaN;',ic);

% %% clean up 
% % First, remove bad eyes
% % MMS1 
% c_eval('e?_Bit1.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN; e?_Tit1.data(:,:) = NaN;',ic);
% 
% % MMS2
% % % % c_eval('e?_Tit5.data(:,:) = NaN; e?_Tit7.data(:,:) = NaN; e?_Tit12.data(:,:) = NaN;',ic);
% % % % c_eval('e?_Bit7.data(:,:) = NaN;',ic);
% % MMS3
% % c_eval('e?_Tit2.data(:,:) = NaN; e?_Tit12.data(:,:) = NaN;',ic);
% % c_eval('e?_Bit2.data(:,:) = NaN; e?_Bit5.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN;',ic);
% 
% % MMS4
% % c_eval('e?_Tit1.data(:,:) = NaN; e?_Tit2.data(:,:) = NaN; e?_Tit7.data(:,:) = NaN;',ic);
% % c_eval('e?_Bit2.data(:,:) = NaN; e?_Bit4.data(:,:) = NaN; e?_Bit5.data(:,:) = NaN;',ic);
% % c_eval('e?_Bit10.data(:,:) = NaN; e?_Bit11.data(:,:) = NaN;',ic);
% 
% % Second, remove bad channels
% % channel 0
% % MMS1
% sensor_id = [2, 5];
% c_eval('e?_Tit!.data(:,1) = NaN;',ic,sensor_id)
% sensor_id = [2, 3, 4, 5, 6, 8, 9, 11, 12];
% c_eval('e?_Bit!.data(:,1) = NaN;',ic, sensor_id)
% % MMS2, 3 & 4
% % % % sensor_id = [1, 2, 3, 4, 5, 9, 10, 11, 12];
% % % % c_eval('e?_Tit!.data(:,1) = NaN;',ic,sensor_id)
% % % % sensor_id = [1, 2, 3, 4, 5, 9, 10, 11, 12];
% % % % c_eval('e?_Bit!.data(:,1) = NaN;',ic, sensor_id)
% % channel 1
% % % MMS1
% c_eval('e?_Tit6.data(:,2) = NaN;',ic);
% c_eval('e?_Bit6.data(:,2) = NaN; e?_Bit11.data(:,2) = NaN;',ic);
% % MMS2
% % c_eval('e?_Tit8.data(:,2) = NaN;',ic);
% % c_eval('e?_Bit7.data(:,2) = NaN; e?_Bit12.data(:,2) = NaN;',ic);
% % % MMS3
% % c_eval('e?_Tit1.data(:,2) = NaN;e?_Tit6.data(:,2) = NaN;e?_Tit7.data(:,2) = NaN;',ic);
% % c_eval('e?_Bit6.data(:,2) = NaN; e?_Bit7.data(:,2) = NaN;',ic);
% % % MMS4
% % c_eval('e?_Tit1.data(:,2) = NaN;e?_Tit6.data(:,2) = NaN;',ic);
% % c_eval('e?_Bit6.data(:,2) = NaN; e?_Bit7.data(:,2) = NaN;e?_Bit8.data(:,2) = NaN; e?_Bit0.data(:,2) = NaN;',ic);



%%
% if isempty(Ne1)
%    c_eval('Ne?(:,2) = NaN;',ic)
%    c_eval('Ni?(:,2) = NaN;',ic)
%    c_eval('gsmVi?(:,2:4)=NaN;',ic)
% end
%% calculate Omni flux
sensors = [1:5, 9:12];
% sensors = [3:5,11:12];
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
sensors   = [1:5,9:12];
nSensors  = length(sensors);
bin_size  = 16.3636; % default
n_pabins  = round(180./bin_size);
pa_bins   = 180*(0:1:n_pabins)./n_pabins;
pa_label  = 180*(0:1:n_pabins-1)/n_pabins + bin_size/2;
dAngResp  = 21.4; % default
delta_pa  = (pa_bins(2)-pa_bins(1))/2.0;

c_eval('dTmp = e?_Bit3;',ic)
c_eval('[sensor?_PA, sensor?_PA_id] = sort(ePA?.data, 2, ''ascend'');',ic)

c_eval([ ...
    'tempOmni = double(omniD?);' ...
    '[tempNt, tempNch, tempNs] = size(tempOmni);' ...
    'tempt = (1:tempNt)'';' ...
    'tempbase = tempt + (sensor?_PA_id-1) * tempNt * tempNch;' ...
    'templin  = permute(tempbase, [1 3 2]) + reshape((0:tempNch-1)*tempNt, 1, tempNch, 1);' ...
    'PAD_Data? = tempOmni(templin);' ...
],ic);

clear tempOmni tempNt tempNch tempNs tempt tempbase templin;

c_eval('PAD_spec_data? = NaN(size(PAD_Data?,1), size(PAD_Data?,2), n_pabins);',ic);

c_eval('PAD0 = PAD_Data?; valid0 = ~isnan(PAD0); PAD0(~valid0) = 0;',ic);
centers   = reshape(pa_label, 1, 1, []);
halfWidth = dAngResp + delta_pa;
c_eval('PA_sort = double(sensor?_PA);',ic);
W = (PA_sort > (centers - halfWidth)) & (PA_sort < (centers + halfWidth));

for ib = 1:n_pabins
    Wb = permute(W(:,:,ib), [1 3 2]);
    c_eval('num = sum(PAD0 .* Wb, 3);',ic);
    c_eval('den = sum(double(valid0) .* Wb, 3);',ic);
    c_eval('tmp = num ./ den; tmp(den==0) = NaN;',ic);
    c_eval('PAD_spec_data?(:,:,ib) = tmp;',ic);
end
clear PAD0 valid0 PA_sort num den tmp;

% find plot range
c_eval('idx_findout = find(eOmni?.time.epochUnix > irf_time(PlotTint?(1),''epochtt>epoch''));', ic)
idx_start = idx_findout(1);
c_eval('idx_findout = find(eOmni?.time.epochUnix < irf_time(PlotTint?(2),''epochtt>epoch''));', ic)
idx_end = idx_findout(end);

anglea = linspace(5,175,18)';
c_eval('PAD_resamp? = NaN(size(PAD_Data?,1), size(PAD_Data?,2), numel(anglea));', ic);

for ix = idx_start:idx_end
    c_eval('ang = double(sensor?_PA(ix,:));', ic); 
    c_eval('flux = double(squeeze(PAD_Data?(ix,:,:)));', ic);

    v = isfinite(ang);
    ang = ang(v);
    flux = flux(:,v);

    [ang_u,~,iu] = unique(ang,'stable');
    if numel(ang_u) < numel(ang)
        flux_u = NaN(size(flux,1), numel(ang_u));
        for k = 1:numel(ang_u)
        cols = (iu==k);
        flux_u(:,k) = mean(flux(:,cols), 2, 'omitnan');
        end
    else
        flux_u = flux;
    end

    c_eval('PAD_resamp?(ix,:,:) = interp1(ang_u, flux_u.'', anglea, ''linear'', NaN).'';', ic);
end
c_eval('PAD_Data? = PAD_resamp?;', ic);

%2 PAD
channel = 14:16;
% c_eval('spe_PAD?_c! = struct(''t'', eOmni?.time.epochUnix);', ic, channel);
% c_eval('spe_PAD?_c!.p = squeeze(PAD_Data?(:,!,:)); spe_PAD?_c!.p_label = {'' ''};', ic, channel);
% c_eval('spe_PAD?_c!.f = anglea; spe_PAD?_c!.f_label = {''PA (deg)''};', ic, channel);

c_eval('spe_PAD?_c! = struct(''t'', eOmni?.time.epochUnix);', ic, channel);
c_eval('spe_PAD?_c!.p = squeeze(PAD_spec_data?(:,!,:)); spe_PAD?_c!.p_label = {'' ''};', ic, channel);
c_eval('spe_PAD?_c!.f = pa_label''; spe_PAD?_c!.f_label = {''PA (deg)''};', ic, channel);


%% Init figure
% ic = 1;
% n=7;
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

%% N plot
% h(i)=irf_subplot(n,1,-i);
% 
% %滤波
% %     irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
% %     irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;
% 
% %非滤波
% % c_eval('irf_plot(h(i),Ne?,''Linewidth'',0.75,''color'',''b'');',ic);
% % % c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
% % c_eval('irf_plot(h(i),Ni?,''Linewidth'',0.75,''color'',''g'');',ic);
% c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;%绘制电子密度数据
% c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold off;
% grid off;
% if ~isnan(Ne1(:,2));
% c_eval("set(gca,'Ylim',[max([0 min([min(Ne?(:,2)) min(Ni?(:,2))])]) max([max(Ne?(:,2)) max(Ni?(:,2))])]);",ic);
% else
% end
% % set(gca,'Ylim',[0.1 0.4])
% % set(gca,'Ylim',[0.15 0.45], 'ytick',[0.1 0.2 0.3 0.4],'fontsize',9);
% % pos1=get(h(1),'pos');
% %  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
% %  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
%   set(gca,'ColorOrder',[[0 0 1];[0 1 0]]); % 图例颜色顺序
%   set(gca,'xtick',[])
%   set(gca,'fontsize',9)
%  irf_legend(gca,{'Ne','Ni'},[0.97 0.92]); % 图例位置
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
% ylabel('N [cm^{-3}]','fontsize',8);
% i=i+1;
% 
% 
% %% Vi plot
% h(i)=irf_subplot(n,1,-i);
% c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% grid off;
% if ~isnan(Ne1(:,2));
% c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+30],'fontsize',9);",ic); % 自动调Y轴显示范围else
% else
% end
% % set(gca,'Ylim',[-100 200], 'ytick',[0 200 400]);
% % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
% set(gca,'xtick',[])
% ylabel('Vi [km/s]','fontsize',8);
% i=i+1;

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
% set(h(i),'Ylim',[200 600]);
% set(h(i),'ytick',[1e4 2e5 4e5 6e5 1e6]);
% set(h(i),'ytick',[2e5 4e5]);
% ytickformat('%g');  % 通用数值格式
% h(i).YRuler.Exponent = 0;
% ytickformat('%.0e'); 
% set(h(i),'ytick',[1e4 1e5]);
set(gca,'fontsize',9)
ylabel(h(i),{'\color{black}E(keV)','\color{black}'},'Interpreter','tex','fontsize',8);
xwidth=0.02; ywidth=0.154;
Pos_c = get(a5,'position'); Pos_s = get(h(i),'position');
set(a5,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a5,'fontsize',8); set(b5,'fontsize',8);
% colormap(h(i),jet)
i=i+1;
% caxis(gca,[0, 3])   

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
Pos_c = get(a6,'position'); Pos_s = get(h(i-1),'position');
set(a6,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b6,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a6,'fontsize',fonts); set(b6,'fontsize',fonts);
colormap(h(i),jet)
% caxis(h(i),[-1, 2])
i=i+1;
%% PAD channel15
h(i)=irf_subplot(n,1,-i);
c_eval('spe_PAD?_c!.p_label = {[''log('' eOmni?.units '')'']};',ic,15);
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
% caxis(h(i),[-1, 2])
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
% caxis(h(i),[-1, 2])
i=i+1;

%%  出图保存部分
title(h(1),sprintf('MMS1 X = %.2f  Y = %.2f  Z = %.2f', markpos1,markpos2, markpos3),'fontsize',10);
set(h(1:n),'fontsize',12);
c_eval('irf_zoom(PlotTint?,''x'',h(1:n));',ic);

% % % c_eval('spec_fpi_e?=struct(''t'',Omni_flux_e?(:,1));',ic)
    % irf_adjust_panel_position;
    % %   irf_plot_axis_align(h)
    irf_plot_axis_align(h)
colormap(jet) % 全局颜色映射设置
set(gca,"XTickLabelRotation",0) % 确保X轴刻度标签方向水平
set(gcf,'render','painters'); % 设置图形矢量渲染
set(gcf,'paperpositionmode','auto') % 打印或导出文件的尺寸与屏幕上显示的图形尺寸完全一致

% str_time = irf_time(flagTime,'epoch>utc');
% fig_dir0 = strcat( OutputDir,'1215MMS4\');
% % fig_dir = strcat(fig_dir0, '\',str_time(1:4),'\',str_time(6:7),'\',str_time(9:10));
% fig_dir = fig_dir0;
% if ~exist( fig_dir, 'dir')
%     mkdir(fig_dir);
% end
% fig_name = strcat(t_center1(1:4),t_center1(6:7),t_center1(9:10),t_center1(12:13),t_center1(15:16),t_center1(18:19));
% fig_path_name = strcat(fig_dir,'\',fig_name,'.png');
% % fig_path_name = strcat(fig_dir,'\',fig_name,'.pdf');
% 
% % set(gcf, 'InvertHardCopy', 'off');
% % set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
% % if(eps)
% % print('-depsc','-painters','-r150',fig_path_name);
% % else
% print('-dpng','-painters','-r300',fig_path_name);
% end
% close all
% catch
%     writematrix([t_center1(1:end-1),'的数据导入4出现问题需要手动下载'],[OutputDir,'errorlog1215MMS4.txt'],...
%         'WriteMode','append','Encoding','UTF-8')
%     continue
% end
% end