%------written by Wending Fu, Oct.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = 'Z:\SPART-WORK\Data\MMS\'; 
DownloadDir = 'Z:\SPART-WORK\Data\MMS\';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
 
ic = 1;
% iic = 1:3;
load([pwd,'/NameTags.mat']);%设置路径

% Date = '2018-01-26T00:00:00.000Z/2019-03-31T00:00:00.000Z';
% Date = '2015-09-11T00:00:00.0Z/2016-09-11T00:00:00.00Z';
Date = '2022-07-11T12:26:33.1Z/2022-07-11T12:30:34.9Z';
splitDate = regexp(Date,'/','split');
startDate = splitDate{1};
endDate   = splitDate{2};
startDT = datetime(startDate, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
endDT   = datetime(endDate,   'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
NameTagsDT = datetime(NameTags, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
mask = (NameTagsDT >= startDT) & (NameTagsDT <= endDT);
NameTags = NameTags(mask);%选择处理时间窗
%% download data
% % % filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% % % filenames2 = SDCFilenames(Date,iic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms');
% % % % filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% % % % filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% % % % filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
% % % filenames = [filenames1,filenames2];
% % % % % % filenames = filenames1;
% % % 
% % % expr = '_[0-9]+\_v';
% % % NameTags = regexp(filenames,expr,'match');
% % % NameTags = unique(cellfun(@cellstr,NameTags));
% % % FileGroups = cell(1,length(NameTags)); 
% % % for j = 1:length(NameTags)
% % %     % FileGroups{j} = [filenames(contains(filenames,NameTags{j})), filenames_srvy(contains(filenames_srvy,NameTags{j}(2:9)))];
% % %     FileGroups{j} = filenames(contains(filenames,NameTags{j}));
% % % end
% % % FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组
%%
%修改文件夹时特别注意SDCFilesDownload需要datamove的文件夹必须是ParentDir，否则需要手动修改
splitDate{1} = erase(splitDate{1},{':','.'});splitDate{2} = erase(splitDate{2},{':','.'});
OutputDir = [ParentDir,'Highelectron_Search/Eev_Search/',splitDate{1},'To',splitDate{2},'/'];
clear mask NameTagsDT
mkdir(OutputDir)
mms.db_init('local_file_db',ParentDir);%输出目录与数据库句柄
%%
units = irf_units;
% parfor_progress(length(NameTags)-1);
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;clear B1 B2 B3 B4 R1 R2 R3 R4 Pos;
PlotFlag = 1;
% parfor_progress;
% disp(['૮₍ ˃ ⤙ ˂ ₎ა',repmat('■',1,round(10*TDT/length(NameTags))),repmat('□',1,10-round(10*TDT/length(NameTags))),'正在光速检索ing...'])
% fprintf(['当前处理时间为:',NameTags{TDT}(2:end-2),'\n'])

% % % if length(FileGroups{TDT})~=4, continue;end
% % % try
% % %     SDCFilesDownload_NAS(FileGroups{TDT},TempDir, 'Threads', 32, 'CheckSize', 0)
% % % catch
% % %     writematrix([NameTags{TDT}(2:end-2),'的Z盘数据导入出现问题1'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % %     continue
% % % end
% % % 
% % % SDCDataMove(TempDir,ParentDir)
% % % 
% % % mms.db_init('local_file_db',ParentDir);
% % % formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];
% % % tempDate = [formatDate(NameTags{TDT}), '/', formatDate(NameTags{TDT+1})];
tempDate = [char(NameTags(TDT)) '/' char(NameTags(TDT+1))];
tempTint=irf.tint(tempDate);

tryTimes = 0;
while tryTimes <= 10
try
    % B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    % tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));  
    mms.db_init('local_file_db',ParentDir);
    Pos = mms.get_data('R_gsm',tempTint,1);
    Pos = irf.ts2mat(Pos);
    tryTimes = 666;
catch
    tryTimes = tryTimes + 1;
    continue
end
end

if tryTimes ~= 666
writematrix([NameTags{TDT}(2:end-2),'的数据导入1出现问题'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
continue
end
%读取R
%% 判据
try
   Pos = mms.get_data('R_gsm',tempTint);
    c_eval('R? = Pos.gsmR?;',ic)
    c_eval('R? = [Pos.time.epochUnix R?(:,4)];',ic)
    RE = units.RE;
    c_eval('R_RE = mean(R?(1,2)/RE*1000,1);',ic) 
%R判据
    if R_RE < 8, continue; end
    tint = irf.tint(irf_time(R1(1,1),'epoch>epochTT'),irf_time(R1(end,1),'epoch>epochTT'));

%FEEPS
eGfact = [1.0, 1.0, 1.0, 1.0];
specie = 'electron';
sensor_id = [1:5,9:12];
c_eval('e?_Tit! = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_!'',tint);', ic, sensor_id);
c_eval('e?_Bit! = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_!'',tint);', ic, sensor_id);
c_eval('Flux_e_feeps = mms.get_data(''Omnifluxelectron_epd_feeps_brst_l2'',tint,?);',ic);
c_eval('energies? = Flux_e_feeps.depend{1, 1}(1,:);',ic);
%calculate Omni flux
sensors = [1:5, 9:12];
nSensors = length(sensors);
c_eval('dTmp= e?_Bit2;', ic)
% omniD = NaN([size(dTmp.data) nSensors*2]);
for iSen = 1:nSensors
    c_eval(['omniD?(:,:,iSen) = ' specie(1) '?_Tit!.data;'...
        'omniD?(:,:,nSensors+iSen) = ' specie(1) '?_Bit!.data;'],ic,sensors(iSen))
end
c_eval([specie(1) 'Omni? = dTmp; ' specie(1) 'Omni?.data =' ...
    'mean(double(omniD?),3,''omitnan'')*' specie(1) 'Gfact(?);'],ic)

c_eval('speOmni? = struct(''t'', eOmni?.time.epochUnix);',ic)
c_eval('speOmni?.p = eOmni?.data; speOmni?.p_label = {[''log('' eOmni?.units '')'']};',ic)
c_eval('speOmni?.f = double(energies?); speOmni?.f_label = {''Energy''};',ic)

%Energy判据
c_eval( 'checkdata? = eOmni?.data;',ic)
c_eval( 'checkdatalast2? = checkdata?(:, end-1:end);',ic)
c_eval( '[maxenergy?, linIdx?] = max(checkdatalast2?(:));',ic)
c_eval( '[rowRel?, colRel?] = ind2sub(size(checkdatalast2?), linIdx?);',ic)
c_eval( 'row? = rowRel?;',ic)                     % 行号（在 A 中相同）
c_eval( 'col?= colRel? + 1;',ic)                  % 列号（映射回 A 的列号：2或3）
c_eval('Flag_StrongEev? = find(maxenergy? >= 10^-25);',ic);
c_eval('Time? = speOmni?.t;',ic);

    if isempty(Flag_StrongEev1) & isempty(Flag_StrongEev2) & isempty(Flag_StrongEev3) & isempty(Flag_StrongEev4) 
        writematrix([NameTags{TDT}(1:end-1) , '无Eev>10^-25事件'],...
            [OutputDir,'Log.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
    end 

    if ~isempty(Flag_StrongEev1)
        writematrix([irf_time(Time1(row1,1),'epoch>utc'),' maxEev1 = ', num2str(maxenergy1)],...
            [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    

%     if ~isempty(Flag_StrongEev2)
%         writematrix([irf_time(Time2(row2,1),'epoch>utc'),' maxEev2 = ', num2str(maxenergy2)],...
%             [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
%     end
%     
%     
%      if ~isempty(Flag_StrongEev3)
%         writematrix([irf_time(Time3(row3,1),'epoch>utc'),' maxEev3 = ', num2str(maxenergy3)],...
%             [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
%      end
%     
%       if ~isempty(Flag_StrongEev4)
%         writematrix([irf_time(Time4(row4,1),'epoch>utc'),' maxEev4 = ', num2str(maxenergy4)],...
%             [OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
%      end
     
catch
    writematrix([NameTags{TDT}(1:end-1),'的数据导入2出现问题'],[OutputDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end

%% 符合判据的继续下载并出图
if PlotFlag == 1
try
%     SDCFilesDownload(FileGroups{TDT},TempDir) 
%     SDCDataMove(TempDir,ParentDir)
    mms.db_init('local_file_db',ParentDir);
    c_eval('flagTime = Time?(row?,1);',ic)
    c_eval('PlotTint?= irf_time([flagTime-30 ,flagTime+30],''epoch>epochTT'');',ic)
%     c_eval('id_flagTime = SDCPlot(PlotTint?,desmoms,ic,NameTags{TDT},flagTime?);',ic)
%% load B
units = irf_units; %单位
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',PlotTint?,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
 % c_eval(['B?=irf_gse2gsm(B?);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);

% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',PlotTint?,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?=irf_gse2gsm(E?_ts);'],ic); % 坐标系转换
c_eval(['E?=irf.ts2mat(E?);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);%电场数据重采样

c_eval(['Bt?_res=irf_resamp(Bt?,Et?);'],ic);

c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);%电场数据转换到以当地磁场为基准的场向FAC坐标系

c_eval('E?_err_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',PlotTint?);',ic);
c_eval('E?_err=irf.ts2mat(E?_err_ts);',ic);%电场误差数据


% load FPI
c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',PlotTint?,?);',ic);
c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',PlotTint?,?);',ic);
% c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
% c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);


c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',PlotTint?,?);',ic); 
% c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],ic);
c_eval(['gsmVi?=irf.ts2mat(gsmVi?_ts);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);
% merge data/time from 2 cdf files
% c_eval('energy_low?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_lowen_brst'',tint);',ic)
% c_eval('energy_mid?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_miden_brst'',tint);',ic)
% c_eval('energy_high?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_highen_brst'',tint);',ic) % 获取电子投掷角分布分析的三个不同能量区间的能谱信息
c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energyspectr_omni_brst'',PlotTint?);',ic)
c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_energyspectr_omni_brst'',PlotTint?);',ic) % 获取电子和离子的全向性能量谱信息



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
c_eval('e?_Tit! = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_!'',PlotTint?);', ic, sensor_id);
c_eval('e?_Bit! = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_!'',PlotTint?);', ic, sensor_id);


c_eval('ePA? = mms.db_get_ts(''mms?_feeps_brst_l2_electron'',''mms?_epd_feeps_brst_l2_electron_pitch_angle'',PlotTint?);', ic);
% % % c_eval('eLow = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_lower_bound'',tint);',ic);
% % % c_eval('eUp = mms.db_get_variable(''mms?_feeps_brst_l2_electron'', ''electron_energy_upper_bound'',tint);',ic);
% % % c_eval('energies? = (eLow.data + eUp.data)/2. + eEcorr(?);',ic)
c_eval('Flux_e_feeps = mms.get_data(''Omnifluxelectron_epd_feeps_brst_l2'',PlotTint?,?);',ic);
c_eval('energies? = Flux_e_feeps.depend{1, 1}(1,:);',ic);

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
    'mean(double(omniD?),3,''omitnan'')*' specie(1) 'Gfact(?);'],ic)

%% Init figure
% ic = 1;
n=5;
i=1;
set(0,'DefaultAxesFontSize',8); %字号
set(0,'DefaultLineLineWidth', 0.5); %线宽
fn=figure(1);clf;
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
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);%在指定位置创建图例
ylabel('B [nT]','fontsize',10); % Y轴标签
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;

%% N plot
h(i)=irf_subplot(n,1,-i);

%滤波
%     irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
%     irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;

%非滤波
c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;%绘制电子密度数据
c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[max([0 min([min(Ne?(:,2)) min(Ni?(:,2))])-0.02]) max([max(Ne?(:,2)) max(Ni?(:,2))])+0.02]);",ic)
% set(gca,'Ylim',[0.1 0.4])
% set(gca,'Ylim',[0.15 0.45], 'ytick',[0.1 0.2 0.3 0.4],'fontsize',9);
% pos1=get(h(1),'pos');
%  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
%  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]); % 图例颜色顺序
  set(gca,'xtick',[])
 irf_legend(gca,{'Ne','Ni'},[0.97 0.92]); % 图例位置
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
ylabel('N [cm^{-3}]','fontsize',8);
i=i+1;


%% Vi plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+30],'fontsize',9);",ic); % 自动调Y轴显示范围
% set(gca,'Ylim',[-100 200], 'ytick',[0 200 400]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
set(gca,'xtick',[])
ylabel('Vi [km/s]','fontsize',8);
i=i+1;

% %% plot low e pad
% %     %0-200eV
% h(i)=irf_subplot(n,1,-i);
% % h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(h(i),jet)
% specrec_p_elow=struct('t',irf_time(energy_low.DEPEND_0.data,'ttns>epoch'));
% specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_elow.p=energy_low.data;%data matrix
% specrec_p_elow.f_label='';
% specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
% 
% [h(i), hcb6]=irf_spectrogram(h(i),specrec_p_elow);
% ylabel('PA low','fontsize',8)
% % set(gca,'yscale','log');
% set(h(i),'ytick',[0 90 180]);
% % caxis(gca,[7 7.7]);
% %irf_legend(h(i),'g',[0.99 0.98],'color','w','fontsize',12);
% poscbar6=get(hcb6,'pos'); % 获取颜色条位置
% poscbar6(3)=poscbar6(3)*0.5; % 颜色条宽度减半
% set(hcb6,'pos',poscbar6); % 宽度信息重新赋给颜色条
% i=i+1;


%% plot FEEPS electron flux spectrom
h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)
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
set(h(i),'ytick',[1e4 1e5]);
ylabel(h(i),{'\color{black}E(eV)','\color{black}'},'Interpreter','tex','fontsize',8);
% xwidth=0.02; ywidth=0.154;
Pos_c = get(a5,'position'); Pos_s = get(h(1),'position');
set(a5,'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);
set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2) 0.02 Pos_c(4)]);
set(a5,'fontsize',8); set(b5,'fontsize',8);
colormap(h(i),jet)
i=i+1;
% caxis(gca,[-0.5, 4.5])   

%% plot e energy spectrom 绘制电子全向能谱
c_eval('Omni_flux_e?(:,1)=irf_time(energy_e?.DEPEND_0.data,''ttns>epoch'');',ic)
c_eval('Omni_flux_e?(:,2:33)=energy_e?.data;;',ic)
c_eval('channel?=transpose(energy_e?.DEPEND_1.data(1,1:32));',ic)

c_eval('spec_fpi_e?=struct(''t'',Omni_flux_e?(:,1));',ic)
c_eval('spec_fpi_e?.f=transpose(energy_e?.DEPEND_1.data(1,1:32));',ic) %energy levels
c_eval('spec_fpi_e?.p=Omni_flux_e?(:,2:33);',ic) %data matrix
c_eval('spec_fpi_e?.f_label=''Energy'';',ic)
c_eval('spec_fpi_e?.p_label={'' '',''(cm^2 s sr keV)^{-1}''};',ic)

c_eval('spec_pl2 = spec_fpi_e?;',ic)

h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)
[a5,b5]=irf_spectrogram(h(i),spec_pl2,'donotfitcolorbarlabel');

% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
grid off;
set(h(i),'yscale','log');  % 对数刻度
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ee(ev)','fontsize',8)
set(gca,'Ylim',[1.5e1 3e4]);
% caxis(gca,[7.1 7.4])
set(gca,'xtick',[])
% irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
Pos_c = get(a5,'position'); Pos_s = get(h(1),'position');
set(a5,'position',[Pos_c(1) Pos_c(2)-0.004 Pos_s(3) Pos_c(4)]);
set(b5,'position',[Pos_c(1)+Pos_s(3)+0.01 Pos_c(2)-0.004 0.02 Pos_c(4)]);
set(a5,'fontsize',8); set(b5,'fontsize',8);
colormap(h(i),jet)

i=i+1;
%%  出图保存部分
colormap(jet) % 全局颜色映射设置
set(gca,"XTickLabelRotation",0) % 确保X轴刻度标签方向水平
set(gcf,'render','painters'); % 设置图形矢量渲染
set(gcf,'paperpositionmode','auto') % 打印或导出文件的尺寸与屏幕上显示的图形尺寸完全一致

% % % cd  /Users/fwd/Documents/Ti~mor~/M/DF&MP/1-Comparison&Implication/Figures/
% % % % rmdir(TempDir,'s'); 
figname = 'mp';
    print(gcf, '-dpng', [figname '.png']); 
catch
 writematrix([irf_time(flagTime,'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
   'WriteMode','append','Encoding','UTF-8')

    continue
end
end
%% 删除文件夹并生成记录文件
% % % try
% % %     cd(OutputDir)
% % %     rmdir(TempDir,'s');    
% % % catch
% % %     writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % % end
end