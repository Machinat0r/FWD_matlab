clear
clc;close all
mms.db_init('local_file_db','/Volumes/FWD-T7Disk/MMS/')
%% Load PDist using mms.make_pdist
%for ic = 1:3; % spacecraft id
ic=1;
% Time='2015-09-19T07:43:31.115Z';
% Time='2019-08-05T16:24:34.500Z';
% time = irf_time(Time,'utc>epochtt');
time = irf.tint('2019-08-05T16:24:30.00Z/2019-08-05T16:24:40.00Z');

%[filepath,filename] = mms.get_filepath('mms2_fpi_brst_l2_des-dist',time);
c_eval('filepath_and_filename = mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',time);',ic);
c_eval('[ePDist?,ePDistError?] = mms.make_pdist(filepath_and_filename);',ic)
%c_eval('[iPDist?,iPDistError?] = mms.make_pdist(filepath_and_filename);',ic)
%% Load supporting data
c_eval('tint = ePDist?.time([1 end]);',ic);
c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('Vdbcs?=mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
%dmpaB1=mms.db_get_ts('mms1_fgm_brst_l2','mms1_fgm_b_dmpa_brst_l2',tint);

%% plots: particle distributions

for ip=1

switch ip
    case 1
Time='2019-08-05T16:24:34.500Z';
% Time='2015-09-19T07:43:31.115Z';
    case 2
Time='2019-08-05T16:24:34.00Z';
end
%for n=0:3  
  %n=3;
%        Time='2015-09-19T07:43:31.403Z';
  %D=Time(18:22);
  %C=str2num(D);
  %C=C+0.03*n;
  %D=num2str(C);
  %Time(18:22)=D;
  % time = irf_time(Time,'utc>epochtt');
  time = irf.tint('2019-08-05T16:24:34.500Z/2019-08-05T16:24:35.000Z');
  
  c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,[17]);',ic);
  c_eval('scpot = scPot?.resample(time);',ic);  
  c_eval('Ve= Vdbcs?.resample(time).data;',ic);
  hatVe = double(irf_norm(Ve));
  c_eval('B0 = dmpaB?.resample(time).data;',ic); 
  hatB0 = double(irf_norm(B0));
  perp1 = mean(cross(cross(hatB0,hatVe),hatB0),1);
  perp2 = mean(cross(hatB0,hatVe),1);
  para = mean(hatB0,1);
  %N=[0.9600, 0.0618, -0.2731];
  %perp1=cross(hatB0,cross(N,hatB0));
  %perp2=cross(N,hatB0);
  %perp1=N;
  %perp2=[0.14,  0.17,  0.98];
  %c_eval('E0 = dslE?.resample(time).data;',ic); 
  %hatE0 = double(irf_norm(E0));
  %hatExB0 = cross(hatE0,hatB0);
  

 % optional input parameters for projection plot
  vlim = 50*1e3; % x and ylim
  elevlim = 20; % angle over plane to include in slice
  strCMap = 'jet'; % colormap
  projclim = [-1.5 0.5]; % colorbar limit
  colormap('jet')

  %x = hatE0;
  %y = hatExB0;
  %z = hatB0;
  x = perp1;
  y = perp2;
  z = para;
  % % % c_eval('tInd = find(abs(ePDist?.time-time)==min(abs(ePDist?.time-time)));',ic);

 % Initialize figure
hhh = figure; 
  nRows = 3; nCols = 1;
    for ii = 1:nRows*nCols;
        h(ii) = subplot(nRows,nCols,ii); 
    end
    isub = 1;

    hca = h(isub); isub = isub + 1;
    xyz=[perp1;perp2;para];
    vlabels={'V_{perp1}';'V_{perp2}';'V_{para}'};
    c_eval('plot_projection_modified_by_fwd(hca,ePDist?.convertto(''s^3/km^6''),''tint'',time,''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''clim'',projclim,''scpot'',scpot,''vlabel'',vlabels);',ic);

    
    %hca = h(isub); isub = isub + 1;     
    %c_eval('mms.plot_skymap(hca,ePDist?,''tint'',time,''energy'',100,''flat'');',ic);

    hca = h(isub); isub = isub + 1;
    xyz=[z;x;y];
    vlabels={'V_{para}';'V_{perp1}';'V_{perp2}'};
    % c_eval('mms.plot_projection(hca,ePDist?.convertto(''s^3/km^6''),''tint'',time,''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''clim'',projclim,''scpot'',scpot,''vlabel'',vlabels);',ic);
    c_eval('plot_projection_modified_by_fwd(hca,ePDist?.convertto(''s^3/km^6''),''tint'',time,''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''clim'',projclim,''scpot'',scpot,''vlabel'',vlabels);',ic);


    %hca = h(isub); isub = isub + 1;      
    %c_eval('mms.plot_skymap(hca,ePDist?,''tint'',time,''energy'',100,''flat'',''log'');',ic);
    %hca.CLim = projclim;

    hca = h(isub); isub = isub +1; 
    xyz=[z;y;x];
    vlabels={'V_{para}';'V_{perp2}';'V_{perp1}'};
    c_eval('plot_projection_modified_by_fwd(hca,ePDist?.convertto(''s^3/km^6''),''tint'',time,''xyz'',xyz,''elevationlim'',elevlim,''vlim'',vlim,''clim'',projclim,''scpot'',scpot,''vlabel'',vlabels);',ic);


    %hca = h(isub); isub = isub + 1;    
    %c_eval('mms.plot_skymap(hca,ePDist?,''tint'',time,''energy'',100,''vectors'',{hatB0,''B''},''log'');',ic);
colormap("jet")

%% save figure
  set(gcf,'render','painters');%矢量图
  % set(gcf,'outerposition',get(0,'screensize'))%%figure 是创建一个绘图窗口,把这个绘图窗口置为全屏显示 
%    c_eval('figname=[''PDist?sc!''];,',tInd,ic);
%   saveas(gcf,figname,'fig');

%   print(gcf, '-dpdf', [figname '.pdf']);
 %end 
 %end
 
% figure('Name','Time');
% print(hhh,'-dpdf','C:\Users\Guo Zhizhong\Desktop\tupian'); 
end
% saveas(gcf,strcat('C:\Users\Guo Zhizhong\Desktop\tupian','pdf'));