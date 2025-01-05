%------modified by Wending Fu, Dec.2024 in Beijing------------
%------written by Bolin Chen, Dec.2024 in Beijing------------
%%
close all
clear;clc

global ParentDir 
ParentDir = 'Z:\/Data/MMS/'; 
DownloadDir = 'Z:\/Data/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);


TT = '2017-05-22T10:42:12.00Z/2017-05-22T10:42:14.00Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

ic = 1:4;
% filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% 
% filenames = [filenames1,filenames2];
% [filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% 
% SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
%% load data
% SDCDataMove(TempDir,ParentDir)
% mms.db_init('local_file_db',ParentDir);

units = irf_units;
c_eval(['B?_ts=mms.get_data(''B_gse_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?_mat=irf.ts2mat(B?_ts);'],ic);
c_eval(['Bt?_mat=irf.ts2mat(Bt?_ts);'],ic);
c_eval('B?=[B?_mat(:,1) B?_mat(:,2:end)* 1e-9];',ic); % T
c_eval('Bt?=[Bt?_mat(:,1) Bt?_mat(:,2:end)* 1e-9];',ic); % T

c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
c_eval(['E?_mat=irf.ts2mat(E?_ts);'],ic);
c_eval('E?=[E?_mat(:,1) E?_mat(:,2:end)* 1e-3];',ic); % V/m

c_eval('Pa_e?_ts = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_prestensor_gse_brst'',tint);',ic);
c_eval('Pa_i?_ts = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_prestensor_gse_brst'',tint);',ic);

c_eval('V_e?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic)
c_eval(['V_e_mat?=irf.ts2mat(V_e?_ts);'],ic);
c_eval('V_e?=[V_e_mat?(:,1) V_e_mat?(:,2:end)* 1e3];',ic); % m/s
c_eval('V_i?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['V_i_mat?=irf.ts2mat(V_i?_ts);'],ic);
c_eval('V_i?=[V_i_mat?(:,1) V_i_mat?(:,2:end)* 1e3];',ic); % m/s

c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
% c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne_mat?=irf.ts2mat(Ne?_ts);'],ic);
c_eval('Ne?=[Ne_mat?(:,1) Ne_mat?(:,2:end)* 1e6];',ic);  % m^-3
% c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni_mat?=irf.ts2mat(Ni?_ts);'],ic);
c_eval('Ni?=[Ni_mat?(:,1) Ni_mat?(:,2:end)* 1e6];',ic); % m^-3

% 将每个 3x3 矩阵展平成行向量,转换为double
c_eval('Pa_e? = double([Pa_e?_ts.time.epochUnix, reshape(Pa_e?_ts.data, size(Pa_e?_ts.data,1), [] )* 1e-9]);', ic); % Pa
c_eval('Pa_i? = double([Pa_i?_ts.time.epochUnix, reshape(Pa_i?_ts.data, size(Pa_i?_ts.data,1), [] )* 1e-9]);', ic); % Pa

Pos = mms.get_data('R_gse',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gseR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)* 1e3];') % m
%% 计算电子和离子标量压力和偏应力张量
for spec = ['e', 'i']
c_eval(['Pa_' spec '?_sca = Pa_' spec '?;'],ic);
c_eval(['Pa_' spec '?_sca(:, [3:5,7:9]) = 0;'],ic);
c_eval(['Pa_' spec '?_sca = [Pa_' spec '?(:,1), 1/3*Pa_' spec '?_sca(:,2:end)];'],ic);
c_eval(['Pa_' spec '?_dev = [Pa_' spec '?(:,1), Pa_' spec '?(:,2:end) - Pa_' spec '?_sca(:,2:end)];'],ic);
end
% Pa_spec_sca = 1/3*P_ii, Pa_spec_dev = Pi_ij
%% 计算电子和离子的应变率
c_eval('R? = irf_resamp(R?,V_e?);')
gradV_e=c_4_grad('R?','V_e?','grad');
    time = gradV_e(:, 1);
    data = gradV_e(:, 2:end);

 Se_data = zeros(size(data));  % S_ij

for i = 1:size(data, 1)
    
    row_data = reshape(data(i, :), 3, 3);  
    Se_row = 1/2 * (row_data + row_data');  
    Se_data(i, :) = Se_row(:)';  
end
Se= [time,Se_data];

c_eval('R? = irf_resamp(R?,V_i?);')
gradV_i=c_4_grad('R?','V_i?','grad');
    time = gradV_i(:, 1);
    data = gradV_i(:, 2:end);

 Si_data = zeros(size(data)); 

for i = 1:size(data, 1)
    
    row_data = reshape(data(i, :), 3, 3);  
    Si_row = 1/2 * (row_data + row_data');  
    Si_data(i, :) = Si_row(:)';  
end
Si= [time,Si_data];



%% 求电子和离子的无迹应变率张量

Se_sca = Se;
Se_sca(:, [3:5,7:9]) = 0; % theta
Se_sca = [Se(:,1), 1/3*Se_sca(:,2:end)];
De = [Se(:,1), Se(:,2:end) - Se_sca(:,2:end)];

Si_sca = Si;
Si_sca(:, [3:5,7:9]) = 0;
Si_sca = [Si(:,1), 1/3*Si_sca(:,2:end)];
Di = [Si(:,1), Si(:,2:end) - Si_sca(:,2:end)];


 %% 求电子和离子的Pi-D
De=De(1:66,:);
c_eval('Pa_e?_dev=Pa_e?_dev(1:66,:);',1:3);
Di=Di(1:13,:);
Pa_i1_dev=Pa_i1_dev(1:13,:);
 c_eval(['Pi_De? = -sum(sum(De(:,2:end).*Pa_e?_dev(:,2:end)));'],ic);
 c_eval(['Pi_Di? = -sum(sum(Di(:,2:end).*Pa_i?_dev(:,2:end)));'],ic);
 
 %% 求平均
 
 c_eval('Ne?=irf_resamp(Ne?,Ne4);',1:3);
Ni1=irf_resamp(Ni1,Ni2);
Ne_average=(Ne1+Ne2+Ne3+Ne4)/4;
Ni_average=(Ni1+Ni2+Ni3+Ni4)/4;

c_eval('md_e?=Ne?(:, 2)  .*units.me;',ic);
c_eval('md_i?=Ni?(:, 2) .*units.mp;',ic);
md_e_average=(md_e1+md_e2+md_e3+md_e4)/4;
md_i_average=(md_i1+md_i2+md_i3+md_i4)/4;

c_eval('Pa_e?=Pa_e?(1:66,:);',1:3);
Pa_i1=Pa_i1(1:13,:);
Pa_e_average=(Pa_e1+Pa_e2+Pa_e3+Pa_e4)/4;
Pa_i_average=(Pa_i1+Pa_i2+Pa_i3+Pa_i4)/4;


%% Favre-filtered (density-weighted-filtered)
%Ve_ff,Vi_ff
c_eval('V_e?=V_e?(1:66,:);',1:3);
V_i1=V_i1(1:13,:);
c_eval('Ve_ff?_date=((Ne?(:,2:end).*V_e?(:,2:end))/4)./Ne_average(:,2:end);',ic);
c_eval('Ve_ff?=[V_e?(:,1), Ve_ff?_date];',ic);
c_eval('Vi_ff?_date=((Ni?(:,2:end).*V_i?(:,2:end))/4)./Ni_average(:,2:end);',ic);
c_eval('Vi_ff?=[V_i?(:,1), Vi_ff?_date];',ic);

%Ve_d_ff,Vi_d_ff
c_eval('Ve_d?_date=V_e?(:,2:end).*V_e?(:,2:end);',ic);
c_eval('Vi_d?_date=V_i?(:,2:end).*V_i?(:,2:end);',ic);
c_eval('Ve_d_ff?_date=((Ne?(:,2:end).*Ve_d?_date)/4)./Ne_average(:,2:end);',ic);
c_eval('Ve_d_ff?=[V_e?(:,1), Ve_ff?_date];',ic);
c_eval('Vi_d_ff?_date=((Ni?(:,2:end).*Vi_d?_date)/4)./Ni_average(:,2:end);',ic);
c_eval('Vi_d_ff?=[V_i?(:,1), Vi_ff?_date];',ic);

%Ve_s_ff,Vi_s_ff
c_eval('Ve_s_ff?_date=Ve_ff?(:,2:end).*Ve_ff?(:,2:end);',ic);
c_eval('Ve_s_ff?=[Ve_ff?(:,1), Ve_s_ff?_date];',ic);
c_eval('Vi_s_ff?_date=Vi_ff?(:,2:end).*Vi_ff?(:,2:end);',ic);
c_eval('Vi_s_ff?=[Vi_ff?(:,1), Vi_s_ff?_date];',ic);

%Be_ff,Bi_ff
c_eval('Be? = irf_resamp(B?,V_e?);',ic);
c_eval('Bi? = irf_resamp(B?,V_i?);',ic);
c_eval('Be_ff?_date=((Ne?(:,2:end).*Be?(:,2:end))/4)./Ne_average(:,2:end);',ic);
c_eval('Be_ff?=[Be?(:,1), Be_ff?_date];',ic);
c_eval('Bi_ff?_date=((Ni?(:,2:end).*Bi?(:,2:end))/4)./Ni_average(:,2:end);',ic);
c_eval('Bi_ff?=[Bi?(:,1), Bi_ff?_date];',ic);

%VexB_ff,VixB_ff
c_eval('VexB? = cross(V_e?(:,2:4),Be?(:,2:4));',ic);
c_eval('VixB? = cross(V_i?(:,2:4),Bi?(:,2:4));',ic);

c_eval('VexB_ff?_date=((Ne?(:,2:end).*VexB?)/4)./Ne_average(:,2:end);',ic);
c_eval('VexB_ff?=[VexB?(:,1), VexB_ff?_date];',ic);
c_eval('VixB_ff?_date=((Ni?(:,2:end).*VixB?)/4)./Ni_average(:,2:end);',ic);
c_eval('VixB_ff?=[VixB?(:,1), VixB_ff?_date];',ic);

%VexB_s_ff,VixB_s_ff
c_eval('VexB_s_ff_date? = cross(Ve_ff?(:,2:4),Be_ff?(:,2:4));',ic);
c_eval('VexB_s_ff?=[VexB?(:,1), VexB_s_ff_date?];',ic);
c_eval('VixB_s_ff_date? = cross(Vi_ff?(:,2:4),Bi_ff?(:,2:4));',ic);
c_eval('VixB_s_ff?=[VixB?(:,1), VixB_s_ff_date?];',ic);

%Ee_ff  Ei_ff

c_eval('Ee? = irf_resamp(E?,Ne?);',ic);
c_eval('Ei? = irf_resamp(E?,Ni?);',ic);
c_eval('Ee_ff?_date=((Ne?(:,2:end).*Ee?(:,2:end))/4)./Ne_average(:,2:end);',ic);
c_eval('Ee_ff?=[Ee?(:,1), Ee_ff?_date];',ic);
c_eval('Ei_ff?_date=((Ni?(:,2:end).*Ei?(:,2:end))/4)./Ni_average(:,2:end);',ic);
c_eval('Ei_ff?=[Ei?(:,1), Ei_ff?_date];',ic);

%Te_u,Ti_u

c_eval('Te?_u_date=Ve_d_ff?_date-Ve_s_ff?_date;',ic);
c_eval('Te?_u=[Ve_ff?(:,1), Te?_u_date];',ic);
c_eval('Ti?_u_date=Vi_d_ff?_date-Vi_s_ff?_date;',ic);
c_eval('Ti?_u=[Vi_ff?(:,1), Ti?_u_date];',ic);

%Te_b,Ti_b
c_eval('Te?_b_date=VexB_ff?_date-VexB_s_ff_date?;',ic);
c_eval('Te?_b=[Ve_ff?(:,1), Te?_b_date];',ic);
c_eval('Ti?_b_date=VixB_ff?_date-VixB_s_ff_date?;',ic);
c_eval('Ti?_b=[Vi_ff?(:,1), Ti?_b_date];',ic);

%%  lable of  energy transfer across scales
%Me，Mi
c_eval('R? = irf_resamp(R?,Ve_ff?);')
gradVe_ff=c_4_grad('R?','Ve_ff?','grad');
gradVe_ff=gradVe_ff(1:66,:);
c_eval('M_ae?_date=-md_e_average.*(Te?_u_date(:,1).*gradVe_ff(:,2)+Te?_u_date(:,2).*gradVe_ff(:,6)+Te?_u_date(:,3).*gradVe_ff(:,10));',ic);
c_eval('R? = irf_resamp(R?,Vi_ff?);')
gradVi_ff=c_4_grad('R?','Vi_ff?','grad');
gradVi_ff=gradVi_ff(1:13,:);
c_eval('M_ai?_date=-md_i_average.*(Ti?_u_date(:,1).*gradVi_ff(:,2)+Ti?_u_date(:,2).*gradVi_ff(:,6)+Ti?_u_date(:,3).*gradVi_ff(:,10));',ic);

c_eval('M_be?_date=-(units.e/units.c).*Ne_average(:,2).*(Te?_b_date(:,1).*Ve_ff?_date(:,1)+Te?_b_date(:,2).*Ve_ff?_date(:,2)+Te?_b_date(:,3).*Ve_ff?_date(:,3));',ic);
c_eval('M_bi?_date=(units.e/units.c).*Ni_average(:,2).*(Ti?_b_date(:,1).*Vi_ff?_date(:,1)+Ti?_b_date(:,2).*Vi_ff?_date(:,2)+Ti?_b_date(:,3).*Vi_ff?_date(:,3));',ic);

c_eval('Me?_date=M_ae?_date-M_be?_date;',ic);
c_eval('Me?=[Ve_ff?(:,1), Me?_date];',ic);
c_eval('Mi?_date=M_ai?_date-M_bi?_date;',ic);
c_eval('Mi?=[Vi_ff?(:,1), Mi?_date];',ic);

%Ae, Ai

c_eval('Ae?_date=-units.e.*Ne_average(:,2).*(Ee_ff?_date(:,1).*Ve_ff?_date(:,1)+Ee_ff?_date(:,2).*Ve_ff?_date(:,2)+Ee_ff?_date(:,3).*Ve_ff?_date(:,3));',ic);
c_eval('Ai?_date=units.e.*Ni_average(:,2).*(Ei_ff?_date(:,1).*Vi_ff?_date(:,1)+Ei_ff?_date(:,2).*Vi_ff?_date(:,2)+Ei_ff?_date(:,3).*Vi_ff?_date(:,3));',ic);
c_eval('Ae?=[Ve_ff?(:,1), Ae?_date];',ic);
c_eval('Ai?=[Vi_ff?(:,1), Ai?_date];',ic);

%Phi_e,Phi_i

Phi_e_x=-(Pa_e_average(:,2).*gradVe_ff(:,2)+Pa_e_average(:,3).*gradVe_ff(:,3)+Pa_e_average(:,4).*gradVe_ff(:,4));
Phi_e_y=-(Pa_e_average(:,5).*gradVe_ff(:,5)+Pa_e_average(:,6).*gradVe_ff(:,6)+Pa_e_average(:,7).*gradVe_ff(:,7));
Phi_e_z=-(Pa_e_average(:,8).*gradVe_ff(:,8)+Pa_e_average(:,9).*gradVe_ff(:,9)+Pa_e_average(:,10).*gradVe_ff(:,10));
Phi_e_date=Phi_e_x+Phi_e_y+Phi_e_z;
Phi_e=[Ve_ff1(:,1), Phi_e_date];

Phi_i_x=-(Pa_i_average(:,2).*gradVi_ff(:,2)+Pa_i_average(:,3).*gradVi_ff(:,3)+Pa_i_average(:,4).*gradVi_ff(:,4));
Phi_i_y=-(Pa_i_average(:,5).*gradVi_ff(:,5)+Pa_i_average(:,6).*gradVi_ff(:,6)+Pa_i_average(:,7).*gradVi_ff(:,7));
Phi_i_z=-(Pa_i_average(:,8).*gradVi_ff(:,8)+Pa_i_average(:,9).*gradVi_ff(:,9)+Pa_i_average(:,10).*gradVi_ff(:,10));
Phi_i_date=Phi_i_x+Phi_i_y+Phi_i_z;
Phi_i=[Vi_ff1(:,1), Phi_i_date];


%% Init figure
n=8;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


%% B plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([B?(:,1) B?(:,2)], 'color','b', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([B?(:,1) B?(:,3)], 'color','g', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([B?(:,1) B?(:,4)], 'color','r', 'Linewidth',0.75);",1); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",1); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",1);
set(gca,'Ylim',[-8e-9 15e-9]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [T]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;

%% Electric field
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",1);hold on;
c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",1);hold on;
c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",1);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",1); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% c_eval("set(gca,'Ylim',[min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))])-0.5 max([max(E?(:,2)) max(E?(:,3)) max(E?(:,4))])+0.5]);",1);
set(gca,'Ylim',[-8e-3 15e-3]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [V/m]','fontsize',8)
i=i+1;
h(i)=irf_subplot(n,1,-i);

%% Me plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Me1(:,1) Me1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Me2(:,1) Me2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Me3(:,1) Me3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Me4(:,1) Me4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Me?(:,1) 0*Me?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(Me?(:,2))]) max(Me?(:,2))]);",ic);
set(gca,'Ylim',[-15e-12 10e-12]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Me ','fontsize',8);
i=i+1;

%% Mi plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Mi1(:,1) Mi1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Mi2(:,1) Mi2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Mi3(:,1) Mi3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Mi4(:,1) Mi4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Mi?(:,1) 0*Mi?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(Mi?(:,2))]) max(Mi?(:,2))]);",ic);
set(gca,'Ylim',[-15e-12 10e-12]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Mi ','fontsize',8);
i=i+1;

%% Ae plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Ae1(:,1) Ae1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Ae2(:,1) Ae2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Ae3(:,1) Ae3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Ae4(:,1) Ae4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Ae?(:,1) 0*Ae?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(Ae?(:,2))]) max(Ae?(:,2))]);",ic);
set(gca,'Ylim',[-15e-12 15e-12]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ae ','fontsize',8);
i=i+1;

%% Ai plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Ai1(:,1) Ai1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Ai2(:,1) Ai2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Ai3(:,1) Ai3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Ai4(:,1) Ai4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Ai?(:,1) 0*Ai?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(Ai?(:,2))]) max(Ai?(:,2))]);",ic);
set(gca,'Ylim',[-15e-12 15e-12]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ai ','fontsize',8);
i=i+1;

%% Phi_e plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Phi_e(:,1) Phi_e(:,2)], 'color','k', 'Linewidth',0.75); hold on;

%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Phi_e(:,1) 0*Phi_e(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% set(gca,'Ylim',[min([min(Phi_e(:,2))])-0.1 max(Phi_e(:,2))+0.1]);
set(gca,'Ylim',[-15e-10 10e-10]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
% set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Phi_e ','fontsize',8);
i=i+1;

%% Phi_i plot
h(i)=irf_subplot(n,1,-i);
irf_plot([Phi_i(:,1) Phi_i(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Phi_i2(:,1) Phi_i2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Phi_i3(:,1) Phi_i3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([Phi_i4(:,1) Phi_i4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([Phi_i(:,1) 0*Phi_i(:,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% set(gca,'Ylim',[min([min(Phi_i(:,2))])-0.1 max(Phi_i(:,2))+0.1]);
set(gca,'Ylim',[-15e-10 10e-10]);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
% set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Phi_i ','fontsize',8);
i=i+1;

irf_zoom(tint,'x',h(1:n));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)

%%  出图保存部分
colormap(jet)
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

% cd  C:\Matlab\bin\新建文件夹\fwd\
% rmdir(TempDir,'s'); 
% figname = 'overview';
%     print(gcf, '-dpdf', [figname '.pdf']);   