clc
clear
close all
mms.db_init('local_file_db','/Volumes/FWD-ResourceDisk/MMS')
Tcomputsta=clock;
%--------------------------------------------
%% load data
tintStr = '19-Sep-2015';
Tsta='2015-09-19T07:43:30.45Z'; 
Tend='2015-09-19T07:43:31.65Z'; 
tint=irf.tint('2015-09-19T07:43:30.12Z/2015-09-19T07:43:31.62Z');
tint2=irf.tint('2015-09-19T07:43:30.00Z/2015-09-19T07:43:30.30Z');
Time='2015-09-19T07:43:31.115Z';
% tintStr = '05-Agu-2019';
% Tsta='2015-09-19T07:43:30.00Z'; 
% Tend='2015-09-19T07:43:32.00Z'; 
% tint=irf.tint('2015-11-20T09:07:34.500Z/2015-11-20T09:07:35.500Z');
% tint2=irf.tint('2015-11-20T09:07:34.900Z/2015-11-20T09:07:35.000Z');
% Time='2015-09-19T07:43:31.115Z';

ic=1:4;
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
c_eval('Bt?=irf.ts2mat(Bt?);',ic);

Pos = mms.get_data('R_gse',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gseR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
%     c_eval('Rgse?=mms.get_data(''R_gse'',tint2,?);',ic);
%     c_eval('R?=irf.ts2mat(Rgse?);',ic);
c_eval(['R?=irf_resamp(R?,B?);']);



%% Null type identification
[Null_info_dis, Null_info_type, Null_info_err] = FOTE_Taylor_Expansion('R?','B?');

% if the raw B data has fluctuations with small time scale, smooth the data
% with specific inputs:
% smooth_span = 2;
% [Null_info_dis, Null_info_type, Null_info_err] = FOTE_Taylor_Expansion('R?','B?',smooth_span);

% smooth_span = 2; smooth_method = 'moving';
% [Null_info_dis, Null_info_type, Null_info_err] = FOTE_Taylor_Expansion('R?','B?',smooth_span, smooth_method);


%% calculate the gradient matrix and SC-null location
% for radial nulls, default e1 in Trans_mat is the eigenvector corresponding to the
% smallest eigenvalues of gradB; for spral nulls, e1 is the eigenvector corresponding to the
% real eigenvalue;

[Null_loc, Trans_mat, dB_null] = FOTE('R?','B?',Time);

%% rotate the reconstructed topology by changing the coordinate system (if needed)
% Relative to the raw data coordinares (GSE/GSM), e1,e2,e3 complete a new
% transform Cartesian coordinates.

% [Null_loc, dB_null] = FOTE_coordupdate('R?','B?',Tnull,e1,e2,e3); 

%% select the points from which you want to trace the magnetic field lines
X_focus = [15:10:35];
Y_focus = [0:10:150];
Z_focus = [0:10:330];
[h,focus_points] = FOTE_fp_select(X_focus, Y_focus, Z_focus, 'Spherical');

%% reconstrcut the field topology with gradient matrix, SC-null location and focus points.
[h2,FOTE_res] = FOTE_reconstruction(dB_null, Null_loc, focus_points, 150, [-37.5,30], 'cline',0,75);

set(gca,'Xlim',[-100 100], 'Ylim',[-100 100], 'Zlim',[-100 100]);
set(gca,'xtick',[-100:50:100], 'ytick',[-100:50:100], 'ztick',[-100:50:100],'fontsize',12);
set(gca,'DataAspectRatio',[1 1.0 1]);
xlabel(gca,'e_{1} [km]');
ylabel(gca,'e_{2} [km]');
zlabel(gca,'e_{3} [km]');
%%
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
% figname = ['C:\Users\fwd\Desktop\Ti~mor~\M\magnetic_monopole\supplementary\' ,'20150919_reconstruction'];
% print(gcf, '-dpdf', [figname '.pdf']);