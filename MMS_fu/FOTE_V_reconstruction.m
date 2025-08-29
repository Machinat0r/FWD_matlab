% gsm
clear;close all
clc
global ParentDir 
ParentDir = '/Users/fwd/Documents/MATLAB/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
%% VORTEX_9图全版
%% load data
iic = 1;
ic=1:4;
hh = 10;
hh1= 5;
thresold = 0.2;  
dspan = 3;
dspan2 = 1;
% Tsta = '2017-07-06T01:44:50.000Z';
% Tend = '2017-07-06T01:45:10.000Z';

% Tsta = '2017-07-12T11:54:33.600Z';
% Tend = '2017-07-12T11:54:35.400Z';
% Tsta = '2018-07-06T12:37:04.500Z';
% Tend = '2018-07-06T12:37:05.250Z';
% Tsta = '2017-05-22T08:23:52.000Z';
% Tend = '2017-05-22T08:23:55.000Z';
% Tsta = '2017-05-25T04:40:29.500Z';
% Tend = '2017-05-25T04:40:33.000Z';
Tsta = '2019-08-05T16:24:00.000Z';
Tend = '2019-08-05T16:25:00.000Z';
tint = irf.tint(Tsta,Tend);
% tint=irf.tint('2017-01-27T12:05:42.50Z/2017-01-27T12:05:43.80Z');
% TT = '2017-07-12T11:54:00.000Z/2017-07-12T11:55:00.000Z';
TT = '2019-08-05T16:24:00.000Z/2019-08-05T16:25:00.000Z';
% % % TT = '2017-07-06T17:31:52.000Z/2017-07-06T17:32:07.000Z';
% TT = '2018-07-06T12:36:30.000Z/2018-07-06T12:38:00.000Z';
% TT = '2017-05-22T08:23:30.000Z/2017-05-22T08:24:30.000Z';
% TT = '2017-05-25T04:40:00.000Z/2017-05-25T04:41:00.000Z';

TT3 = '2019-08-05T16:24:29.000Z/2019-08-05T16:24:33.000Z';
tint2=irf.tint(TT);
tint3 = irf.tint(TT3);

Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

% filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist');
% % filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames = [filenames1, filenames2,filenames4];
% % % % 
% [filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);
% 
% SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 48, 'CheckSize', 0)
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir)
Time='2019-08-05T16:24:31.137Z';
time = irf_time(Time,'utc>epochtt');
c_eval('e_r? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energy_brst'',tint);',ic);

%% load B Fields data
c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
c_eval('B? = irf_resamp(B?, B1);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);

%% smooth
    c_eval('BS?(:,2)=smooth(B?(:,2),hh);',ic);
    c_eval('BS?(:,3)=smooth(B?(:,3),hh);',ic);
    c_eval('BS?(:,4)=smooth(B?(:,4),hh);',ic);
    c_eval('BS?(:,1)=smooth(B?(:,1),hh);',ic);
    c_eval('BS?=irf_abs(BS?);',ic);

c_eval('Ne?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval('ne?=irf.ts2mat(Ne?);',ic);

c_eval('Ni?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval('ni?=irf.ts2mat(Ni?);',ic);
c_eval('Vi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); 
% c_eval('Vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval('vi?=irf.ts2mat(Vi?);',ic);
c_eval('vi? = irf_gse2gsm(vi?);',ic);
%% load Ve Fields data
c_eval('Vexyz?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval('Vexyz? = irf_gse2gsm(Vexyz?);',ic);
% L=[0.05 -0.21 0.98]; % Vi
% M=[-0.23 0.95 0.22];
% N=[-0.97 -0.24 -0.01];
L=[1,0,0];M=[0,1,0];N=[0,0,1];
c_eval('Ve_LMN?=irf_newxyz(Vexyz?,L,M,N);',ic);
c_eval('Ve?=irf.ts2mat(Vexyz?);',ic);

c_eval('B? = irf_resamp(B?, Ve?);',ic);

c_eval('Ve?=irf_resamp(Ve?,B?);',ic);
c_eval('vi?=irf_resamp(vi?,B?);',ic);
c_eval('VeS?=irf.ts2mat(Vexyz?);',ic);
c_eval('VeS_LMN?=irf.ts2mat(Ve_LMN?);',ic);
c_eval('VeS? = irf_resamp(VeS?, VeS1);',ic);
c_eval('VeS_LMN? = irf_resamp(VeS_LMN?, VeS_LMN1);',ic);
%% smooth_step1
kk = size(VeS1,1);
if mod(kk,2) == 1
   le = kk-2; lo = kk-3;%even偶数（这里当它是奇数），odd奇数
else
   le = kk-3; lo = kk-2;%even偶数，odd奇数
end

if e_r1.data(1)>e_r1.data(2)
    for ii = 1:2:le
        c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
        c_eval('VeS_LMN?(ii+1,2:4)=(VeS_LMN?(ii+2,2:4)+VeS_LMN?(ii,2:4))/2;',ic)
    end
else
    for ii = 2:2:lo
        c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
        c_eval('VeS_LMN?(ii+1,2:4)=(VeS_LMN?(ii+2,2:4)+VeS_LMN?(ii,2:4))/2;',ic)
    end
end

%% resample
ic=1:4;
c_eval('VeS?=irf_resamp(VeS?,B?);',ic);
c_eval('VeS?=irf_abs(VeS?);',ic);
Vup=max(VeS2(:,5));
c_eval('VeS_LMN?=irf_resamp(VeS_LMN?,B?);',ic);
c_eval('Vet?=Vexyz?.abs;',ic);

%% smooth_step2
    c_eval('VeSS?(:,1)=VeS?(:,1);',ic);
    c_eval('VeSS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeSS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeSS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeSS?(:,5)=smooth(VeS?(:,5),hh1);',ic);
    
    c_eval('VeSS?(:,1)=VeS?(:,1);',ic);
    c_eval("temp?(:,2)=smooth(VeSS?(:,2),hh1);",ic)
    c_eval('VeSS?(:,2)=temp?(:,2);',ic)
    c_eval("temp?(:,3)=smooth(VeSS?(:,3),hh1);",ic)
    c_eval('VeSS?(:,3)=temp?(:,3);',ic)
    c_eval("temp?(:,4)=smooth(VeSS?(:,4),hh1);",ic)
    c_eval('VeSS?(:,4)=temp?(:,4);',ic)

    c_eval("temp?(:,2)=smooth(VeSS?(:,2),hh1);",ic)
    c_eval('VeSS?(:,2)=temp?(:,2);',ic)
    c_eval("temp?(:,3)=smooth(VeSS?(:,3),hh1);",ic)
    c_eval('VeSS?(:,3)=temp?(:,3);',ic)
    c_eval("temp?(:,4)=smooth(VeSS?(:,4),hh1);",ic)
    c_eval('VeSS?(:,4)=temp?(:,4);',ic)
    % % % 
    % % % c_eval("temp?(:,2)=smooth(VeSS?(:,2),hh1);",ic)
    % % % c_eval('VeSS?(:,2)=temp?(:,2);',ic)
    % % % c_eval("temp?(:,3)=smooth(VeSS?(:,3),hh1);",ic)
    % % % c_eval('VeSS?(:,3)=temp?(:,3);',ic)
    % % % c_eval("temp?(:,4)=smooth(VeSS?(:,4),hh1);",ic)
    % % % c_eval('VeSS?(:,4)=temp?(:,4);',ic)

ic=1:4;    
c_eval('Rgse?=mms.get_data(''R_gsm'',tint2,?);',ic);
c_eval('R?=irf.ts2mat(Rgse?);',ic);

%% FOTE method
VeS2=irf_resamp(VeS2,VeS1);
VeS3=irf_resamp(VeS3,VeS1);
VeS4=irf_resamp(VeS4,VeS1);
VeSS2=irf_resamp(VeSS2,VeSS1);
VeSS3=irf_resamp(VeSS3,VeSS1);
VeSS4=irf_resamp(VeSS4,VeSS1);

tint=[iso2epoch(Tsta) iso2epoch(Tend)];
for ic=1:4
  c_eval(['R?=irf_resamp(R?,VeSS?);'],ic);
  c_eval(['VeSS?=irf_tlim(VeSS?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
end

%% FOTE method
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
% Ve1=irf_resamp(Ve1,B1);
Ve2=irf_resamp(Ve2,Ve1);
Ve3=irf_resamp(Ve3,Ve1);
Ve4=irf_resamp(Ve4,Ve1);
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
for ic=1:4
  c_eval(['R?=irf_resamp(R?,Ve?);'],ic);
  c_eval(['Ve?=irf_tlim(Ve?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
  c_eval(['VemagC?=irf_abs(Ve?);'],ic);
end
gradB=c_4_grad('R?','VeSS?','grad');
d_dot_B=c_4_grad('R?','Ve?','div');
d_cros_B=c_4_grad('R?','Ve?','curl');


%error of curolmeter
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

temp=irf_abs(d_cros_B);
d_cros_B_mag=temp(:,[1 5]);
err_curlmeter=irf_multiply(1,d_dot_B,1,d_cros_B_mag,-1);
err_curlmeter(:,2)=abs(err_curlmeter(:,2))*100;

idxnull=find(gradB(:,1,1)>=iso2epoch(Time)); idxnull=idxnull(1);
dB_null=reshape(gradB(idxnull,2:end),3,3);
[V,D] = eig(dB_null);