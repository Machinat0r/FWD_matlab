%% 找涡旋之零点类型 33.485时刻 平滑 （d）-（i）图数据全部做处理
%% 2/3程序
% gsm
clear;close all
clc
global ParentDir 
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
%% VORTEX_9图全版
%% load data
ic=1:4;
hh = 10;
hh1= 10;
thresold = 0.4;  
% Tsta = '2017-07-06T01:44:50.000Z';
% Tend = '2017-07-06T01:45:10.000Z';

Tsta = '2017-07-12T11:54:33.600Z';
Tend = '2017-07-12T11:54:35.400Z';
% Tsta = '2017-07-06T17:31:55.000Z';
% Tend = '2017-07-06T17:32:05.000Z';
tint = irf.tint(Tsta,Tend);
% tint=irf.tint('2017-01-27T12:05:42.50Z/2017-01-27T12:05:43.80Z');
TT = '2017-07-12T11:54:00.000Z/2017-07-12T11:55:00.000Z';
tint2=irf.tint(TT);

Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames = [filenames1, filenames2];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 48, 'CheckSize', 0)
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db','D:/MMS/')

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
c_eval('Vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval('vi?=irf.ts2mat(Vi?);',ic);
c_eval('vi? = irf_gse2gsm(vi?);',ic);
%% Te
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);%GSE坐标下的B
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);%GSM坐标下的B
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)%GSE坐标下的Te
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)%FAC坐标下的Te
%GSE坐标系下的等离子体温度张量数据旋转到了磁力线坐标系（FAC，Field Aligned Coordinates）
% c_eval('matParTe? = facTe?.xx.resample(gsmB?.time).data;',ic)
% %FAC坐标系下等离子体温度张量的平行（Parallel）分量，.xx是张量的平行分量。
% c_eval('matPerTe?= (facTe?.yy.resample(gsmB?.time).data + facTe?.zz.resample(gsmB?.time).data)/2;',ic)
% %FAC坐标系下等离子体温度张量的垂直（Perpendicular）分量，取yy和zz分量的平均值计算垂直分量。
%温度先不重采样
c_eval('matParTe?(:,2) = facTe?.xx.data;',ic)
%FAC坐标系下等离子体温度张量的平行（Parallel）分量，.xx是张量的平行分量。
c_eval('matPerTe?(:,2)= (facTe?.yy.data + facTe?.zz.data)/2;',ic)
c_eval('matParTe?(:,1)= ne?(:,1);',ic)
c_eval('matPerTe?(:,1)= ne?(:,1);',ic)
c_eval('ne? = irf_resamp(ne?, ne1);',ic);
c_eval('matParTe? = irf_resamp(matParTe?, matParTe1);',ic);
c_eval('matPerTe? = irf_resamp(matPerTe?, matPerTe1);',ic);
%% load Ve Fields data
c_eval('Vexyz?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval('Vexyz? = irf_gse2gsm(Vexyz?);',ic);
% L=[0.05 -0.21 0.98]; % Vi
% M=[-0.23 0.95 0.22];
% N=[-0.97 -0.24 -0.01];
L=[1,0,0];M=[0,1,0];N=[0,0,1];
c_eval('Ve_LMN?=irf_newxyz(Vexyz?,L,M,N);',ic);
c_eval('Ve?=irf.ts2mat(Vexyz?);',ic);
c_eval('Ve?=irf_resamp(Ve?,B?);',ic);
c_eval('vi?=irf_resamp(vi?,B?);',ic);
c_eval('VeS?=irf.ts2mat(Vexyz?);',ic);
c_eval('VeS_LMN?=irf.ts2mat(Ve_LMN?);',ic);
c_eval('VeS? = irf_resamp(VeS?, VeS1);',ic);
c_eval('VeS_LMN? = irf_resamp(VeS_LMN?, VeS_LMN1);',ic);
%% smooth_step1
kk = length (e_r1.data);
if mod(kk,2) == 1
   le = kk-2; lo = kk-3;%even偶数（这里当它是奇数），odd奇数
else
   le = kk-3; lo = kk-2;%even偶数，odd奇数
end
for ic =1:4
    c_eval('e_r = e_r?;',ic)
    if e_r.data(1)>e_r.data(2)
        for ii = 1:2:le
            c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('matParTe?(ii+1,2)=(matParTe?(ii+2,2)+matParTe?(ii,2))/2;',ic)
            c_eval('matPerTe?(ii+1,2)=(matPerTe?(ii+2,2)+matPerTe?(ii,2))/2;',ic)
            c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
            c_eval('VeS_LMN?(ii+1,2:4)=(VeS_LMN?(ii+2,2:4)+VeS_LMN?(ii,2:4))/2;',ic)
        end
    else
        for ii = 2:2:lo
            c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('matParTe?(ii+1,2)=(matParTe?(ii+2,2)+matParTe?(ii,2))/2;',ic)
            c_eval('matPerTe?(ii+1,2)=(matPerTe?(ii+2,2)+matPerTe?(ii,2))/2;',ic)
            c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
            c_eval('VeS_LMN?(ii+1,2:4)=(VeS_LMN?(ii+2,2:4)+VeS_LMN?(ii,2:4))/2;',ic)
        end
    end
end
%% resample
ic=1:4;
c_eval('ne?=irf_resamp(ne?,B?);',ic);
c_eval('ni?=irf_resamp(ni?,B?);',ic);
c_eval('matParTe?=irf_resamp(matParTe?,B?);',ic);
c_eval('matPerTe?=irf_resamp(matPerTe?,B?);',ic);
c_eval('VeS?=irf_resamp(VeS?,B?);',ic);
c_eval('VeS?=irf_abs(VeS?);',ic);
Vup=max(VeS2(:,5));
c_eval('VeS_LMN?=irf_resamp(VeS_LMN?,B?);',ic);
c_eval('Vet?=Vexyz?.abs;',ic);

% for ic =1:4
%     c_eval("temp?(:,2)=smooth(VeS?(:,2),0.018,'rloess');",ic)
%     c_eval('VeS?(:,2)=temp?(:,2);',ic)
%     c_eval("temp?(:,3)=smooth(VeS?(:,3),0.018,'rloess');",ic)
%     c_eval('VeS?(:,3)=temp?(:,3);',ic)
%     c_eval("temp?(:,4)=smooth(VeS?(:,4),0.018,'rloess');",ic)
%     c_eval('VeS?(:,4)=temp?(:,4);',ic)
% end
%% smooth_step2
    c_eval('VeSS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeSS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeSS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeSS?(:,5)=smooth(VeS?(:,5),hh1);',ic);
    
%     c_eval('VeS_LMN?(:,2)=smooth(VeS_LMN?(:,2),hh1);',ic);
%     c_eval('VeS_LMN?(:,3)=smooth(VeS_LMN?(:,3),hh1);',ic);
%     c_eval('VeS_LMN?(:,4)=smooth(VeS_LMN?(:,4),hh1);',ic);
%     c_eval('VeS_LMN?(:,2)=smooth(VeS_LMN?(:,2),hh1);',ic);
%     c_eval('VeS_LMN?(:,3)=smooth(VeS_LMN?(:,3),hh1);',ic);
%     c_eval('VeS_LMN?(:,4)=smooth(VeS_LMN?(:,4),hh1);',ic);
%     c_eval('VeS_LMN?(:,2)=smooth(VeS_LMN?(:,2),hh1);',ic);
%     c_eval('VeS_LMN?(:,3)=smooth(VeS_LMN?(:,3),hh1);',ic);
%     c_eval('VeS_LMN?(:,4)=smooth(VeS_LMN?(:,4),hh1);',ic);
%     c_eval('VeS_LMN?(:,2)=smooth(VeS_LMN?(:,2),hh1);',ic);
%     c_eval('VeS_LMN?(:,3)=smooth(VeS_LMN?(:,3),hh1);',ic);
%     c_eval('VeS_LMN?(:,4)=smooth(VeS_LMN?(:,4),hh1);',ic);
    
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
    c_eval("temp?(:,2)=smooth(VeSS?(:,2),hh1);",ic)
    c_eval('VeSS?(:,2)=temp?(:,2);',ic)
    c_eval("temp?(:,3)=smooth(VeSS?(:,3),hh1);",ic)
    c_eval('VeSS?(:,3)=temp?(:,3);',ic)
    c_eval("temp?(:,4)=smooth(VeSS?(:,4),hh1);",ic)
    c_eval('VeSS?(:,4)=temp?(:,4);',ic)

ic=1:4;    
c_eval('Exyz?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('Exyz? = irf_gse2gsm(Exyz?);',ic);
c_eval('Et?=Exyz?.abs;',ic);
c_eval('E?=irf.ts2mat(Exyz?);',ic);
% c_eval('Vit?=irf.ts2mat(Vi?);',ic);
c_eval('E?=irf_resamp(E?,B?);',ic);

c_eval('Rgse?=mms.get_data(''R_gsm'',tint2,?);',ic);
c_eval('R?=irf.ts2mat(Rgse?);',ic);

c_eval('dsp?(:,1)=B?(:,1);',ic);
c_eval('dsp?(:,2:4)=E?(:,2:4)+cross(Ve?(:,2:4),B?(:,2:4))/1000;',ic);
% c_eval('dsp?=irf_resamp(dsp?,Vit?);',ic);
c_eval('dspt?=irf_abs(dsp?);',ic);

%% FOTE method
VeS2=irf_resamp(VeS2,VeS1);
VeS3=irf_resamp(VeS3,VeS1);
VeS4=irf_resamp(VeS4,VeS1);
VeSS2=irf_resamp(VeSS2,VeSS1);
VeSS3=irf_resamp(VeSS3,VeSS1);
VeSS4=irf_resamp(VeSS4,VeSS1);
ne2=irf_resamp(ne2,ne1);
ne3=irf_resamp(ne3,ne1);
ne4=irf_resamp(ne4,ne1);
ni2=irf_resamp(ni2,ni1);
ni3=irf_resamp(ni3,ni1);
ni4=irf_resamp(ni4,ni1);
vi2=irf_resamp(vi2,vi1);
vi3=irf_resamp(vi3,vi1);
vi4=irf_resamp(vi4,vi1);

tint=[iso2epoch(Tsta) iso2epoch(Tend)];
for ic=1:4
  c_eval(['R?=irf_resamp(R?,VeSS?);'],ic);
  c_eval(['VeSS?=irf_tlim(VeSS?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
  c_eval(['VemagC?=irf_abs(VeSS?);'],ic);
  c_eval('VeS_n?(:,1)=VeSS?(:,1);',ic);
  c_eval('VeS_n?(:,2)=VeSS?(:,2).*ne?(:,2);',ic);
  c_eval('VeS_n?(:,3)=VeSS?(:,3).*ne?(:,2);',ic);
  c_eval('VeS_n?(:,4)=VeSS?(:,4).*ne?(:,2);',ic);
  c_eval('VeS_n?(:,5)=VeSS?(:,5).*ne?(:,2);',ic);
  c_eval('j?(:,1)=VeS?(:,1);',ic);
  c_eval('j?(:,2:4)=(ni?(:,2).*vi?(:,2:4)-ne?(:,2).*VeS?(:,2:4))*1.6/1e10;',ic);
end
n1=length(j1(:,1));
for ij=1:n1
   jmean(ij,1)=j1(ij,1);
   jmean(ij,2)=mean(j1(ij,2)+j2(ij,2)+j3(ij,2)+j4(ij,2));
   jmean(ij,3)=mean(j1(ij,3)+j2(ij,3)+j3(ij,3)+j4(ij,3));
   jmean(ij,4)=mean(j1(ij,4)+j2(ij,4)+j3(ij,4)+j4(ij,4)); 
end
gradVe=c_4_grad('R?','VeSS?','grad');
d_dot_Ve=c_4_grad('R?','VeSS?','div');
d_cros_Ve=c_4_grad('R?','VeSS?','curl');


%error of curolmeter
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','VeS_n?');
temp=irf_abs(j);%t\x\y\z\total
jmag=temp(:,[1 5]);%t\total
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');

% temp=irf_abs(j);%t\x\y\z\total
% jmag=temp(:,[1 5]);%t\total
% err_4C=irf_multiply(1,divB,1,jmag,-1);
% err_4C(:,2)=abs(err_4C(:,2))*100;

temp=irf_abs(d_cros_Ve);
d_cros_B_mag=temp(:,[1 5]);
err_curlmeter=irf_multiply(1,d_dot_Ve,1,d_cros_B_mag,-1);
err_curlmeter(:,2)=abs(err_curlmeter(:,2))*100;


%Null type identification
for ii=1:length(d_dot_Ve(:,1))
    mksizSim(ii)=4;
    
    deltB_null=reshape(gradVe(ii,2:end),3,3);
    [V,D] = eig(deltB_null);
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='w';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            type(ii)='X'; clr(ii)='k'; faceclr(ii)='w';
        end
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='r';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            type(ii)='o'; clr(ii)='k'; faceclr(ii)='w';
        end
    end
    %=========================================================
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='w';
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if min(abs([D(1,1) D(2,2) D(3,3)]))/max(abs([D(1,1) D(2,2) D(3,3)]))<thresold
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
                typeSim(ii)='X'; clrSim(ii)='b'; faceclrSim(ii)='w'; mksizSim(ii)=7;
            end
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='X'; clrSim(ii)='r'; faceclrSim(ii)='w'; mksizSim(ii)=7;
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            typeSim(ii)='X'; clrSim(ii)='k'; faceclrSim(ii)='w'; mksizSim(ii)=7;
        end
        %------------------------------------------------------------------
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='r';
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))/max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) < (thresold*2)
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
                typeSim(ii)='o'; clrSim(ii)='b'; faceclrSim(ii)='w';
            end
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='o'; clrSim(ii)='r'; faceclrSim(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            typeSim(ii)='o'; clrSim(ii)='k'; faceclrSim(ii)='w';
        end
        %------------------------------------------------------------------
    end
    %=========================================================
    
    eigVal_err(ii,2)=abs(real(D(1,1)+D(2,2)+D(3,3)))/max(abs([real(D(1,1)), real(D(2,2)), real(D(3,3))])) * 100;
    sumeigVal(ii,2)=D(1,1)+D(2,2)+D(3,3);
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;
end
eigVal_err(:,1)=d_dot_Ve(:,1);
sumeigVal(:,1)=d_dot_Ve(:,1);
eigVal_err_v2(:,1)=d_dot_Ve(:,1);


%find null position
for ii=1:length(VeSS1(:,1))
    dBeach=reshape(gradVe(ii,2:end),3,3);
    dR1(ii,2:4)=VeSS1(ii,2:4)*inv(dBeach');
    dR2(ii,2:4)=VeSS2(ii,2:4)*inv(dBeach');
    dR3(ii,2:4)=VeSS3(ii,2:4)*inv(dBeach');
    dR4(ii,2:4)=VeSS4(ii,2:4)*inv(dBeach');
end
dR1(:,1)=VeSS1(:,1);
dR2(:,1)=VeSS1(:,1);
dR3(:,1)=VeSS1(:,1);
dR4(:,1)=VeSS1(:,1);

dRmag1=irf_abs(dR1);
dRmag2=irf_abs(dR2);
dRmag3=irf_abs(dR3);
dRmag4=irf_abs(dR4);

dRmin(:,2)=min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1)=dRmag1(:,1);

pause(1)
%% coordinates transform
% L=[0.03 0.24 0.97];% B
% M=[0.1 0.96 -0.25];
% N=[-0.99 0.11 0.01];

% L=[0.05 -0.21 0.98]; % Vi
% M=[-0.23 0.95 0.22];
% N=[-0.97 -0.24 -0.01];
L=[1,0,0];M=[0,1,0];N=[0,0,1];
ic=1:4;
c_eval('B_LMN?=irf_newxyz(Bxyz?,L,M,N);',ic);
c_eval('B_lmn?=irf.ts2mat(B_LMN?);',ic);
c_eval('B_lmn? = irf_resamp(B_lmn?, B_lmn1);',ic);
c_eval('E_LMN?=irf_newxyz(Exyz?,L,M,N);',ic);
c_eval('E_lmn?=irf.ts2mat(E_LMN?);',ic);
c_eval('E_lmn?=irf_resamp(E_lmn?,B?);',ic);
c_eval('Vi_LMN?=irf_newxyz(Vi?,L,M,N);',ic);

[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
d_cros_ve=c_4_grad('R?','VeS?(:,1:4)','curl');
% d_cros_ve=c_4_grad('R?','VeSS?(:,1:4)','curl');
d_cros_ve_LMN=irf_newxyz(d_cros_ve,L,M,N);
n=length(j(:,1)); 
j_total=zeros(n,1);
j_E=zeros(n,1);
j_LMN(:,1)=j(:,1);
j_LMN(:,2:4)=irf_newxyz(j(:,2:4),L,M,N);
E_lmn2=irf_resamp(E_lmn2,j_LMN);


ic=1:4;
for ii=1:n
    j_total(ii,1)= norm(j(ii,2:4));
    d_cros_ve_LMN_abs(ii,1)=norm(d_cros_ve_LMN(ii,2:4));
    c_eval('E0?(ii,2:4)=E?(ii,2:4)+(cross(VeS?(ii,2:4),B?(ii,2:4))/1e3);',ic);
    c_eval('j_dot_E?(ii,1)= dot(E0?(ii,2:4),j?(ii,2:4));',ic);
%     Ve_cross_B(ii,2:4)=cross(VeS2(ii,2:4),B2(ii,2:4));
%     E01(ii,2:4)=E2(ii,2:4)+(Ve_cross_B(ii,2:4)/1e3);
end

for iii=1:n
    c1=dot(d_cros_ve_LMN(iii,2:4),B_lmn1(iii,2:4))/(norm(d_cros_ve_LMN(iii,2:4))*(norm(B_lmn1(iii,2:4))));
    angle_d_cros_ve1(iii,1)=acosd(c1);
end 

%新增阿尔芬速度计算（MMS2）
ne2=irf_resamp(ne2,B2);
ni2=irf_resamp(ni2,B2);
n2=length(ne2(:,1)); 
for i=1:n2
VAL(i,1)= ( B_lmn2(i,2)/1e9 )/(( (4*pi/1e7)*(ne2(i,2)*1e6*9.1/1e31+ni2(i,2)*1e6*1.67/1e27) )^(0.5));%算出来应是m/s
VAL(i,1)=VAL(i,1)/1e3;%（m/s）转（km/s）

VAM(i,1)= ( B_lmn2(i,3)/1e9 )/(( (4*pi/1e7)*(ne2(i,2)*1e6*9.1/1e31+ni2(i,2)*1e6*1.67/1e27) )^(0.5));%算出来应是m/s
VAM(i,1)=VAM(i,1)/1e3;%（m/s）转（km/s）

VAN(i,1)= ( B_lmn2(i,4)/1e9 )/(( (4*pi/1e7)*(ne2(i,2)*1e6*9.1/1e31+ni2(i,2)*1e6*1.67/1e27) )^(0.5));%算出来应是m/s
VAN(i,1)=VAN(i,1)/1e3;%（m/s）转（km/s）
end

for i=1:n
j_dot_E02(i,1)=dot(j_LMN(i,2:4),E_lmn2(i,2:4));
end

%% plot
h=irf_plot(10,'newfigure');
xSize=750; ySize=800;
set(gcf,'Position',[100 100 xSize ySize]);
mmscolors=[0 0 1; 0 1 0; 1 0 0; 0 0 0]; 
set(h,'ColorOrder',mmscolors);
lnwid = 1;

h(1)=irf_panel('delta Blmn');
hold(h(1),'on');
B2Lmean=mean(B_lmn2(:,2));
B2Mmean=mean(B_lmn2(:,3));
B2Nmean=mean(B_lmn2(:,4));
dB_lmn2(:,2)=B_lmn2(:,2)-B2Lmean;
dB_lmn2(:,3)=B_lmn2(:,3)-B2Mmean;
dB_lmn2(:,4)=B_lmn2(:,4)-B2Nmean;
% irf_plot(h(1),Bt2,'Linewidth',lnwid);
% irf_plot(h(1),[B_lmn2(:,1),dB_lmn2(:,2)],'Linewidth',lnwid);
% irf_plot(h(1),[B_lmn2(:,1),dB_lmn2(:,3)],'Linewidth',lnwid);
% irf_plot(h(1),[B_lmn2(:,1),dB_lmn2(:,4)],'Linewidth',lnwid);
irf_plot(h(1),[B_lmn1(:,1),B_lmn1(:,2)],'Linewidth',lnwid);
irf_plot(h(1),[B_lmn1(:,1),B_lmn1(:,3)],'Linewidth',lnwid);
irf_plot(h(1),[B_lmn1(:,1),B_lmn1(:,4)],'Linewidth',lnwid);
irf_plot(h(1),[B_lmn3(:,1) B_lmn3(:,2)*0],'k--', 'Linewidth',0.5);%%铏氱嚎
hold(h(1),'off');
ylabel(h(1),{'B','(nT)'},'Interpreter','tex');
% set(h(1),'ylim',[-7 7],'ytick',[-5:5:5]);
irf_legend(h(1),{'Bx  ','By  ','Bz  '},[0.98 0.03],'fontsize',12)
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',12)
grid(h(1),'off');

h(2)=irf_panel('Elmn');
hold(h(2),'on');
irf_plot(h(2),E_LMN1.x,'Linewidth',lnwid);
irf_plot(h(2),E_LMN1.y,'Linewidth',lnwid);
irf_plot(h(2),E_LMN1.z,'Linewidth',lnwid);
hold(h(2),'off');
ylabel(h(2),{'E','(mV/m)'},'Interpreter','tex');
% set(h(2),'ylim',[-120 70],'ytick',[-100:50:50]);
irf_legend(h(2),{'Ex  ','Ey  ','Ez  '},[0.98 0.03],'fontsize',12)
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12);
grid(h(2),'off');

h(3)=irf_panel('j');
hold(h(3),'on');
% irf_plot(h(3),Ni1,'Linewidth',lnwid,'color','b');
irf_plot(h(3),[j(:,1),j_LMN(:,2)*1e9],'Linewidth',lnwid);
irf_plot(h(3),[j(:,1),j_LMN(:,3)*1e9],'Linewidth',lnwid);
irf_plot(h(3),[j(:,1),j_LMN(:,4)*1e9],'Linewidth',lnwid);
irf_plot(h(3),[B_lmn3(:,1) B_lmn3(:,2)*0],'k--', 'Linewidth',0.5);%%铏氱嚎
hold(h(3),'off');
% set(h(3),'ylim',[-1200 1200],'ytick',[-1000:1000:1000]);
ylabel(h(3),{'j','(nA/m^2)'},'Interpreter','tex');
% irf_legend(h(3),{'n_{i} '},[0.93 0.08],'color','b','fontsize',12)
% irf_legend(h(3),{'n_{e} '},[0.99 0.08],'color','r','fontsize',12)
irf_legend(h(3),{'Jx  ','Jy  ','Jz  '},[0.98 0.03],'fontsize',12)
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',12)
grid(h(3),'off');

h(4)=irf_panel('j_dot_E');
hold(h(4),'on');
irf_plot(h(4),[j(:,1),j_dot_E02(:,1)*1e6],'Linewidth',lnwid,'color','k');
irf_plot(h(4),[B_lmn3(:,1) B_lmn3(:,2)*0],'k--', 'Linewidth',0.5);%%铏氱嚎
hold(h(4),'off');
% set(h(4),'ylim',[-1200 1100],'ytick',[-1000:1000:1000]);
ylabel(h(4),{'j dot E','(nW/m^3)'},'Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)
grid(h(4),'off');

h(5)=irf_panel('temperature');
hold(h(5),'on');
% irf_plot(h(3),Ni1,'Linewidth',lnwid,'color','b');
irf_plot(h(5),[B2(:,1),matParTe1(:,2)],'Linewidth',lnwid,'color','k');
irf_plot(h(5),[B2(:,1),matPerTe1(:,2)],'Linewidth',lnwid,'color','b');
hold(h(5),'off');
% set(h(5),'ylim',[75 125],'ytick',[80:20:120]);
ylabel(h(5),{'T_{e}','(eV)'},'Interpreter','tex');
irf_legend(h(5),{'T_{e//} '},[0.93 0.08],'color','k','fontsize',12)
irf_legend(h(5),{'T_{e⊥} '},[0.99 0.08],'color','b','fontsize',12)
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12)
grid(h(5),'off');

h(6)=irf_panel('VeLMN');
hold(h(6),'on');
irf_plot(h(6),[VeS_LMN2(:,1),VeS_LMN2(:,2)],'Linewidth',lnwid);
irf_plot(h(6),[VeS_LMN2(:,1),VeS_LMN2(:,3)],'Linewidth',lnwid);
irf_plot(h(6),[VeS_LMN2(:,1),VeS_LMN2(:,4)],'Linewidth',lnwid);
irf_plot(h(6),[VeS2(:,1),VeS2(:,5)],'Linewidth',lnwid,'color','k');
irf_plot(h(6),[B_lmn3(:,1) B_lmn3(:,2)*0],'k--', 'Linewidth',0.5);%%铏氱嚎
hold(h(6),'off');
% ylabel(h(2),{'V_{etLMN}','(km/s)'},'Interpreter','tex');
ylabel(h(6),{'V_{e}','(km/s)'},'Interpreter','tex');
% set(h(6),'ylim',[-1100 1100],'ytick',[-1000:1000:1000]);
irf_legend(h(6),{'Vex  ','Vey  ','Vez  '},[0.98 0.03],'fontsize',12)
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)
grid(h(6),'off');

% h(5)=irf_panel('j dot E+VexB');
% hold(h(5),'on');
% irf_plot(h(5),[j(:,1),j_dot_E(:,1)*1e6],'Linewidth',lnwid,'color','k');
% irf_plot(h(5),[B_lmn1(:,1) B_lmn1(:,2)*0],'k--', 'Linewidth',0.5);%%铏氱嚎
% hold(h(5),'off');
% set(h(5),'ylim',[-20 20],'ytick',[-15:10:15]);
% ylabel(h(5),{'j dot(E+V_{e}xB)','(nW/m^3)'},'Interpreter','tex');
% irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(5),'off');

% h(6)=irf_panel('j dot E+VexB');
% hold(h(6),'on');
% irf_plot(h(6),[B3(:,1),j_dot_E2(:,1)*1e6],'Linewidth',lnwid,'color','k');%单个卫星数据算的电流
% irf_plot(h(6),[B_lmn3(:,1) B_lmn3(:,2)*0],'k--', 'Linewidth',0.5);
% hold(h(6),'off');
% set(h(6),'ylim',[-26 26],'ytick',[-20:20:20]);
% ylabel(h(6),{'j\cdot(E+V_{e}xB)','(nW/m^3)'},'Interpreter','tex');
% irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(6),'off');

% h(6)=irf_panel('d_cros_ve');
% mmscolors=[0 0 1; 0 1 0;1 0 0]; %钃濈豢绾?
% set(h(6),'ColorOrder',mmscolors)
% hold(h(6),'on');
% irf_plot(h(6),[d_cros_ve_LMN(:,1),d_cros_ve_LMN_abs(:,1)],'Linewidth',lnwid,'color','k');
% % irf_plot(h(7),[d_cros_ve_LMN(:,1),d_cros_ve_LMN(:,3)],'Linewidth',lnwid);
% % irf_plot(h(7),[d_cros_ve_LMN(:,1),d_cros_ve_LMN(:,4)],'Linewidth',lnwid);
% hold(h(6),'off');
% % ylabel(h(1),{'B_{tLMN}','(nT)'},'Interpreter','tex');
% ylabel(h(6),{'|鈻絰V_{e}|','(s^(-1))'},'Interpreter','tex');
% set(h(6),'ylim',[-20 310],'ytick',[0:100:300]);
% % irf_legend(h(7),{'(鈻絰V_{e})_{L}   ','(鈻絰V_{e})_{M}   ','(鈻絰V_{e})_{N}'},[0.98 0.20],'fontsize',12)
% irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(6),'off');

h(7)=irf_panel('d_cros_ve');
mmscolors=[0 0 1; 0 1 0;1 0 0]; 
set(h(7),'ColorOrder',mmscolors)
hold(h(7),'on');
irf_plot(h(7),[B3(:,1),d_cros_ve_LMN_abs(:,1)],'Linewidth',lnwid,'color','k');
% irf_plot(h(7),[d_cros_ve_LMN(:,1),d_cros_ve_LMN(:,4)],'Linewidth',lnwid);
hold(h(7),'off');
ylabel(h(7),{'|\nabla\timesV_{e}|','(s^-1)'},'Interpreter','tex');
% set(h(7),'ylim',[0 350],'ytick',[0:200:200]);
% irf_legend(h(7),{'(鈻絰V_{e})_{L}   ','(鈻絰V_{e})_{M}   ','(鈻絰V_{e})_{N}'},[0.98 0.20],'fontsize',12)
irf_legend(h(7),'(g)',[0.99 0.98],'color','k','fontsize',12)
grid(h(7),'off');

h(8)=irf_panel('VAM');
hold(h(8),'on');
% irf_plot(h(3),Ni1,'Linewidth',lnwid,'color','b');
irf_plot(h(8),[ne2(:,1),VAM(:,1)/2],'Linewidth',lnwid,'color','k');
irf_plot(h(8),[VeS_LMN2(:,1),VeS_LMN2(:,3)],'Linewidth',lnwid,'color','r');
hold(h(8),'off');
% set(h(9),'ylim',[0 15],'ytick',[5:5:10]);
ylabel(h(8),{'Velocity','(km/s)'},'Interpreter','tex');
irf_legend(h(8),{'|C_{AeM}|/2 '},[0.93 0.08],'color','k','fontsize',12)
irf_legend(h(8),{'V_{eM} '},[0.99 0.08],'color','r','fontsize',12)
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',12)
grid(h(8),'off');

h(9)=irf_panel('distance plot');
hold(h(9),'on');
irf_plot(h(9),[dRmin(:,1),dRmin(:,2)], 'color','b', 'Linewidth',lnwid); hold on;
for ii=1:2:length(dRmin(:,1))
     irf_plot(h(9), dRmin(ii,:), [typeSim(ii) clrSim(ii)], 'MarkerSize',mksizSim(ii),'MarkerFaceColor',faceclrSim(ii), 'Linewidth',lnwid); hold on;
end
% set(h(8), 'Ylim',[0 35],'Ytick',[0:10:30]);
set(h(9), 'Ylim',[0 24],'Ytick',[0:10:20]);
ylabel(h(9),{'|r|','(km)'},'Interpreter','tex');
irf_legend(h(9),'(i)',[0.99 0.25],'color','k','fontsize',12)
grid(h(9),'off');


h(10)=irf_panel('error plot');
hold(h(10),'on');
a1=irf_plot(h(10),[err_4C(:,1) err_4C(:,2)],'k','Linewidth',lnwid);%散度/旋度
fill(h(10), [a1.XData,fliplr(a1.XData)],[zeros(1,length(a1.YData)),fliplr(a1.YData)],[0.8706 0.9216 0.9804]);%浅灰色
irf_plot(h(10),[eigVal_err(:,1) eigVal_err(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;
% a1=irf_plot(h(9),[err_4C(:,1) err_4C(:,2)],'k','Linewidth',lnwid);
% a2=irf_plot(h(9), eigVal_err_v2, 'c','Linewidth',lnwid); 
% 
% % fill(h1(i_plot),[a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
% fill(h(9),[a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
% fill(h(9), [a1.XData,fliplr(a1.XData)],[zeros(1,length(a1.YData)),fliplr(a1.YData)],[0.7 0.7 0.7]);


% irf_plot(h1(i_plot), [err_4C(:,1) err_4C(:,2)*0+40], 'k--','Linewidth',lnwid); 
% irf_plot(h(9),[eigVal_err(:,1) eigVal_err(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;

grid(h(10),'off'); 
% set(h1(i_plot),'Ylim',[0 160]);
%ylabel('\fontcolor{k}|\nabla\cdot\bf{B}|/|\nabla\times\bf{B}|,\fontcolor{b}|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}| [%]');
ylabel(h(10),{'\alpha',' (%)'},'Interpreter','tex');
set(h(10),'Ylim',[0 130],'Ytick',[0 40 100]);
irf_legend(h(10),'(j)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h(10),{'|\nabla\cdot\bf{nVe}|/|\nabla\times\bf{nVe}|'},[0.85 0.9],'color','k','fontsize',10)


%% mark

Tnull='2017-01-27T12:05:31.712Z';
idxnull=find(gradVe(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
mark1=B1(idxnull,1);
for ii = 1:10
    c_eval("irf_pl_mark(h(?),mark1,'k','linewidth',lnwid)",ii);
end

tint3=irf.tint('2017-01-27T12:05:31.67Z/2017-01-27T12:05:31.74Z');
for ii = 1:10
    c_eval('irf_pl_mark(h(?),tint3,[0.7,0.7,0.7])',ii)
end

tint4=irf.tint('2017-01-27T12:05:30.945Z/2017-01-27T12:05:31.045Z');
for ii = 1:10
    c_eval('irf_pl_mark(h(?),tint4,[0.7,0.7,0.7])',ii)
end

Tnull='2017-01-27T12:05:30.99Z';
idxnull=find(gradVe(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
mark1=B1(idxnull,1);
for ii = 1:10
    c_eval("irf_pl_mark(h(?),mark1,'k','linewidth',lnwid)",ii);
end

% xSize = 400; ySize=600;
xSize = 700; ySize=1000;
set(gcf,'position',[200 50 xSize ySize])
set(gcf,'render','painters');
Paper_X = 25; Paper_Y = 17; 
coef=floor(max(xSize/Paper_X,ySize/Paper_Y));
FigSize_X = xSize/coef; FigSize_Y = ySize/coef;
xLeft2 = (Paper_X- FigSize_X)/2;  yTop2 = (Paper_Y- FigSize_Y)/2; 
set(gcf,'PaperSize', [Paper_X Paper_Y]); 
set(gcf,'PaperPosition',[xLeft2 yTop2 FigSize_X FigSize_Y])

irf_plot_ylabels_align(h)
irf_zoom(h,'x',tint);
set(gca,"XTickLabelRotation",0)
% set(gcf,'paperpositionmode','auto')
set(h,'fontsize',12);
% figname = ['vortex_observe_Vt'];
% print(gcf, '-dpdf','-r600',[figname '.pdf']);