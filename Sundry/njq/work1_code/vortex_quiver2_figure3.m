%% 找涡旋之零点类型 33.485时刻 平滑 （d）-（i）图数据全部做处理
%% VeS_LMN2:奇偶平滑后速度，也即要画的速度
clear
clc
mms.db_init('local_file_db','E:\fu\data2')
%% VORTEX_9图全版
%% load data
ic=1:4;
hh = 15;
hh1= 20;
thresold = 0.4;  
tintStr = '27-Jan-2017';
% tint=irf.tint('2017-01-27T12:05:42.50Z/2017-01-27T12:05:43.80Z');
tint =irf.tint('2017-01-27T12:05:30.40Z/2017-01-27T12:05:32.41Z');
tint2 =irf.tint('2017-01-27T12:05:30.00Z/2017-01-27T12:06:00.00Z');
Tsta='2017-01-27T12:05:30.40Z'; 
Tend='2017-01-27T12:05:32.41Z';

c_eval('e_r? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energy_brst'',tint);',ic);

%% load B Fields data
c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
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
%% load Ve Fields data
c_eval('Vexyz?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
L=[0.05 -0.21 0.98]; % Vi
M=[-0.23 0.95 0.22];
N=[-0.97 -0.24 -0.01];
c_eval('Ve_LMN?=irf_newxyz(Vexyz?,L,M,N);',ic);
c_eval('Ve?=irf.ts2mat(Vexyz?);',ic);
c_eval('Ve?=irf_resamp(Ve?,B?);',ic);
c_eval('vi?=irf_resamp(vi?,B?);',ic);
c_eval('VeS?=irf.ts2mat(Vexyz?);',ic);
c_eval('VeS_LMN?=irf.ts2mat(Ve_LMN?);',ic);
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
c_eval('Et?=Exyz?.abs;',ic);
c_eval('E?=irf.ts2mat(Exyz?);',ic);
% c_eval('Vit?=irf.ts2mat(Vi?);',ic);
c_eval('E?=irf_resamp(E?,B?);',ic);

c_eval('Rgse?=mms.get_data(''R_gse'',tint2,?);',ic);
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
    [N,D] = eig(deltB_null);
    
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

L=[0.05 -0.21 0.98]; % Vi
M=[-0.23 0.95 0.22];
N=[-0.97 -0.24 -0.01];
ic=1:4;
c_eval('B_LMN?=irf_newxyz(Bxyz?,L,M,N);',ic);
c_eval('B_lmn?=irf.ts2mat(B_LMN?);',ic);
c_eval('E_LMN?=irf_newxyz(Exyz?,L,M,N);',ic);
c_eval('Vi_LMN?=irf_newxyz(Vi?,L,M,N);',ic);

[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
d_cros_ve=c_4_grad('R?','VeS?(:,1:4)','curl');
% d_cros_ve=c_4_grad('R?','VeSS?(:,1:4)','curl');
d_cros_ve_LMN=irf_newxyz(d_cros_ve,L,M,N);
n=length(j(:,1))-1; 
j_total=zeros(n,1);
j_E=zeros(n,1);
j_LMN=irf_newxyz(j(:,2:4),L,M,N);
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

%% plot
% h=irf_plot(1,'newfigure');
% xSize=300; ySize=300;
% set(gcf,'Position',[100 100 xSize ySize]);
% % mmscolors=[0 0 1; 0 1 0; 1 0 0; 0 0 0]; %蓝绿红黑
% mmscolors=[0 0 1;1 0 0];
% color1=[0 0 1];
% color2=[1 0 0];
% set(h,'ColorOrder',mmscolors)
% lnwid = 1;
% 横纵坐标有问题，都只是xy，不是UV长度  图例单独设置一个矢量，后期挪位置
% h(1)=irf_panel('quiver');
% hold(h(1),'on');
figure;
x = 0:38:950; %%vortex2
x=x';
y = VeS_LMN2(156:181,1)*0;
meanL=mean(VeS_LMN2(156:181,2));
meanM=mean(VeS_LMN2(156:181,3));
meanN=mean(VeS_LMN2(156:181,4));
L_raw = VeS_LMN2(156:181,2);
M_raw = VeS_LMN2(156:181,3);
N_raw = VeS_LMN2(156:181,4);

% x = 0:31:1000; %%vortex1
% x=x';
% y = VeS_LMN2(64:96,1)*0;
% meanL=mean(VeS_LMN2(64:96,2));
% meanM=mean(VeS_LMN2(64:96,3));
% meanN=mean(VeS_LMN2(64:96,4));
% L_raw = VeS_LMN2(64:96,2);
% M_raw = VeS_LMN2(64:96,3);
% N_raw = VeS_LMN2(64:96,4);

count=length(L_raw);
L=zeros(count,1);
M=zeros(count,1);
N=zeros(count,1);
for i=1:count
    L(i,1)=L_raw(i,1)-meanL;
    M(i,1)=M_raw(i,1)-meanM;
    N(i,1)=N_raw(i,1)-meanN;
end

% quiver(x, y, M, N)
quiver(x,y, L, N,0)
hold;
x1=[700,700];
y1=[-400,-400];
standard_x=[100,0];
standard_y=[0,100];
quiver(x1,y1, standard_x, standard_y,0)

xlim([0 1000]);
ylim([-550 320]);

xlabel('L')
ylabel('N')
title('vortex 2')
% legend('Vectors', 'Direction', 'Location', 'best')
% grid(h(1),'off');

% x = 0:0.1:3.2;
% y = VeS_LMN2(64:96,1)*0;
% [U, V] = meshgrid(VeS_LMN2(64:96,2), VeS_LMN2(64:96,4));
% 
% quiver(x, y, U, V)
% plot(0, 0, 'r', 'LineWidth', 2) % 添加代表矢量方向的红色线
% 
% hold(h(1), 'off');
% xlabel('X')
% ylabel('Y')
% title('Vector Field')
% legend('Vectors', 'Direction', 'Location', 'best')
% grid(h(1), 'off');

% xSize = 30; ySize=30;
% set(gcf,'position',[200 50 xSize ySize])
% set(gcf,'render','painters');
% Paper_X = 15; Paper_Y = 20; 
% coef=floor(max(xSize/Paper_X,ySize/Paper_Y));
% FigSize_X = xSize/coef; FigSize_Y = ySize/coef;
% xLeft2 = (Paper_X- FigSize_X)/2;  yTop2 = (Paper_Y- FigSize_Y)/2; 
% set(gcf,'PaperSize', [Paper_X Paper_Y]); 
% set(gcf,'PaperPosition',[xLeft2 yTop2 FigSize_X FigSize_Y])
% set(gca,'unit','centimeters','position',[xLeft2+2 yTop2-2 FigSize_X FigSize_Y])
% irf_plot_ylabels_align(h)
% % irf_zoom(h,'x',tint);
% % set(gcf,'paperpositionmode','auto')
% set(h,'fontsize',12);
% figname = ['quiver_try'];
% print(gcf, '-dpdf','-r600',[figname '.pdf']);