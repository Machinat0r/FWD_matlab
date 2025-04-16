clear
clc

%build database
mms.db_init('local_file_db','Z:/Data/MMS/')
thresold=0.4;  
hh1=10;
%% load data
for ic=1:4
% tintStr = '06-July-2017';
% Tsta='2017-01-27T12:05:30.40Z'; 
% Tend='2017-01-27T12:05:32.41Z'; 
% tint=irf.tint('2017-01-27T12:05:30.40Z/2017-01-27T12:05:32.41Z');


Tsta = '2017-07-06T17:31:58.000Z';
Tend = '2017-07-06T17:31:59.500Z';
tint=irf.tint('2017-07-06T17:31:30.000Z/2017-07-06T17:32:30.000Z');
Time='2017-07-06T17:31:58.700Z';

time = irf_time(Time,'utc>epochtt');
c_eval('Vegse?=mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
% c_eval('Ve?=irf.ts2mat(Vegse?);',ic);
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
c_eval('Bt?=Bxyz?.abs;',ic)
c_eval('Exyz?=mms.get_data(''E_gse_edp_brst_l2'',tint,?);',ic);
c_eval('Et?=Exyz?.abs;',ic)
c_eval('Rgse?=mms.get_data(''R_gse'',tint,?);',ic);
c_eval('R?=irf.ts2mat(Rgse?);',ic);
c_eval('Ne?=mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ne?=irf.ts2mat(Ne?);',ic);
c_eval('e_r? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energy_brst'',tint);',ic);
c_eval('VeS?=irf.ts2mat(Vegse?);',ic);
end
%% smooth_step1
kk = length (e_r1.data);
if mod(kk,2) == 1
   le = kk-2; lo = kk-3;
else
   le = kk-3; lo = kk-2;
end
for ic =1:4
    c_eval('e_r = e_r?;',ic)
    if e_r.data(1)>e_r.data(2)
        for ii = 1:2:le
%             c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
        end
    else
        for ii = 2:2:lo
%             c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
        end
    end
end

VeS1=irf_resamp(VeS1,B1);
VeS2=irf_resamp(VeS2,VeS1);
VeS3=irf_resamp(VeS3,VeS1);
VeS4=irf_resamp(VeS4,VeS1);

for ic =1:4
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('Ve?(:,1:4)=VeS?(:,1:4);',ic);
end
%% coordinate conversion
L=[0.05 -0.21 0.98];
M=[-0.23 0.95 0.22];
N=[-0.97 -0.24 -0.01];
for ic =1:4
c_eval('B_LMN?=irf_newxyz(Bxyz?,L,M,N);',ic);
c_eval('Blmn?=irf.ts2mat(B_LMN?);',ic);
c_eval('E_LMN?=irf_newxyz(Exyz?,L,M,N);',ic);
c_eval('Elmn?=irf.ts2mat(E_LMN?);',ic);
c_eval('V_LMN?=irf_newxyz(Vegse?,L,M,N);',ic);
c_eval('Vlmn?=irf.ts2mat(V_LMN?);',ic);
end


%% FOTE method
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
Ve1=irf_resamp(Ve1,B1);
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
gradB=c_4_grad('R?','Ve?','grad');
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