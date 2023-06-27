clear;
clc;
close all

global ParentDir 
ParentDir = '/Users/fwd/Documents/MATLAB/MMS/'; 
TempDir = '/Users/fwd/Documents/MATLAB/MMS/temp/';mkdir(TempDir);


% Tint=irf.tint('2018-08-22T15:34:32.00Z/2018-08-22T15:34:36.00Z');
% Tsta='2018-06-02T18:28:50Z';   
% Tend='2018-06-02T18:31:50Z';
% Tsta='2017-08-23T15:38:30.00Z';   
% Tend='2017-08-23T15:39:13.00Z';

% Tsta='2020-08-03T02:36:04.00Z';
% Tend='2020-08-03T02:36:30.00Z';   
% Tsta = '2017-08-20T02:02:00.00Z';
% Tend = '2017-08-20T02:03:00.00Z';
% Tsta = '2020-07-05T00:31:00.00Z';
% Tend = '2020-07-05T00:32:00.00Z';
% Tsta='2019-08-05T16:24:00.00Z';   
% Tend='2019-08-05T16:25:00.00Z';
Tsta='2019-08-05T16:24:00.00Z';   
Tend='2019-08-05T16:25:00.00Z';
% Tsta='2021-08-22T06:40:45.000Z';   
% Tend='2021-08-22T06:41:45.000Z';
% Tsta='2021-07-22T12:44:30.00Z';   
% Tend='2021-07-22T12:45:30.00Z';
TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
Tint=irf.tint(Tsta,Tend);
% Tint=irf.tint('2019-07-22T17:09:45.00Z/2019-07-22T17:11:00.00Z');
tint = Tint;
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 2;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filename_scpot = SDCFilenames(Date,iic,'inst','edp','drm','brst','dpt','scpot');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3,filenames4,filename_scpot];

expr1 = '_\d+\_v';
NameTags = regexp(filenames,expr1,'match');
NameTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NameTags)),'UniformOutput',false);

TTlist = strjoin(regexp(TT,'\d+','match'),'');
i = 1;flag = 0;%若flag=0，说明整段时间都在第i-1个Tag里
while str2double(TTlist(17:30)) > NameTags{i}  % 如果时间段刚好仅在某天的最后一个文件里会出bug，可以把时间往前调1ms
    if str2double(TTlist(1:14)) < NameTags{i}
        flag=1; break  %若flag=1，说明时间段的开始在第i-1个Tag里，结束在第i个里
    else
        i=i+1;
    end
end

% '/' for MacOs, '\' for Windows
if flag == 0
    tempTag = num2str(NameTags{i-1});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag(1:4),'/',tempTag(5:6),'/',...
            tempTag(7:8),'/',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
    desmoms1 = desmoms; desmoms2 = desmoms1;
else
    if i == 1
        errordlg('时间起始处无brst数据，请检查时间范围,或使用Overview_srvydownload程序')
    end
    tempTag1 = num2str(NameTags{i-1});
    tempTag2 = num2str(NameTags{i});
    filenames1 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag1)));
    filenames2 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag2)));
    filenames = [filenames1,filenames2];
    desmoms1 = [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag1(1:4),'/',tempTag1(5:6),'/',...
            tempTag1(7:8),'/',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag2(1:4),'/',tempTag2(5:6),'/',...
            tempTag2(7:8),'/',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
end
SDCFilesDownload(filenames,TempDir)
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);
%%
% load B
units = irf_units;
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],iic);
c_eval(['Bt?_ts=B?_ts.abs;'],iic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],iic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],iic);
% lvbo        
c_eval('dfB? = 1/median(diff(B?_ts.time.epochUnix));',iic);
c_eval('Bbf? = B?_ts.filt(0.8,1.1,dfB?,3);',iic);
c_eval(['Bbf?=irf.ts2mat(Bbf?);'],iic);

% c_eval('Bbff? = B?_ts.filt(0,0.8,dfB?,3);',ic);
% c_eval(['Bbff?=irf.ts2mat(Bbff?);'],ic);

%         c_eval('Blmn?=irf_newxyz(Bbf1,L,M,N);',ic);

% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_fast_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

c_eval(['Bt?_res=irf_resamp(Bt?,Et?);'],ic);

c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);

c_eval(['Vexb?=irf_cross(E?,B?);'],ic);
c_eval(['Vexb?(:,2:4)=1e3*Vexb?(:,2:4)./[Bt?_res(:,2).^2 Bt?_res(:,2).^2 Bt?_res(:,2).^2];'],ic);%km/s


c_eval('dfE? =1/median(diff(E?_gsm.time.epochUnix));',ic);
c_eval('Ebf? = E?_gsm.filt(0.8,1.1,dfE?,3);',ic);
c_eval(['Ebf?=irf.ts2mat(Ebf?);'],ic);



% c_eval('Ebff? = E?_gsm.filt(0,0.8,dfE?,3);',ic);
% c_eval(['Ebff?=irf.ts2mat(Ebf?);'],ic);

%         c_eval('Elmn?=irf_newxyz(Ebf1,L,M,N);',ic);


% load FPI
c_eval('Ne?_ts = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('Ni?_ts = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
% c_eval('Ne?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
c_eval(['Ne?=irf.ts2mat(Ne?_ts);'],ic);
% c_eval('Ni?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
c_eval(['Ni?=irf.ts2mat(Ni?_ts);'],ic);
% % % c_eval('dfNe? = 1/median(diff(Ne?_ts.time.epochUnix));',ic);
% % % c_eval('Nebf? = Ne?_ts.filt(0,1,dfNe?,5);',ic);
% % % c_eval(['Nebf?=irf.ts2mat(Nebf?);'],ic);

% % % c_eval('dfNi? = 1/median(diff(Ni?_ts.time.epochUnix));',ic);
% % % c_eval('Nibf? = Ni?_ts.filt(0,1,dfNi?,5);',ic);
% % % c_eval(['Nibf?=irf.ts2mat(Nibf?);'],ic);




c_eval('Te_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_temppara_brst'',tint);',ic);
c_eval(['Te_para?=irf.ts2mat(Te_para?_ts);'],ic);
c_eval('Te_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempperp_brst'',tint);',ic);
c_eval(['Te_perp?=irf.ts2mat(Te_perp?_ts);'],ic);
c_eval(['Te?=[Te_para?(:,1),(Te_para?(:,2)+2*Te_perp?(:,2))/3.0];'],ic);

% c_eval('dfTe_para? = 1/median(diff(Te_para?_ts.time.epochUnix));',ic);
% c_eval('Te_parabf? = Te_para?_ts.filt(0,1.5,dfE?,5);',ic);
% c_eval(['Te_parabf?=irf.ts2mat(Te_parabf?);'],ic);

c_eval('Ti_para?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_temppara_brst'',tint);',ic);
c_eval(['Ti_para?=irf.ts2mat(Ti_para?_ts);'],ic);
c_eval('Ti_perp?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_tempperp_brst'',tint);',ic);
c_eval(['Ti_perp?=irf.ts2mat(Ti_perp?_ts);'],ic);
c_eval(['Ti?=[Ti_para?(:,1),(Ti_para?(:,2)+2*Ti_perp?(:,2))/3.0];'],ic);

c_eval('Ve?_ts = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
% c_eval('Ve?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
c_eval(['Vet?_ts=Ve?_ts.abs;'],ic); 
c_eval(['Ve?=irf.ts2mat(Ve?_ts);'],ic);
c_eval(['gsmVe?_ts=irf_gse2gsm(Ve?_ts);'],ic);
c_eval(['gsmVe?=irf.ts2mat(gsmVe?_ts);'],ic);
c_eval(['Vet?=irf.ts2mat(Vet?_ts);'],ic);

c_eval('dfVe? = 1/median(diff(gsmVe?_ts.time.epochUnix));',ic);
c_eval('Vebf? = gsmVe?_ts.filt(0,1,dfVe?,5);',ic);
c_eval(['Vebf?=irf.ts2mat(Vebf?);'],ic);



c_eval('Vi?_ts = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); 
% c_eval('Vi?_ts=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
c_eval(['Vit?_ts=Vi?_ts.abs;'],ic); 
c_eval(['Vi?=irf.ts2mat(Vi?_ts);'],ic);
c_eval(['gsmVi?_ts=irf_gse2gsm(Vi?_ts);'],ic);
c_eval(['gsmVi?=irf.ts2mat(gsmVi?_ts);'],ic);
c_eval(['Vit?=irf.ts2mat(Vit?_ts);'],ic);

%Vifac
c_eval(['Bt?_resVi=irf_resamp(Bt?,Vi?);'],ic);
c_eval(['Vifac?=irf_convert_fac(Vi?,B?,[1,0,0]);'],ic);

%VxB
c_eval('gsmVi?_res = irf_resamp(gsmVi?,B?);',ic);
c_eval('Evxb? = 1e3*cross(1000*gsmVi?_res(:,2:4),1e-9*B?(:,2:4));',ic); %mV/m
c_eval('Evxb? = Evxb? + E?_resamp(:,2:4);',ic);


% merge data/time from 2 cdf files
c_eval('energy_low?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_lowen_brst'',tint);',ic)
c_eval('energy_mid?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_miden_brst'',tint);',ic)
c_eval('energy_high?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_highen_brst'',tint);',ic)
c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energyspectr_omni_brst'',tint);',ic)
c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_energyspectr_omni_brst'',tint);',ic)
%% pressure & entropy
Pos = mms.get_data('R_gsm',tint);
c_eval('R? = Pos.gsmR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
c_eval('R? = irf_resamp(R?,Bt?);')
mu0 = units.mu0;
kB = units.kB;
 c_eval('Ni?=irf_resamp(Ni?,Bt?);',ic)
 c_eval('Ne?=irf_resamp(Ne?,Bt?);',ic)
 c_eval('Te?=irf_resamp(Te?,Bt?);',ic)
 c_eval('Ti?=irf_resamp(Ti?,Bt?);',ic)
c_eval('Pb?=[Bt?(:,1) ((Bt?(:,2).^2))/(2*mu0)*1e-9];',ic);%nPa
c_eval('Pti? = irf_multiply(11604.505*kB*1e6*1e9,[Ni?(:,1) Ni?(:,2)],1,[Ti?(:,1) Ti?(:,2)],1);',ic);
c_eval('Pte? = irf_multiply(11604.505*kB*1e6*1e9,[Ne?(:,1) Ne?(:,2)],1,[Te?(:,1) Te?(:,2)],1);',ic);
c_eval('Pt? = [Pti?(:,1) Pti?(:,2)+Pte?(:,2)];',ic);
miu0=400*pi;
c_eval('R?(:,2:4) = R?(:,2:4)./units.RE*1000;')
c_eval('VE?=10^0.7368.*sqrt(R?(:,2).^2+R?(:,3).^2).^0.7634.*abs(B?(:,4)).^(-0.3059)./sqrt(B?(:,4).^2+2.*miu0.*Pt?(:,2));',ic);
c_eval('S?=Pt?(:,2).*VE?.^(5/3);',ic);

c_eval('Br? = sqrt(B?(:,2).^2+B?(:,3).^2);')
c_eval('S?(B?(:,4) < Br? | B?(:,4) < 0) = nan;')

c_eval('beta? = [Bt?(:,1) Pt?(:,2)./Pb?(:,2)];',ic);
c_eval('pe_pb? = [Bt?(:,1) Pte?(:,2)./Pb?(:,2)];',ic);
%%
%J
c_eval('Je?_ts = -units.e*Ne?_ts*gsmVe?_ts*1e3*1e6*1e9;',ic);
c_eval('Ji?_ts = units.e*Ne?_ts*gsmVi?_ts.resample(Ne?_ts.time)*1e3*1e6*1e9;',ic);
c_eval('J?_ts = (Je?_ts+Ji?_ts);',ic);
c_eval('J? = irf.ts2mat(J?_ts);',ic);

%JdotE
c_eval('E_resJ = irf_resamp(E?,J?);',ic);
c_eval('JdotEtemp = [E_resJ(:,1) irf_dot(J?(:,2:4),E_resJ(:,2:4))];',ic);
c_eval('JdotE? = irf_resamp(JdotEtemp,Ni?);',ic);


energy_mid1=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_pitchangdist_miden_brst',tint);
energy_high1=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_pitchangdist_highen_brst',tint);
energy_e1=mms.db_get_variable('mms1_fpi_brst_l2_des-moms','mms1_des_energyspectr_omni_brst',tint);
energy_i1=mms.db_get_variable('mms1_fpi_brst_l2_dis-moms','mms1_dis_energyspectr_omni_brst',tint);

energy_low2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_pitchangdist_lowen_brst',tint);
energy_mid2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_pitchangdist_miden_brst',tint);
energy_high2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_pitchangdist_highen_brst',tint);
energy_e2=mms.db_get_variable('mms2_fpi_brst_l2_des-moms','mms2_des_energyspectr_omni_brst',tint);
energy_i2=mms.db_get_variable('mms2_fpi_brst_l2_dis-moms','mms2_dis_energyspectr_omni_brst',tint);

energy_low3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_pitchangdist_lowen_brst',tint);
energy_mid3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_pitchangdist_miden_brst',tint);
energy_high3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_pitchangdist_highen_brst',tint);
energy_e3=mms.db_get_variable('mms3_fpi_brst_l2_des-moms','mms3_des_energyspectr_omni_brst',tint);
energy_i3=mms.db_get_variable('mms3_fpi_brst_l2_dis-moms','mms3_dis_energyspectr_omni_brst',tint);

energy_low4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_pitchangdist_lowen_brst',tint);
energy_mid4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_pitchangdist_miden_brst',tint);
energy_high4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_pitchangdist_highen_brst',tint);
energy_e4=mms.db_get_variable('mms4_fpi_brst_l2_des-moms','mms4_des_energyspectr_omni_brst',tint);
energy_i4=mms.db_get_variable('mms4_fpi_brst_l2_dis-moms','mms4_dis_energyspectr_omni_brst',tint);


% c_eval('dfE? = 1/median(diff(Exyz?.time.epochUnix));',ic);
% c_eval('dfB? = 1/median(diff(Bscm?.time.epochUnix));',ic);
% c_eval('Exyzfachf? = Exyzfac?.filt(9,12,dfE?,5);',ic);
if flag == 1
c_eval(['fpiFilee1 = dataobj(','''',desmoms1,'''',');'],ic);
c_eval('energy_low1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_lowen_brst'');',ic);
c_eval('energy_mid1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_miden_brst'');',ic);
c_eval('energy_high1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_highen_brst'');',ic);
c_eval('energy_e1 = get_variable(fpiFilee1,''mms?_des_energyspectr_omni_brst'');',ic);

c_eval(['fpiFilee2 = dataobj(','''',desmoms2,'''',');'],ic);
c_eval('energy_low2 = get_variable(fpiFilee2,''mms?_des_pitchangdist_lowen_brst'');',ic);
c_eval('energy_mid2 = get_variable(fpiFilee2,''mms?_des_pitchangdist_miden_brst'');',ic);
c_eval('energy_high2 = get_variable(fpiFilee2,''mms?_des_pitchangdist_highen_brst'');',ic);
c_eval('energy_e2 = get_variable(fpiFilee2,''mms?_des_energyspectr_omni_brst'');',ic);
% data merge
data1=energy_low1.data; data2=energy_low2.data; data=[data1;data2];energy_low=energy_low1; energy_low.data=data; energy_low.nrec=energy_low1.nrec+energy_low2.nrec;
data1=energy_mid1.data; data2=energy_mid2.data; data=[data1;data2];energy_mid=energy_mid1; energy_mid.data=data;  energy_mid.nrec=energy_mid1.nrec+energy_mid2.nrec;
data1=energy_high1.data; data2=energy_high2.data; data=[data1;data2];energy_high=energy_high1; energy_high.data=data;  energy_high.nrec=energy_high1.nrec+energy_high2.nrec;
data1=energy_e1.data; data2=energy_e2.data; data=[data1;data2];energy_e=energy_e1; energy_e.data=data;  energy_e.nrec=energy_e1.nrec+energy_e2.nrec;
% time merge
data1=energy_low1.DEPEND_0.data;data2=energy_low2.DEPEND_0.data; data=[data1;data2]; energy_low.DEPEND_0.data=data;
data1=energy_mid1.DEPEND_0.data;data2=energy_mid2.DEPEND_0.data; data=[data1;data2]; energy_mid.DEPEND_0.data=data;
data1=energy_high1.DEPEND_0.data;data2=energy_high2.DEPEND_0.data; data=[data1;data2]; energy_high.DEPEND_0.data=data;
data1=energy_e1.DEPEND_0.data;data2=energy_e2.DEPEND_0.data; data=[data1;data2]; energy_e.DEPEND_0.data=data;
else
c_eval(['fpiFilee1 = dataobj(','''',desmoms1,'''',');'],ic);
c_eval('energy_low1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_lowen_brst'');',ic);
c_eval('energy_mid1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_miden_brst'');',ic);
c_eval('energy_high1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_highen_brst'');',ic);
c_eval('energy_e1 = get_variable(fpiFilee1,''mms?_des_energyspectr_omni_brst'');',ic);

energy_low.DEPEND_0.data=energy_low1.DEPEND_0.data;
energy_mid.DEPEND_0.data=energy_mid1.DEPEND_0.data;
energy_high.DEPEND_0.data=energy_high1.DEPEND_0.data;
energy_low=energy_low1; energy_low.data=energy_low1.data; energy_low.nrec=energy_low1.nrec;
energy_mid=energy_mid1; energy_mid.data=energy_mid1.data;  energy_mid.nrec=energy_mid1.nrec;
energy_high=energy_high1; energy_high.data=energy_high1.data;  energy_high.nrec=energy_high1.nrec;
energy_e=energy_e1; energy_e.data=energy_e1.data;  energy_e.nrec=energy_e1.nrec;
end

% Ne_res = irf_resamp(Ne1,Bt1);
% Va = 1e-9*Bt1(:,2)./sqrt(units.me*units.mu0*1e6*Ne_res(:,2));

%% div and PI
Pos = mms.get_data('R_gsm',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gsmR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
PI(abs(PI)<0.5) = 0;

monopole_index=zeros(length(B1(:,1)),1);
gradB=c_4_grad('R?','B?','grad');
eigVal_err_v2=B1(:,1);
for ii=1:length(B1(:,1))  
    deltB_null=reshape(gradB(ii,2:end),3,3);
    [V,D] = eig(deltB_null);
    % Figure 1o    以最大特征值归一化 |Δ·B|/|ΔxB|
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]); 
    if ~isnan(deltB_null)
    [V,D] = eig(deltB_null);
    % Figure 1o    ξ
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]);  

    if PI(ii)~=0 && isreal(D(1,1)) && isreal(D(2,2)) && isreal(D(3,3)) && D(1,1) > 0 && D(2,2) > 0 && D(3,3) > 0
        monopole_index(ii) = 1;
    elseif PI(ii)~=0 && isreal(D(1,1)) && isreal(D(2,2)) && isreal(D(3,3)) && D(1,1) < 0 && D(2,2) < 0 && D(3,3) < 0
        monopole_index(ii) = -1;
    elseif PI(ii)~=0 && isreal(D(1,1))+isreal(D(2,2))+isreal(D(3,3))==1 && real(D(1,1)) < 0 && real(D(2,2)) < 0 && real(D(3,3)) < 0
        monopole_index(ii) = -0.5;
    elseif PI(ii)~=0 && isreal(D(1,1))+isreal(D(2,2))+isreal(D(3,3))==1 && real(D(1,1)) > 0 && real(D(2,2)) > 0 && real(D(3,3)) > 0
        monopole_index(ii) = 0.5;
    end
    else
    eig  Val_err_v2(ii,2)=nan;
    end
end
%% J dot E'
[J_B,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
J_B(:,2:4) = 1e9*J_B(:,2:4);
c_eval('E_resJ? = irf_resamp(E?,J_B);',ic);
c_eval('JdotE_B? = [E_resJ?(:,1) irf_dot(J_B(:,2:4),E_resJ?(:,2:4))];',ic);

% c_eval('gsmVe? = irf_resamp(gsmVe?,B?);',ic)
% c_eval('Eplus? = [E_resJ?(:,1) E_resJ?(:,2:4) + 1e3*irf_cross(1e3*gsmVe?(:,2:4),1e-9*B?(:,2:4))];',ic)
% c_eval('JdotEplus? = [Eplus?(:,1) irf_dot(J_B(:,2:4),Eplus?(:,2:4))];',ic);
c_eval('B?_resE = irf_resamp(B?,E?);',ic)
c_eval('gsmVe? = irf_resamp(gsmVe?,E?);',ic)
c_eval('J_B = irf_resamp(J_B,E?);',ic)
c_eval('Eplus? = [E?(:,1) E?(:,2:4) + 1e3*irf_cross(1e3*gsmVe?(:,2:4),1e-9*B?_resE(:,2:4))];',ic)
% c_eval('JdotEplus? = [Eplus?(:,1) irf_dot(J_B(:,2:4),Eplus?(:,2:4))];',ic);
c_eval('J? = irf_resamp(J?,E?);',ic);
c_eval('JdotEplus? = [Eplus?(:,1) irf_dot(J?(:,2:4),Eplus?(:,2:4))];',ic);
c_eval('JdotEint = [JdotEplus?(:,1) cumsum(JdotEplus?(:,2))];',ic);
%% Load Data 
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',Tint);',ic);
magB = Bxyz.abs;
B=irf.ts2mat(Bxyz);
Bt=irf.ts2mat(magB);
c_eval('Exyz_gse=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
% E_temp=irf.ts2mat(Exyz_gse);
Exyz=irf_gse2gsm(Exyz_gse);
E=irf.ts2mat(Exyz);
% Exyz = TSeries(Exyz_gse.time,[Exyz_gsm(:,2:4)]);
% Exyz = irf.ts_vec_xyz(Exyz_gse.time,Exyz_gsm(:,2:4));
c_eval('Bscm_ts=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',Tint);',ic);

% c_eval('Epar?_err_ts=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);',ic);
% c_eval('Epar?_err=irf.ts2mat(Epar?_err_ts);',ic);

% Bscm1 = irf_gse2gsm(Bscm_ts);
% Bscm1 = irf.ts2mat(Bscm_ts);

% for ii=1:length(Bscm_cell)
%     c_eval('Bscm?=Bscm_cell{?};',ii); 
%     
%     c_eval('time?=Bscm?.time;',ii); 
% end
% t_Bscm=time1;
% data_Bscm=Bscm1.data;
% for ii=2:length(Bscm_cell)
%     
%     c_eval('t_Bscm.epoch((length(t_Bscm.epoch(:,1))+1):(length(t_Bscm.epoch(:,1))+length(time?.epoch(:,1))),1)=time?.epoch;',ii); 
%     c_eval('data_Bscm((length(data_Bscm(:,1))+1):(length(data_Bscm(:,1))+length(Bscm?.data(:,1))),1:end)=Bscm?.data;',ii); 
% end
% Bscm_gse=irf.ts_vec_xyz(t_Bscm.epoch,data_Bscm);
flag = 1;
try
%     Bscm1=irf_gse2gsm(Bxyz);
    Bscm1=irf_gse2gsm(Bscm_ts);
catch
    Bscm1=irf_gse2gsm(Bscm_ts{1,1});
    Bscm2=irf_gse2gsm(Bscm_ts{1,2});
    if length(Bscm1)~=length(Bscm2)
    flag = 2;
    end
end
% Bscm2=irf_gse2gsm(Bscm_ts{1,2});
% Bscm{1,2}=irf_gse2gsm(Bscm_cell{1,2});

% 
% Bscm_mat1=irf.ts2mat(Bscm_gse1);
% Bscm_mat2=irf.ts2mat(Bscm_gse2);
% Bscm_mat=[Bscm_mat1;Bscm_mat2];
%Bscm=irf_gse2gsm(Bscm_gse);
% Bscm = irf.ts_vec_xyz(irf_time(Bscm_mat(:,1),'epoch>epochtt'),Bscm_mat_gsm(:,2:4));
% Bscm=Bscm{1};            %Bscm是cell
% c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);

% L=-[-0.562 0.719 -0.410];
% M=-[-0.555 -0.695 -0.458];
% N=[0.614 0.030 -0.789];

% B_LMN=irf_newxyz(Bxyz,L,M,N);
% Bxyzmag = TSeries(Bxyz.time,[B_LMN.data magB.data]);
% E_LMN=irf_newxyz(Exyz,L,M,N);
% E_LMN=TSeries(Exyz_gse.time,E_LMN(:,2:4));
%% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac1 = irf_convert_fac(Bscm1,Bxyz,[1 0 0]);
lf=2;
if flag == 2
    Bscmfac2 = irf_convert_fac(Bscm2,Bxyz,[1 0 0]);
    dfB2 = 1/median(diff(Bscm2.time.epochUnix));
    Bscmfachf2 = Bscmfac2.filt(lf,0,dfB2,5);
%     Bscmfachf2 = Bscmfac2.filt(40,300,dfB2,5);
    Bfachf2=irf.ts2mat(Bscmfachf2);
end
%% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB1 = 1/median(diff(Bscm1.time.epochUnix));
% dfE = 128;


% 
Exyzfachf = Exyzfac.filt(lf,0,dfE,5);
Exyzfaclf = Exyzfac.filt(0,lf,dfE,5);
Bscmfachf1 = Bscmfac1.filt(lf,0,dfB1,5);
Bscmfachf1=irf_convert_fac(Bscmfachf1,B,[1,0,0]);
% Bscmfachf1 = Bscmfac1.filt(0.02,0,128,5);
% Bscmfachf1 = Bscmfac1.filt(40,300,dfB1,5);
% % % 
% % % hf = 40;
% % % Exyzfachf = Exyzfac.filt(lf,hf,dfE,5);
% % % Bscmfachf1 = Bscmfac1.filt(lf,hf,dfB1,5);


Efachf=irf.ts2mat(Exyzfachf);
Bfachf1=irf.ts2mat(Bscmfachf1);
%% Wavelet transforms
nf = 100;
% Ewavelet = irf_wavelet(Exyzfachf,'nf',nf,'f',[lf 1]);
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 50000]);
Bwavelet1 = irf_wavelet(Bscmfachf1,'nf',nf,'f',[lf 6000]);
% Bwavelet1 = irf_wavelet(Bscmfachf1,'nf',nf,'f',[lf 1]);
if flag == 2
Bwavelet2 = irf_wavelet(Bscmfachf2,'nf',nf,'f',[lf 6000]);
end

%compress wavelet transform data 10 point average
nc = 20;
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)]
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specperpE=struct('t',Ewavelettimes);
specperpE.f=Ewavelet.f;
specperpE.p=Ewaveletx+Ewavelety;
specperpE.f_label='';
specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};

specparE=struct('t',Ewavelettimes);
specparE.f=Ewavelet.f;
specparE.p=Ewaveletz;
specparE.f_label='';
specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};


idx1 = [nc/2:nc:length(Bwavelet1.t)-nc/2];
Bwavelettimes1 = Bwavelet1.t(idx1);
Bwaveletx1 = zeros(length(idx1),nf);
Bwavelety1 = zeros(length(idx1),nf);
Bwaveletz1 = zeros(length(idx1),nf);
for ii = [1:length(idx1)]
        Bwaveletx1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,1}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
        Bwavelety1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,2}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
        Bwaveletz1(ii,:) = squeeze(irf.nanmean(Bwavelet1.p{1,3}([idx1(ii)-nc/2+1:idx1(ii)+nc/2-1],:),1));
end
specB1=struct('t',Bwavelettimes1);
specB1.f=Bwavelet1.f;
specB1.p=Bwaveletx1+Bwavelety1+Bwaveletz1;
specB1.f_label='';
specB1.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};

if flag ==2
idx2 = [nc/2:nc:length(Bwavelet2.t)-nc/2];
Bwavelettimes2 = Bwavelet2.t(idx2);
Bwaveletx2 = zeros(length(idx2),nf);
Bwavelety2 = zeros(length(idx2),nf);
Bwaveletz2 = zeros(length(idx2),nf);
for ii = [1:length(idx2)]
        Bwaveletx2(ii,:) = squeeze(irf.nanmean(Bwavelet2.p{1,1}([idx2(ii)-nc/2+1:idx2(ii)+nc/2-1],:),1));
        Bwavelety2(ii,:) = squeeze(irf.nanmean(Bwavelet2.p{1,2}([idx2(ii)-nc/2+1:idx2(ii)+nc/2-1],:),1));
        Bwaveletz2(ii,:) = squeeze(irf.nanmean(Bwavelet2.p{1,3}([idx2(ii)-nc/2+1:idx2(ii)+nc/2-1],:),1));
end
specB2=struct('t',Bwavelettimes2);
specB2.f=Bwavelet2.f;
specB2.p=Bwaveletx2+Bwavelety2+Bwaveletz2;
specB2.f_label='';
specB2.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};
end
%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;
B_SI=magB.data*1e-9;
Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fce01=Fce*0.1;
Fce05=Fce*0.5;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
Fcp = irf.ts_scalar(magB.time,Fcp);
Fpe = irf.ts_scalar(magB.time,Fpe);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);
Fce01=irf.ts_scalar(magB.time,Fce01);
Fce05=irf.ts_scalar(magB.time,Fce05);
%% Init figure 2
n_subplots=8;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize = 20; coef=floor(min(600/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([B(:,1) B(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B(:,1) B(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B(:,1) B(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bt(:,1) Bt(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B(:,1) B(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('B [nT]','fontsize',12);
set(h(1),'Ylim',[fix(min([min(B(:,2)) min(B(:,3)) min(B(:,4))])/10)*10-10 fix(max(Bt(:,2))/10)*10+10]);
% set(gca,'Ylim',[-10 25]);
% pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.05 0.06]);
%% E plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([E(:,1) E(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([E(:,1) E(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([E(:,1) E(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Et(:,1) Et(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([E(:,1) E(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('E [mV/m]','fontsize',12);
% set(h(2),'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
set(gca,'Ylim',[-10 10]);
% pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'E_x','E_y','E_z'},[0.05 0.06]);

%% Efachf plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Efachf(:,1) Efachf(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Efachf(:,1) Efachf(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Efachf(:,1) Efachf(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Efachf(:,1) Efachf(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('\deltaE [mV m^{-1}]');
% set(gca,'Ylim',[-10 10], 'ytick',[-9:3:9]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaE_{\perp 1}','\deltaE_{\perp 2}','\deltaE_{||}'},[0.05 0.06]);
% irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.05 0.06]);
%% Especperp
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(gca,jet);
[h(4), hcb]=irf_spectrogram(h(4),specperpE); hold on;
irf_plot(h(4),Fcp,'color','k','LineWidth',1.5);hold on;
irf_plot(h(4),Flh,'color','k','LineWidth',1.5);hold on;
irf_plot(h(4),Fpe,'color','k','LineWidth',1.5);hold on;
irf_plot(h(4),Fce,'color','r','LineWidth',1.5);hold on;
irf_plot(h(4),Fce01,'color','w','LineWidth',1.5);hold on;
irf_plot(h(4),Fce05,'color','c','LineWidth',1.5);hold off;
grid off;
set(gca,'yscale','log');
set(gca,'ytick',[1e1 1e2 1e3 1e4]);
set(gca,'ylim',[8 1e4]);
caxis(h(4),[-6 -2]);
ylabel(h(4),'f (Hz)','fontsize',12);
poscbar=get(hcb,'pos');
poscbar(3)=poscbar(3)*0.5;
set(hcb,'pos',poscbar);
set(hcb,'fontsize',10);
%% Especpar
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(gca,jet);
[h(5), hcb]=irf_spectrogram(h(5),specparE); hold on;
irf_plot(h(5),Fcp,'color','k','LineWidth',1.5);hold on;
irf_plot(h(5),Flh,'color','k','LineWidth',1.5);hold on;
irf_plot(h(5),Fpe,'color','k','LineWidth',1.5);hold on;
irf_plot(h(5),Fce,'color','r','LineWidth',1.5);hold on;
irf_plot(h(5),Fce01,'color','w','LineWidth',1.5);hold on;
irf_plot(h(5),Fce05,'color','c','LineWidth',1.5);hold off;

grid off;
set(gca,'yscale','log');
set(gca,'ytick',[1e1 1e2 1e3 1e4]);
set(gca,'ylim',[8 1e4]);
caxis(h(5),[-6 -2]);
ylabel(h(5),'f (Hz)','fontsize',12);
poscbar=get(hcb,'pos');
poscbar(3)=poscbar(3)*0.5;
set(hcb,'pos',poscbar);
set(hcb,'fontsize',10);
%% Bfachf plot 1
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf1(:,1) Bfachf1(:,2)*0],'k--', 'Linewidth',0.75);hold on;

if flag ==2
irf_plot([Bfachf2(:,1) Bfachf2(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([Bfachf2(:,1) Bfachf2(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([Bfachf2(:,1) Bfachf2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([Bfachf2(:,1) Bfachf2(:,2)*0],'k--', 'Linewidth',0.75);
end
hold off; grid off;
ylabel('\deltaB [nT]');
set(gca,'ylim',[min(min(Bfachf1(:,2:4))) max(max(Bfachf1(:,2:4)))]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\deltaB_{\perp 1}','\deltaB_{\perp 2}','\deltaB_{||}'},[0.05 0.06]);
irf_legend(gca,{'f>',num2str(lf),'Hz'},[0.8 0.12]);
%% Bspec 1
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(gca,jet);
[h(7), hcb]=irf_spectrogram(h(7),specB1); hold on;
if flag ==2
[h(7), hcb]=irf_spectrogram(h(7),specB2); hold on;
end
irf_plot(h(7),Fcp,'color','k','LineWidth',1.5);hold on;
irf_plot(h(7),Flh,'color','k','LineWidth',1.5);hold on;
irf_plot(h(7),Fce,'color','r','LineWidth',1.5);hold on;
irf_plot(h(7),Fce01,'color','w','LineWidth',1.5);hold on;
irf_plot(h(7),Fce05,'color','c','LineWidth',1.5);hold off;
grid off;
set(gca,'yscale','log');
set(gca,'ytick',[1e-1 1e0 1e1 1e2 1e3 1e4]);
set(gca,'ylim',[2 1e4]);
set(gca,'XTickLabelRotation',0)
caxis(h(7),[-7 -2]);
ylabel(h(7),'f (Hz)','fontsize',12);
poscbar=get(hcb,'pos');
poscbar(3)=poscbar(3)*0.5;
set(hcb,'pos',poscbar);
set(hcb,'fontsize',10);
 %% J dot E' plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
c_eval("irf_plot([JdotEplus?(:,1) abs(JdotEplus?(:,2))], 'color','k', 'Linewidth',0.75); hold on;",ic)
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([JdotEplus?(:,1) JdotEplus?(:,2)*0],'k--', 'Linewidth',0.75); hold off;",ic)
grid off;
ylabel('|J\dotE| [pw/m]','fontsize',10);
% set(gca,'Ylim',[fix(min([min(Vi1_gsm(:,2)) min(Vi1_gsm(:,3)) min(Vi1_gsm(:,4))])/10)*10-10 fix(max(Vit1(:,2))/10)*10+10]);
% set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
 %% J dot E' integrate plot
% % % h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% % % c_eval("irf_plot([JdotEint(:,1) JdotEint(:,2)], 'color','k', 'Linewidth',0.75); hold on;",ic)
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([JdotEint(:,1) JdotEint(:,2)*0],'k--', 'Linewidth',0.75); hold off;",ic)
% % % grid off;
% % % ylabel('JdotE [pw/m]','fontsize',10);
% % % % set(gca,'Ylim',[fix(min([min(Vi1_gsm(:,2)) min(Vi1_gsm(:,3)) min(Vi1_gsm(:,4))])/10)*10-10 fix(max(Vit1(:,2))/10)*10+10]);
% % % % set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
 %% 
set(gca,"XTickLabelRotation",0)
colormap(jet);

irf_plot_axis_align(h(1:end));
irf_zoom(h(1:end),'x',Tint);

% irf_pl_mark(h(1:7),[iso2epoch('2015-10-16T13:04:26Z')],'k');
% irf_pl_mark(h(1:7),[iso2epoch('2015-12-14T00:59:04Z')],'k');

set(h(1:end),'fontsize',12);
set(gcf,'paperpositionmode','auto');
% irf_adjust_panel_position;

figname=['/Users/fwd/Documents/lalala/硕士毕业论文/中期/wave&JdotE'];
print(gcf, '-dpdf', [figname '.pdf']);