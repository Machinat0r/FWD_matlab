clear;clc

global ParentDir 
ParentDir = 'D:\MMS\'; 
TempDir = 'D:\MMS\temp\';mkdir(TempDir);
% TT = '2017-08-20T02:01:00.00Z/2017-08-20T02:03:00.00Z';
% TT = '2017-08-07T16:01:00.00Z/2017-08-07T16:01:45.00Z';
TT = '2017-08-07T15:57:10.00Z/2017-08-07T15:57:40.00Z';
% TT = '2017-08-07T16:31:30.00Z/2017-08-07T16:33:00.00Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:4;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3,filenames4];

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

if flag == 0
    tempTag = num2str(NameTags{i-1});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'\fpi\brst\l2\des-moms\',tempTag(1:4),'\',tempTag(5:6),'\',...
            tempTag(7:8),'\',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
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
    desmoms1 = [ParentDir,'mms',num2str(ic),'\fpi\brst\l2\des-moms\',tempTag1(1:4),'\',tempTag1(5:6),'\',...
            tempTag1(7:8),'\',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'\fpi\brst\l2\des-moms\',tempTag2(1:4),'\',tempTag2(5:6),'\',...
            tempTag2(7:8),'\',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
end
SDCFilesDownload(filenames,TempDir)
SDCFilesDownload(filenames_srvy(1),TempDir)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% load data
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

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
c_eval('E? = irf_resamp(E?,E1);',ic)
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval('Et? = irf_resamp(Et?,E1);',ic)
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);
c_eval(['E?_err = mms.db_get_variable(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);'],ic);
c_eval("E?_err.DEPEND_0.data = irf_time(double(E?_err.DEPEND_0.data),'ttns>epoch');",ic)
c_eval('E?_err = [E?_err.DEPEND_0.data,double(E?_err.data(:,1))];',ic);
c_eval('E?_err = irf_resamp(E?_err,Et?);')
c_eval('E?_err(:,2) = abs(E?_err(:,2))./Et?(:,2)*100;')


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
c_eval('Evixb? = 1e3*cross(1000*gsmVi?_res(:,2:4),1e-9*B?(:,2:4));',ic); %mV/m
c_eval('Evixb? = Evixb? + E?_resamp(:,2:4);',ic);

c_eval('gsmVe?_res = irf_resamp(gsmVe?,B?);',ic);
c_eval('Evexb? = 1e3*cross(1000*gsmVe?_res(:,2:4),1e-9*B?(:,2:4));',ic); %mV/m
c_eval('Evexb? = Evexb? + E?_resamp(:,2:4);',ic);

% merge data/time from 2 cdf files
c_eval('energy_low?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_lowen_brst'',tint);',ic)
c_eval('energy_mid?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_miden_brst'',tint);',ic)
c_eval('energy_high?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_highen_brst'',tint);',ic)
c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energyspectr_omni_brst'',tint);',ic)
c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_energyspectr_omni_brst'',tint);',ic)

%Pressure
%Pm
miu0 = 4*pi*10^(-7);
c_eval('Pm = [Bt?(:,1) 10^(-9)*Bt?(:,2).^2 / (2*miu0)];',ic); %nPa
%Pt
c_eval("Pte = [Ne?(:,1) 1.6*10^(-19)*10^(15)*Ne?(:,2).*Te?(:,2)];",ic);
c_eval("Pte_para = [Ne?(:,1) 1.6*10^(-19)*10^(15)*Ne?(:,2).*Te_para?(:,2)];",ic);
c_eval("Pte_perp = [Ne?(:,1) 1.6*10^(-19)*10^(15)*Ne?(:,2).*Te_perp?(:,2)];",ic);
c_eval("Pti = [Ni?(:,1) 1.6*10^(-19)*10^(15)*Ni?(:,2).*Ti?(:,2)];",ic);
c_eval("Pti_para = [Ni?(:,1) 1.6*10^(-19)*10^(15)*Ni?(:,2).*Ti_para?(:,2)];",ic);
c_eval("Pti_perp = [Ni?(:,1) 1.6*10^(-19)*10^(15)*Ni?(:,2).*Ti_perp?(:,2)];",ic);
c_eval('Pte_res = irf_resamp(Pte,Pti);',ic);
c_eval('Pthe = [Pti(:,1) Pte_res(:,2) + Pti(:,2)];',ic);
c_eval('Pthe_res = irf_resamp(Pthe,Pm);',ic);
c_eval('Ptotal = [Pm(:,1) Pthe_res(:,2) + Pm(:,2)];',ic);
c_eval('Beta = [Pm(:,1) Pthe_res(:,2)./Pm(:,2)];',ic);
Pthe_para = irf_resamp(Pte_para,Pti_para)+Pti_para;
Pthe_perp = irf_resamp(Pte_perp,Pti_perp)+Pti_perp;
Pthe_para(:,1) = Pti_para(:,1);Pthe_perp(:,1) = Pti_perp(:,1);

Pm_res = irf_resamp(Pm,Pte);
Pthe_perp = irf_resamp(Pthe_perp,Pte);
c_eval('K = Te_perp?./Te_para?-1-(Pm_res./Pthe_perp);',ic);

%S
% % % miu0=400*pi;
% % % c_eval('Br? = sqrt(B?(:,2).^2+B?(:,3).^2);',ic);
% % % c_eval('BzE? = B?(:,4).*sqrt(1+(Br?.^2./(B?(:,4).^2+2*miu0*Pthe_interp?)));',ic);
% % % % c_eval('BzE? = B?(:,3).*sqrt(1+(Br?.^2./(B?(:,3).^2+2*miu0*Pthe_interp)));',ic);
% % % c_eval('PE? = Pthe_interp?+Br?.^2/(2*miu0);',ic);
% % % C = 0.7368;D = 0.7634;F = -0.3059;
% % % Pos = mms.get_data('R_gsm',tint);
% % % Pos = Pos.gsmR1/6372;
% % % c_eval('VE? = 10^C*(sqrt(Pos(1,1)^2+Pos(1,2)^2)^D).*(BzE?.^F)./sqrt(BzE?.^2+2*miu0*PE?);',ic);
% % % c_eval('S? = Pthe_interp?.*(real(VE?).^(5/3));',ic);
% % % c_eval('S?(S?>=5)=0;',ic);

%u·delta u
c_eval('Vixdiff? = diff(Vi?(:,2));',ic);
c_eval('Vixdiff?(end+1) = Vixdiff?(end);',ic);
c_eval('templen = length(Vixdiff?);',ic);
for temp = 2:templen
    c_eval('Vixdiff?(temp) = mean(Vixdiff?(temp-1:temp));',ic);
end
c_eval('UdU = 1836*units.me*1e12*miu0*gsmVi?(:,2).*Vixdiff?.*Ni?(:,2);',ic);
c_eval('Bdiff? = diff(B?(:,2));',ic);
c_eval('Bdiff?(end+1) = Bdiff?(end);',ic);
c_eval('templen = length(Bdiff?);',ic);
for temp = 2:templen
    c_eval('Bdiff?(temp) = mean(Bdiff?(temp-1:temp));',ic);
end
c_eval('BdB = 1e-18*B?(:,2).*Bdiff?;',ic);
UdU(:,2) = UdU; BdB(:,2) = BdB;
c_eval('UdU(:,1) = Vi?(:,1);',ic); 
c_eval('BdB(:,1) = B?(:,1);',ic);
momentum_equation = irf_resamp(UdU,BdB)-BdB;


%J
c_eval('Je?_ts = -units.e*Ne?_ts*gsmVe?_ts*1e3*1e6*1e9;',ic);
c_eval('Ji?_ts = units.e*Ne?_ts*gsmVi?_ts.resample(Ne?_ts.time)*1e3*1e6*1e9;',ic);
c_eval('J?_ts = (Je?_ts+Ji?_ts);',ic);
c_eval('J? = irf.ts2mat(J?_ts);',ic);

%JdotE
c_eval('E_resJ = irf_resamp(E?,J?);',ic);
c_eval('JdotEtemp = [E_resJ(:,1) irf_dot(J?(:,2:4),E_resJ(:,2:4))];',ic);
c_eval('JdotE? = irf_resamp(JdotEtemp,Ni?);',ic);

%JdotE_B
Pos = mms.get_data('R_gsm',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gsmR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
[J_B,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
J_B(:,2:4) = 1e9*J_B(:,2:4);
c_eval('J_B = irf_resamp(J_B,E?);',1);
c_eval('JdotE_B? = [E?(:,1) irf_dot(J_B(:,2:4),E?(:,2:4))];',ic);

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
% c_eval(['fpiFilee1 = dataobj(','''',desmoms1,'''',');'],ic);
% c_eval('energy_low1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_lowen_brst'');',ic);
% c_eval('energy_mid1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_miden_brst'');',ic);
% c_eval('energy_high1 = get_variable(fpiFilee1,''mms?_des_pitchangdist_highen_brst'');',ic);
% c_eval('energy_e1 = get_variable(fpiFilee1,''mms?_des_energyspectr_omni_brst'');',ic);
% 
% energy_low.DEPEND_0.data=energy_low1.DEPEND_0.data;
% energy_mid.DEPEND_0.data=energy_mid1.DEPEND_0.data;
% energy_high.DEPEND_0.data=energy_high1.DEPEND_0.data;
% energy_low=energy_low1; energy_low.data=energy_low1.data; energy_low.nrec=energy_low1.nrec;
% energy_mid=energy_mid1; energy_mid.data=energy_mid1.data;  energy_mid.nrec=energy_mid1.nrec;
% energy_high=energy_high1; energy_high.data=energy_high1.data;  energy_high.nrec=energy_high1.nrec;
% energy_e=energy_e1; energy_e.data=energy_e1.data;  energy_e.nrec=energy_e1.nrec;
end

%% Init figure
n=5;
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
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([B?(:,1) B?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([B?(:,1) B?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([B?(:,1) B?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
% % % ylabel('B [nT]','fontsize',10);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Bx & Bz plot
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([B2(:,1) B2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([B3(:,1) B3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([B4(:,1) B4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,4))])-1 max([max(B?(:,2)),max(B?(:,4))])+1]);",iic);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
irf_legend(gca,{'Bx','Bz'},[0.97 0.92]);
ylabel('Bx Bz [nT]','fontsize',12);
i=i+1;
%% By plot
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([B2(:,1) B2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([B3(:,1) B3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([B4(:,1) B4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(B?(:,3))])-1 max(B?(:,3))+1]);",iic);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[0 1 0]);
irf_legend(gca,{'By'},[0.97 0.92]);
ylabel('By [nT]','fontsize',12);
i=i+1;
%% Efac
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % ylabel('Efac [mV/m]','fontsize',12)
% % % % set(gca,'Ylim',[-10 10], 'ytick',[-50 -25 0 25 50 75 100]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[fix(min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))]))-1 fix(max(Et?(:,2)))+1]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
% % % irf_legend(gca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.1 0.12]);
% % % 
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % i=i+1;
%% Vi plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)], 'color','b', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,3)], 'color','g', 'Linewidth',0.75);",1); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,4)], 'color','r', 'Linewidth',0.75);",1); hold on;
% c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([gsmVi?(:,1) gsmVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+10],'fontsize',9);",1);
%set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',12);
i=i+1;
%% E+VixB field
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([B?(:,1) Evixb?(:,1)], 'color','b', 'Linewidth',0.75);",1);hold on;
c_eval("irf_plot([B?(:,1) Evixb?(:,2)], 'color','g', 'Linewidth',0.75);",1);hold on;
c_eval("irf_plot([B?(:,1) Evixb?(:,3)], 'color','r', 'Linewidth',0.75);",1);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",1); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(Evixb?(:,1)) min(Evixb?(:,2)) min(Evixb?(:,3))])-2 max([max(Evixb?(:,1)) max(Evixb?(:,2)) max(Evixb?(:,3))])+2]);",1);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E+VixB_x','E+VixB_y','E+VixB_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',12)
i=i+1;
%% Eerr
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([E?_err(:,1) E?_err(:,2)], 'color','k', 'Linewidth',0.75);",1); hold on;
% c_eval("irf_plot([Efac?(:,1) Efac?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Efac?(:,1) Efac?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Efac?(:,1) Efac?(:,2)*0],'k--', 'Linewidth',0.75);",1); hold off;
grid off;
ylabel('E_e_r_r [%]','fontsize',12)
set(gca,'Ylim',[0 100], 'ytick',[-50 -25 0 25 50 75 100]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% c_eval("set(gca,'Ylim',[fix(min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))]))-1 fix(max(Et?(:,2)))+1]);",ic);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
% irf_legend(gca,{'E_e_r_r'},[0.1 0.12]);

pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
i=i+1;
%% JdotE
% % % h(i)=irf_subplot(n,1,-i);
% % % cor = 'bgrk';
% % % c_eval("irf_plot([JdotE_B?(:,1) JdotE_B?(:,2)], 'color',cor(?), 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([JdotE_B?(:,1) JdotE_B?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % c_eval("set(gca,'Ylim',[fix(min(JdotE?)/10)*10-10 fix(max(JdotE?)/10)*10+10],'fontsize',9);",ic);
% % % %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[0 0 0]);
% % % % irf_legend(gca,{'JdotE'},[0.97 0.92]);
% % % ylabel('J\cdotE [pW/m^3] ','fontsize',10);
% % % i=i+1;
%% E+VexB field
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([B?(:,1) Evexb?(:,1)], 'color','b', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([B?(:,1) Evexb?(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([B?(:,1) Evexb?(:,3)], 'color','r', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% % % % set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[min([min(Evexb?(:,1)) min(Evexb?(:,2)) min(Evexb?(:,3))])-2 max([max(Evexb?(:,1)) max(Evexb?(:,2)) max(Evexb?(:,3))])+2]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'E+VexB_x','E+VexB_y','E+VexB_z'},[0.97 0.92]);
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('E [mV/m]','fontsize',12)
% % % i=i+1;

%   set(h(1:n),'fontsize',8);
% %   irf_zoom(tint,'x',h(1:n));
irf_zoom(tint,'x',h(1:n));
% irf_adjust_panel_position;
% %   irf_plot_axis_align(h)
irf_plot_axis_align(h)

%   irf_pl_mark(h(1:n),[iso2epoch('2017-07-24T12:56:52.00Z')],'g');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.159Z')],'k');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.399Z')],'k');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.589Z')],'b');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:29.789Z')],'k'); 
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30Z')],'g');
%   irf_pl_mark(h(1:8),[iso2epoch('2015-10-16T13:04:30.209Z')],'k');
%  add_position(gca,gseR1), xlabel(gca,'')
%  irf_zoom(tintlmn,'x',h(4:7))

%%  出图保存部分
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
cd  C:\Matlab\bin\新建文件夹\fwd\
rmdir(TempDir,'s'); 
Dir = 'C:\Users\fwd\Desktop\Ti~mor~\M\Electron Fermi acceleration by flow vortex\2-Figures\Figure4\';
figname = [Dir TTlist(1:8) 'T' TTlist(9:14)];
    print(gcf, '-dpdf', [figname '.pdf']);    