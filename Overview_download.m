close all
clear;clc

global ParentDir 
ParentDir = '/Users/fwd/Documents/MATLAB/MMS/'; 
TempDir = '/Users/fwd/Documents/MATLAB/MMS/temp/';mkdir(TempDir);
% TT = '2021-08-15T03:35:15.00Z/2021-08-15T03:35:30.00Z';
% TT = '2021-08-22T06:39:30.00Z/2021-08-22T06:43:00.00';
% TT = '2018-02-06T13:29:00.00Z/2018-02-06T13:30:30.00Z';
TT = '2019-08-05T16:24:00.00Z/2019-08-05T16:25:00.00Z';
% TT = '2015-11-04T04:34:00.00Z/2015-11-04T04:37:00.00Z';
% TT = '2018-08-27T12:15:30.00Z/2018-08-27T12:16:30.00Z';
% TT = '2019-08-16T01:03:33.00Z/2019-08-16T01:05:13.00Z';
% TT='2017-08-23T15:38:30.00Z/2017-08-23T15:39:15.00Z';
% TT = '2021-07-10T12:41:23.00Z/2021-07-10T12:42:23.00Z';
% TT = '2018-07-03T15:50:00.00Z/2018-07-03T15:51:00.00Z';
% TT = '2017-08-20T02:01:30.00Z/2017-08-20T02:03:00.00Z';
% TT = '2017-08-07T16:01:00.00Z/2017-08-07T16:02:00.00Z';
% TT = '2021-07-21T12:46:20.00Z/2021-07-21T12:46:40.00Z';
% TT = '2017-05-05T20:06:30.00Z/2017-05-05T20:07:10.00';
% TT = '2022-08-18T23:53:00.00Z/2022-08-18T23:54:00.00Z';
% TT = '2022-08-19T01:13:40.00Z/2022-08-19T01:14:40.00Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1;
iic = 1:4;
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
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
    if i > length(NameTags)
        break
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
% SDCFilesDownload(filenames_srvy(1),TempDir)
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
%% pressure & entropy (old)
%Pressure
%Pm
units = irf_units;
% % % c_eval('Pm = [Bt?(:,1) 10^(-9)*Bt?(:,2).^2 / (2*units.mu0)];',ic); %nPa
% % % Pm = irf_resamp(Pm,Ni1);
% % % Pte_para = [Ne1(:,1) units.e*10^(15)*Ne1(:,2).*Te_para1(:,2)];
% % % Pte_perp = [Ne1(:,1) units.e*10^(15)*Ne1(:,2).*Te_perp1(:,2)];
% % % Pti_para = [Ni1(:,1) units.e*10^(15)*Ni1(:,2).*Ti_para1(:,2)];
% % % Pti_perp = [Ni1(:,1) units.e*10^(15)*Ni1(:,2).*Ti_perp1(:,2)];
% % % Pthe_para = irf_resamp(Pte_para,Pti_para)+Pti_para;
% % % Pthe_perp = irf_resamp(Pte_perp,Pti_perp)+Pti_perp;
% % % Beta_para = [Ni1(:,1) Pthe_para(:,2)./Pm(:,2)];
% % % Beta_perp = [Ni1(:,1) Pthe_perp(:,2)./Pm(:,2)];
%Pt
% % % c_eval("Pte = [Ne?(:,1) 1.6*10^(-19)*10^(15)*Ne?(:,2).*Te?(:,2)];",ic);
% % % c_eval("Pte_para = [Ne?(:,1) 1.6*10^(-19)*10^(15)*Ne?(:,2).*Te_para?(:,2)];",ic);
% % % c_eval("Pte_perp = [Ne?(:,1) 1.6*10^(-19)*10^(15)*Ne?(:,2).*Te_perp?(:,2)];",ic);
% % % c_eval("Pti = [Ni?(:,1) 1.6*10^(-19)*10^(15)*Ni?(:,2).*Ti?(:,2)];",ic);
% % % c_eval("Pti_para = [Ni?(:,1) 1.6*10^(-19)*10^(15)*Ni?(:,2).*Ti_para?(:,2)];",ic);
% % % c_eval("Pti_perp = [Ni?(:,1) 1.6*10^(-19)*10^(15)*Ni?(:,2).*Ti_perp?(:,2)];",ic);
% % % c_eval('Pte_res = irf_resamp(Pte,Pti);',ic);
% % % c_eval('Pthe = [Pti(:,1) Pte_res(:,2) + Pti(:,2)];',ic);
% % % c_eval('Pthe_res = irf_resamp(Pthe,Pm);',ic);
% % % c_eval('Ptotal = [Pm(:,1) Pthe_res(:,2) + Pm(:,2)];',ic);
% % % c_eval('Beta= [Pm(:,1) Pthe_res(:,2)./Pm(:,2)];',ic);
% % % Pthe_para = irf_resamp(Pte_para,Pti_para)+Pti_para;
% % % Pthe_perp = irf_resamp(Pte_perp,Pti_perp)+Pti_perp;
% % % 
% % % Pthe_para(:,1) = Pti_para(:,1);Pthe_perp(:,1) = Pti_perp(:,1);
% % % 
% % % Pm_res = irf_resamp(Pm,Pte);
% % % Pthe_perp = irf_resamp(Pthe_perp,Pte);
% % % c_eval('K = Te_perp?./Te_para?-1-(Pm_res./Pthe_perp);',ic);

%S
% % % miu0=400*pi;
% % % % % % c_eval('Br? = sqrt(B?(:,2).^2+B?(:,3).^2);',ic);
% % % % % % c_eval('BzE? = B?(:,4).*sqrt(1+(Br?.^2./(B?(:,4).^2+2*miu0*Pthe_interp?)));',ic);
% % % % % % % c_eval('BzE? = B?(:,3).*sqrt(1+(Br?.^2./(B?(:,3).^2+2*miu0*Pthe_interp)));',ic);
% % % % % % c_eval('PE? = Pthe_interp?+Br?.^2/(2*miu0);',ic);
% % % % % % C = 0.7368;D = 0.7634;F = -0.3059;
% % % Pos = mms.get_data('R_gsm',tint);
% % % Pos = Pos.gsmR1/units.RE/1000;
% % % c_eval('v?=10^0.7368.*sqrt(R?(:,2).^2+R?(:,3).^2).^0.7634.*B?(:,4).^(-0.3059)./sqrt(B?(:,4).^2+2.*miu0.*Pthe?(:,2));',ic);
% % % c_eval('etr?=Pt?(:,2).*v?.^(5/3);',ic);
% % % c_eval('VE? = 10^C*(sqrt(Pos(1,1)^2+Pos(1,2)^2)^D).*(BzE?.^F)./sqrt(BzE?.^2+2*miu0*PE?);',ic);
% % % c_eval('S? = Pthe_interp?.*(real(VE?).^(5/3));',ic);
% % % c_eval('S?(S?>=5)=0;',ic);

%u·delta u
% % % c_eval('Vixdiff? = diff(Vi?(:,2));',ic);
% % % c_eval('Vixdiff?(end+1) = Vixdiff?(end);',ic);
% % % c_eval('templen = length(Vixdiff?);',ic);
% % % for temp = 2:templen
% % %     c_eval('Vixdiff?(temp) = mean(Vixdiff?(temp-1:temp));',ic);
% % % end
% % % c_eval('UdU = 1836*units.me*1e12*miu0*gsmVi?(:,2).*Vixdiff?.*Ni?(:,2);',ic);
% % % c_eval('Bdiff? = diff(B?(:,2));',ic);
% % % c_eval('Bdiff?(end+1) = Bdiff?(end);',ic);
% % % c_eval('templen = length(Bdiff?);',ic);
% % % for temp = 2:templen
% % %     c_eval('Bdiff?(temp) = mean(Bdiff?(temp-1:temp));',ic);
% % % end
% % % c_eval('BdB = 1e-18*B?(:,2).*Bdiff?;',ic);
% % % UdU(:,2) = UdU; BdB(:,2) = BdB;
% % % c_eval('UdU(:,1) = Vi?(:,1);',ic); 
% % % c_eval('BdB(:,1) = B?(:,1);',ic);
% % % momentum_equation = irf_resamp(UdU,BdB)-BdB;
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
if length(B1) ~= length(B2) || length(B1) ~= length(B3) || length(B1) ~= length(B4)
c_eval('B? = irf_resamp(B?,B1);',2:4);
end
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
[J_B,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
J_B(:,2:4) = 1e9*J_B(:,2:4);
c_eval('E_resJ = irf_resamp(E?,J_B);',ic);
c_eval('JdotE_B = [E_resJ(:,1) irf_dot(J_B(:,2:4),E_resJ(:,2:4))];',ic);

%% lmn
% irf_minvar_gui(B1);
% % % L=[0.69 0.69 -0.24];%最大变化方向
% % % M=[0.72 -0.69 0.08];%N x L
% % % N=[0.10 0.23 0.97];%外法向

L=[0.90 -0.43 0.07];%最大变化方向
N=[0.03 0.23 0.97];%外法向
M=cross(N,L);

%65，21

% L=[0 0 1];
% M=[0 1 0];
% N=[1 0 0];
% for ic=1:1,
c_eval('Blmn?=irf_newxyz(B?,L,M,N);',ic);

% % % c_eval(['Elmn?=irf_newxyz(E?,N,M,L);'],ic);
% % % 
% % % c_eval(['lmnVe?_ts=irf_newxyz(gsmVe?_ts,N,M,L);'],ic);
% % % c_eval("lmnVe? = irf.ts2mat(lmnVe?_ts);",ic);
% % % c_eval(['lmnVi?_ts=irf_newxyz(gsmVi?_ts,N,M,L);'],ic);
% % % c_eval("lmnVi? = irf.ts2mat(lmnVi?_ts);",ic);
% % % c_eval('lmnJe?_ts = -units.e*Ne?_ts*lmnVe?_ts*1e3*1e6*1e9;',ic);
% % % c_eval('lmnJi?_ts = units.e*Ne?_ts*lmnVi?_ts.resample(Ne?_ts.time)*1e3*1e6*1e9;',ic);
% % % c_eval('lmnJ?_ts = (lmnJe?_ts+lmnJi?_ts);',ic);
% % % c_eval('lmnJ? = irf.ts2mat(lmnJ?_ts);',ic);
% end
%% Init figure
n=12;
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
c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([B?(:,1) B?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(B?(:,2)) min(B?(:,3)) min(B?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
% set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z','|B|'},[0.97 0.92]);
ylabel('B [nT]','fontsize',10);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bx plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) B1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) B2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) B3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) B4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,2))])-0.1 max(B?(:,2))+0.1]);",iic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('Bx [nT]','fontsize',8);
% % % i=i+1;
%% By plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) B1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) B2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) B3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) B4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,3))])-0.1 max(B?(:,3))+0.1]);",iic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('By [nT]','fontsize',8);
% % % i=i+1;
%% Bz plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) B1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) B2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) B3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) B4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(B?(:,4))])-0.1 max(B?(:,4))+0.1]);",iic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('Bz [nT]','fontsize',8);
% % % i=i+1;
%% Btotal plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([Bt2(:,1) Bt2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([Bt3(:,1) Bt3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Bt4(:,1) Bt4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Bt?(:,1) 0*Bt?(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min([min(Bt?(:,2))])-0.1 max(Bt?(:,2))+0.1]);",iic);
% % % % set(gca,'Ylim',[-10 20], 'ytick',[-30 -20 -10 0 10 20 30],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('Btotal [nT]','fontsize',8);
% % % i=i+1;
%% diffB plot
% h(i)=irf_subplot(n,1,-i);
% % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% B1(:,5)=diff(B1(:,4));
% 
% irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,5)], 'color','g', 'Linewidth',0.75); hold on
% % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
% % irf_plot([Bt1(:,1) Bt1], 'color','k', 'Linewidth',0.75); hold on;
% % irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% grid off;
% ylabel('B [nT]','fontsize',10);
% % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% set(gca,'Ylim',[-5 15], 'ytick',[ -5 0 5 10]);
% pos1=get(gca,'pos');
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% i=i+1;

%% Bsmall plot
% h(i)=irf_subplot(n,1,-i);
% irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
% irf_plot([Bt1(:,1) Bt1], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% grid off;
% ylabel('B [nT]','fontsize',10);
% % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% set(gca,'Ylim',[5 10], 'ytick',[  0 5 6 7 8 9 10]);
% pos1=get(gca,'pos');
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% i=i+1;

%% Bbf plot
% h(i)=irf_subplot(n,1,-i);
% % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Bbf1(:,1) Bbf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([Bbf1(:,1) Bbf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([Bbf1(:,1) Bbf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
% % irf_plot([Bbft1(:,1) Bbft1], 'color','k', 'Linewidth',0.75); hold on;
% % irf_plot([Bbf1(:,1) Bbf1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% grid off;
% ylabel('Bbf [nT]','fontsize',10);
% % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% set(gca,'Ylim',[-0.5 0.5], 'ytick',[0 5 6 7 8 9 10]);
% pos1=get(gca,'pos');
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% i=i+1;

%% Blmn plot
% % % h(i)=irf_subplot(n,1,-i);
% % % % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([Bt?(:,1) Bt?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Blmn?(:,1) Blmn?(:,4)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Blmn?(:,1) Blmn?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Blmn?(:,1) Blmn?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Blmn?(:,1) Blmn?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Blmn?(:,1) Blmn?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Blmn?(:,1) Blmn?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
% % % % irf_plot([Bbft1(:,1) Bbft1], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([Blmn1(:,1) Blmn1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% % % grid off;
% % % ylabel('Blmn [nT]','fontsize',10);
% % % c_eval("set(gca,'Ylim',[min([min(Blmn?(:,2)) min(Blmn?(:,3)) min(Blmn?(:,4))])-1 max(Bt?(:,2))+1]);",ic);
% % % % set(gca,'Ylim',[-5 4], 'ytick',[-0.5 0 0.5]);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'B_l','B_m','B_n'},[0.05 0.92]);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;


%% dfN plot
% h(i)=irf_subplot(n,1,-i);
% irf_plot([Ne1(:,1) Ne1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
% irf_plot([Ni1(:,1) Ni1(:,2)], 'color','g', 'Linewidth',0.75); hold off;
% grid off;
% ylabel('N [cm^{-3}]','fontsize',10);
% %set(gca,'Ylim',[fix(min([min(Ne1(:,2)) min(Ni1(:,2))]))-2 fix(max([max(Ne1(:,2)) max(Ni1(:,2))]))+2]);
% set(gca,'Ylim',[0.3 1], 'ytick',[0.2 0.5 0.8]);
% % pos1=get(h(1),'pos');
%  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
%  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
% i=i+1;

%% N plot
h(i)=irf_subplot(n,1,-i);

%滤波
%     irf_plot([Nebf1(:,1) Nebf1(:,2)], 'color','b', 'Linewidth',0.75);hold on;
%     irf_plot([Nibf1(:,1) Nibf1(:,2)], 'color','g', 'Linewidth',0.75); hold off;

%非滤波
c_eval("irf_plot([Ne?(:,1) Ne?(:,2)], 'color','b', 'Linewidth',0.75);",ic);hold on;
c_eval("irf_plot([Ni?(:,1) Ni?(:,2)], 'color','g', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[max([0 min([min(Ne?(:,2)) min(Ni?(:,2))])-0.02]) max([max(Ne?(:,2)) max(Ni?(:,2))])+0.02]);",ic)
%     set(gca,'Ylim',[0.15 0.45], 'ytick',[0.1 0.2 0.3 0.4],'fontsize',9);
% pos1=get(h(1),'pos');
%  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
%  irf_legend(gca,{'Ne','Ni'},[0.1 0.12]);
  set(gca,'ColorOrder',[[0 0 1];[0 1 0]]);
 irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',12)
ylabel('N [cm^{-3}]','fontsize',8);
i=i+1;



%% deElectric field
% h(i)=irf_subplot(n,1,-i);
% irf_plot([Ebf1(:,1) Ebf1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([Ebf1(:,1) Ebf1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([Ebf1(:,1) Ebf1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([Ebf1(:,1) Ebf1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% grid off;
% ylabel('Ebf [mV/m]','fontsize',10)
% % set(gca,'Ylim',[-1 3], 'ytick',[-0.5 0 1 2  ]);
% % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% %set(gca,'Ylim',[fix(min([min(E1(:,2)) min(E1(:,3)) min(E1(:,4))])/10)*10-10 fix(max(Et1(:,2))/10)*10+10]);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'dfE_x','dfE_y','dfE_z'},[0.1 0.12]);
% pos3=get(gca,'pos');
% set(gca,'ColorOrder',[[0 1 0]]);
% %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% i=i+1;
%% BandE plot
% h(i)=irf_subplot(n,1,-i);
% B1(:,5)=(Bbf1(:,2).^2+Bbf1(:,3).^2+Bbf1(:,4).^2).^0.5/9
% B1(:,6)=(Ebf1(:,2).^2+Ebf1(:,3).^2+Ebf1(:,4).^2).^0.5
% % irf_plot([Bt1(:,1) Bt1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,5)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([B1(:,1) B1(:,6)], 'color','g', 'Linewidth',0.75); hold on;
% 
% 
% grid off;
% ylabel('B [nT]','fontsize',10);
% % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% % set(gca,'Ylim',[-5 15], 'ytick',[ -5 0 5 10]);
% pos1=get(gca,'pos');
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% i=i+1;
%% Efac
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Efac?(:,1) Efac?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % ylabel('Efac [mV/m]','fontsize',8)
% % % % set(gca,'Ylim',[-10 10], 'ytick',[-50 -25 0 25 50 75 100]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[min([min(Efac?(:,2)) min(Efac?(:,3)) min(Efac?(:,4))])-1 max(Et?(:,2))+1]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
% % % irf_legend(gca,{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.1 0.12]);
% % % 
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % i=i+1;



%% Elmn plot
% % % h(i)=irf_subplot(n,1,-i);
% % % % E_fac = irf_convert_fac(E1,B1,[1 0 0]);
% % % c_eval("irf_plot([Elmn?(:,1) Elmn?(:,4)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Elmn?(:,1) Elmn?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Elmn?(:,1) Elmn?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Et1(:,1) Et1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([E1(:,1) E1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% % % grid off;
% % % ylabel('Elmn [mV/m]','fontsize',10)
% % % %set(gca,'Ylim',[-5 5], 'ytick',[-3 0 3]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[fix(min([min(Elmn?(:,2)) min(Elmn?(:,3)) min(Elmn?(:,4))]))-1 fix(max(Et?(:,2)))+1]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'E_L','E_M','E_N'},[0.05 0.92]);
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % % irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % i=i+1;



%% Ve plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vebf?(:,1) Vebf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vebf?(:,1) Vebf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% c_eval("irf_plot([Vebf?(:,1) Vebf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% c_eval("quiver(gsmVe?(:,1),0*gsmVe?(:,1),gsmVe?(:,2),gsmVe?(:,3));",ic);hold on;
% irf_plot([Vet1(:,1) Vet1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([gsmVe?(:,1) gsmVe?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;

grid off;
ylabel('Ve [km/s]','fontsize',8);
c_eval("set(gca,'Ylim',[fix(min([min(gsmVe?(:,2)) min(gsmVe?(:,3)) min(gsmVe?(:,4))])/10)*10-10 fix(max(Vet?(:,2))/10)*10+10]);",ic);

% c_eval("set(gca,'Ylim',[-700 700]);",ic);
% set(gca,'Ylim',[-1000 1000], 'ytick',[-600 -400 -200 0 200 400 600]);
%irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Ve_N','Ve_M','Ve_L','|Ve|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Ve_x','Ve_y','Ve_z'},[0.05 0.92]);
i=i+1;

%% Ve_lmn plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([lmnVe?(:,1) lmnVe?(:,4)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([lmnVe?(:,1) lmnVe?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([lmnVe?(:,1) lmnVe?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % 
% % % % irf_plot([Vet1(:,1) Vet1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([lmnVe?(:,1) lmnVe?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % 
% % % grid off;
% % % ylabel('Ve [km/s]','fontsize',10);
% % % c_eval("set(gca,'Ylim',[fix(min([min(lmnVe?(:,2)) min(lmnVe?(:,3)) min(lmnVe?(:,4))])/10)*10-10 fix(max(Vet?(:,2))/10)*10+10]);",ic);
% % % % set(gca,'Ylim',[-800 800], 'ytick',[-600 -400 -200 0 200 400 600]);
% % % %irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Ve_N','Ve_M','Ve_L','|Ve|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Ve_L','Ve_M','Ve_N'},[0.05 0.92]);
% % % i=i+1;
%% Vexb plot
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([Vexb1(:,1) Vexb1(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% % % irf_plot([Vexb1(:,1) Vexb1(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([Vexb1(:,1) Vexb1(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % irf_plot([Vexb1(:,1) Vexb1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% % % grid off;
% % % ylabel('Vexb [km/s]','fontsize',10);
% % % set(gca,'Ylim',[fix(min([min(Vexb1(:,2)) min(Vexb1(:,3)) min(Vexb1(:,4))]))-1 fix(max(max(Vexb1(:,2:4))))+1]);
% % % % set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Vexb_x','Vexb_y','Vexb_z'},[0.05 0.92]);
% % % i=i+1;

%% Electric field
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% % % % set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))])-0.5 max([max(E?(:,2)) max(E?(:,3)) max(E?(:,4))])+0.5]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('E [mV/m]','fontsize',8)
% % % i=i+1;
%% E+VxB field
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([B?(:,1) Evxb?(:,1)], 'color','b', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([B?(:,1) Evxb?(:,2)], 'color','g', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([B?(:,1) Evxb?(:,3)], 'color','r', 'Linewidth',0.75);",ic);hold on;
% % % c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% % % % set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% % % % irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
% % % c_eval("set(gca,'Ylim',[min([min(Evxb?(:,1)) min(Evxb?(:,2)) min(Evxb?(:,3))])-2 max([max(Evxb?(:,1)) max(Evxb?(:,2)) max(Evxb?(:,3))])+2]);",ic);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'E+VixB_x','E+VixB_y','E+VixB_z'},[0.97 0.92]);
% % % pos3=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('E [mV/m]','fontsize',12)
% % % i=i+1;
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
c_eval("set(gca,'Ylim',[fix(min([min(gsmVi?(:,2)) min(gsmVi?(:,3)) min(gsmVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+10],'fontsize',9);",ic);
%set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Vi_x','Vi_y','Vi_z'},[0.97 0.92]);
ylabel('Vi [km/s]','fontsize',8);
i=i+1;

%% Vifac plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Vifac?(:,1) Vifac?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vifac?(:,1) Vifac?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Vifac?(:,1) Vifac?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % c_eval("irf_plot([Bt?(:,2) Vn], 'color','r', 'Linewidth',0.75)",ic);
% % % % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % % % % c_eval("irf_plot([Vibf?(:,1) Vibf?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([Vifac?(:,1) Vifac?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[fix(min(min(Vifac?(:,2:4)))/10)*10-10 fix(max(Vit?(:,2))/10)*10+10],'fontsize',9);",ic);
% % % %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Vi_{\perp 1}','Vi_{\perp 2}','Vi_{||}'},[0.97 0.92]);
% % % ylabel('Vi_fac [km/s]','fontsize',8);
% % % i=i+1;
%% Vi_lmn plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([lmnVi?(:,1) lmnVi?(:,4)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([lmnVi?(:,1) lmnVi?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([lmnVi?(:,1) lmnVi?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([lmnVi?(:,1) lmnVi?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[fix(min([min(lmnVi?(:,2)) min(lmnVi?(:,3)) min(lmnVi?(:,4))])/10)*10-10 fix(max(Vit?(:,2))/10)*10+10],'fontsize',9);",ic);
% % % %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'Vi_L','Vi_M','Vi_N'},[0.97 0.92]);
% % % ylabel('Vi [km/s]','fontsize',12);
% % % i=i+1;
%% momentum equation
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([B?(:,1) momentum_equation(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([B?(:,1) B?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % c_eval("set(gca,'Ylim',[fix(min(momentum_equation(:,2))/10)*10 fix(max(momentum_equation(:,2))/10)*10],'fontsize',9);",ic);
% % % %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % % irf_legend(gca,{'J_x','J_y','J_z'},[0.97 0.92]);
% % % % ylabel('J [nA/m^2]','fontsize',10);
% % % i=i+1;
%% Jplasma plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([J?(:,1) J?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([J?(:,1) J?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([J?(:,1) J?(:,4)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("Jtotal? = [J?(:,1) sqrt(J?(:,2).^2+J?(:,3).^2+J?(:,4).^2)];",ic); hold on;
% % % c_eval("irf_plot([Jtotal?(:,1) Jtotal?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([J?(:,1) J?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[fix(min([min(J?(:,2)) min(J?(:,3)) min(J?(:,4))])/10)*10-10 fix(max(max(Jtotal?(:,2)))/10)*10+10],'fontsize',9);",ic);
% % % %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'J_x','J_y','J_z'},[0.97 0.92]);
% % % ylabel('J [nA/m^2]','fontsize',8);
% % % i=i+1;
%% J_B plot
h(i)=irf_subplot(n,1,-i);
irf_plot([J_B(:,1) J_B(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([J_B(:,1) J_B(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([J_B(:,1) J_B(:,4)], 'color','r', 'Linewidth',0.75); hold on;
J_B = irf_abs(J_B);
% irf_plot([J_B(:,1) J_B(:,5)], 'color','r', 'Linewidth',0.75); hold on;
% c_eval("Jtotal_B = [J_B(:,1) sqrt(J_B(:,2).^2+J_B(:,3).^2+J_B(:,4).^2)];",ic); hold on;
% % % c_eval("irf_plot([Jtotal_B(:,1) Jtotal_B(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
c_eval("irf_plot([J_B(:,1) J_B(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min([min(J_B(:,2)) min(J_B(:,3)) min(J_B(:,4))])-1 max([max(J_B(:,2)) max(J_B(:,3)) max(J_B(:,4))])+1],'fontsize',9);",ic);
% set(gca,'Ylim',[-8000 8000]);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'J_x','J_y','J_z'},[0.97 0.92]);
ylabel('J [nA/m^2]','fontsize',8);
i=i+1;
%% JdotE
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([JdotE_B(:,1) Ne1(:,2) JdotE_B(:,2)],'yy',1);
% % % % irf_plot(JdotE_B,['yy',2], 'color','k', 'Linewidth',0.75);
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([JdotE?(:,1) JdotE?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
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
%% Jplasma_lmn plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([lmnJ?(:,1) lmnJ?(:,4)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([lmnJ?(:,1) lmnJ?(:,3)], 'color','g', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([lmnJ?(:,1) lmnJ?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
% % % % irf_plot([Vit1(:,1) Vit1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([Vexbt1(:,1) Vexbt1(:,2)*1e-3], 'color',[1 0 1], 'Linewidth',0.75); hold on;
% % % c_eval("irf_plot([J?(:,1) J?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[fix(min([min(lmnJ?(:,2)) min(lmnJ?(:,3)) min(lmnJ?(:,4))])/10)*10-10 fix(max([max(lmnJ?(:,2)) max(lmnJ?(:,3)) max(lmnJ?(:,4))])/10)*10+10],'fontsize',9);",ic);
% % % %set(gca,'Ylim',[-200 400], 'ytick',[-100 0 300]);
% % % % irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',12);
% % % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0];[1 0 1]]);
% % % % irf_legend(gca,{'Vi_N','Vi_M','Vi_L','|Vi|','|Vexb|'},[0.1 0.12]);
% % % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % % irf_legend(gca,{'J_L','J_M','J_N'},[0.97 0.92]);
% % % ylabel('J [nA/m^2]','fontsize',12);
% % % i=i+1;
%% Te plot
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Te_para?(:,1) (Te_para?(:,2)+2*Te_perp?(:,2))/3], 'color','k', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Te_para?(:,1) Te_para?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Te_perp?(:,1) Te_perp?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[fix(min([min(Te_para?(:,2)) min(Te_para?(:,2)) min(Te_perp?(:,2))])/10)*10-10 fix(max([max(Te_para?(:,2)) max(Te_para?(:,2)) max(Te_perp?(:,2))])/10)*10+10],'fontsize',9);",ic);
% set(gca,'Ylim',[500 2500]);
% irf_legend(gca,'e',[0.99 0.98],'color','k','fontsize',12);
set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
irf_legend(gca,{'Te','T_/_/','T_⊥'},[0.97 0.92]);
ylabel('Te [eV]','fontsize',8);
i=i+1;
%% Ti plot
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Ti_para?(:,1) (Ti_para?(:,2)+2*Ti_perp?(:,2))/3], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Ti_para?(:,1) Ti_para?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Ti_perp?(:,1) Ti_perp?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % ylabel('Ti [eV]','fontsize',8);
% % % c_eval("set(gca,'Ylim',[fix(min([min(Ti_para?(:,2)) min(Ti_para?(:,2)) min(Ti_perp?(:,2))])/10)*10-10 fix(max([max(Ti_para?(:,2)) max(Ti_para?(:,2)) max(Ti_perp?(:,2))])/10)*10+10]);",ic);
% % % %set(gca,'Ylim',[-100 300]);
% % % % irf_legend(gca,'e',[0.99 0.98],'color','k','fontsize',12);
% % % set(gca,'ColorOrder',[[0 0 0];[0 0 1];[1 0 0]]);
% % % irf_legend(gca,{'Ti','Tipara','Tiperp'},[0.97 0.92]);
% % % i=i+1;
%% Pressure
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([Pb?(:,1) Pb?(:,2)], 'color','b', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Pt?(:,1) Pt?(:,2)], 'color','r', 'Linewidth',0.75);",ic); hold on;
c_eval("irf_plot([Pt?(:,1) Pt?(:,2)+Pb?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% irf_plot([Pthe_para(:,1) Pthe_para(:,2)], 'color','g', 'Linewidth',0.75);hold on;
% irf_plot([Pthe_perp(:,1) Pthe_perp(:,2)], 'color','y', 'Linewidth',0.75);hold on;
grid off;
% set(h(i),'yscale','log');
% set(h(i),'ytick',[0 0.25 0.5],'fontsize',9);
% set(gca,'Ylim',[0 0.5]);
c_eval("set(gca,'Ylim',[0 max(Pt?(:,2)+Pb?(:,2))+0.01]);",ic);
set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
irf_legend(gca,{'Pm','Pthe','Ptotal'},[0.97 0.92]);
% pos3=get(gca,'pos');
% set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('P [nPa]','fontsize',8)
i=i+1; 
%% K
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Pte(:,1) K(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % c_eval("irf_plot([Pte(:,1) 0*K(:,2)],'k--', 'Linewidth',0.75);",ic); hold on;
% % % 
% % % 
% % % grid off;
% % % % set(h(i),'yscale','log');
% % % % set(h(i),'ytick',[0 0.25 0.5],'fontsize',9);
% % % % set(gca,'Ylim',[0 0.5]);
% % % c_eval("set(gca,'Ylim',[min(K(:,2))-0.01 max(K(:,2))+0.01]);",ic);
% % % % set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
% % % % irf_legend(gca,{'Pm','Pthe','Ptotal'},[0.97 0.92]);
% % % % pos3=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('K','fontsize',8)
% % % i=i+1;
%% Pte
% % % h(i)=irf_subplot(n,1,-i);
% % % c_eval("irf_plot([Pte?(:,1) Pte?(:,2)], 'color','k', 'Linewidth',0.75);",ic); hold on;
% % % 
% % % grid off;
% % % % set(h(i),'yscale','log');
% % % % set(h(i),'ytick',[0 0.25 0.5],'fontsize',9);
% % % % set(gca,'Ylim',[0 0.5]);
% % % c_eval("set(gca,'Ylim',[0 max(Pte?(:,2))+0.05]);",ic);
% % % % set(gca,'ColorOrder',[[0 0 1];[1 0 0];[0 0 0]]);
% % % % irf_legend(gca,{'Pm','Pthe','Ptotal'},[0.97 0.92]);
% % % % pos3=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % ylabel('Pte [nPa]','fontsize',10)
% % % i=i+1; 
%% β/S
% % % h(i)=irf_subplot(n,1,-i);
% % % % c_eval("irf_plot([beta?(:,1) beta?(:,2)], 'color','k', 'Linewidth',0.75);",1); hold on;
% % % c_eval("irf_plot([beta?(:,1) S?], 'color','b', 'Linewidth',0.75); hold on;",1);
% % % grid off;
% % % % set(h(i),'yscale','log');
% % % % set(h(i),'ytick',[1 2 3 4],'fontsize',9);
% % % set(gca,'Ylim',[0 1]);
% % % % set(gca,'Ylim',[round(min(S1))-1 round(max(S1))]);
% % % set(gca,'ColorOrder',[0 0 1]);
% % % % irf_legend(gca,{'/beta'},[0.97 0.92]);
% % % % pos3=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 1 0]]);
% % % %irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
% % % % ylabel('\beta','fontsize',12)
% % % ylabel('Entropy','fontsize',12)
% % % i=i+1; 
%% plot low e pad
% % % %     %0-200eV
h(i)=irf_subplot(n,1,-i);
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(i),jet)
specrec_p_elow=struct('t',irf_time(energy_low.DEPEND_0.data,'ttns>epoch'));
specrec_p_elow.f=transpose(energy_low.DEPEND_1.data(1,1:30));%energy levels
specrec_p_elow.p=energy_low.data;%data matrix
specrec_p_elow.f_label='';
specrec_p_elow.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb6]=irf_spectrogram(h(i),specrec_p_elow);
ylabel('PA low','fontsize',8)
% set(gca,'yscale','log');
set(h(i),'ytick',[0 90 180]);
% caxis(gca,[7 7.7]);
%irf_legend(h(i),'g',[0.99 0.98],'color','w','fontsize',12);
poscbar6=get(hcb6,'pos');
poscbar6(3)=poscbar6(3)*0.5;
set(hcb6,'pos',poscbar6);
i=i+1;
%% plot mid e pad
% % % %     %200-2000eV
h(i)=irf_subplot(n,1,-i);
%h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(i),jet)

specrec_p_emid=struct('t',irf_time(energy_mid.DEPEND_0.data,'ttns>epoch'));
specrec_p_emid.f=transpose(energy_mid.DEPEND_1.data(1,1:30));%energy levels
specrec_p_emid.p=energy_mid.data;%data matrix
specrec_p_emid.f_label='';
specrec_p_emid.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb7]=irf_spectrogram(h(i),specrec_p_emid);
ylabel('PA mid','fontsize',8)
%set(gca,'yscale','log');
set(h(i),'ytick',[0 90 180]);
caxis(gca,[5.8 6.8]);
%irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
poscbar7=get(hcb7,'pos');
poscbar7(3)=poscbar7(3)*0.5;
set(hcb7,'pos',poscbar7);
i=i+1;
%% plot high e pad
% % % %2k-30keV
h(i)=irf_subplot(n,1,-i);
%h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
colormap(h(i),jet)

specrec_p_ehigh=struct('t',irf_time(energy_high.DEPEND_0.data,'ttns>epoch'));
specrec_p_ehigh.f=transpose(energy_high.DEPEND_1.data(1,1:30));%energy levels
specrec_p_ehigh.p=energy_high.data;%data matrix
specrec_p_ehigh.f_label='';
specrec_p_ehigh.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb6]=irf_spectrogram(h(i),specrec_p_ehigh);
ylabel('PA high','fontsize',8)

set(h(i),'ytick',[0 90 180]);
caxis(gca,[7 7.8]);
%irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
poscbar6=get(hcb6,'pos');
poscbar6(3)=poscbar6(3)*0.5;
set(hcb6,'pos',poscbar6);
i=i+1;

%% plot high e pad2
% h(i)=irf_subplot(n,1,-i);
% %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(h(i),jet)
% 
% specrec_p_ehigh2=struct('t',irf_time(energy_high2.DEPEND_0.data,'ttns>epoch'));
% specrec_p_ehigh2.f=transpose(energy_high2.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_ehigh2.p=energy_high2.data;%data matrix
% specrec_p_ehigh2.f_label='';
% specrec_p_ehigh2.p_label={' ','keV/(cm^2 s sr keV)'};
% [h(i), hcb7]=irf_spectrogram(h(i),specrec_p_ehigh2);
% ylabel('PA high','fontsize',10)
% 
% set(h(i),'ytick',[0 90 180]);
%  caxis(gca,[6.3 7.1]);
% %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
% poscbar7=get(hcb7,'pos');
% poscbar7(3)=poscbar7(3)*0.5;
% set(hcb7,'pos',poscbar7);
% i=i+1;

%% plot high e pad3
% h(i)=irf_subplot(n,1,-i);
% %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(h(i),jet)
% 
% specrec_p_ehigh3=struct('t',irf_time(energy_high3.DEPEND_0.data,'ttns>epoch'));
% specrec_p_ehigh3.f=transpose(energy_high3.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_ehigh3.p=energy_high3.data;%data matrix
% specrec_p_ehigh3.f_label='';
% specrec_p_ehigh3.p_label={' ','keV/(cm^2 s sr keV)'};
% [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_ehigh3);
% ylabel('PA high','fontsize',10)
% 
% set(h(i),'ytick',[0 90 180]);
%  caxis(gca,[6.3 7.1]);
% %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
% poscbar8=get(hcb8,'pos');
% poscbar8(3)=poscbar8(3)*0.5;
% set(hcb8,'pos',poscbar8);
% i=i+1;

%% plot high e pad4
% h(i)=irf_subplot(n,1,-i);
% %h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% colormap(h(i),jet)
% 
% specrec_p_ehigh4=struct('t',irf_time(energy_high4.DEPEND_0.data,'ttns>epoch'));
% specrec_p_ehigh4.f=transpose(energy_high4.DEPEND_1.data(1,1:30));%energy levels
% specrec_p_ehigh4.p=energy_high4.data;%data matrix
% specrec_p_ehigh4.f_label='';
% specrec_p_ehigh4.p_label={' ','keV/(cm^2 s sr keV)'};
% [h(i), hcb9]=irf_spectrogram(h(i),specrec_p_ehigh4);
% ylabel('PA high','fontsize',10)
% 
% set(h(i),'ytick',[0 90 180]);
%  caxis(gca,[6.3 7.1]);
% %irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
% poscbar9=get(hcb9,'pos');
% poscbar9(3)=poscbar9(3)*0.5;
% set(hcb9,'pos',poscbar9);
% i=i+1;

%% pad pingjun plot
% h(i)=irf_subplot(n,1,-i);
% % specrec_p_ehigh1.p(:,32)=specrec_p_ehigh1.t(:,1)
% % average_perp2_data=(specrec_p_ehigh2.p(:,14)+specrec_p_ehigh2.p(:,15)+specrec_p_ehigh2.p(:,16)+specrec_p_ehigh2.p(:,17))/4;
% % specrec_p_ehigh3.p(:,31)=(specrec_p_ehigh3.p(:,14)+specrec_p_ehigh3.p(:,15)+specrec_p_ehigh3.p(:,16)+specrec_p_ehigh3.p(:,17))/4;
% % specrec_p_ehigh4.p(:,31)=(specrec_p_ehigh4.p(:,14)+specrec_p_ehigh4.p(:,15)+specrec_p_ehigh4.p(:,16)+specrec_p_ehigh4.p(:,17))/4;
% average_perp1_data=(specrec_p_ehigh1.p(:,14)+specrec_p_ehigh1.p(:,15)+specrec_p_ehigh1.p(:,16)+specrec_p_ehigh1.p(:,17))/4;
% average_perp1=[irf_time(energy_high1.DEPEND_0.data,'ttns>epoch') average_perp1_data];
% % average_perp1=irf_tlim(average_perp1,[iso2epoch('2017-06-11T01:59:34Z') iso2epoch('2017-06-11T01:59:46Z')]);
% 
% irf_plot([average_perp1(:,1) average_perp1(:,2)], 'color','k', 'Linewidth',0.75);
% % irf_plot([specrec_p_ehigh1.t(:,1) specrec_p_ehigh1.p(:,31)], 'color','k', 'Linewidth',0.75); hold off;
% % plot(specrec_p_ehigh1.p(:,32),specrec_p_ehigh1.p(:,31),'color','k','Linewidth',0.75);hold on;
% % irf_plot([specrec_p_ehigh2.p(:,32) specrec_p_ehigh2.p(:,31)], 'color','b', 'Linewidth',0.75); hold on;
% % irf_plot([specrec_p_ehigh3.t(:,1) specrec_p_ehigh3.p(:,31)], 'color','g', 'Linewidth',0.75); hold on;
% % irf_plot([specrec_p_ehigh3.t(:,1) specrec_p_ehigh4.p(:,31)], 'color','r', 'Linewidth',0.75); hold on;
% % B(:,:)=sqrt(B1(:,2).^2+B1(:,3).^2+B1(:,4).^2); hold on;
% % irf_plot([Bt1(:,1) Bt1], 'color','k', 'Linewidth',0.75); hold on;
% % irf_plot([B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.75); hold off;
% grid off;
% ylabel('f [keV/(cm^2 s sr keV]','fontsize',10);
% % set(gca,'Ylim',[fix(min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])/10)*10-10 fix(max(Bt1(:,2))/10)*10+10]);
% % set(h(3),'Ylim',[3 8], 'ytick',[5 6]);
% % set(h(9),'Yscale','log');
% % pos1=get(gca,'pos');
% % set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
% % irf_legend(gca,{'B_x','B_y','B_z','B_t'},[0.05 0.92]);
% % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% i=i+1;
%% plot e energy spectrom
h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)

specrec_p_e=struct('t',irf_time(energy_e.DEPEND_0.data,'ttns>epoch'));
specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
specrec_p_e.p=energy_e.data;%data matrix
specrec_p_e.f_label='';
specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
grid off;
set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ee(ev)','fontsize',8)
set(gca,'Ylim',[1e1 3e4]);
caxis(gca,[6.4 7.4])

% irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
poscbar8=get(hcb8,'pos');
poscbar8(3)=poscbar8(3)*0.5;
%poscbar6(1)=poscbar6(1)*0.5;
set(hcb8,'pos',poscbar8);
i=i+1;
%% plot ION energy spectrom
h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)

c_eval("specrec_p_i=struct('t',irf_time(energy_i?.DEPEND_0.data,'ttns>epoch'));",ic);
c_eval("specrec_p_i.f=transpose(energy_i?.DEPEND_1.data(1,1:32));",ic);%energy levels
c_eval("specrec_p_i.p=energy_i?.data;",ic);%data matrix
specrec_p_i.f_label='';
specrec_p_i.p_label={' ','keV/(cm^2 s sr keV)'};
[h(i), hcb7]=irf_spectrogram(h(i),specrec_p_i);
% hold on;
% irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
grid off;
set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
ylabel('Ei(ev)','fontsize',8)
set(gca,'Ylim',[1e1 3e4]);
caxis(gca,[5 6])

% irf_legend(gca,'f',[0.99 0.98],'colo6.4r','k','fontsize',12);
poscbar7=get(hcb7,'pos');
poscbar7(3)=poscbar7(3)*0.5;
set(hcb7,'pos',poscbar7);
i=i+1;
%% plot e energy spectrom
% h(i)=irf_subplot(n,1,-i);
% colormap(h(i),jet)
% 
% specrec_p_e=struct('t',irf_time(energy_e.DEPEND_0.data,'ttns>epoch'));
% specrec_p_e.f=transpose(energy_e.DEPEND_1.data(1,1:32));%energy levels
% specrec_p_e.p=energy_e.data;%data matrix
% specrec_p_e.f_label='';
% specrec_p_e.p_label={' ','keV/(cm^2 s sr keV)'};
% [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);
% % hold on;
% % irf_plot([Energy_exb1(:,1) Energy_exb1(:,2)], 'color','k', 'Linewidth',0.75); hold off;
% grid off;
% set(h(i),'yscale','log');
% set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);
% ylabel('Ee(ev)','fontsize',12)
% set(gca,'Ylim',[1e1 3e4]);
% caxis(gca,[6.4 7.4])
% % irf_legend(gca,'f',[0.99 0.98],'color','k','fontsize',12);
% poscbar8=get(hcb8,'pos');
% poscbar8(3)=poscbar8(3)*0.5;
% %poscbar6(1)=poscbar6(1)*0.5;
% set(hcb8,'pos',poscbar8);
% i=i+1;

%% plot waves
% % % c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
% % % c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
% % % c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
% % % % Bscm=Bscm{1};            %Bscm??cell
% % % c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
% % % magB = Bxyz.abs;
% % % 
% % % %gse2gsm
% % % c_eval(['Egse=irf.ts2mat(Exyz);'],ic);
% % % c_eval(['Egsm=irf_gse2gsm(Egse);'],ic);
% % % Exyz.data=Egsm(:,2:4);
% % % try
% % % c_eval(['Bscmgse=irf.ts2mat(Bscm);'],ic);
% % % c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
% % % Bscm.data=Bscmgsm(:,2:4);
% % % 
% % % % Rotate E and B into field-aligned coordinates
% % % Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
% % % Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);
% % % % Bandpass filter E and B waveforms
% % % dfE = 1/median(diff(Exyz.time.epochUnix));
% % % dfB = 1/median(diff(Bscm.time.epochUnix));
% % % Exyzfachf = Exyzfac.filt(10,0,dfE,5);
% % % Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
% % % Bscmfachf = Bscmfac.filt(10,0,dfB,5);
% % % catch
% % % % % 当Bscm发生bug时其会变为{1,2}的cell，点进去发现两部分是一样的，有时候重启matlab会好使有时候不好使就用下面这部分（到wave transforms之前）
% % % c_eval(['Bscmgse=irf.ts2mat(Bscm{1,1});'],ic);
% % % c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
% % % Bscm{1,1}.data=Bscmgsm(:,2:4);
% % % 
% % % % Rotate E and B into field-aligned coordinates
% % % Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
% % % Bscmfac = irf_convert_fac(Bscm{1,1},Bxyz,[1 0 0]);
% % % % Bandpass filter E and B waveforms
% % % dfE = 1/median(diff(Exyz.time.epochUnix));
% % % dfB = 1/median(diff(Bscm{1,1}.time.epochUnix));
% % % Exyzfachf = Exyzfac.filt(10,0,dfE,5);
% % % Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
% % % Bscmfachf = Bscmfac.filt(10,0,dfB,5);
% % % end
% % % 
% % % % Wavelet transforms
% % % nf = 100;
% % % Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 4000]);
% % % Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 50000]);
% % % Bwavelet = irf_wavelet(Bscmfac,'nf',nf,'f',[5 4000]);
% % % 
% % % %compress wavelet transform data 10 point average
% % % nc = 20;
% % % idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
% % % Ewavelettimes = Ewavelet.t(idx);
% % % Ewaveletx = zeros(length(idx),nf);
% % % Ewavelety = zeros(length(idx),nf);
% % % Ewaveletz = zeros(length(idx),nf);
% % % for ii = [1:length(idx)];
% % %         Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % % end
% % % specperpE=struct('t',Ewavelettimes);
% % % specperpE.f=Ewavelet.f;
% % % specperpE.p=Ewaveletx+Ewavelety;
% % % specperpE.f_label='';
% % % specperpE.p_label={'log_{10} E_{\perp}^2','mV^2 m^{-2} Hz^{-1}'};
% % % 
% % % specparE=struct('t',Ewavelettimes);
% % % specparE.f=Ewavelet.f;
% % % specparE.p=Ewaveletz;
% % % specparE.f_label='';
% % % specparE.p_label={'log_{10} E_{||}^2','mV^2 m^{-2} Hz^{-1}'};
% % % 
% % % specE=struct('t',Ewavelettimes);
% % % specE.f=Ewavelet.f;
% % % specE.p=Ewaveletx+Ewavelety+Ewaveletz;
% % % specE.f_label='';
% % % specE.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};
% % % 
% % % 
% % % idx = [nc/2:nc:length(Bwavelet.t)-nc/2];
% % % Bwavelettimes = Bwavelet.t(idx);
% % % Bwaveletx = zeros(length(idx),nf);
% % % Bwavelety = zeros(length(idx),nf);
% % % Bwaveletz = zeros(length(idx),nf);
% % % for ii = [1:length(idx)];
% % %         Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % %         Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
% % % end
% % % specB=struct('t',Bwavelettimes);
% % % specB.f=Bwavelet.f;
% % % specB.p=Bwaveletx+Bwavelety+Bwaveletz;
% % % specB.f_label='';
% % % specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};
% % % 
% % % 
% % % % Compute characteristic frequencies
% % % Units=irf_units; % read in standard units
% % % Me=Units.me;
% % % Mp=Units.mp;
% % % e=Units.e;
% % % epso=Units.eps0;
% % % mu0=Units.mu0;
% % % Mp_Me = Mp/Me;
% % % B_SI=magB.data*1e-9;
% % % Wpe = sqrt(ne.resample(Bxyz).data*1e6*e^2/Me/epso);
% % % Wce = e*B_SI/Me;
% % % Wpp = sqrt(ne.resample(Bxyz).data*1e6*e^2/Mp/epso);
% % % Fce = Wce/2/pi;
% % % Fce01=Fce*0.1;
% % % Fce05=Fce*0.5;
% % % Fpe = Wpe/2/pi;
% % % Fcp = Fce/Mp_Me;
% % % Fpp = Wpp/2/pi;
% % % Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);
% % % Fpe = irf.ts_scalar(magB.time,Fpe);
% % % Fce = irf.ts_scalar(magB.time,Fce);
% % % Flh = irf.ts_scalar(magB.time,Flh);
% % % Fpp = irf.ts_scalar(magB.time,Fpp);
% % % Fce01=irf.ts_scalar(magB.time,Fce01);
% % % Fce05=irf.ts_scalar(magB.time,Fce05);
% % % 
% % % h(i)=irf_subplot(n,1,-i);
% % % colormap(h(i),jet)
% % % [a8,b8]=irf_spectrogram(h(i),specE,'log');
% % % 
% % % hold(h(i),'on');
% % % irf_plot(h(i),Fpe,'color','k','LineWidth',1)
% % % irf_plot(h(i),Flh,'color','k','LineWidth',1)
% % % irf_plot(h(i),Fce,'color','r','LineWidth',1)
% % % irf_plot(h(i),Fce01,'color','w','LineWidth',1)
% % % irf_plot(h(i),Fce05,'color','c','LineWidth',1)
% % % hold(h(i),'off');
% % % 
% % % % irf_legend(h(i),'(h)',[0.99 0.97],'color','w','fontsize',12)
% % % caxis(h(i),[-12 0]);
% % % set(h(i),'yscale','log');
% % % set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
% % % ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
% % % ylabel(b8,{'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',10);
% % % grid(h(i),'off');
% % % poscbar8=get(b8,'pos');
% % % poscbar8(3)=poscbar8(3)*0.5;
% % % set(b8,'pos',poscbar8);
% % % i=i+1;
% % % 
% % % h(i)=irf_subplot(n,1,-i);
% % % colormap(h(i),jet)
% % % [a9,b9]=irf_spectrogram(h(i),specB,'log');
% % % %[h(i),b9]=irf_spectrogram(h(i),specB,'log');
% % % 
% % % hold(h(i),'on');
% % % irf_plot(h(i),Flh,'color','k','LineWidth',1)
% % % irf_plot(h(i),Fce,'color','r','LineWidth',1)
% % % irf_plot(h(i),Fce01,'color','w','LineWidth',1)
% % % irf_plot(h(i),Fce05,'color','c','LineWidth',1)
% % % hold(h(i),'off');
% % % 
% % % % irf_legend(h(i),'(zhaomingjie)',[0.99 0.97],'color','w','fontsize',12)
% % % caxis(h(i),[-12 0]);
% % % set(h(i),'yscale','log');
% % % set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
% % % ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
% % % ylabel(b9,{'log_{10} B^2','nT^2 Hz^{-1}'},'fontsize',10);
% % % grid(h(i),'off');
% % % poscbar9=get(b9,'pos');
% % % poscbar9(3)=poscbar9(3)*0.5;
% % % set(b9,'pos',poscbar9);
% % % i=i+1;
% % % % 
%   set(h(1:n),'fontsize',8);
% %   irf_zoom(tint,'x',h(1:n));

% c_eval("irf_pl_mark(h(?),B1(tempidx_B1,1),'k');",1:n)
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
%   irf_pl_marak(h(1:8),[iso2epoch('2015-10-16T13:04:30.209Z')],'k');
%  add_position(gca,gseR1), xlabel(gca,'')
%  irf_zoom(tintlmn,'x',h(4:7))

%%  出图保存部分
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
% cd  C:\Matlab\bin\新建文件夹\fwd\
% rmdir(TempDir,'s'); 
% figname = 'overview';
%     print(gcf, '-dpdf', [figname '.pdf']);    