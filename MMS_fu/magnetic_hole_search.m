%gsm
clear;
clc;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       南无电子阿弥陀佛驱散仿生bug
%                                _ooOoo_
%                               o8888888o
%                               88" . "88
%                               (| -_- |)
%                               O\  =  /O
%                            ____/`---'\____
%                          .'  \\|     |//  `.
%                         /  \\|||  :  |||//  \
%                        /  _||||| -:- |||||-  \
%                        |   | \\\  -  /// |   |
%                        | \_|  ''\-/''  |   |
%                        \  .-\__  `-`  ___/-. /
%                      ___`. .'  /-.-\  `. . __
%                   ."" '<  `.___\_<|>_/___.'  >'"".
%                  | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%                  \  \ `-.   \_ __\ /__ _/   .-` /  /
%             ======`-.____`-.___\_____/___.-`____.-'======
% 	                   `=-='
%                 天地玄宗，万气本根。广修亿劫，证吾神通。
%                 三界内外，惟道独尊。体有金光，覆映吾身。
%                 视之不见，听之不闻。包罗天地，养育群生。
%                 受持万遍，身有光明。三界侍卫，五帝司迎。
%                 万神朝礼，役使雷霆。鬼妖丧胆，精怪忘形。
%                 内有霹雳，雷神隐名。洞慧交彻，五炁腾腾。
%                金光速现，覆护真人。急急如律令，bug全去除！
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Date = '2018-06-01/2018-09-01';
% Date = '2017-04-30/2017-05-01';
% Date = '2017-01-01/2021-01-01';
% Date = '2022-06-09/2022-07-02';

splitDate = regexp(Date,'/','split');
ic = 1;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames3 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %To get loaction
filenames = filenames1;

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
NameTags = unique(cellfun(@cellstr,NameTags));
FileGroups = cell(1,length(NameTags)); 
for j = 1:length(NameTags)
    FileGroups{j} = filenames(contains(filenames,NameTags{j}));
end
FileGroups = cellfun(@cellstr,FileGroups,'UniformOutput',false);%按时间分类整理后的文件名组

global OutputDir ParentDir
ParentDir = '/Volumes/FWD-T7Disk/MMS/'; 
TempDir = [ParentDir, 'temp/'];mkdir(TempDir);
%The dir of "SDCFilesDownload" to "datamove" must be the ParentDir!
OutputDir = [ParentDir,'MH_Search/',splitDate{1},'To',splitDate{2},'/'];
if ~isfolder([OutputDir,'OverviewFig/'])
    mkdir([OutputDir,'OverviewFig/']);
end
%%
units = irf_units;
NameTags{end+1} = ['_' strrep(splitDate{2},'-','') '235959_v'];
for TDT = 1:length(NameTags)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
tempDir = [OutputDir,NameTags{TDT}(2:end-2),'/'];
clc
fprintf(['Current Date:',NameTags{TDT}(2:end-2),'\n'])

tempDate = [NameTags{TDT}(2:5),'-',NameTags{TDT}(6:7),'-',NameTags{TDT}(8:9),'T',...
    NameTags{TDT}(10:11),':',NameTags{TDT}(12:13),':',NameTags{TDT}(14:15),'.000Z/',...
    NameTags{TDT+1}(2:5),'-',NameTags{TDT+1}(6:7),'-',NameTags{TDT+1}(8:9),'T',...
    NameTags{TDT+1}(10:11),':',NameTags{TDT+1}(12:13),':',NameTags{TDT+1}(14:15),'.000Z'];
tempTint=irf.tint(tempDate);

%% Poincare Index  
% srvyIdx = find(contains(filenames_srvy,NameTags{TDT}(2:9))==1);
flag = 0;flag_t = 0;
try
    SDCFilesDownload(FileGroups{TDT},tempDir);
    SDCDataMove(tempDir,ParentDir); mms.db_init('local_file_db',ParentDir);
%% Data Load
    B1_ts=mms.get_data('B_gsm_brst',tempTint,1);%先导入一个文件看看文件中包含的时间段
    tint = irf.tint(B1_ts.time.epoch(1),B1_ts.time.epoch(end));
    c_eval("B?_ts=mms.get_data('B_gsm_brst',tint,?);"),ic;
    c_eval('B?_gsm = irf.ts2mat(B?_ts);',ic); 
    c_eval('B? = irf_abs(B?_gsm);',ic);

    Pos = mms.get_data('R_gsm',tint);
    c_eval('R? = Pos.gsmR?;',ic)
    c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];',ic)

%% Data Filter
    RE = units.RE;
    c_eval('R_RE = mean(R?(1,2:4)/RE*1000,1);',ic) 
    if R_RE(1) <= -10 && abs(R_RE(2)) <= 15 && abs(R_RE(3)) <= 10
%% MH Choose
        t_range = 20; % 30s for ambient magnetic field
        t_start = t_range/2*128+1; t_end = size(B1,1)-128*t_range/2-1;
        for t = t_start:5:t_end
            clc;
            disp(['Current Date:' NameTags{TDT}(2:end-2)])
            disp(['Current Calculate:',repmat('■',1,round(10*(t-t_range/2*128)/length(t_start:t_end)))...
                ,repmat('□',1,10-round(10*(t-t_range/2*128)/length(t_start:t_end))),' ᕙ(`▿´)ᕗ'...
                ,num2str(t-t_range/2*128),'/',num2str(length(t_start:t_end))]);

            mag_amb = mean(B1(t-128*t_range/2:t+128*t_range/2,:),1);

            if abs(mag_amb(4)) >= abs(mag_amb(3)) && abs(mag_amb(4) >= abs(mag_amb(2))) %偶极场
                mag_fore = mean(B1(t-128*t_range/2:t-1,:),1); mag_post = mean(B1(t+1:t+128*t_range/2,:),1);
                if flag == 1, break;end
                if B1(t,5) <= 0.75 * mag_amb(5) && acosd(dot(mag_fore(2:4),mag_post(2:4))...
                        / (norm(mag_fore(5)) * norm(mag_post(5)))) <= 15 % B < 0.75*Bambient & <Bfore,Bpost>  <=  15
                    flag = 1;
                    flag_t = B1(t);
                    writematrix(['MH find at: ',datestr(datenum(1970,1,1,0,0,0)+flag_t/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [OutputDir,'MH_case.txt'],'WriteMode','append','Encoding','UTF-8')
                else
                    DeleteFolder(tempDir)
                    continue
                end
            else
                DeleteFolder(tempDir)
                continue
            end
        end
    else
        DeleteFolder(tempDir)
        continue
    end  
catch
    writematrix(['Data Read in ',NameTags{TDT}(2:end-2),' Error'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    DeleteFolder(tempDir)
    continue
end

%% Plot
if flag == 1
fprintf(['Case find at: ',datestr(datenum(1970,1,1,0,0,0)+flag_t/86400,'yyyymmdd HH:MM:SS.FFF'), ' φ(≧ω≦*)♪\n'])
fprintf('Ploting...')
try
%% Data Download
TT1 = datestr(datenum(1970,1,1,0,0,0)+B1(t-128*t_range/2,1)/86400,'yyyy-mm-ddTHH:MM:SS.FFFZ');
TT2 = datestr(datenum(1970,1,1,0,0,0)+B1(t+128*t_range/2,1)/86400,'yyyy-mm-ddTHH:MM:SS.FFFZ');
TT = [TT1,'/',TT2];
tint = irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = [filenames1,filenames2,filenames3,filenames4];

expr1 = '_\d+\_v';
NamesTags = regexp(filenames,expr1,'match');
NamesTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NamesTags)),'UniformOutput',false);

TTlist = strjoin(regexp(TT,'\d+','match'),'');
i = 1;flag = 0;%若flag=0，说明整段时间都在第i-1个Tag里
while str2double(TTlist(18:31)) > NamesTags{i}  % 如果时间段刚好仅在某天的最后一个文件里会出bug，可以把时间往前调1ms
    if str2double(TTlist(1:14)) < NamesTags{i}
        flag=1; break  %若flag=1，说明时间段的开始在第i-1个Tag里，结束在第i个里
    else
        i=i+1;
    end
    if i > length(NamesTags)
        break
    end
end

% '/' for MacOs, '\' for Windows
if flag == 0
    tempTag = num2str(NamesTags{i-1});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag(1:4),'/',tempTag(5:6),'/',...
            tempTag(7:8),'/',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
    desmoms1 = desmoms; desmoms2 = desmoms1;
else
    if i == 1
        errordlg('时间起始处无brst数据，请检查时间范围,或使用Overview_srvydownload程序')
    end
    tempTag1 = num2str(NamesTags{i-1});
    tempTag2 = num2str(NamesTags{i});
    filenames1 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag1)));
    filenames2 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag2)));
    filenames = [filenames1,filenames2];
    desmoms1 = [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag1(1:4),'/',tempTag1(5:6),'/',...
            tempTag1(7:8),'/',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag2(1:4),'/',tempTag2(5:6),'/',...
            tempTag2(7:8),'/',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
end
SDCFilesDownload(filenames,tempDir)
%% load data
SDCDataMove(tempDir,ParentDir)
mms.db_init('local_file_db',ParentDir);

% load B
units = irf_units;
c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
%  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);
c_eval(['Bt?=irf.ts2mat(Bt?_ts);'],ic);


% load E
c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
%%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_fast_l2'',tint,?);'],ic);
c_eval(['Et?_ts=E?_ts.abs;'],ic); 
c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
c_eval(['Et?=irf.ts2mat(Et?_ts);'],ic);
c_eval(['E?_resamp=irf_resamp(E?,B?);'],ic);

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

% merge data/time from 2 cdf files
c_eval('energy_low?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_lowen_brst'',tint);',ic)
c_eval('energy_mid?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_miden_brst'',tint);',ic)
c_eval('energy_high?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_pitchangdist_highen_brst'',tint);',ic)
c_eval('energy_e?=mms.db_get_variable(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energyspectr_omni_brst'',tint);',ic)
c_eval('energy_i?=mms.db_get_variable(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_energyspectr_omni_brst'',tint);',ic)
%%
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
%% Electric field
h(i)=irf_subplot(n,1,-i);
c_eval("irf_plot([E?(:,1) E?(:,2)], 'color','b', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,3)], 'color','g', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,4)], 'color','r', 'Linewidth',0.75); ",ic);hold on;
c_eval("irf_plot([E?(:,1) E?(:,2)*0],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[-8 8], 'ytick',[-10:4:10],'fontsize',9);
% set(gca,'Ylim',[-40 50], 'ytick',[-60 -40 -20 0 20 40 60]);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',12);
c_eval("set(gca,'Ylim',[min([min(E?(:,2)) min(E?(:,3)) min(E?(:,4))])-0.5 max([max(E?(:,2)) max(E?(:,3)) max(E?(:,4))])+0.5]);",ic);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
irf_legend(gca,{'E_x','E_y','E_z'},[0.97 0.92]);
pos3=get(gca,'pos');
set(gca,'ColorOrder',[[0 1 0]]);
%irf_legend(gca,{'MMS3'},[pos3(1)+1.15*pos3(3),pos3(2)]);
ylabel('E [mV/m]','fontsize',8)
i=i+1;
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
% caxis(gca,[5.8 6.8]);
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
% caxis(gca,[7 7.8]);
%irf_legend(h(i),'h',[0.99 0.98],'color','w','fontsize',12);
poscbar6=get(hcb6,'pos');
poscbar6(3)=poscbar6(3)*0.5;
set(hcb6,'pos',poscbar6);
i=i+1;
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

%% plot waves
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('Bscm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
% Bscm=Bscm{1};            %Bscm??cell
c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
magB = Bxyz.abs;

%gse2gsm
c_eval(['Egse=irf.ts2mat(Exyz);'],ic);
c_eval(['Egsm=irf_gse2gsm(Egse);'],ic);
Exyz.data=Egsm(:,2:4);
try
c_eval(['Bscmgse=irf.ts2mat(Bscm);'],ic);
c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
Bscm.data=Bscmgsm(:,2:4);

% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);
% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB = 1/median(diff(Bscm.time.epochUnix));
Exyzfachf = Exyzfac.filt(10,0,dfE,5);
Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
Bscmfachf = Bscmfac.filt(10,0,dfB,5);
catch
% % 当Bscm发生bug时其会变为{1,2}的cell，点进去发现两部分是一样的，有时候重启matlab会好使有时候不好使就用下面这部分（到wave transforms之前）
c_eval(['Bscmgse=irf.ts2mat(Bscm{1,1});'],ic);
c_eval(['Bscmgsm=irf_gse2gsm(Bscmgse);'],ic);
Bscm{1,1}.data=Bscmgsm(:,2:4);

% Rotate E and B into field-aligned coordinates
Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac = irf_convert_fac(Bscm{1,1},Bxyz,[1 0 0]);
% Bandpass filter E and B waveforms
dfE = 1/median(diff(Exyz.time.epochUnix));
dfB = 1/median(diff(Bscm{1,1}.time.epochUnix));
Exyzfachf = Exyzfac.filt(10,0,dfE,5);
Exyzfaclf = Exyzfac.filt(0,10,dfE,5);
Bscmfachf = Bscmfac.filt(10,0,dfB,5);
end

% Wavelet transforms
nf = 100;
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 4000]);
Ewavelet = irf_wavelet(Exyzfac,'nf',nf,'f',[5 50000]);
Bwavelet = irf_wavelet(Bscmfac,'nf',nf,'f',[5 4000]);

%compress wavelet transform data 10 point average
nc = 20;
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
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

specE=struct('t',Ewavelettimes);
specE.f=Ewavelet.f;
specE.p=Ewaveletx+Ewavelety+Ewaveletz;
specE.f_label='';
specE.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};


idx = [nc/2:nc:length(Bwavelet.t)-nc/2];
Bwavelettimes = Bwavelet.t(idx);
Bwaveletx = zeros(length(idx),nf);
Bwavelety = zeros(length(idx),nf);
Bwaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
        Bwaveletx(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Bwavelety(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Bwaveletz(ii,:) = squeeze(irf.nanmean(Bwavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specB=struct('t',Bwavelettimes);
specB.f=Bwavelet.f;
specB.p=Bwaveletx+Bwavelety+Bwaveletz;
specB.f_label='';
specB.p_label={'log_{10} B^2','nT^2 Hz^{-1}'};


% Compute characteristic frequencies
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
Fpe = irf.ts_scalar(magB.time,Fpe);
Fce = irf.ts_scalar(magB.time,Fce);
Flh = irf.ts_scalar(magB.time,Flh);
Fpp = irf.ts_scalar(magB.time,Fpp);
Fce01=irf.ts_scalar(magB.time,Fce01);
Fce05=irf.ts_scalar(magB.time,Fce05);

h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)
[a8,b8]=irf_spectrogram(h(i),specE,'log');

hold(h(i),'on');
irf_plot(h(i),Fpe,'color','k','LineWidth',1)
irf_plot(h(i),Flh,'color','k','LineWidth',1)
irf_plot(h(i),Fce,'color','r','LineWidth',1)
irf_plot(h(i),Fce01,'color','w','LineWidth',1)
irf_plot(h(i),Fce05,'color','c','LineWidth',1)
hold(h(i),'off');

% irf_legend(h(i),'(h)',[0.99 0.97],'color','w','fontsize',12)
caxis(h(i),[-12 0]);
set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
ylabel(b8,{'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'},'fontsize',10);
grid(h(i),'off');
poscbar8=get(b8,'pos');
poscbar8(3)=poscbar8(3)*0.5;
set(b8,'pos',poscbar8);
i=i+1;

h(i)=irf_subplot(n,1,-i);
colormap(h(i),jet)
[a9,b9]=irf_spectrogram(h(i),specB,'log');
%[h(i),b9]=irf_spectrogram(h(i),specB,'log');

hold(h(i),'on');
irf_plot(h(i),Flh,'color','k','LineWidth',1)
irf_plot(h(i),Fce,'color','r','LineWidth',1)
irf_plot(h(i),Fce01,'color','w','LineWidth',1)
irf_plot(h(i),Fce05,'color','c','LineWidth',1)
hold(h(i),'off');

% irf_legend(h(i),'(zhaomingjie)',[0.99 0.97],'color','w','fontsize',12)
caxis(h(i),[-12 0]);
set(h(i),'yscale','log');
set(h(i),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(i),{'f (Hz)'},'fontsize',12,'Interpreter','tex');
ylabel(b9,{'log_{10} B^2','nT^2 Hz^{-1}'},'fontsize',10);
grid(h(i),'off');
poscbar9=get(b9,'pos');
poscbar9(3)=poscbar9(3)*0.5;
set(b9,'pos',poscbar9);
i=i+1;
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
colormap(jet)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname = [OutputDir,'OverviewFig/',datestr(datenum(1970,1,1,0,0,0)+B1(round(length(B1)/2),1)/86400,'yyyy-mm-ddTHH-MM-SS-FFF'),'-Overview'];
print(gcf, '-dpng', [figname '.png']);
clf
catch
writematrix(['Plot in ',datestr(datenum(1970,1,1,0,0,0)+B1(round(length(B1)/2),1)/86400,'yyyy-mm-ddTHH:MM:SS.FFF'),' Error, BUT Case In'],[OutputDir,'errorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
end
end
DeleteFolder(tempDir)
end
%% Delete Folder
function DeleteFolder(tempDir)
try
    global OutputDir
    cd(OutputDir)
    rmdir(tempDir,'s');
    fclose all;
catch
    fprintf(['Delete Folder ',tempDir,' False\n'])
end
end