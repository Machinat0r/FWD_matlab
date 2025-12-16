%----------written by Wending Fu, Dec.2025 in Beijing------------
clear;clc
close all
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
%% load data
% % % global ParentDir 
% % % ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
% % % DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
% % % TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
cd /Volumes/SPART-WORK/Data/Cluster/
% cd /Volumes/172.17.190.41/Data/Cluster/
% cd D:\Cluster\
ic=1:4;
% 
% Tsta='2002-02-18T08:41:00.00Z';
% Tend='2002-02-18T08:42:00.00Z';%1,114
% Tsta='2003-12-30T18:55:00.00Z';
% Tend='2003-12-30T18:56:00.00Z';%7,268
Tsta = '2004-01-17T19:31:00.000Z';
Tend = '2004-01-17T19:32:00.000Z';
% % % Tsta = '2004-01-21T08:56:00.000Z';
% % % Tend = '2004-01-21T08:58:00.000Z';
% Tsta='2002-08-21T08:18:30.000Z';
% Tend='2002-08-21T08:21:00.000Z';

% Tsta = '2000-01-01T00:00:00.000Z';
% Tend = '2020-01-01T00:00:00.000Z';


% Tsta='2004-01-01T11:10:00.00Z';
% Tend='2004-01-01T11:11:00.00Z';
% Tsta='2004-05-30T20:22:00.00Z';
% Tend='2004-05-30T20:23:00.00Z';
% Tsta='2003-12-29T13:40:00.00Z';
% Tend='2003-12-29T13:41:00.00Z';

% Tsta='2004-05-03T05:02:00.00Z';
% Tend='2004-05-03T05:03:00.00Z';


tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH
% 
% Tsta2='2002-02-18T08:41:07.000Z';
% Tend2='2002-02-18T08:41:09.000Z';
% Tend2='2002-02-18T08:41:10.000Z';
% Tsta2='2002-06-13T18:41:30.00Z';
% Tend2='2002-06-13T18:41:50.00Z';
% Tsta2='2003-12-17T01:20:05.00Z';
% Tend2='2003-12-17T01:20:09.00Z';
% Tsta2='2003-12-30T18:55:08.560Z';
% Tend2='2003-12-30T18:55:08.740Z';
% Tsta2='2003-12-30T18:55:06.800Z';
% Tend2='2003-12-30T18:55:09.800Z';
% Tsta2 = '2004-01-17T19:31:45.765Z';
% Tend2 = '2004-01-17T19:31:45.925Z';
Tsta2 = '2004-01-17T19:31:44.900Z';
Tend2 = '2004-01-17T19:31:46.500Z';
% Tsta2='2004-05-03T05:02:28.00Z';
% Tend2='2004-05-03T05:02:30.00Z';

% Tsta2='2004-01-01T11:10:03.500Z';
% Tend2='2004-01-01T11:10:04.500Z';
% % % Tsta2='2004-01-21T08:56:58.200Z';
% % % Tend2='2004-01-21T08:57:00.200Z';
% Tsta2='2004-05-30T20:22:30.00Z';
% Tend2='2004-05-30T20:22:32.00Z';
% Tsta2='2003-12-29T13:40:25.00Z';
% Tend2='2003-12-29T13:40:26.600Z';
% Tsta2='2003-12-29T13:40:26.150Z';
% Tend2='2003-12-29T13:40:26.320Z';
% % % Tsta2='2004-01-21T08:56:59.275Z';
% % % Tend2='2004-01-21T08:56:59.435Z';
% % % Tsta2='2004-01-21T08:56:57.750Z';
% % % Tend2='2004-01-21T08:57:00.500Z';
% Tsta2='2002-08-21T08:19:30.000Z';
% Tend2='2002-08-21T08:20:00.000Z';
tint2=[iso2epoch(Tsta2) iso2epoch(Tend2)];
% tint2=tint;

% caa_download(tint,'C?_CP_EFW_L3_E3D_INERT');
% caa_download(tint,'C?_CP_EFW_L?_E');
try
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
%     c_eval("caa_load('C?_CP_FGM_FULL',Tsta,Tend);",ic);
%     c_eval("caa_load('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
catch

    %    Magnetic fields
c_eval("caa_download(tint,'C?_CP_FGM_FULL')",ic);
c_eval("caa_download(tint,'C?_CP_AUX_POSGSE_1M')",ic);  % position & velocity for each sc
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
end
% caa_download(tint,'CL_SP_AUX')% position,attitude.. for all sc
% caa_download(tint,'C2_CP_FGM_FULL');
% caa_download(tint,'C4_CP_FGM_FULL');

% caa_load C  %load data from datebase form C

%background magnetic field

% % % dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
% % % varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
%smooth
% for rec = 2:4
% c_eval('B?_gse(:,rec) = smooth(B?_gse(:,rec),9);')
% end

% c_eval('B?_gse = irf_gse2gsm(B?_gse);',ic);
c_eval('Blong?=irf_abs(B?_gse);',ic);
c_eval('Blong? = irf_resamp(Blong?,Blong1);',2:4);

c_eval('R?_gse = c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');',ic);
c_eval('Rlong? = R?_gse;')
% c_eval('R? = irf_gse2gsm(R?_gse);',ic);
c_eval('Rlong? = irf_resamp(Rlong?,Blong1);')
CenterPoint = (Rlong1(:,2:4)+Rlong2(:,2:4)+Rlong3(:,2:4)+Rlong4(:,2:4))/4;
c_eval('Rlong?(:,2:4) = Rlong?(:,2:4)-CenterPoint;');

try
    [~,Tsta_id] = sort(abs(Blong1(:,1)-iso2epoch(Tsta2)));
    Tsta_id = Tsta_id(1);
    [~,Tend_id] = sort(abs(Blong1(:,1)-iso2epoch(Tend2)));
    Tend_id = Tend_id(1);
    c_eval('B? = Blong?(Tsta_id:Tend_id,:);');
    c_eval('B?_gse = B?;');
    c_eval('R? = Rlong?(Tsta_id:Tend_id,:);');
    tint = [iso2epoch(Tsta2) iso2epoch(Tend2)];
catch
    disp('no tint2,calculate all over tint')
end
%% PI
% PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
% [PI,~,~,~,~]=poincare_index_and_solid_angle(B1_gse(2:4),B2_gse(2:4),B3_gse(2:4),B4_gse(2:4));
PI(abs(PI)<0.5) = 0;
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
%% divB
gradB=c_4_grad('R?','B?(:,1:4)','grad');
divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度
PI(isnan(divB(:,2)))=0;
% tint = [iso2epoch(Tsta2) iso2epoch(Tend2)];
%% FOTE 误差
LocPoint = zeros(length(PI),3)*nan;
LocRes = cell(length(PI),1);
Q = zeros(length(PI),1)*nan;
resQ = cell(length(PI),1);

Qerror = ones(length(PI),1)*1000;
Locerror = ones(length(PI),1)*200;

PI_id = find(PI~=0)';
noise = estimate_noise_strength_4sc(Blong1(:,2:4), Blong2(:,2:4), Blong3(:,2:4), Blong4(:,2:4));
BGnoise = noise.sigma_scalar;

for i = 23
clc;
res = scan_confidence_sigma([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)], ...
    [B1(i,2:4);B2(i,2:4);B3(i,2:4);B4(i,2:4)], 1e4, 100, BGnoise);
end