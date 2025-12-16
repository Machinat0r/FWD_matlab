clear;clc;close all
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

global ParentDir 
ParentDir = 'Z:\SPART-WORK\Data\MMS/'; 
TempDir = 'Z:\SPART-WORK\Data\MMS/temp/';mkdir(TempDir);

% TT = '2019-01-16T04:09:50.00Z/2019-01-16T04:10:00.00Z';
TT = '2019-01-16T04:09:55.220Z/2019-01-16T04:09:56.000Z'; %no boundary, 10,78-81
% TT = '2019-01-16T04:09:55.420Z/2019-01-16T04:09:55.800Z'; %no boundary, 10,78-81
% TT = '2018-08-27T12:15:30.00Z/2018-08-27T12:15:50.00Z';

% TT = '2016-01-06T00:33:07.00Z/2016-01-06T00:33:07.200Z';
% TT = '2016-01-27T03:03:24.000Z/2016-01-27T03:03:30.000Z';
% TT = '2017-02-20T04:43:57.00Z/2017-02-20T04:43:58.00Z';
% TT = '2021-03-30T07:58:01.00Z/2021-03-30T07:58:02.00Z';
% TT = '2020-07-27T11:30:23.00Z/2020-07-27T11:30:24.00Z';
% TT = '2016-01-10T09:09:38.000Z/2016-01-10T09:09:45.000Z';
% TT = '2017-07-06T22:13:05.00Z/2017-07-06T22:13:20.00Z';
% % % TT = '2016-01-06T00:32:36.000Z/2016-01-06T00:32:37.000Z';
% TT = '2016-01-08T03:25:20.000Z/2016-01-08T03:25:30.000Z';
% TT = '2016-03-07T00:38:52.000Z/2016-03-07T00:38:52.500Z';
TTlong = '2019-01-16T04:00:00.000Z/2019-01-16T04:20:00.000Z'; %no boundary, 10,78-81

tintlong = irf.tint(TTlong);
tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

ic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
[filenames,~,~] = findFilenames(TT,filenames1,'brst',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
%% Poincare Index  
SDCDataMove(TempDir,ParentDir); 
mms.db_init('local_file_db',ParentDir);

c_eval("B?_ts=mms.get_data('B_gse_brst',tint,?);");
c_eval('B?_gse = irf.ts2mat(B?_ts);'); 
c_eval('B? = irf_abs(B?_gse);');
c_eval('B? = irf_resamp(B?,B1);',2:4);

c_eval("Blong?_ts=mms.get_data('B_gse_brst',tintlong,?);");
c_eval('Blong?_gse = irf.ts2mat(Blong?_ts);'); 
c_eval('Blong? = irf_abs(Blong?_gse);');
c_eval('Blong? = irf_resamp(Blong?,Blong1);',2:4);

% PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
PI(abs(PI)<0.5) = 0;
%% div     
Pos = mms.get_data('R_gse',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gseR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
c_eval('R? = irf_resamp(R?,B1);')
CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

% monopole_index=zeros(length(B1(:,1)),1);
gradB=c_4_grad('R?','B?','grad');
divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度

%% solve monopole
units = irf_units;

LocPoint = zeros(size(B1,1),3);
LocRes = cell(size(B1,1),1);
Q = zeros(size(B1,1),1);
resQ = cell(size(B1,1),1);
dLoc = zeros(size(B1,1),15);

RR12 = irf_abs(R1-R2); RR13 = irf_abs(R1-R3); RR14 = irf_abs(R1-R4); 
RR23 = irf_abs(R2-R3); RR24 = irf_abs(R2-R4); RR34 = irf_abs(R3-R4); 
RR_mean = (RR12(:,5) + RR13(:,5) + RR14(:,5) + RR23(:,5) + RR24(:,5) + RR34(:,5))/6;
id = nchoosek(1:6,2);

noise = estimate_noise_strength_4sc(Blong1(:,2:4), Blong2(:,2:4), Blong3(:,2:4), Blong4(:,2:4));
BGnoise = noise.sigma_scalar;

for i = 50
clc;
% res = scan_confidence_sigma([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)], ...
%     [B1(i,2:4);B2(i,2:4);B3(i,2:4);B4(i,2:4)], 1e4, 100, BGnoise);
res = scan_confidence_sigma(randn(4,3), ...
    randn(4,3), 1e4, 100, BGnoise);
end