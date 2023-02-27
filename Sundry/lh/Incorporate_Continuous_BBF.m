%% bug驱散令
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
%% prepare data & folder
load('C:\Matlab\bin\新建文件夹\fwd\Sundry\lh\AA.mat');

global ParentDir 
ParentDir = 'C:\THEMIS\'; 
TempDir = 'C:\THEMIS\temp\'; 
if ~isfolder(TempDir)
    mkdir(TempDir);
end
%%
flag_diff = find(diff(AA(:,1))>5);

for i = 210:211
% for i = 1:length(flag_diff)
%% 
if i == length(AA)
    disp('The breakpoint is at the end of data!')
    break
end
%%
clc
disp(['check point:',num2str(i),'/',num2str(length(flag_diff))])
%%
tint = [AA(flag_diff(i),1),AA(flag_diff(i)+1,1)];
ic = {'b'};
TT = epoch2iso(tint(1));
% c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','fgm',TempDir)",ic);
% c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','efi',TempDir)",ic);
c_eval("THEMISDownload(strrep(TT(1:10),'-',''),'th?','esa',TempDir)",ic);
THEMISDataMove(TempDir,ParentDir)

% c_eval("B_? = th_read_l2_change_by_fwd('th?_fgl_gsm',tint);",ic);
c_eval("Vi? = th_read_l2_change_by_fwd('th?_peir_velocity_gsm',tint);",ic);

if isempty(find(abs(Vib(:,2))<50,1))
    flag_diff(i) = nan;
end
end

flag_diff(isnan(flag_diff)) = [];