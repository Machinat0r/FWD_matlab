clear;clc;close all
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

global ParentDir 
ParentDir = '/Volumes/FWD-WorkDisk/MMS/'; 
TempDir = '/Volumes/FWD-WorkDisk/MMS/temp'; mkdir(TempDir);
% TT = '2019-01-16T04:09:50.00Z/2019-01-16T04:10:00.00Z';
TT = '2019-01-16T04:09:55.220Z/2019-01-16T04:09:56.000Z'; %no boundary, 10,78-81
% TT = '2019-01-16T04:09:55.420Z/2019-01-16T04:09:55.800Z'; %no boundary, 10,78-81
% TT = '2018-08-27T12:15:30.00Z/2018-08-27T12:15:50.00Z';

% TT = '2016-01-06T00:33:07.00Z/2016-01-06T00:33:07.200Z';
% TT = '2015-09-19T07:43:28.000Z/2015-09-19T07:43:33.000Z';
% TT = '2017-02-20T04:43:57.00Z/2017-02-20T04:43:58.00Z';
% TT = '2021-03-30T07:58:01.00Z/2021-03-30T07:58:02.00Z';
% TT = '2020-07-27T11:30:23.00Z/2020-07-27T11:30:24.00Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:4;
iic = 1:4;
try
filenames1 = SDCFilenames(Date,iic,'inst','fgm','drm','brst');
% filenames_srvy = SDCFilenames(Date,iic,'inst','fgm','drm','srvy'); %为了知道坐标
filenames = filenames1;

expr1 = '_\d+\_v';
NameTags = regexp(filenames,expr1,'match');
NameTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NameTags)),'UniformOutput',false);

TTlist = regexp(TT,'\d+','match');
i = 1;flag = 0;%若flag=0，说明整段时间都在第i-1个Tag里
while i<=length(NameTags) && str2double(strjoin(TTlist(8:13),'')) > NameTags{i}
    if str2double(strjoin(TTlist(1:6),'')) < NameTags{i}
        flag=1; break  %若flag=1，说明时间段的开始在第i-1个Tag里，结束在第i个里
    else,i=i+1;
    end
end
tempTag = num2str(NameTags{i-1});
filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));

SDCFilesDownload(filenames,TempDir)
catch
    disp('Download Files Failed!')
end
% SDCFilesDownload(filenames_srvy(1:2:8),TempDir)
% % % id_flagTime = OverView_download(tint,desmoms,IC,Name,flagTime)
%% Poincare Index  
SDCDataMove(TempDir,ParentDir); mms.db_init('local_file_db',ParentDir);

c_eval("B?_ts=mms.get_data('B_gse_brst',tint,?);");
c_eval('B?_gse = irf.ts2mat(B?_ts);'); 
c_eval('B? = irf_abs(B?_gse);');
c_eval('B? = irf_resamp(B?,B1);',2:4);

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

for i = 1:length(PI)
clc;
disp(['current calculate:',num2str(i),'/',num2str(length(PI))]);
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
    'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
% if PI(i)~=0
    [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?',i,i*sign(divB(i,2)),10,1);
id = nchoosek(1:6,2);
c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
tempd = [];
c_eval('tempd = [tempd,tempd?(4)/RR_mean];',1:15);
dLoc(i,:) = tempd;

% else
%     Q(i) = nan; resQ{i} = nan; LocPoint(i,:) = [nan,nan,nan]; LocRes{i} = nan;
% end
end

%% calculate error
Qerror = zeros(length(PI),1);
Locerror = zeros(length(PI),1);

for i = 1:length(PI)
if ~isnan(resQ{i})
    Qerror(i) = abs(100*std(resQ{i})/Q(i));
else
    Qerror(i) = 1001;
end

if ~isnan(LocRes{i})
tri_a = delaunayTriangulation([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)]);
[~,volume_a] = convexHull(tri_a);
tri = delaunayTriangulation(LocRes{i});%%delaunay三角剖分
    if size(tri.Points,1)==1
        volume = 0;
    else
    [~,volume] = convexHull(tri);%%计算多面体体积
    end
Locerror(i) = 100*volume/volume_a;
else
Locerror(i) = 200;
end
end
%% Init figure 1
figure(1)
n=8;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Btotal
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2(:,1) B2(:,5)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3(:,1) B3(:,5)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,5)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% c_eval("irf_plot([B?_gse(:,1) 0*B?_gse(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[min([min(B1(:,5)) min(B2(:,5)) min(B3(:,5)) min(B4(:,5))])-10 ...
%     max([max(B1(:,5)) max(B2(:,5)) max(B3(:,5)) max(B4(:,5))])+10]);
% c_eval("set(gca,'Ylim',[min([min(B?_gse(:,2))])-100 max([max(B?_gse(:,2))])+100]);",ic);
set(gca,'Ylim',[0 20], 'ytick',[0:10:20],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('|B| [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bx
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gse(:,1) B1_gse(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gse(:,1) B2_gse(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gse(:,1) B3_gse(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gse(:,1) B4_gse(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% c_eval("irf_plot([B?_gse(:,1) 0*B?_gse(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% set(gca,'Ylim',[min([min(B1_gse(:,2)) min(B2_gse(:,2)) min(B3_gse(:,2)) min(B4_gse(:,2))])-10 ...
%     max([max(B1_gse(:,2)) max(B2_gse(:,2)) max(B3_gse(:,2)) max(B4_gse(:,2))])+10]);
% c_eval("set(gca,'Ylim',[min([min(B?_gse(:,2))])-100 max([max(B?_gse(:,2))])+100]);",ic);
set(gca,'Ylim',[-15 15], 'ytick',[-10 0 10],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bx [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% By
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gse(:,1) B1_gse(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gse(:,1) B2_gse(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gse(:,1) B3_gse(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gse(:,1) B4_gse(:,3)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% c_eval("irf_plot([B?_gse(:,1) 0*B?_gse(:,3)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?_gse(:,3))])-100 max([max(B?_gse(:,3))])+100]);",ic);
set(gca,'Ylim',[min([min(B1_gse(:,3)) min(B2_gse(:,3)) min(B3_gse(:,3)) min(B4_gse(:,3))])-10 ...
    max([max(B1_gse(:,3)) max(B2_gse(:,3)) max(B3_gse(:,3)) max(B4_gse(:,3))])+10]);
% set(gca,'Ylim',[-18 20], 'ytick',[-10:10:20],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('By [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bz
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gse(:,1) B1_gse(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gse(:,1) B2_gse(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gse(:,1) B3_gse(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gse(:,1) B4_gse(:,4)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% c_eval("irf_plot([B?_gse(:,1) 0*B?_gse(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?_gse(:,4))])-100 max([max(B?_gse(:,4))])+100]);",ic);
set(gca,'Ylim',[min([min(B1_gse(:,4)) min(B2_gse(:,4)) min(B3_gse(:,4)) min(B4_gse(:,4))])-10 ...
    max([max(B1_gse(:,4)) max(B2_gse(:,4)) max(B3_gse(:,4)) max(B4_gse(:,4))])+15]);
% set(gca,'Ylim',[-20 15], 'ytick',[-20:20:10],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bz [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% η
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([err_4C(:,1) err_4C(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % % % c_eval("set(gca,'Ylim',[min([min(B?_gse(:,4))])-3 max([max(B?_gse(:,4))])+3]);",ic);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % c_eval("set(gca,'Ylim',[min(err_4C(:,2))-0.1 max(err_4C(:,2))+0.1]);",ic);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('η','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% ξ
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([eigVal_err_v2(:,1) eigVal_err_v2(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(eigVal_err_v2(:,2))-0.1 max(eigVal_err_v2(:,2))+0.1]);",ic);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('ξ','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% PI
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) PI], 'color','k', 'Linewidth',0.75); hold on;

c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('PI ','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% divB
h(i)=irf_subplot(n,1,-i);
irf_plot([divB(:,1) divB(:,2)], 'color','k', 'Linewidth',0.75); hold on;

c_eval("irf_plot([divB(:,1) 0*divB(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min(divB(:,2))-0.05 max(divB(:,2))+0.05]);",ic);
set(gca,'Ylim',[-0.5 1.5], 'ytick',[-0.5 0 0.5 1 1.5],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('divB [nT/km]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% monopole_index
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) monopole_index], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
% % % set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('MI','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Location Error
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) Locerror], 'color','k', 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca,'yscale','log')
% % % % set(gca,'Ylim',[1e7 1e14], 'ytick',[1e7,1e10,1e13],'fontsize',9);
% % % set(gca,'Ylim',[0,250], 'ytick',[0 100 200],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Loc Err [%]','fontsize',12);
% % % i=i+1;
%% dLoc Error
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) max(dLoc,[],2)], 'color','k', 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca,'Ylim',[0 max(Locerror)]);
% % % set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Loc Err','fontsize',12);
% % % i=i+1;
%% meandLoc Error
h(i)=irf_subplot(n,1,-i);

irf_plot([B1(:,1) mean(dLoc,2)], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Locerror)]);
set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('mean','fontsize',12);
i=i+1;
%% stddLoc Error
% % % h(i)=irf_subplot(n,1,-i);
% % % 
% % % irf_plot([B1(:,1) std(dLoc,0,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca,'Ylim',[0 max(Locerror)]);
% % % set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('std','fontsize',12);
% % % i=i+1;
%% Q Error
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) Qerror], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 1100]);
set(gca,'Ylim',[0,200], 'ytick',[0 100 200],'fontsize',9);
pos1=get(gca,'pos');
ylabel('Q Err [%]','fontsize',12);
i=i+1;
%% Q
% h(i)=irf_subplot(n,1,-i);
% irf_plot([B1(:,1) Q], 'color','k', 'Linewidth',0.75); hold on;
% grid off;
% set(gca,'Ylim',[min(Q) max(Q)]);
% % set(gca,'Ylim',[0 max(Q)]);
% % set(gca,'Ylim',[0 100], 'ytick',[0 500 1000],'fontsize',9);
% pos1=get(gca,'pos');
% ylabel('Q [nT·km^2]','fontsize',12);
% set(gca,'Yscale','log')
% i=i+1;
%% Adjust the position
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
% figname = [OutputDir,'OverviewFig\',NameTags{TDT}(2:end-2)];    
% print(gcf, '-dpng', [figname '.png']);

%% Init Figure 2
figure(2)
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Index id
length(find(Locerror<50 & Qerror<=1000))
% tempidx_B1 = find(PI~=0);
% tempidx_B1 = tempidx_B1(mean(dLoc(tempidx_B1,:),2) == min(mean(dLoc(tempidx_B1,:),2)));
% tempidx_B1 = find(mean(dLoc,2) == min(mean(dLoc,2)));
tempidx_B1 = 398;
% tempidx_B1 = tempidx_B1(1);
%     [~,tempidx_B] = max(abs(divB(:,2)));
c_eval('[~,tempidx_B?] = sort(abs(B?_gse(:,1)-B1_gse(tempidx_B1,1)));',2:4);
c_eval('tempidx_B? = tempidx_B?(1);')
[~,tempidx_R] = sort(abs(R1(:,1)-B1(tempidx_B1,1)));
tempidx_R = tempidx_R(1);
[Null_loc, Trans_mat, dB_null] = FOTE('R?','B?',epoch2iso(B1(tempidx_B1,1)));
%% disp
% Qid = find(Qerror == min(Qerror));
% % % Qid = tempidx_B1;
% % % disp(['Q:',num2str(Q(Qid))])
% % % disp(['Loc:',num2str((LocPoint(Qid)+CenterPoint(Qid,:)))])
% % % Qeb = [min(resQ{Qid}),max(resQ{Qid})];
% % % Loceb = [min(LocRes{Qid});max(LocRes{Qid})];
% % % disp(['QerrorBoundary:',num2str(Qeb(1)),' ',num2str(Qeb(2))])
% % % disp(['LocerrorBoundary:',num2str(Loceb(1,:)),' ',num2str(Loceb(2,:))])
%% Coordinate
% c_eval('R?(:,2:4) = R?(:,2:4)-LocPoint(tempidx_B1,:);');
% LocRes{tempidx_B1} = LocRes{tempidx_B1}-LocPoint(tempidx_B1,:);
% LocPoint = LocPoint - LocPoint(tempidx_B1,:);
Null_loc = [Null_loc.Rbc(1),Null_loc.Rbc(2),Null_loc.Rbc(3)];
c_eval('R?(:,2:4) = R?(:,2:4)-Null_loc;');
Null_loc = Null_loc - Null_loc;
%% Location 
cor = 'krgb';
c_eval("plot3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),'s' ,'color',cor(?),'linewidth',5);hold on;")

%% Tetrahedron configuration
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(tempidx_R,2),R',num2str(ii),'(tempidx_R,3),R',num2str(ii),'(tempidx_R,4);',...
    'R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
plot3(RR12(:,1),RR12(:,2),RR12(:,3),'--k');hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'--k');hold on;  
plot3(RR14(:,1),RR14(:,2),RR14(:,3),'--k');hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'--k');hold on;  
plot3(RR34(:,1),RR34(:,2),RR34(:,3),'--k');hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'--k');hold on;  
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
xlabel('e_1 [km]','fontsize',12);
ylabel('e_2 [km]','fontsize',12);
zlabel('e_3 [km]','fontsize',12);

%% Quiver
c_eval('B?_gse = irf_abs(B?_gse);');
maxB = max([B1_gse(tempidx_B1,5),B2_gse(tempidx_B1,5),B3_gse(tempidx_B1,5),B4_gse(tempidx_B1,5)]);
c_eval("quiver3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),RR_mean*B?_gse(tempidx_B1,2)/maxB,RR_mean*B?_gse(tempidx_B1,3)/maxB,RR_mean*B?_gse(tempidx_B1,4)/maxB,'color',cor(?));hold on;")

%% Loc res
% plotPolyhedron(LocRes{tempidx_B1}(:,1),LocRes{tempidx_B1}(:,2),LocRes{tempidx_B1}(:,3),'#FFB8CE',0.3);
% % % id = nchoosek(1:6,2);
% % % for i = 1:length(LocRes)
% % %     c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
% % %     ttd = [];
% % %     c_eval('ttd = [ttd,tempd?(4)/RR_mean(i)];',1:15);
% % %     dLoc(i,:) = ttd;
% % % end

% idx = [78,80,81];
% plot3(LocPoint(idx,1),LocPoint(idx,2),LocPoint(idx,3),'*','color','#FFBAF1');
% line(LocPoint(78:81,1),LocPoint(78:81,2),LocPoint(78:81,3),'color','#FFBAF1');
plot3(Null_loc,Null_loc,Null_loc,'*','color','m');
% plot3(LocPoint(tempidx_B1,1),LocPoint(tempidx_B1,2),LocPoint(tempidx_B1,3),'*','color','m');
% plot3(LocRes{tempidx_B1}(:,1),LocRes{tempidx_B1}(:,2),LocRes{tempidx_B1}(:,3),'*','color','#FFB8CE');
%% adjust the position
axis equal 
view(75,55)
% axis([-2*RR_mean 2*RR_mean -2*RR_mean 2*RR_mean -2*RR_mean 2*RR_mean])

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
