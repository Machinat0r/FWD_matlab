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
ParentDir = 'C:\MMS\'; 
TempDir = 'C:\MMS\temp\';mkdir(TempDir);
% TT = '2017-10-21T22:26:28.00Z/2017-10-21T22:26:32.00Z';

TT =  '2019-01-16T04:09:55.605Z/2019-01-16T04:09:55.630Z';
% TT = '2021-03-30T07:58:45.00Z/2021-03-30T07:58:55.00Z';%近似平面的构型，做反例
% TT = '2017-10-08T03:44:05.00Z/2017-10-08T03:44:15.00Z';


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

TTlist = strjoin(regexp(TT,'\d+','match'),'');
i = 1;flag = 0;%若flag=0，说明整段时间都在第i-1个Tag里
while i<=length(NameTags) && str2double(TTlist(17:30)) > NameTags{i}
    if str2double(TTlist(1:14)) < NameTags{i}
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

c_eval("B?_ts=mms.get_data('B_gsm_brst',tint,?);");
c_eval('B?_gsm = irf.ts2mat(B?_ts);'); 
c_eval('B? = irf_abs(B?_gsm);');
c_eval('B? = irf_resamp(B?,B1);',2:4);

% PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
PI(abs(PI)<0.5) = 0;
%% Monopole Index        
Pos = mms.get_data('R_gsm',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gsmR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
c_eval('R? = irf_resamp(R?,B1);')
CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

% monopole_index=zeros(length(B1(:,1)),1);
gradB=c_4_grad('R?','B?_gsm','grad');
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
% if PI(i) ~= 0
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
    Qerror(i) = 100*std(resQ{i})/Q(i);
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
%% Index id
tempidx_B1 = find(mean(dLoc,2) == min(mean(dLoc,2)));
% tempidx_B1 = 79;
%     [~,tempidx_B] = max(abs(divB(:,2)));
c_eval('[~,tempidx_B?] = sort(abs(B?_gsm(:,1)-B1_gsm(tempidx_B1,1)));',2:4);
c_eval('tempidx_B? = tempidx_B?(1);')
[~,tempidx_R] = sort(abs(R1(:,1)-B1(tempidx_B1,1)));
tempidx_R = tempidx_R(1);
%% Init Figure 2
figure(2)
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Coordinate
c_eval('R?(:,2:4) = R?(:,2:4)-LocPoint(tempidx_B1,:);');
LocRes{tempidx_B1} = LocRes{tempidx_B1}-LocPoint(tempidx_B1,:);
LocPoint = LocPoint - LocPoint(tempidx_B1,:);
%% Location 
cor = 'krgb';
c_eval("plot3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),'s' ,'color',cor(?),'linewidth',5);hold on;")

%% Tetrahedron configuration
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
    'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4); 
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
c_eval(['plot3(RR',num2str(ii),'?(:,1),RR',num2str(ii),'?(:,2),RR',num2str(ii),'?(:,3),''--w '');hold on;'],ii+1:4)
end
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
xlabel('e_1 [km]','fontsize',12);
ylabel('e_2 [km]','fontsize',12);
zlabel('e_3 [km]','fontsize',12);

%% Quiver
c_eval('B?_gsm = irf_abs(B?_gsm);');
maxB = max([B1_gsm(tempidx_B1,5),B2_gsm(tempidx_B1,5),B3_gsm(tempidx_B1,5),B4_gsm(tempidx_B1,5)]);
% c_eval("quiver3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),RR_mean*B?_gsm(tempidx_B1,2)/maxB,RR_mean*B?_gsm(tempidx_B1,3)/maxB,RR_mean*B?_gsm(tempidx_B1,4)/maxB,'color',cor(?));hold on;")

%% adjust the position
axis equal 
%% Reconstruction
% RR_max = 1.2*max(abs([R1(tempidx_R,2:4),R2(tempidx_R,2:4),R3(tempidx_R,2:4),R4(tempidx_R,2:4)]));
RR_max = 25;
ss = 200;
step = 100/ss;
RR_s = 25;
[gridX,gridY,gridZ] = meshgrid([-RR_s:step:RR_s],[-RR_s:step:RR_s],[-RR_s:step:RR_s]);
gridB = cell(size(gridX));gridBt = zeros(size(gridX));gridBX = gridBt;gridBY = gridBX; gridBZ = gridBX;
for i = 1:numel(gridX)
    [gridBX(i),gridBY(i),gridBZ(i),gridBt(i)] = calculateB(gridX(i),gridY(i),gridZ(i),Q(tempidx_B1),LocPoint(tempidx_B1,:));
    gridB{i} = [gridBX(i),gridBY(i),gridBZ(i)];
end

% [stX,stY,stZ] = meshgrid([-5,-3,-1.5,-0.5,0,0.5,1.5,3,5],[-5,-3,-1.5,-0.5,0,0.5,1.5,3,5],[-5,-3,-1.5,-0.5,0,0.5,1.5,3,5]);
ls = [-1,0,1];
[stX,stY,stZ] = meshgrid(ls,ls,ls);
lines = stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ,[1,100]);
for i = 1:length(lines)
    lines_len = length(lines{i});
    Bt = ones(lines_len,1);
    for j = 1:lines_len
        [~,~,~,Bt(j)] = calculateB(lines{i}(j,1),lines{i}(j,2),lines{i}(j,3),Q(tempidx_B1),LocPoint(tempidx_B1,:));
    end
    try
    cline(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),Bt,0,2,'Blues9',log(Bt)/log(5)); hold on;
    catch
        disp(['error at lines ',num2str(i)])
    end
end

ls = [-2:1:2];
[stX,stY,stZ] = meshgrid(ls,ls,ls);
lines = stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ,[1,100]);
for i = 1:length(lines)
    lines_len = length(lines{i});
    Bt = ones(lines_len,1);
    for j = 1:lines_len
        [~,~,~,Bt(j)] = calculateB(lines{i}(j,1),lines{i}(j,2),lines{i}(j,3),Q(tempidx_B1),LocPoint(tempidx_B1,:));
    end
    try
    cline(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),Bt,0,2,'Blues9',0.5); hold on;
    catch
        disp(['error at lines ',num2str(i)])
    end
end

cmap = othercolor('Blues9');
cmap = cmap(32:256,:);
colormap(cmap);
c = colorbar;
caxis([0,2])
set(c,'YTick',[0,1,2]);
set(c,'YTickLabel',{'10^0','10^1','10^2'}) ;
% set(c,'Color','w')
set(get(c,'ylabel'),'string','|B| [nT]','fontsize',12);
%% adjust the position
axis equal 
axis([-RR_max RR_max -RR_max RR_max -RR_max RR_max])
view(40,42)
% axis off

box on;
set(gca,'Color','k')
set(gca,'XColor','w');set(gca,'YColor','w');set(gca,'ZColor','w')
set(gcf,'color','k')
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname = ['C:\Users\fwd\Desktop\Ti~mor~\M\magnetic_monopole\Figure\new\reconstruction\' ,'20190116_reconstruction-79'];
% print(gcf, '-dpdf', [figname '.pdf']);