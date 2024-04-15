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
ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

% TT = '2019-01-16T04:09:50.00Z/2019-01-16T04:10:00.00Z';
TT = '2019-01-16T04:09:55.605Z/2019-01-16T04:09:55.630Z'; %no boundary, 10,78-81
% TT = '2019-01-16T04:09:55.220Z/2019-01-16T04:09:56.000Z'; %no boundary, 10,78-81
% TT = '2019-01-16T04:09:55.420Z/2019-01-16T04:09:55.800Z'; %no boundary, 10,78-81
% TT = '2018-08-27T12:15:30.00Z/2018-08-27T12:15:50.00Z';

% TT = '2020-03-07T12:38:25.00Z/2020-03-07T12:38:30.00Z';
% TT = '2020-05-03T01:58:25.00Z/2020-05-03T01:58:35.00Z';
% TT = '2017-10-04T08:11:45.00Z/2017-10-04T08:11:50.00Z';
% TT = '2021-03-30T07:58:45.00Z/2021-03-30T07:58:55.00Z';%近似平面的构型，做反例
% TT = '2017-10-08T03:44:05.00Z/2017-10-08T03:44:15.00Z';
% TT = '2017-10-06T22:07:40.00Z/2017-10-06T22:07:50.00Z';

% TT = '2015-09-19T07:43:28.000Z/2015-09-19T07:43:33.000Z';

tint=irf.tint(TT);
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];
ic = 1:4;
filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
[filenames,~,~] = findFilenames(TT,filenames1,'brst',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 32, 'CheckSize', 0)
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
    [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?',i,i*sign(divB(i,2)),10 ,1);
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
%% Coordinate
% % % c_eval('R?(:,2:4) = R?(:,2:4)-LocPoint;');

maxR = nan*zeros(1,3); minR = nan*zeros(1,3);
c_eval('maxR = max(maxR,max(R?(:,2:4),[],1));',1:4);
c_eval('minR = min(minR,min(R?(:,2:4),[],1));',1:4);
GridX = linspace(minR(1),maxR(1),20);
GridY = linspace(minR(2),maxR(2),20);
GridZ = linspace(minR(3),maxR(3),20);
[MGX, MGY, MGZ] = meshgrid(GridX, GridY, GridZ);
X = MGX(:);Y = MGY(:);Z = MGZ(:);
points = [X, Y, Z];

c_eval('V? = [R1(?,2:4); R2(?,2:4); R3(?,2:4); R4(?,2:4)];',1:size(B1,1));
c_eval('Tri? = delaunayTriangulation(V?);',1:size(B1,1));
c_eval('points? = points(pointLocation(Tri?,points) == 1,:);',1:size(B1,1));

for i = 1:length(PI)
    LocRes{i} = LocRes{i}-LocPoint(i,:);
end
% % % LocPoint = LocPoint - LocPoint;
%% VSH
c_eval('R?_grid = cell(size(points1,1));', 1:4);
coeff_grid = cell(size(points1,1)); res_grid = cell(size(points1,1));
c_eval('rate1?_grid = zeros(size(B?,1),size(points1,1));', 1:4);
c_eval('rate2?_grid = zeros(size(B?,1),size(points1,1));', 1:4);
for i_grid = 1:size(points1,1)
c_eval('R?_grid{i_grid} = [R?(:,1), R?(:,2:4) - points1(i_grid,:)];',1:4);
[coeff, res] = VSH_Expand('R?_grid{i_grid}', 'B?(:,1:4)', [0,0,0], 1);
coeff_grid{i_grid} = coeff; res_grid{i_grid} = res;


% % % [coeff, res] = VSH_Expand('R?', 'B?(:,1:4)', [0,0,0], 1);
% coeff0 = coeff([1,5,9,13,14,15]);
% coeff1 = coeff([2,3,4,6,7,8,10,11,12]);
alpha_00 = coeff(:,1); beta_00 = coeff(:,5); gamma_00 = coeff(:,9);
alpha_12 = coeff(:,2); beta_12 = coeff(:,6); gamma_12 = coeff(:,10);
alpha_10 = coeff(:,3); beta_10 = coeff(:,7); gamma_10 = coeff(:,11);
alpha_11 = coeff(:,4); beta_11 = coeff(:,8); gamma_11 = coeff(:,12);
c_eval('[r?, theta?, phi?] = Coor_Trans(R?(:,2:4));', 1:4);
[r1, r2, r3, r4, ~] = Normal_R(r1, r2, r3, r4);
for temp_ic = 1:4
temp_ic = num2str(temp_ic);
eval(['Br0',temp_ic,' = Cal_Br0(r',temp_ic,', theta',temp_ic,', alpha_00);']);
eval(['Br1',temp_ic,' = Cal_Br1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', alpha_12, alpha_10, alpha_11);']);
eval(['Bt1',temp_ic,' = Cal_Bt1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', beta_12, beta_10, beta_11, gamma_12, gamma_11);']);
eval(['Bp1',temp_ic,' = Cal_Bp1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', beta_12, beta_11, gamma_12, gamma_10, gamma_11);']);
end
c_eval('[~, B?_sph] = Coor_Trans(R?(:,2:4), B?(:,2:4));', 1:4);
c_eval('B0? = irf_abs([B?(:,1), Br0?, zeros(size(Br0?)), zeros(size(Br0?))]);', 1:4);
c_eval('B1? = irf_abs([B?(:,1), Br1?, Bt1?, Bp1?]);', 1:4);
c_eval('B2? = irf_abs([B?(:,1), B?_sph(:,1) - Br0? - Br1?, B?_sph(:,2) - Bt1?, B?_sph(:,3) - Bp1?]);', 1:4);
c_eval('rate1? = B0?(:,5)./B1?(:,5);', 1:4);
c_eval('rate2? = B0?(:,5)./B2?(:,5);', 1:4);

c_eval('rate1?_grid(:,i_grid) = rate1?;',1:4);
c_eval('rate2?_grid(:,i_grid) = rate2?;',1:4);
end
%% cal B
for i = 1:length(PI)
c_eval('[tempBx?,tempBy?,tempBz?,tempBt?]=calculateB(R?(i,2),R?(i,3),R?(i,4),Q(i),LocPoint(i,:));')
c_eval('Bm?(i,:) = [B1(i,1),tempBx?,tempBy?,tempBz?,tempBt?];')
end
gradBm=c_4_grad('R?','Bm?(:,1:4)','grad');
divBm=[gradBm(:,1) sum([gradBm(:,2) gradBm(:,6) gradBm(:,10)],2)]; 
divErr = [divB(:,1) 100*abs(divB(:,2)-divBm(:,2))./abs(divBm(:,2))]; %|deltB-deltBm|/deltB
%% B error
c_eval('R? = irf_abs(R?);')
c_eval('Bcross? = irf_dot(B?(:,2:4),Bm?(:,2:4));')
c_eval('Btheta? = acosd(Bcross?./(B?(:,5).*Bm?(:,5)))./RR_mean.*R?(:,5);')
c_eval('R? = irf_abs(R?);')
c_eval('theta? = acos(Bcross?./(B?(:,5).*Bm?(:,5)));')
Rres = cellfun(@irf_abs,LocRes,'UniformOutput',false);
Rres = cellfun(@(x) [max(x([1,4,5],4)),max(x([1,2,6],4)),max(x([2,3,4],4)),max(x([3,5,6],4)),mean(x(:,4))],Rres,'UniformOutput',false);
Rres = cell2mat(Rres);
% c_eval('Btheta? = 100*theta?./RR_mean.*Rres(:,?);')
% % c_eval('Btheta? = 100*theta?./RR_mean.*R?(:,5);')
% c_eval('Btheta? = asin(sqrt(1-(Bcross?./(B?(:,5).*Bm?(:,5))).^2))*180/pi;')
% c_eval('Blength? = abs(B?(:,5)-Bm?(:,5))./Bm?(:,5)./RR_mean.*Rres(:,?);')
c_eval('Blength? = abs(B?(:,5)-Bm?(:,5))./Bm?(:,5)./RR_mean.*R?(:,5);')
% % % c_eval('Blength? = abs(B?(:,5)-Bm?(:,5))./Bm?(:,5)./RR_mean.*R?(:,5);')
Blength = 0.25*(Blength1+Blength2+Blength3+Blength4);
c_eval('BErr? = 100*(1-(1-Btheta?).*(1-Blength?));')
BErr = mean([BErr1,BErr2,BErr3,BErr4],2);
%% Init figure 3
h = figure(3);
set(gcf,'PaperUnits','centimeters')
xSize = 80; ySize = 150; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% plot distance
h(1) = subplot(7,1,1);
for i = 1:length(resQ)
    RR_temp = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
    'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
c_eval(['RR_temp=RR_temp+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean(i) = RR_temp(4)/6;
end

id = nchoosek(1:6,2);
for i = 1:length(LocRes)
    c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
    ttd = [];
    c_eval('ttd = [ttd,tempd?(4)/RR_mean(i)];',1:15);
    dLoc(i,:) = ttd;
end
    meand = mean(dLoc,2);
    
% 15 errorbars
resd = dLoc-meand;
resdpos = resd; resdneg = resd;
resdpos(resdpos<0) = 0;resdneg(resdneg>=0) = 0;
c_eval('d?neg = errorbar(B1(:,1),meand,resdneg(:,?),resdpos(:,?));hold on;',1:15);
line(B1(:,1),meand,'color','#BEBEBE','LineWidth',1);
scatter(B1(:,1),meand,'d','filled','MarkerFaceColor','#EDB120')
% c_eval('d?neg = errorbar(B1(:,1),meand,resdneg(:,?),zeros(size(resdneg(:,?))));hold on;',1:6);
% c_eval('d?pos = errorbar(B1(:,1),meand,zeros(size(resdpos(:,?))),resdpos(:,?));hold on;',1:6);
% c_eval("d?pos.Marker = 'o';",1:6)
% c_eval("d?pos.Color = 'k';",1:6)
% c_eval("d?neg.Marker = 'o';",1:6)
c_eval("d?neg.Color = '#696969';",1:15)
c_eval("d?neg.LineWidth = 0.8;",1:15)
c_eval("d?neg.CapSize = 10;",1:15)

% errorbar
% % % resd = dLoc-meand;
% % % resdpos = max(resd,[],2); resdneg = min(resd,[],2);
% % % e1 = errorbar(B1(:,1),meand,resdneg,resdpos);hold on;
% % % e1.LineWidth = 1;
% % % e1.Color = '#696969';
% % % quart = quantile(dLoc,[.25,.75],2);
% % % c_eval("line([B1(?,1),B1(?,1)],[quart(?,1),quart(?,2)],'color','#BEBEBE','linewidth',5);hold on;",1:length(Q))
% % % line(B1(:,1),meand,'linewidth',2,'color','#BEBEBE');hold on;
% % % scatter(B1(:,1),meand,'d','filled','MarkerFaceColor','#EDB120')

ax = gca;
ax.YLim = [0,2];
% irf_timeaxis(ax,'date')
ylabel('C_R ','fontsize',12);
ax.XTickLabel = '';
%% subplot Q
h(2) = subplot(7,1,2);

% QErrorbar(B1(:,1),Q,resQ)

% 6 errorbars
for i = 1:length(resQ)
c_eval('resQ?(i) = resQ{i}(?)-Q(i);',1:6)
end
c_eval('resQ?pos = resQ?;',1:6);c_eval('resQ?neg = resQ?;',1:6);
c_eval('resQ?pos(resQ?<0) = 0;',1:6);
c_eval('resQ?neg(resQ?>=0) = 0;',1:6);

c_eval('e?pos = errorbar(B1(:,1),Q,resQ?neg,resQ?pos);hold on;',1:6)
% c_eval('e?pos = errorbar(B1(:,1),Q,resQ?neg,zeros(size(resQ?pos)));hold on;',1:6)
% c_eval('e?neg = errorbar(B1(:,1),Q,zeros(size(resQ?neg)),resQ?pos);hold on;',1:6)
% c_eval("e?pos.LineWidth = 1.5;",1:6)
scatter(B1(:,1),Q,'d','filled','MarkerFaceColor','#EDB120')
c_eval("e?pos.Color = 'k';",1:6)
% c_eval("e?neg.Marker = 'o';",1:6)
% c_eval("e?neg.Color = 'k';",1:6)

% errorbar
% % % c_eval('resQmat(?,:) = resQ{?}-Q(?);',1:length(Q))
% % % 
% % % resQpos = max(resQmat,[],2); resQneg = min(resQmat,[],2);
% % % e2 = errorbar(B1(:,1),Q,resQneg,resQpos);hold on;
% % % e2.LineWidth = 1;
% % % e2.Color = '#696969';
% % % e2.CapSize = 0.01;
% % % quart = quantile(transpose(reshape(cell2mat(resQ),6,[])),[.25,.75],2);
% % % c_eval("line([B1(?,1),B1(?,1)],[quart(?,1),quart(?,2)],'color','#BEBEBE','linewidth',5);hold on;",1:length(Q))
% % % line(B1(:,1),Q,'linewidth',2,'color','#BEBEBE');hold on;
% % % scatter(B1(:,1),Q,'d','filled','MarkerFaceColor','#EDB120')
% % % scatter(B1(:,1),resQpos+Q,'filled','MarkerFaceColor','k')

ylabel('Q [nT·km^2] ','fontsize',12);
ax = gca;
ax.YScale = 'log';
ax.YLim = [1e3,1e6];
% irf_timeaxis(ax,'date')
xticks('');
%% plot div err
h(3) = subplot(7,1,3);

line(B1(:,1),divErr(:,2),'color','#8470FF');hold on;
s1 = scatter(B1(:,1),divErr(:,2),100,'p','filled','MarkerFaceColor','#7E2F8E');hold on;

ax = gca;
ax.YLim = [0,100];
% irf_timeaxis(ax,'date')
box('on')
ylabel('\delta [%]','fontsize',12);
ax.XTickLabel = '';
%% Model direction Error
h(4) = subplot(7,1,4);
cor = {'#000000','#D95319','#77AC30','#0072BD'};
c_eval("line(B?(:,1),Btheta?,'color',cor{?});hold on;")
s1 = scatter(B1(:,1),Btheta1,50,'k','filled');hold on;
s2 = scatter(B2(:,1),Btheta2,50,'filled');hold on;
s2.MarkerFaceColor = '#D95319';
s3 = scatter(B3(:,1),Btheta3,50,'filled');hold on;
s3.MarkerFaceColor = '#77AC30';
s4 = scatter(B4(:,1),Btheta4,50,'filled');hold on;
s4.MarkerFaceColor = '#0072BD';


ax = gca;
ax.YLim = [0,100];
box('on')
% irf_timeaxis(ax,'date')
ylabel('\kappa_\theta [%]','fontsize',12);
ax.XTickLabel = '';
%% Model strength Error
h(5) = subplot(7,1,5);
c_eval("line(B?(:,1),100*Blength?,'color',cor{?});hold on;")
% cor = {'#000000','#D95319','#77AC30','#0072BD'};
s1 = scatter(B1(:,1),100*Blength1,50,'k','filled');hold on;
s2 = scatter(B2(:,1),100*Blength2,50,'filled');hold on;
s2.MarkerFaceColor = '#D95319';
s3 = scatter(B3(:,1),100*Blength3,50,'filled');hold on;
s3.MarkerFaceColor = '#77AC30';
s4 = scatter(B4(:,1),100*Blength4,50,'filled');hold on;
s4.MarkerFaceColor = '#0072BD';


ax = gca;
ax.YLim = [0,100];
box('on')
% irf_timeaxis(ax,'date')
ylabel('\kappa_\rho [%]','fontsize',12);
ax.XTickLabel = '';
%% SPH rate 1: l0/l1
h(6)=subplot(7,1,6);
c_eval("line(B?(:,1), rate1?, 'color', cor{?}); hold on;");
line(B1(:,1), ones(size(B1(:,1))),'LineStyle','--', 'Linewidth',0.25); hold on;

s1 = scatter(B1(:,1),rate11,60,'k','filled','hexagram');hold on;
s2 = scatter(B2(:,1),rate12,60,'filled','hexagram');hold on;
s2.MarkerFaceColor = '#D95319';
s3 = scatter(B3(:,1),rate13,60,'filled','hexagram');hold on;
s3.MarkerFaceColor = '#77AC30';
s4 = scatter(B4(:,1),rate14,60,'filled','hexagram');hold on;
s4.MarkerFaceColor = '#0072BD';
grid off;

% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);

ylabel('\Re_1','fontsize',12);
ax = gca;
ax.YScale = 'log';
ax.YLim = [1e0,1e2];
ax.XTickLabel = '';
% ax.YTick = [0.1, 5, 10];
box('on')
% irf_timeaxis(ax,'date')
%% SPH rate 2: l0/l2
h(7)=subplot(7,1,7);
c_eval("line(B?(:,1), rate2?, 'color', cor{?}); hold on;");
line(B1(:,1), ones(size(B1(:,1))),'LineStyle','--', 'Linewidth',0.25); hold on;
s1 = scatter(B1(:,1),rate21,60,'k','filled','hexagram');hold on;
s2 = scatter(B2(:,1),rate22,60,'filled','hexagram');hold on;
s2.MarkerFaceColor = '#D95319';
s3 = scatter(B3(:,1),rate23,60,'filled','hexagram');hold on;
s3.MarkerFaceColor = '#77AC30';
s4 = scatter(B4(:,1),rate24,60,'filled','hexagram');hold on;
s4.MarkerFaceColor = '#0072BD';
grid off;


% set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('\Re_2','fontsize',12);
ax = gca;
ax.YScale = 'log';
ax.YLim = [1e0,1e2];
% ax.YTick = [0.1, 5, 10];
box('on')
irf_timeaxis(ax,'date')
set(ax,'XTickLabelRotation',0);
%% Adjust the position
set(gcf,'render','painters');
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align(h);
% ax.YLimMode = 'auto';
set(gcf,'paperpositionmode','auto')



%% two-step component function
function Br0 = Cal_Br0(r, theta, alpha_00)
Br0 = r.^-2 .* alpha_00 .* Cal_SPH(0, 0, theta)';
end

function Bt0 = Cal_Bt0(r, theta, beta_00)
Bt0 = r.^-2 .* beta_00 .* Cal_dPlm(0, 0, theta)';
end

function Bp0 = Cal_Bp0(r, theta, gamma_00)
Bp0 = r.^-2 .* gamma_00 .* Cal_dPlm(0, 0, theta)';
end

function Br1 = Cal_Br1(r, theta, phi, alpha_12, alpha_10, alpha_11)
Br1 = r.^-3 .* (alpha_12 .* cos(-phi) .* Cal_SPH(1, -1, theta)'...
    + alpha_10 .* Cal_SPH(1, 0, theta)' + ...
    alpha_11 .* cos(phi) .* Cal_SPH(1, 1, theta)');
end

function Bt1 = Cal_Bt1(r, theta, phi, beta_12, beta_10, beta_11, gamma_12, gamma_11)
Bt1 = r.^-3 .* (gamma_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta)' + ...
    beta_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta)' + ...
    beta_10 ./ 2 .* Cal_dPlm(1, 0, theta)' + ...
    gamma_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta)' + ...
    beta_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta)');
end

function Bp1 = Cal_Bp1(r, theta, phi, beta_12, beta_11, gamma_12, gamma_10, gamma_11)
Bp1 = r.^-3 .* (gamma_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta)' -...
    beta_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta)' + ...
    gamma_10 ./ 2 .* Cal_dPlm(1, 0, theta)' + ...
    gamma_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta)' -...
    beta_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta)');
end

%% Calculate Spherical Harmonic
function Ylm = Cal_SPH(l, m, theta)
% theta, phi should be rad
% without exp(i*m*phi), it is placed in the magnetic field calculation

if abs(m) > l, Ylm = 0; return; end

% unnormal legendre funciton
Plm = legendre(l, cos(theta));
if m < 0, Plm = Plm(-m + 1,:); else, Plm = Plm(m + 1,:); end

% orthonormalized coefficient 
a = (2*l+1)*factorial(l-m);
b = 4*pi*factorial(l+m);
C = sqrt(a/b);

Ylm = sqrt(2) .* C .* Plm;
end

%% Calculate Difference of Legendre Function
function dPlm = Cal_dPlm(l, m, theta)
Clm1 = 0.5 * sqrt((l + m) * (l - m + 1));
Clm2 = 0.5 * sqrt((l + m + 1) * (l - m));
Plm1 = Cal_SPH(l, m-1, theta);
Plm2 = Cal_SPH(l, m+1, theta);

dPlm = Clm1 .* Plm1 - Clm2 .* Plm2;
end
%% R normalization
function [R1, R2, R3, R4, maxR] = Normal_R(R1, R2, R3, R4)
maxR = 0;
c_eval('maxR = max(maxR,R?(:,1));',1:4);
c_eval('R? = R?./maxR;',1:4);
end