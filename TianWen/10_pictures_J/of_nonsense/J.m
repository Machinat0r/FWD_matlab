clear; clc; close all
% this script should be in the same folder with BRcart_year.mat
%------written by Wending Fu, Jun.2024 in Beijing------------
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
% ha=tight_subplot(1,3,[.1 .1],[.05 .05],[.1 .1]);
% LB={'All','Day','Night'};

coor = 'aeh';
gridNuma = 100; gridNume = 200; gridNumh = 4; % 200-300、300-500、500-700、700-900
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%day%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input data (sph)
inputFile = 'E:\Martian\programs\2_Maven_Download\10_pictures_J/T_BRsphere_za_year_2015_day.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year_day', 'var'), warning('wrong input data'); end
Time_day = T_BRsphere_za_year_day(:,1);
Bxyz_day = T_BRsphere_za_year_day(:,2:4);
Rsph_day = T_BRsphere_za_year_day(:,5:7); % lon, lat, alt = Azimuth, Elevation, High = a, e, h
%%
units = irf_units;
RM = units.Mars.radius; % RM就是火星半径
Rsph_day(:,3) = Rsph_day(:,3) + RM/1e3; % h = h + RM
Rsph_day(:,2) = Rsph_day(:,2) * pi / 180; % 化为rad，且e = [-pi/2, pi/2]
Rsph_day(:,1) = (Rsph_day(:,1) - 180) * pi / 180; % 化为rad，且a = [-pi, pi]
Bsph_day = CoorTrans(Bxyz_day, Rsph_day); % 进行同样的变化，phi, the ,r = a, e, h
%% divide grid
%coor = 'aeh';
minRa_day = min(Rsph_day(:,1)); maxRa_day = max(Rsph_day(:,1));
minRe_day = min(Rsph_day(:,2)); maxRe_day = max(Rsph_day(:,2));
minRh_day = min(Rsph_day(:,3)); maxRh_day = max(Rsph_day(:,3));
%gridNuma = 100; gridNume = 200; gridNumh = 4; 

% coor,gridNum?公用
% 若try中语句错误，则跳转到catch
try
load(['Ba_day',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat']);
catch
c_eval('grid?_day = linspace(minR?_day, maxR?_day, gridNum?);', coor);% 让grid?等于从min到max，均分共gridNum的数字
[aGrid_day, eGrid_day, hGrid_day] = meshgrid(grida_day, gride_day, gridh_day); 
c_eval('?Grid_day = permute(?Grid_day, [2, 1, 3]);',coor); % 交换第一个和第二个维度

[~, ~, aIndex_day] = histcounts(Rsph_day(:,1), grida_day);
[~, ~, eIndex_day] = histcounts(Rsph_day(:,2), gride_day);
[~, ~, hIndex_day] = histcounts(Rsph_day(:,3), gridh_day);

Bsph_day(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex_day, eIndex_day, hIndex_day);
% hIndex_day中最小值得大于0，用其他分法。

%% average B in grid
c_eval('B?Grid_day = nan(size(?Grid_day));', coor);

% parpool(12)
parfor iGrid_day = 1:numel(aGrid_day)
idx_day = Bsph_day(:,4) == iGrid_day;
if any(idx_day)
tempB_day = mean(Bsph_day(idx_day, 1:3), 1);
BaGrid_day(iGrid_day) = tempB_day(1);
BeGrid_day(iGrid_day) = tempB_day(2);
BhGrid_day(iGrid_day) = tempB_day(3);
end
end

save(['Ba_day',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat'],...
    "aGrid_day", "eGrid_day", "hGrid_day", "BaGrid_day", "BeGrid_day", "BhGrid_day");
end

%% calculate J
hGrid_day = hGrid_day .* 1e3; % m
c_eval('B?Grid_day = B?Grid_day .* 1e-9;', coor); % T
[curlBa_day, curlBe_day, curlBh_day] = CurlSph(BaGrid_day, BeGrid_day, BhGrid_day, aGrid_day, eGrid_day, hGrid_day); % A/m^2
c_eval('J?_day = 1e9 * curlB?_day / units.mu0;', coor);% nA/m^2
hGrid_day = hGrid_day / 1e3; %km
aGrid_day = aGrid_day * 180 / pi; %deg
eGrid_day = eGrid_day * 180 / pi; %deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%day%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nig%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input data (sph)
inputFile = 'E:\Martian\programs\2_Maven_Download\10_pictures_J/T_BRsphere_za_year_2015_nig.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year_nig', 'var'), warning('wrong input data'); end
Time_nig = T_BRsphere_za_year_nig(:,1);
Bxyz_nig = T_BRsphere_za_year_nig(:,2:4);
Rsph_nig = T_BRsphere_za_year_nig(:,5:7); % lon, lat, alt = Azimuth, Elevation, High = a, e, h
%%
units = irf_units;
RM = units.Mars.radius; % RM就是火星半径
Rsph_nig(:,3) = Rsph_nig(:,3) + RM/1e3; % h = h + RM
Rsph_nig(:,2) = Rsph_nig(:,2) * pi / 180; % 化为rad，且e = [-pi/2, pi/2]
Rsph_nig(:,1) = (Rsph_nig(:,1) - 180) * pi / 180; % 化为rad，且a = [-pi, pi]
Bsph_nig = CoorTrans(Bxyz_nig, Rsph_nig); % 进行同样的变化，phi, the ,r = a, e, h
%% divide grid
%coor = 'aeh';
minRa_nig = min(Rsph_nig(:,1)); maxRa_nig = max(Rsph_nig(:,1));
minRe_nig = min(Rsph_nig(:,2)); maxRe_nig = max(Rsph_nig(:,2));
minRh_nig = RM/1000 + 100; maxRh_nig = RM/1000 + 900;
%gridNuma = 100; gridNume = 200; gridNumh = 4; 

% coor,gridNum?公用
% 若try中语句错误，则跳转到catch
try
load(['Ba_nig',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat']);
catch
c_eval('grid?_nig = linspace(minR?_nig, maxR?_nig, gridNum?);', coor);% 让grid?等于从min到max，均分共gridNum的数字
[aGrid_nig, eGrid_nig, hGrid_nig] = meshgrid(grida_nig, gride_nig, gridh_nig); 
c_eval('?Grid_nig = permute(?Grid_nig, [2, 1, 3]);',coor); % 交换第一个和第二个维度

[~, ~, aIndex_nig] = histcounts(Rsph_nig(:,1), grida_nig);
[~, ~, eIndex_nig] = histcounts(Rsph_nig(:,2), gride_nig);
[~, ~, hIndex_nig] = histcounts(Rsph_nig(:,3), gridh_nig);

Bsph_nig(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex_nig, eIndex_nig, hIndex_nig);

%% average B in grid
c_eval('B?Grid_nig = nan(size(?Grid_nig));', coor);

% parpool(12)
parfor iGrid_nig = 1:numel(aGrid_nig)
idx_nig = Bsph_nig(:,4) == iGrid_nig;
if any(idx_nig)
tempB_nig = mean(Bsph_nig(idx_nig, 1:3), 1);
BaGrid_nig(iGrid_nig) = tempB_nig(1);
BeGrid_nig(iGrid_nig) = tempB_nig(2);
BhGrid_nig(iGrid_nig) = tempB_nig(3);
end
end

save(['Ba_nig',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat'],...
    "aGrid_nig", "eGrid_nig", "hGrid_nig", "BaGrid_nig", "BeGrid_nig", "BhGrid_nig");
end

%% calculate J
hGrid_nig = hGrid_nig .* 1e3; % m
c_eval('B?Grid_nig = B?Grid_nig .* 1e-9;', coor); % T
[curlBa_nig, curlBe_nig, curlBh_nig] = CurlSph(BaGrid_nig, BeGrid_nig, BhGrid_nig, aGrid_nig, eGrid_nig, hGrid_nig); % A/m^2
c_eval('J?_nig = 1e9 * curlB?_nig / units.mu0;', coor);% nA/m^2
hGrid_nig = hGrid_nig / 1e3; %km
aGrid_nig = aGrid_nig * 180 / pi; %deg
eGrid_nig = eGrid_nig * 180 / pi; %deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nig%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%all%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input data (sph)
inputFile = 'E:\Martian\programs\2_Maven_Download\10_pictures_J/T_BRsphere_za_year_2015_all.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year_all', 'var'), warning('wrong input data'); end
Time_all = T_BRsphere_za_year_all(:,1);
Bxyz_all = T_BRsphere_za_year_all(:,2:4);
Rsph_all = T_BRsphere_za_year_all(:,5:7); % lon, lat, alt = Azimuth, Elevation, High = a, e, h
%%
units = irf_units;
RM = units.Mars.radius; % RM就是火星半径
Rsph_all(:,3) = Rsph_all(:,3) + RM/1e3; % h = h + RM
Rsph_all(:,2) = Rsph_all(:,2) * pi / 180; % 化为rad，且e = [-pi/2, pi/2]
Rsph_all(:,1) = (Rsph_all(:,1) - 180) * pi / 180; % 化为rad，且a = [-pi, pi]
Bsph_all = CoorTrans(Bxyz_all, Rsph_all); % 进行同样的变化，phi, the ,r = a, e, h
%% divide grid
%coor = 'aeh';
minRa_all = min(Rsph_all(:,1)); maxRa_all = max(Rsph_all(:,1));
minRe_all = min(Rsph_all(:,2)); maxRe_all = max(Rsph_all(:,2));
minRh_all = RM/1000 + 100; maxRh_all = RM/1000 + 900;
%gridNuma = 100; gridNume = 200; gridNumh = 4; 

% coor,gridNum?公用
% 若try中语句错误，则跳转到catch
try
load(['Ba_all',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat']);
catch
c_eval('grid?_all = linspace(minR?_all, maxR?_all, gridNum?);', coor);% 让grid?等于从min到max，均分共gridNum的数字
[aGrid_all, eGrid_all, hGrid_all] = meshgrid(grida_all, gride_all, gridh_all); 
c_eval('?Grid_all = permute(?Grid_all, [2, 1, 3]);',coor); % 交换第一个和第二个维度

[~, ~, aIndex_all] = histcounts(Rsph_all(:,1), grida_all);
[~, ~, eIndex_all] = histcounts(Rsph_all(:,2), gride_all);
[~, ~, hIndex_all] = histcounts(Rsph_all(:,3), gridh_all);

Bsph_all(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex_all, eIndex_all, hIndex_all);

%% average B in grid
c_eval('B?Grid_all = nan(size(?Grid_all));', coor);

% parpool(12)
parfor iGrid_all = 1:numel(aGrid_all)
idx_all = Bsph_all(:,4) == iGrid_all;
if any(idx_all)
tempB_all = mean(Bsph_all(idx_all, 1:3), 1);
BaGrid_all(iGrid_all) = tempB_all(1);
BeGrid_all(iGrid_all) = tempB_all(2);
BhGrid_all(iGrid_all) = tempB_all(3);
end
end

save(['Ba_all',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat'],...
    "aGrid_all", "eGrid_all", "hGrid_all", "BaGrid_all", "BeGrid_all", "BhGrid_all");
end

%% calculate J
hGrid_all = hGrid_all .* 1e3; % m
c_eval('B?Grid_all = B?Grid_all .* 1e-9;', coor); % T
[curlBa_all, curlBe_all, curlBh_all] = CurlSph(BaGrid_all, BeGrid_all, BhGrid_all, aGrid_all, eGrid_all, hGrid_all); % A/m^2
c_eval('J?_all = 1e9 * curlB?_all / units.mu0;', coor);% nA/m^2
hGrid_all = hGrid_all / 1e3; %km
aGrid_all = aGrid_all * 180 / pi; %deg
eGrid_all = eGrid_all * 180 / pi; %deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%all%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%需要重新布局图形
%需要分高度：100~300
%需要分日夜侧


%% Init figure
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 200; ySize = 70; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% plot slide
slide = 5;
Jt = sqrt(Ja(:,:,slide).^2 + Je(:,:,slide).^2 + Jh(:,:,slide).^2);
Bt = sqrt(BaGrid(:,:,slide).^2 + BeGrid(:,:,slide).^2 + BhGrid(:,:,slide).^2);
% pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Jt);
pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Jh(:,:,slide));
%画图，注意J?需要更改！
%后续只用区分日夜侧和高度即可
shading interp

% colormap('jet'); 
colormap(othercolor('BuDRd_12'))
c = colorbar;
c.Label.String = 'J_{elevation} [nA/m^2]';% 记得改：phi、theta、r
xlabel('Lon [\circ]'); ylabel('Lat [\circ]');
hold on;

% quiver(aGrid(:,:,slide),eGrid(:,:,slide),Ja(:,:,slide),Je(:,:,slide),'w')

%%
set(gca,'FontName','Times New Roman','FontSize',10);

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gcf,'color','w')

%% coor Trans函数
function Bsph = CoorTrans(B, Rsph)
theta = pi/2 - Rsph(:,2); phi = Rsph(:,1);
Br = B(:,1).* sin(theta).* cos(phi) + B(:,2).* sin(theta).* sin(phi) + B(:,3).* cos(theta);
Bt = B(:,1).* cos(theta).* cos(phi) + B(:,2).* cos(theta).* sin(phi) - B(:,3).* sin(theta);
Bp = -B(:,1).* sin(phi) + B(:,2).* cos(phi);
Bsph = [Bp, Bt, Br];
end