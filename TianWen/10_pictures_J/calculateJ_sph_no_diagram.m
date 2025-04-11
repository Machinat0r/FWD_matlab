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
%% input data (sph)
inputFile = 'E:\Martian\programs\2_Maven_Download\10_pictures_J/T_BRsphere_za_year_2015.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year', 'var'), warning('wrong input data'); end
Time = T_BRsphere_za_year(:,1);
Bxyz = T_BRsphere_za_year(:,2:4);
Rsph = T_BRsphere_za_year(:,5:7); 

% lon, lat, alt = phi, theta, r = Azimuth, Elevation, High = a, e, h
%Time Bxyz Rsph

%%
units = irf_units;
RM = units.Mars.radius; % RM就是火星半径
Rsph(:,3) = Rsph(:,3) + RM/1e3; % h = h + RM
Rsph(:,2) = Rsph(:,2) * pi / 180; % 化为rad，且e = [-pi/2, pi/2]
Rsph(:,1) = (Rsph(:,1) - 180) * pi / 180; % 化为rad，且a = [-pi, pi]
Bsph = CoorTrans(Bxyz, Rsph); % 进行同样的变化，phi, the ,r = a, e, h

%将Rsph变换成h a e，Bxyz依照Rsph的规律转换成Bsph

%% divide grid
coor = 'aeh';
minRa = min(Rsph(:,1)); maxRa = max(Rsph(:,1));
minRe = min(Rsph(:,2)); maxRe = max(Rsph(:,2));
minRh = min(Rsph(:,3)); maxRh = max(Rsph(:,3));

gridNuma = 100; gridNume = 200; gridNumh = 5; 
% 若try中语句错误，则跳转到catch
try
load(['E:\Martian\programs\2_Maven_Download\10_pictures_J\B\', ...
    'Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat']);
catch
c_eval('grid? = linspace(minR?, maxR?, gridNum?);', coor);
% 即grida = linspace(minRa, maxRa, gridNuma);

[aGrid, eGrid, hGrid] = meshgrid(grida, gride, gridh); 
% [X, Y, Z] = meshgrid(x, y, z)生成网格，大小为length(y)*length(x)*length(z)

c_eval('?Grid = permute(?Grid, [2, 1, 3]);',coor); 
% aGrid = permute(aGrid, [2, 1, 3]);交换1、2维度

[~, ~, aIndex] = histcounts(Rsph(:,1), grida);
%以grida为边界，将Ra装进bin里，aIndex(i)表示Rsph(i, 1)落入了第几个bin里
[~, ~, eIndex] = histcounts(Rsph(:,2), gride);
[~, ~, hIndex] = histcounts(Rsph(:,3), gridh);

Bsph(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex, eIndex, hIndex);
% (aIndex, eIndex, hIndex)为每一行R对应网格矩阵的坐标，将坐标转换成序号值得到Bsph(:,4)

%% average B in grid
c_eval('B?Grid = nan(size(?Grid));', coor);
% BaGrid = nan(size(aGrid));将B网格初始值令为nan，大小同网格矩阵aGrid

% parpool(12)
parfor iGrid = 1:numel(aGrid)
% numel用于获取aGrid元素数量，iGrid从1开始遍历aGrid的所有序号
idx = Bsph(:,4) == iGrid;
% 找到被放在aGrid(iGrid)这一bin里的所有B
if any(idx)
tempB = mean(Bsph(idx, 1:3), 1);
% mean(Bsph(idx, 1:3), 1),求列均值，此处为求格子内B的均值
BaGrid(iGrid) = tempB(1);
BeGrid(iGrid) = tempB(2);
BhGrid(iGrid) = tempB(3);
% 令BaGrid(iGrid)为均值
end
end
% (aIndex, eIndex, hIndex)都是相对于grida、gride、gridh的空位而言的，
% 所以它们的最大值分别为gridNuma-1，gridNume-1，gridNumh-1。
% 所以在每一个矩阵的最后一行（列、页），都是不存在任何数据的。
% 现在更改函数，将矩阵缩小再进行尝试。

save(['E:\Martian\programs\2_Maven_Download\10_pictures_J\B\', ...
    'Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat'],...
    "aGrid", "eGrid", "hGrid", "BaGrid", "BeGrid", "BhGrid");
% 储存变量aGgrid eGrid hGrid（网格矩阵），BaGrid BeGrid BhGrid（对应的值）
end

%% calculate J
hGrid = hGrid .* 1e3; % m
c_eval('B?Grid = B?Grid .* 1e-9;', coor); % T
% 转换单位

[curlBa, curlBe, curlBh] = CurlSph(BaGrid, BeGrid, BhGrid, aGrid, eGrid, hGrid); % A/m^2
c_eval('J? = 1e9 * curlB? / units.mu0;', coor);% nA/m^2
%Ja = 1e9 * curlBa / units.mu0;求J

hGrid = hGrid / 1e3; %km
aGrid = aGrid * 180 / pi; %deg
eGrid = eGrid * 180 / pi; %deg
% 转换单位

savename_J=['E:\Martian\programs\2_Maven_Download\10_pictures_J\J\',...
    'Ja',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat'];
save(savename_J, "Ja", "Je", "Jh");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input data (sph)
inputFile = 'E:\Martian\programs\2_Maven_Download\10_pictures_J/T_BRsphere_za_year_2015_day.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year_day', 'var'), warning('wrong input data'); end
Time_day = T_BRsphere_za_year_day(:,1);
Bxyz_day = T_BRsphere_za_year_day(:,2:4);
Rsph_day = T_BRsphere_za_year_day(:,5:7); 

% lon, lat, alt = phi, theta, r = Azimuth, Elevation, High = a, e, h
%Time Bxyz Rsph

%%
units = irf_units;
RM = units.Mars.radius; % RM就是火星半径
Rsph_day(:,3) = Rsph_day(:,3) + RM/1e3; % h = h + RM
Rsph_day(:,2) = Rsph_day(:,2) * pi / 180; % 化为rad，且e = [-pi/2, pi/2]
Rsph_day(:,1) = (Rsph_day(:,1) - 180) * pi / 180; % 化为rad，且a = [-pi, pi]
Bsph_day = CoorTrans(Bxyz_day, Rsph_day); % 进行同样的变化，phi, the ,r = a, e, h

%将Rsph变换成h a e，Bxyz依照Rsph的规律转换成Bsph

%% divide grid
% coor = 'aeh';
minRa_day = min(Rsph_day(:,1)); maxRa_day = max(Rsph_day(:,1));
minRe_day = min(Rsph_day(:,2)); maxRe_day = max(Rsph_day(:,2));
minRh_day = min(Rsph_day(:,3)); maxRh_day = max(Rsph_day(:,3));

% gridNuma = 100; gridNume = 200; gridNumh = 8; 
% 若try中语句错误，则跳转到catch
try
load(['E:\Martian\programs\2_Maven_Download\10_pictures_J\B\', ...
    'Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'_day','.mat']);
catch
c_eval('grid? = linspace(minR?_day, maxR?_day, gridNum?);', coor);
% 即grida = linspace(minRa, maxRa, gridNuma);

[aGrid, eGrid, hGrid] = meshgrid(grida, gride, gridh); 
% [X, Y, Z] = meshgrid(x, y, z)生成网格，大小为length(y)*length(x)*length(z)

c_eval('?Grid = permute(?Grid, [2, 1, 3]);',coor); 
% aGrid = permute(aGrid, [2, 1, 3]);交换1、2维度

[~, ~, aIndex_day] = histcounts(Rsph_day(:,1), grida);
[~, ~, eIndex_day] = histcounts(Rsph_day(:,2), gride);
[~, ~, hIndex_day] = histcounts(Rsph_day(:,3), gridh);
%以grida为边界，将Ra装进bin里，aIndex(i)表示Rsph(i, 1)落入了第几个bin里

Bsph_day(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex_day, eIndex_day, hIndex_day);
% (aIndex, eIndex, hIndex)为每一行R对应网格矩阵的坐标，将坐标转换成序号值得到Bsph(:,4)

%% average B in grid
c_eval('B?Grid_day = nan(size(?Grid));', coor);
% BaGrid = nan(size(aGrid));将B网格初始值令为nan，大小同网格矩阵aGrid

% parpool(12)
parfor iGrid = 1:numel(aGrid)
% numel用于获取aGrid元素数量，iGrid从1开始遍历aGrid的所有序号
idx_day = Bsph_day(:,4) == iGrid;
% 找到被放在aGrid(iGrid)这一bin里的所有B
if any(idx_day)
tempB_day = mean(Bsph_day(idx_day, 1:3), 1);
% mean(Bsph(idx, 1:3), 1),求列均值，此处为求格子内B的均值
BaGrid_day(iGrid) = tempB_day(1);
BeGrid_day(iGrid) = tempB_day(2);
BhGrid_day(iGrid) = tempB_day(3);
% 令BaGrid(iGrid)为均值
end
end

save(['E:\Martian\programs\2_Maven_Download\10_pictures_J\B\', ...
    'Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh), '_day','.mat'],...
    "aGrid", "eGrid", "hGrid", "BaGrid_day", "BeGrid_day", "BhGrid_day");
% 储存变量aGgrid eGrid hGrid（网格矩阵），BaGrid BeGrid BhGrid（对应的值）
end

%% calculate J
hGrid = hGrid .* 1e3; % m
c_eval('B?Grid_day = B?Grid_day .* 1e-9;', coor); % T
% 转换单位

[curlBa_day, curlBe_day, curlBh_day] = CurlSph(BaGrid_day, BeGrid_day, BhGrid_day, aGrid, eGrid, hGrid); % A/m^2
c_eval('J?_day = 1e9 * curlB?_day / units.mu0;', coor);% nA/m^2
%Ja = 1e9 * curlBa / units.mu0;求J

hGrid = hGrid / 1e3; %km
aGrid = aGrid * 180 / pi; %deg
eGrid = eGrid * 180 / pi; %deg
% 转换单位

savename_J=['E:\Martian\programs\2_Maven_Download\10_pictures_J\J\',...
    'Ja',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh), '_day','.mat'];
save(savename_J, "Ja_day", "Je_day", "Jh_day");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input data (sph)
inputFile = 'E:\Martian\programs\2_Maven_Download\10_pictures_J/T_BRsphere_za_year_2015_nig.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year_nig', 'var'), warning('wrong input data'); end
Time_nig = T_BRsphere_za_year_nig(:,1);
Bxyz_nig = T_BRsphere_za_year_nig(:,2:4);
Rsph_nig = T_BRsphere_za_year_nig(:,5:7); 

% lon, lat, alt = phi, theta, r = Azimuth, Elevation, High = a, e, h
%Time Bxyz Rsph

%%
units = irf_units;
RM = units.Mars.radius; % RM就是火星半径
Rsph_nig(:,3) = Rsph_nig(:,3) + RM/1e3; % h = h + RM
Rsph_nig(:,2) = Rsph_nig(:,2) * pi / 180; % 化为rad，且e = [-pi/2, pi/2]
Rsph_nig(:,1) = (Rsph_nig(:,1) - 180) * pi / 180; % 化为rad，且a = [-pi, pi]
Bsph_nig = CoorTrans(Bxyz_nig, Rsph_nig); % 进行同样的变化，phi, the ,r = a, e, h

%将Rsph变换成h a e，Bxyz依照Rsph的规律转换成Bsph

%% divide grid
% coor = 'aeh';
minRa_nig = min(Rsph_nig(:,1)); maxRa_nig = max(Rsph_nig(:,1));
minRe_nig = min(Rsph_nig(:,2)); maxRe_nig = max(Rsph_nig(:,2));
minRh_nig = min(Rsph_nig(:,3)); maxRh_nig = max(Rsph_nig(:,3));

% gridNuma = 100; gridNume = 200; gridNumh = 8; 
% 若try中语句错误，则跳转到catch
try
load(['Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'_nig','.mat']);
catch
c_eval('grid? = linspace(minR?_nig, maxR?_nig, gridNum?);', coor);
% 即grida = linspace(minRa, maxRa, gridNuma);

% [aGrid, eGrid, hGrid] = meshgrid(grida, gride, gridh); 
% % [X, Y, Z] = meshgrid(x, y, z)生成网格，大小为length(y)*length(x)*length(z)

% c_eval('?Grid = permute(?Grid, [2, 1, 3]);',coor); 
% % aGrid = permute(aGrid, [2, 1, 3]);交换1、2维度

[~, ~, aIndex_nig] = histcounts(Rsph_nig(:,1), grida);
[~, ~, eIndex_nig] = histcounts(Rsph_nig(:,2), gride);
[~, ~, hIndex_nig] = histcounts(Rsph_nig(:,3), gridh);
%以grida为边界，将Ra装进bin里，aIndex(i)表示Rsph(i, 1)落入了第几个bin里

Bsph_nig(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex_nig, eIndex_nig, hIndex_nig);
% (aIndex, eIndex, hIndex)为每一行R对应网格矩阵的坐标，将坐标转换成序号值得到Bsph(:,4)

%% average B in grid
c_eval('B?Grid_nig = nan(size(?Grid));', coor);
% BaGrid = nan(size(aGrid));将B网格初始值令为nan，大小同网格矩阵aGrid

% parpool(12)
parfor iGrid = 1:numel(aGrid)
% numel用于获取aGrid元素数量，iGrid从1开始遍历aGrid的所有序号
idx_nig = Bsph_nig(:,4) == iGrid;
% 找到被放在aGrid(iGrid)这一bin里的所有B
if any(idx_nig)
tempB_nig = mean(Bsph_nig(idx_nig, 1:3), 1);
% mean(Bsph(idx, 1:3), 1),求列均值，此处为求格子内B的均值
BaGrid_nig(iGrid) = tempB_nig(1);
BeGrid_nig(iGrid) = tempB_nig(2);
BhGrid_nig(iGrid) = tempB_nig(3);
% 令BaGrid(iGrid)为均值
end
end

save(['E:\Martian\programs\2_Maven_Download\10_pictures_J\B\', ...
    'Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh), '_nig','.mat'],...
    "aGrid", "eGrid", "hGrid", "BaGrid_nig", "BeGrid_nig", "BhGrid_nig");
% 储存变量aGgrid eGrid hGrid（网格矩阵），BaGrid BeGrid BhGrid（对应的值）
end

%% calculate J
hGrid = hGrid .* 1e3; % m
c_eval('B?Grid_nig = B?Grid_nig .* 1e-9;', coor); % T
% 转换单位

[curlBa_nig, curlBe_nig, curlBh_nig] = CurlSph(BaGrid_nig, BeGrid_nig, BhGrid_nig, aGrid, eGrid, hGrid); % A/m^2
c_eval('J?_nig = 1e9 * curlB?_nig / units.mu0;', coor);% nA/m^2
%Ja = 1e9 * curlBa / units.mu0;求J

hGrid = hGrid / 1e3; %km
aGrid = aGrid * 180 / pi; %deg
eGrid = eGrid * 180 / pi; %deg
% 转换单位

savename_J=['E:\Martian\programs\2_Maven_Download\10_pictures_J\J\',...
    'Ja',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh), '_nig','.mat'];
save(savename_J, "Ja_nig", "Je_nig", "Jh_nig");

% %% Init figure
% set(0,'DefaultAxesFontSize',8);
% set(0,'DefaultLineLineWidth', 0.5);
% fn=figure(1);clf;
% set(gcf,'PaperUnits','centimeters')
% xSize = 200; ySize = 70; coef=floor(min(800/xSize,800/ySize));
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% 
% %% plot slide
% slide = 1;
% Jt_nig = sqrt(Ja_nig(:,:,slide).^2 + Je_nig(:,:,slide).^2 + Jh_nig(:,:,slide).^2);
% Bt_nig = sqrt(BaGrid_nig(:,:,slide).^2 + BeGrid_nig(:,:,slide).^2 + BhGrid_nig(:,:,slide).^2);
% % pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Jt);
% pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Ja_nig(:,:,slide));
% shading interp
% % colormap('jet'); 
% colormap(othercolor('BuDRd_12'))
% c = colorbar;
% c.Label.String = 'J_{h} [nA/m^2]';
% xlabel('Lon [\circ]'); ylabel('Lat [\circ]');
% hold on;
% 
% % quiver(aGrid(:,:,slide),eGrid(:,:,slide),Ja(:,:,slide),Je(:,:,slide),'w')
% 
% %%
% set(gca,'FontName','Times New Roman','FontSize',10);
% 
% set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
% set(gcf,'color','w')

%% coor Trans
function Bsph = CoorTrans(B, Rsph)
theta = pi/2 - Rsph(:,2); phi = Rsph(:,1);
Br = B(:,1).* sin(theta).* cos(phi) + B(:,2).* sin(theta).* sin(phi) + B(:,3).* cos(theta);
Bt = B(:,1).* cos(theta).* cos(phi) + B(:,2).* cos(theta).* sin(phi) - B(:,3).* sin(theta);
Bp = -B(:,1).* sin(phi) + B(:,2).* cos(phi);
Bsph = [Bp, Bt, Br];
end
