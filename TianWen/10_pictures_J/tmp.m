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

save(['Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh), '_nig','.mat'],...
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

savename_J=['E:\Martian\programs\2_Maven_Download\10_pictures_J\',...
    'Ja',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh), '_nig','.mat'];
save(savename_J, "Ja_nig", "Je_nig", "Jh_nig");
