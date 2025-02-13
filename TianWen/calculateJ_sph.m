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
inputFile = '/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/TianWen/T_BRsphere_za_year_2015.mat';
load(inputFile);
if ~exist('T_BRsphere_za_year', 'var'), warning('wrong input data'); end
Time = T_BRsphere_za_year(:,1);
Bxyz = T_BRsphere_za_year(:,2:4);
Rsph = T_BRsphere_za_year(:,5:7); % lon, lat, alt = Azimuth, Elevation, High = a, e, h
%%
units = irf_units;
RM = units.Mars.radius; 
Rsph(:,3) = Rsph(:,3) + RM/1e3; % h = h + RM
Rsph(:,2) = Rsph(:,2) * pi / 180; % e = [-pi/2, pi/2]
Rsph(:,1) = (Rsph(:,1) - 180) * pi / 180; % a = [-pi, pi]
Bsph = CoorTrans(Bxyz, Rsph); % phi, the ,r = a, e, h
%% divide grid
coor = 'aeh';
minRa = min(Rsph(:,1)); maxRa = max(Rsph(:,1));
minRe = min(Rsph(:,2)); maxRe = max(Rsph(:,2));
minRh = min(Rsph(:,3)); maxRh = max(Rsph(:,3));

gridNuma = 50; gridNume = 100; gridNumh = 8; 
try
load(['Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat']);
catch
c_eval('grid? = linspace(minR?, maxR?, gridNum?);', coor);
[aGrid, eGrid, hGrid] = meshgrid(grida, gride, gridh); 
c_eval('?Grid = permute(?Grid, [2, 1, 3]);',coor); 

[~, ~, aIndex] = histcounts(Rsph(:,1), grida);
[~, ~, eIndex] = histcounts(Rsph(:,2), gride);
[~, ~, hIndex] = histcounts(Rsph(:,3), gridh);

Bsph(:,4) = sub2ind([gridNuma, gridNume, gridNumh], aIndex, eIndex, hIndex);

%% average B in grid
c_eval('B?Grid = nan(size(?Grid));', coor);

% parpool(12)
parfor iGrid = 1:numel(aGrid)
idx = Bsph(:,4) == iGrid;
if any(idx)
tempB = mean(Bsph(idx, 1:3), 1);
BaGrid(iGrid) = tempB(1);
BeGrid(iGrid) = tempB(2);
BhGrid(iGrid) = tempB(3);
end
end

save(['Ba',num2str(gridNuma),'e',num2str(gridNume),'h',num2str(gridNumh),'.mat'],...
    "aGrid", "eGrid", "hGrid", "BaGrid", "BeGrid", "BhGrid", "gridh");
end

%% calculate J
hGrid = hGrid .* 1e3; % m
c_eval('B?Grid = B?Grid .* 1e-9;', coor); % T
[curlBa, curlBe, curlBh] = CurlSph(BaGrid, BeGrid, BhGrid, aGrid, eGrid, hGrid); % A/m^2
c_eval('J? = 1e9 * curlB? / units.mu0;', coor);% nA/m^2
hGrid = hGrid / 1e3; %km
aGrid = aGrid * 180 / pi; %deg
eGrid = eGrid * 180 / pi; %deg

%% plot slide
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 200; ySize = 90; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


slide = 1;
Jt = sqrt(Ja(:,:,slide).^2 + Je(:,:,slide).^2 + Jh(:,:,slide).^2);
Bt = sqrt(BaGrid(:,:,slide).^2 + BeGrid(:,:,slide).^2 + BhGrid(:,:,slide).^2);
% pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Jt);
pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Je(:,:,slide));
shading interp
clim([-800,800]);
colormap('jet'); 
c = colorbar;
c.Label.String = 'J_{elevation} [nA/m^2]';
xlabel('Lon [\circ]'); ylabel('Lat [\circ]');
hold on;

% quiver(aGrid(:,:,slide),eGrid(:,:,slide),Ja(:,:,slide),Je(:,:,slide),'w')
%% plot stack
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 100; ySize = 200; coef=floor(min(1000/xSize,1000/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])


Jt = nan(size(Je));
for slide = 1:gridNumh
hGrid(:,:,slide) = hGrid(:,:,slide) - RM/1e3;
Jt(:,:,slide) = sqrt(Ja(:,:,slide).^2 + Je(:,:,slide).^2 + Jh(:,:,slide).^2);
Bt = sqrt(BaGrid(:,:,slide).^2 + BeGrid(:,:,slide).^2 + BhGrid(:,:,slide).^2);
% pcolor(aGrid(:,:,slide),eGrid(:,:,slide),Jt);
s = surf(aGrid(:,:,slide),eGrid(:,:,slide),hGrid(:,:,slide),Jh(:,:,slide));hold on
s.EdgeColor = 'none';
shading interp
colormap('jet'); 
clim([-800,800]);
c = colorbar;
c.Label.String = 'J_{elevation} [nA/m^2]';
xlabel('Lon [\circ]'); ylabel('Lat [\circ]'); zlabel('Alt [km]')
hold on;

% quiver(aGrid(:,:,slide),eGrid(:,:,slide),Ja(:,:,slide),Je(:,:,slide),'w')
end
%%
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gcf,'color','w')

%% coor Trans
function Bsph = CoorTrans(B, Rsph)
theta = pi/2 - Rsph(:,2); phi = Rsph(:,1);
Br = B(:,1).* sin(theta).* cos(phi) + B(:,2).* sin(theta).* sin(phi) + B(:,3).* cos(theta);
Bt = B(:,1).* cos(theta).* cos(phi) + B(:,2).* cos(theta).* sin(phi) - B(:,3).* sin(theta);
Bp = -B(:,1).* sin(phi) + B(:,2).* cos(phi);
Bsph = [Bp, Bt, Br];
end