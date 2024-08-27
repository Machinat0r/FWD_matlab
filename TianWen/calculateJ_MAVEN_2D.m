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
%% input data (sph, abolish)
% % % inputFile = '/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/TianWen/T_BRsphere_za_year_2015.mat';
% % % load(inputFile);
% % % if ~exist('T_BRsphere_za_year', 'var'), warning('wrong input data'); end
% % % Time = T_BRsphere_za_year(:,1);
% % % Bxyz = T_BRsphere_za_year(:,2:4);
% % % Rsph = T_BRsphere_za_year(:,5:7);
% % % 
% % % Rsph(:,3) = Rsph(:,3) + 3396;
% % % Rsph(:,1:2) = Rsph(:,1:2) / 180 * pi;
% % % [Rx, Ry, Rz] = sph2cart(Rsph(:,1),Rsph(:,2),Rsph(:,3)); % azi=lon, ele=lat, r=alt+Rm
% % % Rxyz = [Rx, Ry, Rz];

%% input data (cart, replace the path)
try
load('BRcart_year.mat')
catch
warning('No whole-year data, combining and saving data in local folder ( つ•̀ω•́)つ');
inputDirPath = '/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/TianWen/2015/';
inputDir = dir(inputDirPath);
BRcart_year = [];
for i = 1:length(inputDir)
if contains(inputDir(i).name, '.mat')
    tempFilePath = [inputDirPath, inputDir(i).name];
    tempFile = load(tempFilePath);
    BRcart_year = [BRcart_year; tempFile.T_B_R_za];
end
end

if isempty(BRcart_year), warning('Please replace the path! (ง •̀_•́)ง‼');end
save('BRcart_year.mat')
end

Time = BRcart_year(:,1);
Bxyz = BRcart_year(:,2:4);
Rxyz = BRcart_year(:,5:7);

%% divide grid
coor = 'xyz';
minRx = min(Rxyz(:,1)); maxRx = max(Rxyz(:,1));
minRy = min(Rxyz(:,2)); maxRy = max(Rxyz(:,2));
minRz = min(Rxyz(:,3)); maxRz = max(Rxyz(:,3));

gridNumx = 100; gridNumy = 100; gridNumz = 10; 
try
load(['Bx',num2str(gridNumx),'y',num2str(gridNumy),'z',num2str(gridNumz),'.mat']);
catch
c_eval('grid? = linspace(minR?, maxR?, gridNum?);', coor);
[xGrid, yGrid, zGrid] = meshgrid(gridx, gridy, gridz);

[~, ~, xIndex] = histcounts(Rxyz(:,1), gridx);
[~, ~, yIndex] = histcounts(Rxyz(:,2), gridy);
[~, ~, zIndex] = histcounts(Rxyz(:,3), gridz);

Bxyz(:,4) = sub2ind([gridNumx, gridNumy, gridNumz], xIndex, yIndex, zIndex);

%% average B in grid
c_eval('B?Grid = nan(size(?Grid));', coor);

% parpool(12)
parfor iGrid = 1:numel(xGrid)
idx = Bxyz(:,4) == iGrid;
if any(idx)
tempB = mean(Bxyz(idx, 1:3), 1);
BxGrid(iGrid) = tempB(1);
ByGrid(iGrid) = tempB(2);
BzGrid(iGrid) = tempB(3);
end
end

save(['Bx',num2str(gridNumx),'y',num2str(gridNumy),'z',num2str(gridNumz),'.mat'],...
    "xGrid", "yGrid", "zGrid", "BxGrid", "ByGrid", "BzGrid");
end
%% calculate J
units = irf_units;
RM = units.Mars.radius; 
c_eval('?Grid = ?Grid .* RM;', coor); % m
c_eval('B?Grid = B?Grid .* 1e-9;', coor); % T
[curlx,curly,curlz,cav] = curl(xGrid,yGrid,zGrid,BxGrid,ByGrid,BzGrid); % A/m^2
c_eval('J? = 1e9 * curl? / units.mu0;', coor);% nA/m^2
c_eval('?Grid = ?Grid / 1e3;', coor); %km
%% plot slide
slide = 1;
Jt = sqrt(Jx(:,:,slide).^2+Jy(:,:,slide).^2);
pcolor(xGrid(:,:,slide),yGrid(:,:,slide),Jt);
shading interp
colormap('jet'); 
c = colorbar;
c.Label.String = 'J [nA/m^2]';
xlabel('x [km]'); ylabel('y [km]');
hold on;

quiver(xGrid(:,:,slide),yGrid(:,:,slide),Jx(:,:,slide),Jy(:,:,slide),'w')
