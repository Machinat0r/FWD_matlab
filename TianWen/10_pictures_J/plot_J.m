clear; clc; close all;clf;

coor = 'a';% aeh
LA = {'', '_day', '_nig'};
LB = {'All', 'Day', 'Night'};
ha = tight_subplot(1, 3, [.1 .1], [.05 .05], [.1 .14]);

%% Set the paramiters
gridNuma = 720; gridNume = 300; gridNumh = 5; 
ParentDir = '/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/TianWen/10_pictures_J/';
load([ParentDir, 'J/Ja' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '.mat']);
load([ParentDir, 'J/Ja' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_day', '.mat']);
load([ParentDir, 'J/Ja' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_nig', '.mat']);

load([ParentDir, 'B/Ba' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '.mat']);
load([ParentDir, 'B/Ba' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_day', '.mat']);
load([ParentDir, 'B/Ba' ...
    , num2str(gridNuma), 'e', num2str(gridNume), 'h', num2str(gridNumh), '_nig', '.mat']);

n_slice = size(Ja, 3);

%%
units = irf_units;
RM = units.Mars.radius;
hGrid = hGrid - RM/1e3;
mina = min(aGrid, [], 'all'); maxa = max(aGrid, [], 'all');
mine = min(eGrid, [], 'all'); maxe = max(eGrid, [], 'all');
minH = round(min(hGrid, [], 'all')); maxH = round(max(hGrid, [], 'all'));
dH = linspace(minH, maxH, n_slice);

%% Plot
for k = 1 : 3
    axes(ha(k))
    
    for slice = 1 : n_slice - 1
        H = hGrid(1, 1, slice);
        z = H * ones(size(eGrid(:, :, 1)));

        % eval(['b = surf(aGrid(:, :, 1), eGrid(:, :, 1), z, B', coor, 'Grid', cell2mat(LA(k)), '(:,:,slice));'])
        eval(['b = surf(aGrid(:, :, 1), eGrid(:, :, 1), z, J', coor, cell2mat(LA(k)), '(:,:,slice));'])
        % b=surf(aGrid(:, :, 1), eGrid(:, :, 1), z, Ja(:, :, slice));

        grid off; box on;

        % set plot
        clim([-1000 1000])
        set(b,'linestyle','none');
        view(-30,28)
        hold on
        xl=xlabel('Longitude (°)','rotation',18);
        yl=ylabel(['La' ...
            'titude (°)'],'rotation',-39);
        zl=zlabel('Altitude (km)');
        set(gca,'FontName','Times New Roman','FontSize',5);
        xlim([mina maxa]); 
        ylim([mine maxe]); 
        zlim([minH dH(n_slice - 1)]); 
        zticks(dH);
        pbaspect([2 1 2])  % 长宽高之比
        set(gca,'linewidth',0.6)
    end
    tt=title(LB{k},'FontWeight','bold');
    tt.FontName='Times New Roman';
    hold on
end
set(gcf,'color','w');%
set(gcf,'render','painters');%
hbr=colorbar;
hbr.Label.String={'J[nA/m^2]'};
% hbr.Label.String={'Dgree'};
hbr.Label.FontName='Times New Roman';
hbr.Label.FontSize=8;
set(hbr,'position',[0.92,0.34,0.01,0.35]);
colormap(othercolor('BuDRd_12'))
% colormap('jet')

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