clear;clc;close all
cd '/Volumes/SPART-NAS/Data/Cluster/'
ParentDir = '/Volumes/SPART-NAS/Data/Cluster/';
ic=2;

% TT = '2002-03-17\2003-01-01';
TT = '2001-01-01\2006-01-01';
Datelist = regexp(TT,'\d+-\d+-\d+','match');
TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'/']; mkdir(TaskDir)

load([TaskDir,'BRdata.mat'], 'all_MLT','all_MLat','all_Bphi','all_L');

%% Counts
L_edges    = [1.5 2.5 3.5 4.5 5.5 6.5];
MLT_edges  = linspace(0.5,24.5,25);
MLat_edges = linspace(-60,60,31);

bin_L    = discretize(all_L,    L_edges, 'IncludedEdge','right');
bin_MLT  = discretize(all_MLT,  MLT_edges, 'IncludedEdge','right');
bin_MLat = discretize(all_MLat, MLat_edges, 'IncludedEdge','right');

valid_MLT  = ~isnan(bin_L)  & ~isnan(bin_MLT);
valid_MLat = ~isnan(bin_L)  & ~isnan(bin_MLat);

counts_MLT = accumarray(...
    [ bin_L(valid_MLT), bin_MLT(valid_MLT) ], ...
    1, ...
    [ numel(L_edges)-1, numel(MLT_edges)-1 ] ...
);

counts_MLat = accumarray(...
    [ bin_L(valid_MLat), bin_MLat(valid_MLat) ], ...
    1, ...
    [ numel(L_edges)-1, numel(MLat_edges)-1 ] ...
);

L_centers    = (L_edges(1:end-1)   + L_edges(2:end))   / 2;
MLT_centers  = (MLT_edges(1:end-1) + MLT_edges(2:end)) / 2;
MLat_centers = (MLat_edges(1:end-1)+ MLat_edges(2:end))/ 2;

%% Bphi
valid_Bphi_MLT  = ~isnan(bin_L) & ~isnan(bin_MLT) & ~isnan(all_Bphi);
valid_Bphi_MLat = ~isnan(bin_L) & ~isnan(bin_MLat) & ~isnan(all_Bphi);

sum_Bphi_MLT = accumarray(...
    [ bin_L(valid_Bphi_MLT), bin_MLT(valid_Bphi_MLT) ], ...
    all_Bphi(valid_Bphi_MLT), ...
    [ numel(L_edges)-1, numel(MLT_edges)-1 ], ...
    @sum, NaN);

count_Bphi_MLT = accumarray(...
    [ bin_L(valid_Bphi_MLT), bin_MLT(valid_Bphi_MLT) ], ...
    1, ...
    [ numel(L_edges)-1, numel(MLT_edges)-1 ], ...
    @sum, NaN);

mean_Bphi_MLT = sum_Bphi_MLT ./ count_Bphi_MLT;

sum_Bphi_MLat = accumarray(...
    [ bin_L(valid_Bphi_MLat), bin_MLat(valid_Bphi_MLat) ], ...
    all_Bphi(valid_Bphi_MLat), ...
    [ numel(L_edges)-1, numel(MLat_edges)-1 ], ...
    @sum, NaN);

count_Bphi_MLat = accumarray(...
    [ bin_L(valid_Bphi_MLat), bin_MLat(valid_Bphi_MLat) ], ...
    1, ...
    [ numel(L_edges)-1, numel(MLat_edges)-1 ], ...
    @sum, NaN);

mean_Bphi_MLat = sum_Bphi_MLat ./ count_Bphi_MLat;

mean_Bphi_MLT = [mean_Bphi_MLT, mean_Bphi_MLT(:,1)];
counts_MLT = [counts_MLT, counts_MLT(:,1)];
MLT_centers   = [MLT_centers, 25];

%% counts plot
max_counts = max(counts_MLat(:));
figure('Position',[100 100 1200 500]);
colormap(jet); 

% === 左：MLT–L 极坐标柱状图 ===
subplot(1,2,1);

[TH, R] = meshgrid(MLT_centers/24*2*pi, L_centers);
X = R .* cos(TH);
Y = R .* sin(TH);

pcolor(X', Y', counts_MLT');  
shading flat; 
axis equal off;
caxis([0 max_counts]);
% title('MLT–L');


hold on;
theta_line = linspace(0,2*pi,25);
for r = 1:6
    plot(r*cos(theta_line), r*sin(theta_line), 'k--', 'LineWidth',0.5);
end

text(0, L_edges(end),  '12','HorizontalAlignment','center');
text(L_edges(end), 0,   '06','HorizontalAlignment','center');
text(0, -L_edges(end),'00','HorizontalAlignment','center');
text(-L_edges(end),0,  '18','HorizontalAlignment','center');
% colorbar('Location','eastoutside');


% === 右：MLat–L 楔形柱状图 ===
subplot(1,2,2);
[TH2, R2] = meshgrid(deg2rad(MLat_centers), L_centers);
X2 = R2 .* cos(TH2);
Y2 = R2 .* sin(TH2);
pcolor(X2', Y2', counts_MLat');
shading flat;
axis equal off;
caxis([0 max_counts]);
% title('MLat–L');


hold on;
theta_wedge = deg2rad(linspace(-55,55,31));
for r = 1:6
    plot(r*cos(theta_wedge), r*sin(theta_wedge), 'k-', 'LineWidth',1);
end
plot([0 L_edges(end)*cos(deg2rad(-60))], [0 L_edges(end)*sin(deg2rad(-55))], 'k-', 'LineWidth',1);
plot([0 L_edges(end)*cos(deg2rad(60))],  [0 L_edges(end)*sin(deg2rad(55))],  'k-', 'LineWidth',1);
colorbar('Location','eastoutside');

%% Bphi plot
max_Bphi = 100;

figure('Position',[100 100 1200 500]);
n_half = 128;
neg = [linspace(0,1,n_half)' linspace(0,1,n_half)' ones(n_half,1)];
pos = [ones(n_half,1) linspace(1,0,n_half)' linspace(1,0,n_half)'];
cmap = [neg; pos];
colormap(cmap);

% === 左：MLT–L 极坐标图（Bphi） ===
subplot(1,2,1);
[TH, R] = meshgrid(MLT_centers/24*2*pi, L_centers);
X = R .* cos(TH);
Y = R .* sin(TH);

pcolor(X', Y', mean_Bphi_MLT');  
shading flat; 
axis equal off;
% title('B_\phi in MLT–L');
caxis([-max_Bphi max_Bphi]); % 对称颜色条
% colorbar('Location','eastoutside');

hold on;
theta_line = linspace(0,2*pi,25);
for r = 1:6
    plot(r*cos(theta_line), r*sin(theta_line), 'k--', 'LineWidth',0.5);
end
text(0, L_edges(end),  '12','HorizontalAlignment','center');
text(L_edges(end), 0,   '06','HorizontalAlignment','center');
text(0, -L_edges(end),'00','HorizontalAlignment','center');
text(-L_edges(end),0,  '18','HorizontalAlignment','center');


% === 右：MLat–L 楔形图（Bphi） ===
subplot(1,2,2);
[TH2, R2] = meshgrid(deg2rad(MLat_centers), L_centers);
X2 = R2 .* cos(TH2);
Y2 = R2 .* sin(TH2);
pcolor(X2', Y2', mean_Bphi_MLat');
shading flat;
axis equal off;
% title('B_\phi in MLat–L');
caxis([-max_Bphi max_Bphi]);
colorbar('Location','eastoutside');
colorbar('Ticks',[-max_Bphi 0 max_Bphi], 'TickLabels',{num2str(-max_Bphi),'0',num2str(max_Bphi)});

hold on;
theta_wedge = deg2rad(linspace(-55,55,31));
for r = 1:6
    plot(r*cos(theta_wedge), r*sin(theta_wedge), 'k-', 'LineWidth',1);
end
plot([0 L_edges(end)*cos(deg2rad(-60))], [0 L_edges(end)*sin(deg2rad(-45))], 'k-', 'LineWidth',1);
plot([0 L_edges(end)*cos(deg2rad(60))],  [0 L_edges(end)*sin(deg2rad(45))],  'k-', 'LineWidth',1);

