%------written by Wending Fu, Nov.2023 in Beijing------------
clear;clc;close all
% 18.0969
% file_path = '/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/aurora/PS.APL_V0105S027CB0006_SC.U_DI.A_GP.F18-SSUSI_PA.APL-EDR-AURORA_DD.20190805_SN.027894-00_DF.NC';%1808
% file_path = '/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/aurora/PS.APL_V0105S027CB0006_SC.U_DI.A_GP.F18-SSUSI_PA.APL-EDR-AURORA_DD.20190805_SN.027893-00_DF.NC';
% file_path = '/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/aurora/dmspf17_ssusi_edr-aurora_2019217T150651-2019217T164842-REV65793_vA8.2.0r000.nc';%%1600
file_path = '/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/aurora/dmspf18_ssusi_edr-aurora_2019217T153253-2019217T171445-REV27893_vA8.2.0r000.nc'; %%1624
% file_path = '/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/aurora/dmspf17_ssusi_edr-aurora_2019217T164847-2019217T182424-REV65794_vA8.2.0r000.nc';%%1800
info = ncinfo(file_path);
ncid = netcdf.open(file_path, 'NOWRITE');

for i=1:length(info.Variables)
variable_name = info.Variables(i).Name;
% eval([variable_name,' = netcdf.inqVar(',num2str(ncid),',', num2str(i-1),');']);
eval([variable_name,'.data = ncread(file_path,variable_name);'])

% dimensions = netcdf.inqVar(ncid, variable_data);
% eval([variable_name, '.dimensions = ncreadatt(file_path, variable_name,
% ''_FillValue'');']);i
eval([variable_name, '.attributes = ncinfo(file_path, variable_name);']);
end
netcdf.close(ncid);
%% data
lon_NORTH_mesh = LONGITUDE_GEOMAGNETIC_NORTH_GRID_MAP.data;
lon_SOUTH_mesh = LONGITUDE_GEOMAGNETIC_SOUTH_GRID_MAP.data;
lat_mesh = LATITUDE_GEOMAGNETIC_GRID_MAP.data;

NORTH_DQI = DISK_DQI_NORTH.data;
SOUTH_DQI = DISK_DQI_SOUTH.data;
NORTH_DQI(NORTH_DQI ~= 2048 | 4096) = 0; NORTH_DQI(NORTH_DQI == 2048 | 4096) = 1;
SOUTH_DQI(SOUTH_DQI ~= 2048 | 4096) = 0; SOUTH_DQI(SOUTH_DQI == 2048 | 4096) = 1;

ARC_NORTH_data = DISK_RADIANCEDATA_INTENSITY_NORTH.data;
ARC_NORTH_data = ARC_NORTH_data .* NORTH_DQI;
ARC_NORTH_data = ARC_NORTH_data(:,:,5);

ARC_SOUTH_data = DISK_RADIANCEDATA_INTENSITY_SOUTH.data;
ARC_SOUTH_data = ARC_SOUTH_data .* SOUTH_DQI;
ARC_SOUTH_data = ARC_SOUTH_data(:,:,5);
% ARC_NORTH_data = sum(ARC_NORTH_data,3);

% lon_NORTH_mesh = mod(lon_NORTH_mesh, 360);
% lon_SOUTH_mesh = mod(lon_SOUTH_mesh, 360);

% lon_index = find(lon_NORTH_mesh(1, :) >= 60 & lon_NORTH_mesh(1, :) <= 120);
% lon_NORTH_mesh = lon_NORTH_mesh(:, lon_index);
% lon_SOUTH_mesh = lon_SOUTH_mesh(:, lon_index);
% ARC_NORTH_data = ARC_NORTH_data(:, lon_index);
% ARC_SOUTH_data = ARC_SOUTH_data(:, lon_index);
% lat_mesh = lat_mesh(:, lon_index);

Energy_flux_south_data = ENERGY_FLUX_SOUTH_MAP.data;
Mean_energy_south_data = ELECTRON_MEAN_NORTH_ENERGY_MAP.data;
Mean_energy_south_data(Energy_flux_south_data<80) = nan;
%% plot
tt = UT_S.data;
tt(tt==0) = nan;
tt_0 = ones(size(tt));
tt = datetime(tt_0*2017,tt_0*7,tt_0*16,floor(tt),floor(100*(tt-floor(tt))),floor(100*(100*tt-floor(100*tt))));
lon_SOUTH_mesh(lon_SOUTH_mesh==0) = nan;
mlt = nan*zeros(size(tt));
for row = 1:size(tt,1)
    for pan = 1:size(tt,2)
        [mlt(row,pan),~,~] = magneticLocalTime(tt(row,pan),lon_SOUTH_mesh(row,pan));
    end
end

theta = mlt*15/180*pi;
rho = 90-lat_mesh;

figure;
polarscatter(theta(:), rho(:), 15, log10(abs(ARC_SOUTH_data(:))), 'filled'); hold on;
% polarscatter(theta(:), rho(:), 15, abs(Mean_energy_south_data(:)), 'filled');
% polarscatter(theta(:), rho(:), 15, abs(Energy_flux_south_data(:)), 'filled');
polarscatter(40.33*pi/180,90-74.21,150,'w','pentagram','filled')
ax=gca;
% ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'bottom';
hold on;


% % % aaa = [189,233,238]/255;
% % % bbb = [23,129,171]/255;
% % % ccc = [60,70,108]/255; % 
% % % cm = zeros(256,3);
% % % cm(:,1) = [linspace(ccc(1),bbb(1),128),linspace(bbb(1),aaa(1),128)];
% % % cm(:,2) = [linspace(ccc(2),bbb(2),128),linspace(bbb(2),aaa(2),128)];
% % % cm(:,3) = [linspace(ccc(3),bbb(3),128),linspace(bbb(3),aaa(3),128)];
% % % colormap(cm)

load('/Users/fwd/Documents/MATLAB/Code/zwz/color_rwb.mat');
colormap(color_rwb)

c = colorbar;
c.Label.String = 'SSUSI LBHS Log10(Intensity) [R]';
clim([2.5 4]);
rlim([0 30]);
rticks = {''};

title('Aurora Intensity in Polar Coordinates');
%%
% % % figure;
% % % subplot(2, 1, 1, polaraxes);
% % % pcolor(lon_NORTH_mesh, lat_mesh, ARC_NORTH_data);
% % % % polarplot(deg2rad(lon_NORTH_mesh), ARC_NORTH_data, 'LineWidth', 2);
% % % shading interp;
% % % xlabel('经度 (磁地方时)');
% % % ylabel('纬度');
% % % clim([0, 1500]);
% % % colormap('hot');
% % % title('北极区极光强度分布图');
% % % colorbar;
% % % 
% % % % 绘制南极区图
% % % subplot(2, 1, 2, polaraxes);
% % % pcolor(deg2rad(lon_SOUTH_mesh), deg2rad(lat_mesh), ARC_SOUTH_data);
% % % % polarplot(deg2rad(lon_SOUTH_mesh), ARC_SOUTH_data, 'LineWidth', 2);
% % % shading interp;
% % % xlabel('经度 (磁地方时)');
% % % ylabel('纬度');
% % % clim([0, 500]);
% % % colormap('hot');
% % % title('南极区极光强度分布图');
% % % colorbar;
