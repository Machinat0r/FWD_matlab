clear; clc; close all
%%
% 路径按你的习惯
cd /Volumes/SPART-WORK/Data/Cluster/

ic = 1:4;

%% ====== 时间======
Tsta = '2001-03-31T14:00:00.000Z';
Tend = '2001-03-31T17:00:00.000Z';

% Tsta = '2002-10-04T01:03:00.000Z';
% Tend = '2002-10-04T01:13:00.000Z';
tint = [iso2epoch(Tsta) iso2epoch(Tend)];

%% ====== 1) 读取/下载数据======
try
    c_eval("caa_load_skip_interact('C?_CP_FGM_FULL',Tsta,Tend);", ic);
    c_eval("caa_load_skip_interact('C?_CP_AUX_POSGSE_1M',Tsta,Tend);", ic);

    % c_eval("caa_load_skip_interact('C?_PP_CIS',Tsta,Tend);", ic);
    % c_eval("caa_load_skip_interact('C?_PP_PEA',Tsta,Tend);", ic);


    caa_load_skip_interact('C1_CP_RAP_ESPCT6', Tsta, Tend);   % RAPID electrons
    caa_load_skip_interact('C1_CP_RAP_HSPCT',  Tsta, Tend);   % RAPID protons (high-energy ions proxy)

catch
    c_eval("caa_download(tint,'C?_CP_FGM_FULL');", ic);
    c_eval("caa_download(tint,'C?_CP_AUX_POSGSE_1M');", ic);

    c_eval("caa_download(tint,'C?_PP_CIS');", ic);
    c_eval("caa_download(tint,'C?_PP_PEA');", ic);

    caa_download(tint,'C1_CP_PEA_PITCH_SPIN_DEFlux');
    caa_download(tint,'C1_CP_CIS_HIA_HS_1D_PEF');


    caa_download(tint,'C1_CP_RAP_ESPCT6');
    caa_download(tint,'C1_CP_RAP_HSPCT');

    c_eval("caa_load_skip_interact('C?_CP_FGM_FULL',Tsta,Tend);", ic);
    c_eval("caa_load_skip_interact('C?_CP_AUX_POSGSE_1M',Tsta,Tend);", ic);
end


%% ====== 2) 取出变量 ======

% ---- B (FGM FULL, GSE) ----
c_eval("B?_gse = c_caa_var_get('B_vec_xyz_gse__C?_CP_FGM_FULL','mat');", ic);
c_eval("B?     = irf_abs(B?_gse);", ic);     % 生成第5列 |B|
c_eval("B?     = irf_resamp(B?, B1);", 2:4); % 与 C1 对齐采样
c_eval("B?_gse = B?;", ic);                  % 统一成 [t Bx By Bz |B|]

% ---- Density: Ne (PEACE PP) / Ni (CIS PP, proton density) ----
c_eval("Ne? = c_caa_var_get('N_e_den__C?_PP_PEA','mat');", 1); % [t Ne]
c_eval("Ni? = c_caa_var_get('N_p__C?_PP_CIS','mat');", 1);     % [t Ni]


% % 重采样到 B1 时间轴（若某星为空则跳过）
% c_eval("if ~isempty(Ne?), Ne? = irf_resamp(Ne?, B1); end", ic);
% c_eval("if ~isempty(Ni?), Ni? = irf_resamp(Ni?, B1); end", ic);

% ---- Ion bulk velocity Vi (CIS PP) ----
c_eval("Vi?_gse = c_caa_var_get('V_p_xyz_gse__C?_PP_CIS','mat');", 1); % [t Vx Vy Vz]
% c_eval("if ~isempty(Vi?_gse), Vi?_gse = irf_resamp(Vi?_gse, B1); end", 1);

% ---- Spectrograms (C1) ----
% 电子通量谱
eSpec_ts = dataobj('/Volumes/SPART-WORK/Data/Cluster/CAA/C1_CP_PEA_PITCH_SPIN_DEFlux/C1_CP_PEA_PITCH_SPIN_DEFlux__20010331_140000_20010331_170000_V170609.cdf');
% eSpec = get_variable(eSpec_ts, 'Electron_Dif_flux__C1_CP_RAP_ESPCT6');
% % % cisDir = '/Volumes/SPART-WORK/Data/Cluster/CAA/C1_CP_CIS-HIA_HS_1D_PEF';  % 你的本地目录
% % % cisCdf = caa_pick_cdf_by_tint(cisDir, tint);  % 自动选文件（不交互）
% % % iSpec_ts = dataobj(cisCdf);
eSpec_data = get_variable(eSpec_ts, 'Data__C1_CP_PEA_PITCH_SPIN_DEFlux');
DEF_tE = pea_collapse_to_tE(eSpec_data.data, 44);
DEF_tE(DEF_tE<0) = 0;

eSpec_t = get_variable(eSpec_ts,'time_tags__C1_CP_PEA_PITCH_SPIN_DEFlux');
eSpec_dimension = eSpec_data.DEPEND_2;
eSpec.DEPEND_0.data = eSpec_t.data;
eSpec.DEPEND_1.data = eSpec_dimension.data(1,:);
eSpec.data = DEF_tE;

% % 离子通量谱
% iSpec = c_caa_var_get('flux__C1_CP_CIS_HIA_HS_1D_PEF');

iSpec_t = c_caa_var_get('time_tags__C1_CP_CIS_HIA_HS_1D_PEF');
iSpec_data = c_caa_var_get('flux__C1_CP_CIS_HIA_HS_1D_PEF');
iSpec_dimension = c_caa_var_get('energy_table__C1_CP_CIS_HIA_HS_1D_PEF');
iSpec.DEPEND_0.data = iSpec_t.data;
iSpec.DEPEND_1.data = iSpec_dimension.data;
iSpec.data = iSpec_data.data;

%% ====== RAPID spectrograms (C1) ======
rESpec = [];  % RAPID electron
rPSpec = [];  % RAPID proton (high-energy ion proxy)

try
    % --- RAPID Electron: ESPCT6 ---
    rapE_struct = c_caa_var_get('Electron_Dif_flux__C1_CP_RAP_ESPCT6');
    rapE_time = c_caa_var_get('Time_tags__C1_CP_RAP_ESPCT6');
    E_e = c_caa_var_get('Dimension_E__C1_CP_RAP_ESPCT6');

    rESpec.DEPEND_0.data = rapE_time.data;
    rESpec.DEPEND_1.data = E_e.data;
    rESpec.data          = rapE_struct.data;
catch
    rESpec = [];
end

try
    % --- RAPID Proton: HSPCT ---
    rapP_struct = c_caa_var_get('Proton_Dif_flux__C1_CP_RAP_HSPCT');
    rapP_time = c_caa_var_get('Time_tags__C1_CP_RAP_HSPCT');
    E_p = c_caa_var_get('Dimension_E__C1_CP_RAP_HSPCT');
    
    temprapP = double(rapP_struct.data);
    temprapP(temprapP<0) = 0;
    rPSpec.DEPEND_0.data = rapP_time.data;
    rPSpec.DEPEND_1.data = double(E_p.data);
    rPSpec.data          = temprapP;

catch
    rPSpec = [];
end

%% ====== 3) 作图：8 panels ======
n = 10;
i = 1;

set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth',0.75);

figure(1); clf;

cor = {'k','r','g','b'};

% ---------- Panel 1: Btotal ----------
h(i) = irf_subplot(n,1,-i);
for sc = 1:4
    eval(sprintf("irf_plot([B%d(:,1) B%d(:,5)], 'color','%s','Linewidth',0.75); hold on;", sc, sc, cor{sc}));
end
hold off; grid off;
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
ylabel('|B| (nT)','fontsize',12);
i = i+1;

% ---------- Panel 2: Bx ----------
h(i) = irf_subplot(n,1,-i);
for sc = 1:4
    eval(sprintf("irf_plot([B%d_gse(:,1) B%d_gse(:,2)], 'color','%s','Linewidth',0.75); hold on;", sc, sc, cor{sc}));
end
hold off; grid off;
ylabel('Bx (nT)','fontsize',12);
i = i+1;

% ---------- Panel 3: By ----------
h(i) = irf_subplot(n,1,-i);
for sc = 1:4
    eval(sprintf("irf_plot([B%d_gse(:,1) B%d_gse(:,3)], 'color','%s','Linewidth',0.75); hold on;", sc, sc, cor{sc}));
end
hold off; grid off;
ylabel('By (nT)','fontsize',12);
i = i+1;

% ---------- Panel 4: Bz ----------
h(i) = irf_subplot(n,1,-i);
for sc = 1:4
    eval(sprintf("irf_plot([B%d_gse(:,1) B%d_gse(:,4)], 'color','%s','Linewidth',0.75); hold on;", sc, sc, cor{sc}));
end
hold off; grid off;
ylabel('Bz (nT)','fontsize',12);
i = i+1;

% ---------- Panel 5: Density (Ne solid, Ni dashed) ----------
h(i) = irf_subplot(n,1,-i);
for sc = 1
    % Ne
    eval(sprintf("if ~isempty(Ne%d), irf_plot([Ne%d(:,1) Ne%d(:,2)], 'color','%s','Linewidth',0.75); hold on; end", sc, sc, sc, 'b'));
end
for sc = 1
    % Ni dashed
    eval(sprintf("if ~isempty(Ni%d), irf_plot([Ni%d(:,1) Ni%d(:,2)], 'color','%s','Linewidth',0.75,'LineStyle','--'); hold on; end", sc, sc, sc, 'g'));
end
hold off; grid off;
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
ylabel('n (cm^{-3})','fontsize',12);
irf_legend(gca,{'Ne','Ni'},[0.97 0.92]);
i = i+1;

% ---------- Panel 6: Vi (Vx -, Vy --, Vz :) ----------
h(i) = irf_subplot(n,1,-i);
for sc = 1
    eval(sprintf("if ~isempty(Vi%d_gse), irf_plot([Vi%d_gse(:,1) Vi%d_gse(:,2)], 'color','%s','Linewidth',0.75); hold on; end", sc, sc, sc, 'b'));
    eval(sprintf("if ~isempty(Vi%d_gse), irf_plot([Vi%d_gse(:,1) Vi%d_gse(:,3)], 'color','%s','Linewidth',0.75); hold on; end", sc, sc, sc, 'g'));
    eval(sprintf("if ~isempty(Vi%d_gse), irf_plot([Vi%d_gse(:,1) Vi%d_gse(:,4)], 'color','%s','Linewidth',0.75); hold on; end", sc, sc, sc, 'r'));
end
hold off; grid off;
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0];[0 0 0]]);
ylabel('V_i (km/s)','fontsize',12);
irf_legend(gca,{'Vx','Vy','Vz'},[0.97 0.92]);
i = i+1;
ylim([-200 300]);

% ---------- Panel 7: Electron flux spectrogram (C1) ----------
h(i) = irf_subplot(n,1,-i);
if ~isempty(eSpec)
    try
        % irf_spectrogram(gca, eSpec, 'log');

        specrec_p_e=struct('t',eSpec.DEPEND_0.data);
        specrec_p_e.f=transpose(eSpec.DEPEND_1.data);%energy levels
        specrec_p_e.p=eSpec.data;%data matrix
        specrec_p_e.f_label='';
        specrec_p_e.p_label={' ','log10(keV/(cm^2 s sr keV))'};
        [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_e);

       set(h(i),'yscale','log');
        set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);

        ylabel('E_e (eV)','fontsize',10);
    catch
        text(0.05,0.5,'Electron spectrogram plot failed: check eSpec format','Units','normalized');
    end
else
    text(0.05,0.5,'No Electron_Dif_flux__C1_CP_RAP_ESPCT6','Units','normalized');
end
grid off;
i = i+1;

% ---------- Panel 8: Ion flux spectrogram (C1) ----------
h(i) = irf_subplot(n,1,-i);
if ~isempty(iSpec)
    try
        
        specrec_p_i=struct('t',iSpec.DEPEND_0.data);
        specrec_p_i.f=iSpec.DEPEND_1.data;%energy levels
        specrec_p_i.p=iSpec.data;%data matrix
        specrec_p_i.f_label='';
        specrec_p_i.p_label={' ','log10(keV/(cm^2 s sr keV))'};
        [h(i), hcb8]=irf_spectrogram(h(i),specrec_p_i);

        set(h(i),'yscale','log');
        set(h(i),'ytick',[1e1 1e2 1e3 1e4],'fontsize',9);

        ylabel('E_i (eV)','fontsize',10);
    catch
        text(0.05,0.5,'Ion spectrogram plot failed: check iSpec format','Units','normalized');
    end
else
    text(0.05,0.5,'No flux__C1_CP_CIS_HIA_HS_1D_PEF','Units','normalized');
end
grid off;

i = i+1;

% ---------- Panel 9: RAPID Electron flux spectrogram (C1) ----------
h(i) = irf_subplot(n,1,-i);
if ~isempty(rESpec)
    try
        specrec_re = struct('t', rESpec.DEPEND_0.data);
        specrec_re.f = transpose(rESpec.DEPEND_1.data);  % energy
        specrec_re.p = rESpec.data;
        specrec_re.f_label = '';
        specrec_re.p_label = {' ','log10(DifFlux)'};

        [h(i), hcb9] = irf_spectrogram(h(i), specrec_re);
        set(h(i),'yscale','log');
        set(h(i),'ytick',[1e1 1e2 1e3],'fontsize',9);
        ylabel('E_{e,RAPID} (keV)','fontsize',10);
    catch
        text(0.05,0.5,'RAPID electron spectrogram plot failed','Units','normalized');
    end
else
    text(0.05,0.5,'No RAPID ESPCT6 (Electron_Dif_flux)','Units','normalized');
end
grid off;
i = i+1;

% ---------- Panel 10: RAPID Proton flux spectrogram (C1) ----------
h(i) = irf_subplot(n,1,-i);
if ~isempty(rPSpec)
    try
        specrec_rp = struct('t', rPSpec.DEPEND_0.data);
        specrec_rp.f = transpose(rPSpec.DEPEND_1.data);
        specrec_rp.p = rPSpec.data;
        specrec_rp.f_label = '';
        specrec_rp.p_label = {' ','log10(DifFlux)'};

        [h(i), hcb10] = irf_spectrogram(h(i), specrec_rp);
        set(h(i),'yscale','log');
        set(h(i),'ytick',[1e1 1e2 1e3],'fontsize',9);
        ylabel('E_{p,RAPID} (keV)','fontsize',10);
    catch
        text(0.05,0.5,'RAPID proton spectrogram plot failed','Units','normalized');
    end
else
    text(0.05,0.5,'No RAPID HSPCT (Proton_Dif_flux)','Units','normalized');
end
grid off;


%% ====== 4) 统一时间轴、对齐、配色 ======
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;

set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)


%%
function DEF_tE = pea_collapse_to_tE(DEF, nE)
% DEF: 从 c_caa_var_get 读出的通量数据数组（可能是 3-D 或 4-D）
% nE : 能量通道数（length(energy)）
%
% 输出 DEF_tE: [Nt x nE]，角度维度已平均（忽略 NaN）

    sz = size(DEF);

    % 找能量维：与 nE 相等的那个维度
    eDim = find(sz == nE, 1, 'first');
    if isempty(eDim)
        error('找不到能量维：size(DEF) 中没有等于 nE=%d 的维度。请检查 nE 或 DEF 的维度。', nE);
    end

    % 把能量维 permute 到第2维，时间默认取第1维（若不是请你在外部先调整）
    nd = ndims(DEF);
    perm = 1:nd;
    perm([2 eDim]) = perm([eDim 2]);   % swap dim-2 and eDim
    DEFp = permute(DEF, perm);

    % 现在 DEFp 的维度应类似 [Nt, nE, ...angles...]
    % 对第3维及之后的所有维度做平均
    DEF_tE = DEFp;
    for k = nd:-1:3
        DEF_tE = nanmean(DEF_tE, k);
    end

    % squeeze 成 [Nt x nE]
    DEF_tE = squeeze(DEF_tE);
end
