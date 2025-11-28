%% 初始化 
clear; clc; close all

% 初始化 MMS 本地数据库路径（按你自己的路径改）
mms.db_init('local_file_db','/Volumes/SPART-WORK/Data/MMS/');

ic = 1:4;   % MMS1-4
units = irf_units; 

speed_file = '/Users/fwd/Documents/Ti~mor~/M/ESW/Searching/DRAFTCaseList.txt';   % 改成你的 txt 文件名

fid = fopen(speed_file,'r','n','UTF-8');
C = textscan(fid,'%s%f', ...
    'Delimiter',{'的最大速度为','\n'}, ...
    'MultipleDelimsAsOne',true);
fclose(fid);

eventTimeStr = C{1};  % 字符串时间
vmax         = C{2};  % 最大速度 (数值)

% 按速度降序排序
[v_sorted, idx_sorted] = sort(vmax,'descend');
eventTimeStr_sorted    = eventTimeStr(idx_sorted);

Nevents = min(100, numel(v_sorted));   % 只取前 10 个

%%
for iev = 1:Nevents
    clear low_FGM1 low_FGM2 low_FGM3 low_FGM4 high_SCM1 high_SCM2 high_SCM3 high_SCM4...
        B_brst1 B_brst2 B_brst3 B_brst4
    fprintf('绘制事件 %d / %d: %s, vmax = %.1f\n', ...
        iev, Nevents, eventTimeStr_sorted{iev}, v_sorted(iev));

    %% 2.1 事件时间 ±1 s 的时间段
    t0_dt = datetime(eventTimeStr_sorted{iev}, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS''Z''', ...
        'TimeZone','UTC');

    t_start = t0_dt - seconds(0.2);
    t_end   = t0_dt + seconds(0.2);

    % irf.tint 接受的时间字符串（到毫秒就够用）
    t_start_str = datestr(t_start,'yyyy-mm-ddTHH:MM:SS.FFFZ');
    t_end_str   = datestr(t_end,  'yyyy-mm-ddTHH:MM:SS.FFFZ');

    tint = irf.tint([t_start_str '/' t_end_str]);

    %% 2.2 读取 4 颗卫星的 B_gsm 和 E_gsm
    try
        % 磁场：FGM GSM 三分量
        c_eval('Bxyz?  = mms.db_get_ts(''mms?_fgm_brst_l2'', ''mms?_fgm_b_gsm_brst_l2'', tint);', ic);
        c_eval('B_brst0? = irf.ts2mat(Bxyz?);', ic);   % -> [t, Bx, By, Bz]

        % 电场：EDP GSE -> GSM
        c_eval('Exyz?  = mms.db_get_ts(''mms?_edp_brst_l2_dce'', ''mms?_edp_dce_gse_brst_l2'', tint);', ic);
        c_eval('E_brst? = irf.ts2mat(Exyz?);', ic);   % [t, Ex, Ey, Ez] (GSE)
        c_eval('E_brst? = irf_gse2gsm(E_brst?);', ic);% -> GSM

        c_eval('Bscm_gse?=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',ic);
        try   %scm经典Bug,可能读到cell格式
            c_eval('B_scm_gse?=irf.ts2mat(Bscm_gse?);',ic);
        catch
            c_eval('B_scm1_gse?=irf.ts2mat(Bscm_gse?{1});',ic);
            c_eval('B_scm2_gse?=irf.ts2mat(Bscm_gse?{2});',ic);
            c_eval('B_scm_gse?=[B_scm1_gse?;B_scm2_gse?];',ic);
        end
        c_eval('B_scm?=irf_gse2gsm(B_scm_gse?);',ic);

        Pos = mms.get_data('R_gsm',tint);
        R1 = Pos.gsmR1;R2 = Pos.gsmR2; R3 = Pos.gsmR3;R4 = Pos.gsmR4;
        RR = [R1(1,1:3);R2(1,1:3);R3(1,1:3);R4(1,1:3)];
        R = R1(1,1:3)./units.RE*1e3;
        RR_mean = mean(pdist(RR));


    catch ME
        warning('事件 %d (%s) 数据读取失败: %s', ...
            iev, eventTimeStr_sorted{iev}, ME.message);
        continue;
    end

    %% 处理磁场数据

    % filter FGM <10Hz
    Fs=128;    % 采样频率 (Hz)
    Fc=10;  % 截止频率 (Hz)
    filterOrder = 5;    % 滤波器阶数
    [b2, a2] = butter(filterOrder, Fc/(Fs/2), 'low');  % 设计滤波器
    c_eval('low_FGM?(:,1) = B_brst0?(:,1);',ic);
    c_eval('low_FGM?(:,2) = filtfilt(b2, a2, B_brst0?(:,2));',ic);
    c_eval('low_FGM?(:,3) = filtfilt(b2, a2, B_brst0?(:,3));',ic);
    c_eval('low_FGM?(:,4) = filtfilt(b2, a2, B_brst0?(:,4));',ic);

    % filter SCM >10Hz
    Fs=8192;    % 采样频率 (Hz)
    Fc=10;  % 截止频率 (Hz)
    filterOrder = 5;    % 滤波器阶数
    [b2, a2] = butter(filterOrder, Fc/(Fs/2), 'high');  % 设计滤波器
    c_eval('high_SCM?(:,1) = B_scm?(:,1);',ic);
    c_eval('high_SCM?(:,2) = filtfilt(b2, a2, B_scm?(:,2));',ic);
    c_eval('high_SCM?(:,3) = filtfilt(b2, a2, B_scm?(:,3));',ic);
    c_eval('high_SCM?(:,4) = filtfilt(b2, a2, B_scm?(:,4));',ic);


    c_eval('B_re?=irf_resamp(low_FGM?,high_SCM?);',ic);
    c_eval('B_brst?(:,1)=B_re?(:,1);',ic);
    c_eval('B_brst?(:,2:4)=B_re?(:,2:4)+high_SCM?(:,2:4);',ic);
    %% 2.3 画图：6 个 panel (Bx, By, Bz, Ex, Ey, Ez)
    nplot = 6;
    h = irf_plot(nplot,'newfigure');

    set(gcf,'PaperUnits','centimeters')
    xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

    xwidth = 0.8;
    ywidth = (0.95 - 0.15)/nplot;
    mmscolors = [0 0 0; 1 0 0; 0 1 0; 0 0 1];  % MMS1-4: 黑红绿蓝
    set(h,'ColorOrder',mmscolors);
    c_eval('set(h(?), ''position'', [0.15 0.960-?*(ywidth+0.005) xwidth ywidth]);', 1:nplot);
    lnwid = 1.2;

    % --------- Bx ----------
    i = 1;
    h(i) = irf_panel('B_GSM_x');
    hold(h(i),'on');
    irf_plot(h(i), [B_brst1(:,1) B_brst1(:,2)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst2(:,1) B_brst2(:,2)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst3(:,1) B_brst3(:,2)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst4(:,1) B_brst4(:,2)], 'LineWidth', lnwid);
    hold(h(i),'off');
    ylabel(h(i), {'B_{x} (nT)'}, 'Interpreter','tex','FontSize',10);
    irf_legend(h(i), {'MMS1','MMS2','MMS3','MMS4'}, [0.02 0.98], 'FontSize',10);
    set(h(i),'xgrid','off','ygrid','off');

    % --------- By ----------
    i = i+1;
    h(i) = irf_panel('B_GSM_y');
    hold(h(i),'on');
    irf_plot(h(i), [B_brst1(:,1) B_brst1(:,3)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst2(:,1) B_brst2(:,3)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst3(:,1) B_brst3(:,3)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst4(:,1) B_brst4(:,3)], 'LineWidth', lnwid);
    hold(h(i),'off');
    ylabel(h(i), {'B_{y} (nT)'}, 'Interpreter','tex','FontSize',10);
    set(h(i),'xgrid','off','ygrid','off');

    % --------- Bz ----------
    i = i+1;
    h(i) = irf_panel('B_GSM_z');
    hold(h(i),'on');
    irf_plot(h(i), [B_brst1(:,1) B_brst1(:,4)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst2(:,1) B_brst2(:,4)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst3(:,1) B_brst3(:,4)], 'LineWidth', lnwid);
    irf_plot(h(i), [B_brst4(:,1) B_brst4(:,4)], 'LineWidth', lnwid);
    hold(h(i),'off');
    ylabel(h(i), {'B_{z} (nT)'}, 'Interpreter','tex','FontSize',10);
    set(h(i),'xgrid','off','ygrid','off');

    % --------- Ex ----------
    i = i+1;
    h(i) = irf_panel('E_GSM_x');
    hold(h(i),'on');
    irf_plot(h(i), [E_brst1(:,1) E_brst1(:,2)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst2(:,1) E_brst2(:,2)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst3(:,1) E_brst3(:,2)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst4(:,1) E_brst4(:,2)], 'LineWidth', lnwid);
    hold(h(i),'off');
    ylabel(h(i), {'E_{x} (mV/m)'}, 'Interpreter','tex','FontSize',10);
    set(h(i),'xgrid','off','ygrid','off');

    % --------- Ey ----------
    i = i+1;
    h(i) = irf_panel('E_GSM_y');
    hold(h(i),'on');
    irf_plot(h(i), [E_brst1(:,1) E_brst1(:,3)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst2(:,1) E_brst2(:,3)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst3(:,1) E_brst3(:,3)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst4(:,1) E_brst4(:,3)], 'LineWidth', lnwid);
    hold(h(i),'off');
    ylabel(h(i), {'E_{y} (mV/m)'}, 'Interpreter','tex','FontSize',10);
    set(h(i),'xgrid','off','ygrid','off');

    % --------- Ez ----------
    i = i+1;
    h(i) = irf_panel('E_GSM_z');
    hold(h(i),'on');
    irf_plot(h(i), [E_brst1(:,1) E_brst1(:,4)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst2(:,1) E_brst2(:,4)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst3(:,1) E_brst3(:,4)], 'LineWidth', lnwid);
    irf_plot(h(i), [E_brst4(:,1) E_brst4(:,4)], 'LineWidth', lnwid);
    hold(h(i),'off');
    ylabel(h(i), {'E_{z} (mV/m)'}, 'Interpreter','tex','FontSize',10);
    set(h(i),'xgrid','off','ygrid','off');

    %% 2.4 对齐坐标轴 & 缩放到 ±1 s
    Tint_plot = tint;
    irf_plot_axis_align(h(1:nplot));
    irf_zoom(h(1:nplot), 'x', Tint_plot);
    set(h(1:nplot),'FontSize',10);
    grid(h(1:nplot),'off');

    % 给第一行加个标题（可以按需打开/修改）
    title(h(1), sprintf('Top %d  Vmax=%.1f km/s (= %.1f c)  @ [%.1f %.1f %.1f] Re  Separation: %.1f km', ...
        iev, v_sorted(iev)/1e3, v_sorted(iev)/units.c, R(1), R(2), R(3), RR_mean), 'Interpreter','none');

    % 底部 x 轴刻度不旋转
    set(gca,"XTickLabelRotation",0);
    set(gcf,'render','painters');
    set(gcf,'paperpositionmode','auto')

    figname = ['/Users/fwd/Documents/Ti~mor~/M/ESW/Searching/Figures/', strrep(char(t0_dt),':','')];
    print(gcf, '-dpdf', [figname '.pdf']); 

    close all
end
