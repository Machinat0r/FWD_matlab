%% ======================= 0. 基本初始化 ==========================
clear; clc;

% 初始化 MMS 本地数据库路径（按你自己的路径改）
mms.db_init('local_file_db','Z:\SPART-WORK\Data\MMS\');

ic = 1:4;   % MMS1-4
units = irf_units; %#ok<NASGU>

%% ======================= 1. 读入并排序事件 ======================
% txt 每行格式类似：
% 2018-08-13T22:44:24.007000000Z的最大速度为28467097.3975
speed_file = 'C:\MMS\ESWSearch\2018-01-01To2019-01-01\DRAFTCaseList.txt';   % 改成你的 txt 文件名

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

Nevents = min(10, numel(v_sorted));   % 只取前 10 个

%% ======================= 2. 循环画 10 个事件 =====================
for iev = 1:Nevents
    fprintf('绘制事件 %d / %d: %s, vmax = %.1f\n', ...
        iev, Nevents, eventTimeStr_sorted{iev}, v_sorted(iev));

    %% 2.1 事件时间 ±1 s 的时间段
    t0_dt = datetime(eventTimeStr_sorted{iev}, ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSSS''Z''', ...
        'TimeZone','UTC');

    t_start = t0_dt - seconds(1);
    t_end   = t0_dt + seconds(1);

    % irf.tint 接受的时间字符串（到毫秒就够用）
    t_start_str = datestr(t_start,'yyyy-mm-ddTHH:MM:SS.FFFZ');
    t_end_str   = datestr(t_end,  'yyyy-mm-ddTHH:MM:SS.FFFZ');

    tint = irf.tint([t_start_str '/' t_end_str]);

    %% 2.2 读取 4 颗卫星的 B_gsm 和 E_gsm
    try
        % 磁场：FGM GSM 三分量
        c_eval('Bxyz?  = mms.db_get_ts(''mms?_fgm_brst_l2'', ''mms?_fgm_b_gsm_brst_l2'', tint);', ic);
        c_eval('B_brst? = irf.ts2mat(Bxyz?);', ic);   % -> [t, Bx, By, Bz]

        % 电场：EDP GSE -> GSM
        c_eval('Exyz?  = mms.db_get_ts(''mms?_edp_brst_l2_dce'', ''mms?_edp_dce_gse_brst_l2'', tint);', ic);
        c_eval('E_brst? = irf.ts2mat(Exyz?);', ic);   % [t, Ex, Ey, Ez] (GSE)
        c_eval('E_brst? = irf_gse2gsm(E_brst?);', ic);% -> GSM
    catch ME
        warning('事件 %d (%s) 数据读取失败: %s', ...
            iev, eventTimeStr_sorted{iev}, ME.message);
        continue;
    end

    %% 2.3 画图：6 个 panel (Bx, By, Bz, Ex, Ey, Ez)
    nplot = 6;
    h = irf_plot(nplot,'newfigure');
    set(gcf,'Position',[50 50 400 620],'color','w');

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
    title(h(1), sprintf('Top %d  v_{max}=%.1f km/s  @ %s', ...
        iev, v_sorted(iev)/1e3, char(t0_dt)), 'Interpreter','none');

    % 底部 x 轴刻度不旋转
    set(gca,"XTickLabelRotation",0);
end
