function displayAudit = Plot_V1_Interstellar_Daily_Overview(result,visible,includeS4)
% Display approved daily statistics through 2021-12-16 inclusive.
% Standalone runs read original CDFs again. Log axes mask display values
% only; original daily values stay unchanged. NaNs break thin lines.
if nargin<2
    visible = true;
end
if nargin<3
    includeS4 = false;
end
root=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(root),'Case1_PPT_VerticalLine_Events_7d_fu'));
if nargin<1 || isempty(result)
    if includeS4
        result = Run_V1_Interstellar_Daily_WithS4(visible);
    else
        result = Run_V1_Interstellar_Daily_Overview('SyncArchive',false,'Visible',visible);
    end
    displayAudit = result.DisplayAudit;
    return
end
%% 时间范围与面板数据（来自本次原始 CDF 读取）
out = char(result.OutputFolder);
startUTC=datetime(2012,8,25,'TimeZone','UTC');
stopUTC=datetime(2021,12,17,'TimeZone','UTC');
daily=result.Daily;
use=daily.EpochUTC>=startUTC & daily.EpochUTC<stopUTC;
daily=daily(use,:);
sectors = [1 2 3 5 6 7];
sumField = 'P1SixSectorSum';
panelDLabel = 'P1 sum: S1,S2,S3,S5,S6,S7';
fileStem = 'V1_interstellar_daily_4panels';
if includeS4
    sectors = 1:7;
    sumField = 'P1SevenSectorSum';
    panelDLabel = 'P1 sum: S1-S7 (including S4)';
    fileStem = 'V1_interstellar_daily_4panels_with_S4';
end
v = [daily.BMean_nT daily.P1Mean daily.P1Median daily.(sumField)];
shown = v;
shown(:,2:4) = maskLog(shown(:,2:4));
%% 绘图
visibility = 'off';
if visible
    visibility = 'on';
end
fig=figure('Color','w','Position',[80 80 1600 1050],'Visible',visibility);
layout=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
labels={{'|B| daily mean','(nT)'},{'P1 daily mean','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'P1 daily median','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {panelDLabel,'(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}};
colors=[0.15 0.15 0.15;0.48 0.12 0.62;0.48 0.12 0.62;0.08 0.36 0.62];
ax = gobjects(4,1);
lines = gobjects(4,1);
for k=1:4
    ax(k)=nexttile(layout);
    lines(k)=plot(ax(k),daily.EpochUTC,shown(:,k),'.-', ...
        'Color',colors(k,:),'MarkerSize',3,'LineWidth',0.4);
    if k>1
        set(ax(k),'YScale','log');
    end
    ylabel(ax(k),labels{k});
    set(ax(k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',0.12);
    xlim(ax(k),[startUTC stopUTC]);
    text(ax(k),0.008,0.89,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
    if k<4
        ax(k).XTickLabel = [];
    end
end
linkaxes(ax,'x');
xlabel(ax(4),'UTC');
xtickformat(ax(4),'yyyy');
title(layout,'Voyager 1 beyond the heliopause | 2012-08-25 to 2021-12-16 | P1 0.57-1.78 MeV', ...
    'FontWeight','bold','FontSize',15);
%% 检查图形对象并保存结果
for k=1:4
    assert(isequaln(lines(k).YData(:),shown(:,k)));
    assert(strcmp(lines(k).LineStyle,'-') && lines(k).LineWidth==0.4);
    assert(isequal(ax(k).XLim,[startUTC stopUTC]));
    if k>1
        assert(strcmp(ax(k).YScale,'log'));
    end
end
files = string(fullfile(out,{[fileStem,'.png'],[fileStem,'.pdf'],[fileStem,'.fig']}));
exportgraphics(fig,files(1),'Resolution',220);
exportgraphics(fig,files(2),'ContentType','vector');
savefig(fig,files(3));
displayAudit=struct('CreatedUTC',datetime('now','TimeZone','UTC'), ...
    'StartUTC',startUTC,'EndUTCExclusive',stopUTC,'DailyRows',height(daily), ...
    'PanelDSectors',sectors,'PanelDSumField',sumField, ...
    'LineWidth',0.4,'AxesScales',{{'linear','log','log','log'}}, ...
    'NonpositiveHiddenByPanel',sum(isfinite(v)&v<=0,1).*[0 1 1 1], ...
    'MissingDisplayRule','NaN breaks lines; no connection across missing or nonpositive log-axis values.', ...
    'SourceAudit',string(fullfile(out,'V1_daily_overview_audit.mat')), ...
    'PlotCodeSHA256',string(Case1_File_SHA256(mfilename('fullpath')+".m")), ...
    'OutputFiles',files,'GraphicsChecksPassed',true);
save(fullfile(out,'V1_daily_display_audit.mat'),'displayAudit');
disp(displayAudit);
end

function x=maskLog(x)
x(x<=0)=NaN;
end



