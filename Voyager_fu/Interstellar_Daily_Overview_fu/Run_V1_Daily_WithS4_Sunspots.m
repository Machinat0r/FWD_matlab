function result = Run_V1_Daily_WithS4_Sunspots(visible)
% 原始 SILSO 日黑子数 + 原始 Voyager CDF 日统计，五面板。
% 黑子数 -1 转 NaN，零保留；无平滑、插值或传播时间平移。
if nargin<1, visible=true; end
%% 路径
CodeDir = 'C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/';
DataDir = 'Z:/SPART-WORK/Data/Voyager/';
SolarFile = 'Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
OutputDir = 'C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/with_S4_sunspots/';
addpath([CodeDir,'Interstellar_Daily_Overview_fu']);
addpath([CodeDir,'Case1_PPT_VerticalLine_Events_7d_fu']);
if ~isfolder(OutputDir), mkdir(OutputDir); end
startUTC=datetime(2012,8,25,'TimeZone','UTC');
stopUTC=datetime(2021,12,17,'TimeZone','UTC');
%% 读取原始发布日黑子数（CSV 为 SILSO 原始发布格式）
x=readmatrix(SolarFile,'Delimiter',';','FileType','text');
t=datetime(x(:,1),x(:,2),x(:,3),'TimeZone','UTC')+hours(12);
use=t>=startUTC & t<stopUTC;
x=x(use,:); t=t(use); sourceRows=find(use);
assert(numel(unique(t))==numel(t),'Duplicate SILSO dates.');
sn=x(:,5); sn(sn==-1)=NaN;
assert(all(isnan(sn)|sn>=0));
%% Voyager 直接从原始 CDF 读取，保持已有日统计和 L1 优先规则
source=Run_V1_Interstellar_Daily_Overview('DataRoot',DataDir, ...
    'OutputRoot',fullfile(OutputDir,'source_recompute'), ...
    'SyncArchive',false,'Visible',false,'MakePlot',false);
daily=source.Daily;
daily=daily(daily.EpochUTC>=startUTC & daily.EpochUTC<stopUTC,:);
seven=sum(daily.SectorDailyMean(:,1:7),2);
seven(~all(isfinite(daily.SectorDailyMean(:,1:7)),2))=NaN;
[matched,idx]=ismember(t,daily.EpochUTC);
assert(all(matched));
sunspot=NaN(height(daily),1); sunspot(idx)=sn;
values=[sunspot daily.BMean_nT daily.P1Mean daily.P1Median seven];
shown=values; q=shown(:,3:5); q(q<=0)=NaN; shown(:,3:5)=q;
%% 五面板图
visibility='off'; if visible, visibility='on'; end
fig=figure('Color','w','Position',[80 40 1600 1300],'Visible',visibility);
layout=tiledlayout(fig,5,1,'TileSpacing','compact','Padding','compact');
labels={{'Daily sunspot number'}, {'|B| daily mean','(nT)'}, ...
    {'P1 daily mean','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'P1 daily median','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'P1 sum: S1-S7 (including S4)','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}};
colors=[0.70 0.30 0.06;0.15 0.15 0.15;0.48 0.12 0.62;0.48 0.12 0.62;0.08 0.36 0.62];
ax=gobjects(5,1);
for k=1:5
    ax(k)=nexttile(layout);
    h=plot(ax(k),daily.EpochUTC,shown(:,k),'.-', ...
        'Color',colors(k,:),'MarkerSize',3,'LineWidth',0.4);
    if k>=3, set(ax(k),'YScale','log'); end
    ylabel(ax(k),labels{k});
    set(ax(k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',0.12);
    xlim(ax(k),[startUTC stopUTC]);
    text(ax(k),0.008,0.89,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
    if k<5, ax(k).XTickLabel=[]; end
    assert(isequaln(h.YData(:),shown(:,k)));
end
% Match the supplied daily figure's P1 intensity limits.
ylim(ax(3),[1e-2 1e1]); ylim(ax(4),[1e-2 1e1]);
linkaxes(ax,'x'); xlabel(ax(5),'UTC'); xtickformat(ax(5),'yyyy');
title(layout,'Voyager 1 beyond the heliopause | 2012-08-25 to 2021-12-16 | P1 0.57-1.78 MeV', ...
    'FontWeight','bold','FontSize',15);
fileStem=fullfile(OutputDir,'V1_daily_5panels_with_S4_sunspots');
exportgraphics(fig,[fileStem,'.png'],'Resolution',220);
exportgraphics(fig,[fileStem,'.pdf'],'ContentType','vector');
savefig(fig,[fileStem,'.fig']);
%% 来源、数值和缺测审计
result=struct;
result.CreatedUTC=datetime('now','TimeZone','UTC');
result.SolarSourceFile=SolarFile;
result.SolarSourceSHA256=Case1_File_SHA256(SolarFile);
result.SolarSourceURL='https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv';
result.SolarCredit='WDC-SILSO, Royal Observatory of Belgium, daily total sunspot number V2.0';
result.SolarSourceRows=sourceRows; result.SolarOriginalColumns=x;
result.SolarMethod='Original daily total sunspot number column 5; -1 mapped to NaN; zero retained. UTC noon display anchor. No averaging, smoothing, filling or time shift.';
result.VoyagerSourceAudit=fullfile(source.OutputFolder,'V1_daily_overview_audit.mat');
result.VoyagerMethod=source.Method;
result.PanelESectors=1:7;
result.TimeRange=[startUTC stopUTC];
result.ValidCount=sum(isfinite(values),1);
result.MissingCount=sum(~isfinite(values),1);
result.SolarZeroCount=nnz(sunspot==0);
result.CodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
result.OutputFiles=string({[fileStem,'.png'],[fileStem,'.pdf'],[fileStem,'.fig']});
result.Daily=table(daily.EpochUTC,sunspot,daily.BMean_nT,daily.P1Mean,daily.P1Median,seven, ...
    'VariableNames',{'EpochUTC','SunspotNumber','BMean_nT','P1Mean','P1Median','P1SevenSectorSum'});
writetable(result.Daily,fullfile(OutputDir,'daily_5panels_values.csv'));
save(fullfile(OutputDir,'daily_5panels_audit.mat'),'result','-v7.3');
disp(result.ValidCount); disp(result.MissingCount); fprintf('Output: %s\n',OutputDir);
if ~visible, close(fig); end
end
