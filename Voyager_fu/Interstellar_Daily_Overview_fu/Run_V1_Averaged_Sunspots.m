function result = Run_V1_Averaged_Sunspots(visible)
% 三日和自然月太阳黑子数平均；Voyager 原始CDF重算。
if nargin<1, visible=true; end
%% 路径
CodeDir='C:/Users/Administrator/Documents/FWD_matlab/Voyager_fu/Interstellar_Daily_Overview_fu/';
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/averaged_no_S4_sunspots/';
SolarFile='Z:/SPART-WORK/Data/Solar_Indices/raw/sunspot/SN_d_tot_V2.0.csv';
addpath(CodeDir);
if ~isfolder(OutputDir), mkdir(OutputDir); end
%% 复用原始CDF统计入口，独立输出，保留原四面板图
source=Run_V1_ThreeDay_Monthly_Overview(false,'both',fullfile(OutputDir,'source_recompute'),false);
x=readmatrix(SolarFile,'Delimiter',';','FileType','text');
t=datetime(x(:,1),x(:,2),x(:,3),'TimeZone','UTC')+hours(12);
use=t>=source.Method.StartUTC & t<source.Method.EndUTCExclusive;
x=x(use,:); t=t(use); sourceRows=find(use);
assert(numel(unique(t))==numel(t));
sn=x(:,5); sn(sn==-1)=NaN;
assert(all(isnan(sn)|sn>=0));
%% 原有窗口内对有效每日黑子数等权算术平均
result=struct;
result.SourceAudit=fullfile(OutputDir,'source_recompute','three_day_monthly_audit.mat');
result.SolarSource=SolarFile;
result.SolarSHA256=Case1_File_SHA256(SolarFile);
result.SolarSourceRows=sourceRows; result.SolarOriginalColumns=x;
result.Credit='WDC-SILSO, Royal Observatory of Belgium, Daily total sunspot number V2.0';
result.Method='Use original Voyager bin boundaries; mean of finite daily SILSO values. -1 to NaN, zero retained. No smoothing, interpolation, lag or coverage threshold. Partial endpoint windows retained.';
result.VoyagerMethod=source.Method;
result.CodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
modes={'three_day','monthly'};
for m=1:2
    tag=modes{m}; w=source.(tag).Windows;
    avg=NaN(height(w),1); count=zeros(height(w),1); bins=zeros(numel(t),1);
    for k=1:height(w)
        rows=t>=w.StartUTC(k) & t<w.EndUTCExclusive(k);
        assert(all(bins(rows)==0)); bins(rows)=k;
        v=sn(rows & isfinite(sn)); count(k)=numel(v);
        if ~isempty(v), avg(k)=mean(v); end
    end
    assert(all(bins>0));
    good=isfinite(sn);
    check=accumarray(bins(good),sn(good),[height(w) 1],@mean,NaN);
    assert(isequaln(avg,check) && sum(count)==nnz(good));
    cName='MeanOfDailyP1Medians';
    if strcmp(tag,'three_day'), cName='ThreeDayP1Median'; end
    values=[avg w.BMean_nT w.MeanOfDailyP1Means w.(cName) w.MeanOfDailySixSectorSums];
    w.SunspotMean=avg; w.SunspotValidDays=count;
    files=drawFive(w.EpochUTC,values,source.Method,tag,OutputDir,visible);
    result.(tag)=struct('Windows',w,'SolarBinIndex',bins,'OutputFiles',files);
    writetable(w,fullfile(OutputDir,[tag,'_values.csv']));
    fprintf('%s: %d windows, %d valid solar windows, %d daily samples\n',tag,height(w),nnz(isfinite(avg)),sum(count));
end
result.CreatedUTC=datetime('now','TimeZone','UTC');
save(fullfile(OutputDir,'averaged_sunspots_audit.mat'),'result','-v7.3');
end

function files=drawFive(t,values,method,tag,out,visible)
%% 图形标签；计算说明仅保存在README与审计
vis='off'; if visible, vis='on'; end
fig=figure('Color','w','Position',[80 40 1600 1300],'Visible',vis);
layout=tiledlayout(fig,5,1,'TileSpacing','compact','Padding','compact');
period='Monthly'; medianLabel='Mean of daily P1 medians';
if strcmp(tag,'three_day'), period='3-day'; medianLabel='3-day P1 median'; end
labels={{[period,' mean sunspot number']},{'Mean of daily |B|','(nT)'}, ...
    {'Mean of daily P1 means','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {medianLabel,'(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}, ...
    {'Mean P1 sum: S1,S2,S3,S5,S6,S7','(cm^{-2} s^{-1} sr^{-1} MeV^{-1})'}};
colors=[0.70 0.30 0.06;0.15 0.15 0.15;0.48 0.12 0.62;0.48 0.12 0.62;0.08 0.36 0.62];
ax=gobjects(5,1);
for k=1:5
    ax(k)=nexttile(layout); y=values(:,k);
    if k>=3, y(y<=0)=NaN; end
    h=plot(ax(k),t,y,'.-','Color',colors(k,:),'MarkerSize',3,'LineWidth',0.4);
    assert(isequaln(h.YData(:),y));
    if k>=3, set(ax(k),'YScale','log'); end
    ylabel(ax(k),labels{k});
    set(ax(k),'FontSize',11,'TickDir','out','Box','on','XGrid','on','GridAlpha',0.12);
    xlim(ax(k),[method.StartUTC method.EndUTCExclusive]);
    text(ax(k),0.008,0.89,sprintf('(%c)',96+k),'Units','normalized','FontWeight','bold');
    if k<5, ax(k).XTickLabel=[]; end
end
linkaxes(ax,'x'); xlabel(ax(5),'UTC'); xtickformat(ax(5),'yyyy');
title(layout,[period,' | Voyager 1 | 2012-08-25 to 2021-12-16 | P1 0.57-1.78 MeV'], ...
    'FontWeight','bold','FontSize',15);
stem=fullfile(out,['V1_',tag,'_5panels_no_S4_sunspots']);
files=string({[stem,'.png'],[stem,'.pdf'],[stem,'.fig']});
exportgraphics(fig,files(1),'Resolution',220);
exportgraphics(fig,files(2),'ContentType','vector'); savefig(fig,files(3));
if ~visible, close(fig); end
end
