function checks=Validate_Voyager_From1990_Overviews
% Output-only audit for the requested 1990 start; never a science input.
%% Paths
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/';
ReviewDir=fullfile(OutputDir,'start1990_update_20260911');
startUTC=datetime(1990,1,1,'TimeZone','UTC');
modes={'daily','three_day','monthly'};
checks=struct;
%% Compare unchanged daily/monthly values and original 2008 products
for sc=1:2
    folder=sprintf('V%d',sc); out=fullfile(OutputDir,folder);
    previous=fullfile(ReviewDir,folder);
    for tag={'daily','monthly'}
        name=['all_',tag{1},'_values.csv'];
        old=readtable(fullfile(previous,name)); new=readtable(fullfile(out,name));
        old=old(year(old.StartUTC)>=1990,:);
        assert(isequaln(old,new),'Daily/monthly overlap changed: V%d %s',sc,tag{1});
    end
    a=load(fullfile(out,'extended_overviews_audit.mat'),'audit'); a=a.audit;
    before=load(fullfile(previous,'extended_overviews_audit.mat'),'audit');
    for k=1:3
        name=['from2008_',modes{k}];
        assert(isequaln(a.(name),before.audit.(name)),'2008 audit changed.');
        old=readtable(fullfile(previous,[name,'_values.csv']));
        new=readtable(fullfile(out,[name,'_values.csv'])); assert(isequaln(old,new));
    end
    daily=a.all_daily.Windows;
    w=a.all_three_day.Windows;
    assert(w.StartUTC(1)==startUTC && all(diff(w.StartUTC)==days(3)));
    bins=discretize(daily.EpochUTC,[w.StartUTC;w.EndUTCExclusive(end)]);
    columns={'SunspotNumber','B_nT','P1Mean'};
    for p=1:numel(columns)
        x=daily.(columns{p}); valid=isfinite(x);
        expected=accumarray(bins(valid),x(valid),[height(w),1],@mean,NaN);
        assert(isequaln(expected,w.(columns{p})),'3-day means mismatch.');
    end
    raw=load(fullfile(out,'source_recompute','V1_daily_overview_audit.mat'),'raw'); raw=raw.raw;
    valid=raw.EpochUTC>=startUTC & raw.EpochUTC<w.EndUTCExclusive(end) & isfinite(raw.P1);
    bins=discretize(raw.EpochUTC(valid),[w.StartUTC;w.EndUTCExclusive(end)]);
    expected=accumarray(bins,raw.P1(valid),[height(w),1],@median,NaN);
    count=accumarray(bins,1,[height(w),1],@sum,0);
    assert(isequaln(expected,w.P1MedianPanel) && isequal(count,w.P1MedianHourlySamples));
    %% Inspect all six delivered MATLAB figures
    for m=1:3
        item=a.(['all_',modes{m}]); w=item.Windows;
        assert(w.StartUTC(1)==startUTC);
        assert(all(w.StartUTC(2:end)==w.EndUTCExclusive(1:end-1)));
        for j=1:3, assert(isfile(item.OutputFiles(j))); end
        f=openfig(item.OutputFiles(3),'invisible'); ax=findall(f,'Type','axes'); assert(numel(ax)==5);
        data=[w.SunspotNumber w.B_nT w.P1Mean w.P1MedianPanel w.SectorSum];
        for j=1:5
            k=ax(j).Layout.Tile; h=findall(ax(j),'Type','line');
            assert(numel(h)==1 && h.LineWidth==0.4 && strcmp(h.LineStyle,'-'));
            assert(isequal(ax(j).XLim,[startUTC w.EndUTCExclusive(end)]));
            y=data(:,k);
            if k>=3, y(y<=0)=NaN; assert(strcmp(ax(j).YScale,'log'));
            else, assert(strcmp(ax(j).YScale,'linear')); end
            assert(isequaln(h.YData(:),y));
        end
        close(f);
    end
    checks.(folder)=struct('StartUTC',startUTC,'EndUTCExclusive',daily.EndUTCExclusive(end), ...
        'DailyPoints',height(a.all_daily.Windows),'ThreeDayPoints',height(a.all_three_day.Windows), ...
        'MonthlyPoints',height(a.all_monthly.Windows),'DailyMonthlyOverlapUnchanged',true, ...
        'ThreeDayHourlyMedianPassed',true,'Existing2008AuditsUnchanged',true,'ThreeFiguresPassed',true);
end
%% Retain the previously verified V2 magnetic supplementation
s=load(fullfile(OutputDir,'V2','source_recompute','daily_with_MAG_supplement_audit.mat'),'result');
d=s.result.Daily; supplied=s.result.MAGSupplement.Daily{:,6};
assert(nnz(supplied)==340 && nnz(supplied & year(d.EpochUTC)==2021)==329);
checks.V2SuppliedMagneticDays=nnz(supplied);
checks.CreatedUTC=datetime('now','TimeZone','UTC');
save(fullfile(ReviewDir,'validation_audit.mat'),'checks'); disp(checks);
%% Refresh the output-only twelve-figure contact sheet
f=figure('Visible','off','Color','w','Position',[0 0 1800 2100]);
t=tiledlayout(f,4,3,'TileSpacing','none','Padding','none');
for sc=1:2
    for range={'from2008','all'}
        for m=1:3
            path=fullfile(OutputDir,sprintf('V%d',sc),sprintf('V%d_%s_%s_5panels.png',sc,range{1},modes{m}));
            ax=nexttile(t); image(ax,imread(path)); axis(ax,'image'); axis(ax,'off');
        end
    end
end
exportgraphics(f,fullfile(OutputDir,'twelve_figures_contact_sheet.png'),'Resolution',110);
exportgraphics(f,fullfile(OutputDir,'twelve_figures_contact_sheet.jpg'),'Resolution',110); close(f);
end