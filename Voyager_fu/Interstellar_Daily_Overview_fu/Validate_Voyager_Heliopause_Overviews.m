function checks=Validate_Voyager_Heliopause_Overviews
% Output-only audit: no scientific calculations use derived MAT/CSV inputs.
root='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Heliopause_Overviews/';
checks=struct;
for sc=1:2
    folder=fullfile(root,sprintf('V%d',sc));
    v=load(fullfile(folder,'extended_overviews_audit.mat'),'audit'); a=v.audit;
    names={'heliopause_daily','heliopause_three_day','heliopause_monthly'};
    for j=1:3
        item=a.(names{j}); w=item.Windows;
        assert(all(w.StartUTC(2:end)==w.EndUTCExclusive(1:end-1)));
        assert(all(w.CalendarDays>0) && all(all(w.ValidDays<=w.CalendarDays)));
        for k=1:3, assert(isfile(item.OutputFiles(k))); end
        f=openfig(item.OutputFiles(3),'invisible'); ax=findall(f,'Type','axes'); assert(numel(ax)==5);
        for k=1:5
            tile=ax(k).Layout.Tile; h=findall(ax(k),'Type','line'); assert(numel(h)==1 && h.LineWidth==0.4);
            if tile>=3, assert(strcmp(ax(k).YScale,'log')); else, assert(strcmp(ax(k).YScale,'linear')); end
            y=h.YData; y=y(isfinite(y)); assert(all(y>=ax(k).YLim(1)&y<=ax(k).YLim(2)));
        end
        close(f);
    end
    checks.(sprintf('V%d',sc))='Three figures: bins, counts, axes, line width and plotted ranges passed.';
end
v=load(fullfile(root,'V2','source_recompute','daily_with_MAG_supplement_audit.mat'),'result'); s=v.result;
d=s.Daily; a=s.MAGSupplement.Daily; before=a{:,2}; used=a{:,6};
assert(isequaln(d.BMean_nT(isfinite(before)),before(isfinite(before))));
assert(all(~isfinite(before(used))) && all(isfinite(d.BMean_nT(used))));
assert(nnz(used)==338 && nnz(used & year(d.EpochUTC)==2021)==329);
assert(nnz(year(d.EpochUTC)==2021 & ~isfinite(d.BMean_nT))==36);
base=load(fullfile(root,'V2','source_recompute','V1_daily_overview_audit.mat'),'result');
assert(isequaln(d.P1Mean,base.result.Daily.P1Mean) && isequaln(d.P1Median,base.result.Daily.P1Median));
assert(isequaln(d.SectorDailyMean,base.result.Daily.SectorDailyMean));
% Direct two-stage mean cross-check on three supplemented 2021 dates.
r=s.MAGSupplement.RawSelectedRecords;
for day=datetime(2021,[1 6 11],[1 15 26],'TimeZone','UTC')
    rows=r.EpochUTC>=day&r.EpochUTC<day+days(1)&isfinite(r.F1_nT);
    t=r.EpochUTC(rows); b=r.F1_nT(rows); hr=dateshift(t,'start','hour'); u=unique(hr); avg=NaN(numel(u),1);
    for k=1:numel(u), avg(k)=mean(b(hr==u(k))); end
    expected=mean(avg); actual=d.BMean_nT(d.EpochUTC==day+hours(12));
    assert(abs(actual-expected)<1e-12);
end
checks.MAGRealDataDaysSupplied=338; checks.MAG2021DaysSupplied=329;
checks.OriginalFiniteCOHOAndAllParticleValuesUnchanged=true;
checks.DirectMAGHourlyThenDailyMeansPassed=true;
checks.CreatedUTC=datetime('now','TimeZone','UTC');
save(fullfile(root,'validation_audit.mat'),'checks'); disp(checks);
end
