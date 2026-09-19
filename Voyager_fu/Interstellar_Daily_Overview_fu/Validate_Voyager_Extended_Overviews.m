function checks=Validate_Voyager_Extended_Overviews
% 只读输出作验证；不作为科学输入或默认重画入口。
out='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/';
checks=struct;
for sc=1:2
    folder=fullfile(out,sprintf('V%d',sc));
    loaded=load(fullfile(folder,'extended_overviews_audit.mat'),'audit'); a=loaded.audit;
    names={'from2008_daily','from2008_three_day','from2008_monthly','all_daily','all_three_day','all_monthly'};
    for n=1:numel(names)
        item=a.(names{n}); w=item.Windows;
        assert(all(w.EndUTCExclusive(1:end-1)==w.StartUTC(2:end)));
        assert(all(w.CalendarDays>0) && all(all(w.ValidDays<=w.CalendarDays)));
        for k=1:3, assert(isfile(item.OutputFiles(k))); end
        f=openfig(item.OutputFiles(3),'invisible'); ax=findall(f,'Type','axes'); assert(numel(ax)==5);
        for k=1:numel(ax)
            tile=ax(k).Layout.Tile; h=findall(ax(k),'Type','line'); assert(numel(h)==1);
            assert(h.LineWidth==0.4 && strcmp(h.LineStyle,'-'));
            if tile>=3, assert(strcmp(ax(k).YScale,'log')); else, assert(strcmp(ax(k).YScale,'linear')); end
            y=h.YData; good=isfinite(y); assert(all(y(good)>=ax(k).YLim(1) & y(good)<=ax(k).YLim(2)));
        end
        close(f);
    end
    checks.(sprintf('V%d',sc))='Six figure groups, counts, shared boundaries, linear/log scales, thin lines and plotted limits passed.';
end
% V1 原日统计重叠部分逐值比较，包括太阳黑子数与七扇区求和。
old=readtable('C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V1_Interstellar_Daily_Overview/with_S4_sunspots/daily_5panels_values.csv');
old.EpochUTC.TimeZone='UTC';
a=load(fullfile(out,'V1','extended_overviews_audit.mat'),'audit'); w=a.audit.from2008_daily.Windows;
[found,idx]=ismember(old.EpochUTC,w.EpochUTC); assert(all(found));
oldValues=[old.SunspotNumber old.BMean_nT old.P1Mean old.P1Median old.P1SevenSectorSum];
newValues=[w.SunspotNumber(idx) w.B_nT(idx) w.P1Mean(idx) w.P1MedianPanel(idx) w.SectorSum(idx)];
assert(isequal(isnan(oldValues),isnan(newValues)),'Old/new missing masks differ.');
difference=abs(oldValues-newValues); tolerance=32*eps(max(1,max(abs(oldValues),abs(newValues))));
fprintf('Max old CSV/new MAT differences by panel:\n'); disp(max(difference,[],1,'omitnan'));
assert(all(difference(isfinite(difference))<=tolerance(isfinite(difference))),'Changes exceed CSV floating-point roundoff.');
checks.PreviousV1DailyValuesAgreeWithinCSVRoundoff=true;
checks.MaxPreviousDailyDifference=max(difference,[],1,'omitnan');
checks.CreatedUTC=datetime('now','TimeZone','UTC');
save(fullfile(out,'validation_audit.mat'),'checks'); disp(checks);
end
