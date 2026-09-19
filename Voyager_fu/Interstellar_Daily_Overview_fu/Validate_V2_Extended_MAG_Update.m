function checks=Validate_V2_Extended_MAG_Update
% Output-only validation against pre-update statistics and original COHO.
root='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Extended_Overviews/V2/';
before=fullfile(root,'mag_update_20260910','before_values');
files=dir(fullfile(before,'*_values.csv')); assert(numel(files)==6);
checks=struct;
for k=1:numel(files)
    old=readtable(fullfile(before,files(k).name)); new=readtable(fullfile(root,files(k).name));
    names=setdiff(old.Properties.VariableNames,{'B_nT','ValidDays_2'},'stable');
    assert(isequaln(old(:,names),new(:,names)),'Nonmagnetic data or window changed.');
    if contains(files(k).name,'daily')
        good=isfinite(old.B_nT); assert(isequaln(old.B_nT(good),new.B_nT(good)),'Finite daily COHO changed.');
        assert(all(isfinite(new.B_nT(good))));
    end
end
v=load(fullfile(root,'source_recompute','daily_with_MAG_supplement_audit.mat'),'result'); s=v.result;
a=s.MAGSupplement.Daily; used=a{:,6}; oldB=a{:,2}; d=s.Daily;
assert(all(~isfinite(oldB(used))) && all(isfinite(d.BMean_nT(used))));
assert(isequaln(oldB(isfinite(oldB)),d.BMean_nT(isfinite(oldB))));
checks.TotalDaysSupplied=nnz(used);
checks.DaysSuppliedFrom2008=nnz(used & year(d.EpochUTC)>=2008);
checks.DaysSupplied2021=nnz(used & year(d.EpochUTC)==2021);
checks.Missing2021=nnz(year(d.EpochUTC)==2021 & ~isfinite(d.BMean_nT));
checks.OtherPanelsAndWindowsUnchanged=true;
checks.FiniteOriginalDailyBUnchanged=true;
% Compare magnetic daily values in the shared interval to verified heliopause plot.
hp=readtable('C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/Voyager_Heliopause_Overviews/V2/heliopause_daily_values.csv');
new=readtable(fullfile(root,'all_daily_values.csv'));
[yes,idx]=ismember(hp.EpochUTC,new.EpochUTC); assert(all(yes)); assert(isequaln(hp.B_nT,new.B_nT(idx)));
checks.HeliopauseMagneticValuesIdentical=true;
v=load(fullfile(root,'extended_overviews_audit.mat'),'audit'); a=v.audit;
keys={'from2008_daily','from2008_three_day','from2008_monthly','all_daily','all_three_day','all_monthly'};
for j=1:6
    item=a.(keys{j}); w=item.Windows; assert(all(w.StartUTC(2:end)==w.EndUTCExclusive(1:end-1)));
    for k=1:3, assert(isfile(item.OutputFiles(k))); end
    f=openfig(item.OutputFiles(3),'invisible'); ax=findall(f,'Type','axes'); assert(numel(ax)==5);
    for k=1:5
        h=findall(ax(k),'Type','line'); assert(numel(h)==1 && h.LineWidth==0.4 && strcmp(h.LineStyle,'-'));
        if ax(k).Layout.Tile>=3, assert(strcmp(ax(k).YScale,'log')); else, assert(strcmp(ax(k).YScale,'linear')); end
        y=h.YData; y=y(isfinite(y)); assert(all(y>=ax(k).YLim(1)&y<=ax(k).YLim(2)));
    end
    close(f);
end
checks.SixFigureGroupsPassed=true;
checks.CreatedUTC=datetime('now','TimeZone','UTC');
save(fullfile(root,'mag_update_20260910','validation_audit.mat'),'checks'); disp(checks);
end
