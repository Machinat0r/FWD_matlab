function summary = Validate_V2_Plasma_Overview
% 独立核对交付统计、CDF抽样和已保存图形，不作为重画输入。
%% 路径
CodeDir=fileparts(mfilename('fullpath')); ProjectCode=fileparts(CodeDir);
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V2_Plasma_Overview/';
addpath(fullfile(ProjectCode,'Case1_PPT_VerticalLine_Events_7d_fu'));
Case1_Add_IRFU_Path('C:/Users/Administrator/Documents/irfu-matlab-master');
data=load(fullfile(OutputDir,'V2_plasma_overview_audit.mat'),'result');
r=data.result; d=r.Daily;
assert(height(d)==12784);
assert(d.EpochUTC(1)==datetime(1990,1,1,12,0,0,'TimeZone','UTC'));
assert(d.EpochUTC(end)==datetime(2024,12,31,12,0,0,'TimeZone','UTC'));
assert(isequal(r.ProtectedBefore,r.ProtectedAfter));
for k=1:height(r.ProtectedBefore)
    assert(strcmp(Case1_File_SHA256(r.ProtectedBefore.File(k)),r.ProtectedBefore.SHA256(k)));
end
%% 从原始 CDF 重新抽查跨时期的日均，而非以派生 CSV 计算
checkDate=datetime([1990 1995 2000 2007 2007 2008 2015 2018 2018 2019 2024], ...
    [1 1 1 8 9 1 1 11 11 1 12],[2 15 1 30 1 15 15 4 5 1 31],'TimeZone','UTC').';
variable=["protonFlux1_LECP","protonDensity","protonTemp","V"];
dailyName=["P1Mean","ProtonDensity_cm3","ProtonTemperature_K","ProtonSpeed_kms"];
countName=["P1SampleCount","DensitySampleCount","TemperatureSampleCount","SpeedSampleCount"];
for k=1:numel(checkDate)
    token=datestr(checkDate(k),'yyyymm');
    paths=string(cellfun(@(s) char(s.File),r.COHOSources,'UniformOutput',false));
    selected=contains(paths,['_',token,'01_']); assert(nnz(selected)==1);
    q=Voyager_Read_CDF_Product(paths(selected),'coho');
    use=dateshift(q.Epoch,'start','day')==checkDate(k);
    row=dateshift(d.EpochUTC,'start','day')==checkDate(k);assert(nnz(row)==1);
    for j=1:numel(variable)
        values=q.(variable(j))(use);values=values(isfinite(values));
        actual=d.(dailyName(j))(row);
        if isempty(values),assert(isnan(actual));
        else,assert(abs(mean(values)-actual)<=1e-12*max(1,abs(actual)));end
        assert(d.(countName(j))(row)==numel(values));
    end
end
rawNames=["P1","ProtonDensity_cm3","ProtonTemperature_K","ProtonSpeed_kms"];
for j=1:4
    assert(sum(d.(countName(j)))==nnz(isfinite(r.RawCOHO.(rawNames(j)))));
    assert(isequal(d.(countName(j))==0,isnan(d.(dailyName(j)))));
end
lastPLS=datetime(2018,11,4,12,0,0,'TimeZone','UTC');
for j=2:4
    assert(max(d.EpochUTC(isfinite(d.(dailyName(j)))))==lastPLS);
end
assert(all(r.Preservation.MaximumRelativeDifference<1e-12));
%% 图中每条曲线、坐标范围、对数轴及 panel 标签
f=openfig(r.OutputFiles(3),'invisible'); cleanup=onCleanup(@() close(f));
values=d{:,2:7}; logarithmic=[false false true true true false];
for k=1:6
    ax=findall(f,'Type','axes','Tag',sprintf('panel_%d',k));assert(numel(ax)==1);
    h=findall(ax,'Type','line','Tag',d.Properties.VariableNames{k+1});assert(numel(h)==1);
    y=values(:,k);if logarithmic(k),y(y<=0)=NaN;end
    assert(isequaln(h.YData(:),y));
    assert(isequal(h.XData(:),datenum(d.EpochUTC)));
    assert(isequal(ax.XLim,datenum([datetime(1990,1,1) datetime(2025,1,1)])));
    assert(strcmp(ax.YScale,'log')==logarithmic(k));
    assert(strcmp(h.LineStyle,'-') && h.LineWidth==.4);
end
assert(numel(findall(f,'Type','axes'))==6);
panelC=findall(f,'Type','axes','Tag','panel_3');
assert(contains(string(panelC.YLabel.String{1}),'LECP P1'));
assert(contains(r.Method.EnergyLabel,'0.52-1.45'));
%% 保存验证结果
summary=struct('Passed',true,'CheckedUTC',char(datetime('now','TimeZone','UTC')), ...
    'DailyRows',height(d),'FigurePanels',6,'CDFSpotCheckDays',numel(checkDate), ...
    'OriginalABCMatch',true,'PreviousFigureHashesUnchanged',true, ...
    'RawCDFsVerified',r.OfficialInventory.SelectedCDFs,'NewCDFDownloads',r.OfficialInventory.DownloadedCDFs, ...
    'MAGDaysSupplemented',nnz(d.MAGDailySource=="MAG48s_hourly_then_daily"), ...
    'LastValidPLSDate','2018-11-04','P1EnergyLabel','0.52-1.45 MeV');
summary.ValidationCodeSHA256=Case1_File_SHA256([mfilename('fullpath'),'.m']);
file=fopen(fullfile(OutputDir,'final_validation.json'),'w','n','UTF-8');
assert(file~=-1);fprintf(file,'%s\n',jsonencode(summary,'PrettyPrint',true));fclose(file);
disp(summary);
end
