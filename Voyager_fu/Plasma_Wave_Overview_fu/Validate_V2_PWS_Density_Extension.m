function summary = Validate_V2_PWS_Density_Extension
% 检验已交付图形，逐点对照官方网站原生 CSV。
%% 原始产品与本次结果
CodeDir=fileparts(mfilename('fullpath'));ProjectCode=fileparts(CodeDir);
addpath(fullfile(ProjectCode,'Case1_PPT_VerticalLine_Events_7d_fu'));
OutputDir='C:/Users/Administrator/Documents/Recovery-Work-Voyager_betatron/V2_Plasma_Overview/With_Official_PWS_Density/';
DataDir='Z:/SPART-WORK/Data/Voyager/voyager2/';
q=load(fullfile(OutputDir,'V2_plasma_overview_audit.mat'),'result');r=q.result;d=r.Daily;
[native,source]=V2_Read_Official_PWS_Density(DataDir);
assert(isequaln(native,r.PWSDensity));
assert(height(native)==5 && source.RecordCount==5);
assert(height(d)==13028);
assert(r.Method.EndUTCExclusive==datetime(2025,9,2,'TimeZone','UTC'));
assert(all(r.Preservation.MaximumRelativeDifference<1e-12));
assert(isequal(r.ProtectedBefore,r.ProtectedAfter));
for k=1:height(r.ProtectedBefore)
    assert(strcmp(Case1_File_SHA256(r.ProtectedBefore.File(k)),r.ProtectedBefore.SHA256(k)));
end
% PWS 电子密度不写入 PLS 质子密度、温度或速度数组。
afterPLS=d.EpochUTC>=datetime(2018,11,5,'TimeZone','UTC');
assert(all(all(isnan(d{afterPLS,5:7}))));
afterCOHO=d.EpochUTC>=datetime(2025,1,1,'TimeZone','UTC');
assert(all(all(isnan(d{afterCOHO,3:4}))));
assert(all(d.SunspotNumber(afterCOHO)>=0 | isnan(d.SunspotNumber(afterCOHO))));
rawNames=["P1","ProtonDensity_cm3","ProtonTemperature_K","ProtonSpeed_kms"];
countNames=["P1SampleCount","DensitySampleCount","TemperatureSampleCount","SpeedSampleCount"];
for k=1:4,assert(nnz(isfinite(r.RawCOHO.(rawNames(k))))==sum(d.(countNames(k))));end
%% 六个 panel 坐标和原始曲线
f=openfig(r.OutputFiles(3),'invisible'); cleaner=onCleanup(@()close(f)); %#ok<NASGU>
assert(numel(findall(f,'Type','axes'))==6);
isLog=[false false true true true false];
for k=1:6
    ax=findall(f,'Type','axes','Tag',sprintf('panel_%d',k));assert(numel(ax)==1);
    h=findall(ax,'Type','line','Tag',d.Properties.VariableNames{k+1});assert(numel(h)==1);
    value=d{:,k+1};if isLog(k),value(value<=0)=NaN;end
    assert(isequaln(h.YData(:),value));
    assert(isequal(h.XData(:),datenum(d.EpochUTC)));
    assert(strcmp(ax.YScale,'log')==isLog(k));
    assert(strcmp(h.LineStyle,'-') && h.LineWidth==.4);
end
ax=findall(f,'Type','axes','Tag','panel_4');
h=findall(ax,'Type','line','Tag','PWS_ElectronDensity_native');assert(numel(h)==1);
assert(isequal(h.XData(:),datenum(native.EpochUTC)));
assert(isequal(h.YData(:),native.ElectronDensity_cm3));
assert(strcmp(h.LineStyle,'none'));
assert(all(h.YData>ax.YLim(1) & h.YData<ax.YLim(2)));
assert(all(h.XData>ax.XLim(1) & h.XData<ax.XLim(2)));
assert(contains(string(ax.YLabel.String{1}),'Density'));
assert(~contains(string(ax.YLabel.String{1}),'proton'));
assert(isempty(findall(f,'Type','legend')));
summary=struct('Passed',true,'CheckedUTC',char(datetime('now','TimeZone','UTC')), ...
    'SourceMD5',source.MD5,'PWSRecordsProvided',source.RecordCount,'PWSRecordsPlotted',numel(h.YData), ...
    'ExactNativeTimesAndDensities',true,'NoPWSValuesWrittenToPLS',true, ...
    'AllPWSPointsInsideAxes',true,'OldFiguresUnchanged',true, ...
    'FirstPWSTimeUTC',char(native.EpochUTC(1)),'LastPWSTimeUTC',char(native.EpochUTC(end)), ...
    'PLSTemperatureAndSpeedEnd','2018-11-04','DailyGridRows',height(d), ...
    'ValidationCodeSHA256',Case1_File_SHA256([mfilename('fullpath'),'.m']));
fid=fopen(fullfile(OutputDir,'final_validation.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));fclose(fid);disp(summary);
end
