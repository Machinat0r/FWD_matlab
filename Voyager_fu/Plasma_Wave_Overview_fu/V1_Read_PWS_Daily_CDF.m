function [wave,sourceAudit] = V1_Read_PWS_Daily_CDF(files,firstUTC,lastUTC,useParallel)
% 原始 PWS CDF -> 逐频道 UTC 日算术均值；无插值、频率合并或背景扣除。
% LASTUTC 为不包含的右端点。所有源记录必须落在其文件标称 UTC 日内。
if nargin<4, useParallel=true; end
day=(firstUTC:caldays(1):lastUTC-caldays(1)).';
files=string(files); files=files(:); n=numel(files); assert(n>0,'No original PWS CDF files.');
sourceAudit=cell(n,1);
if useParallel && license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if isempty(pool), pool=parpool('Processes',4); end %#ok<NASGU>
    queue=parallel.pool.DataQueue; completed=0;
    afterEach(queue,@reportProgress);
    parfor k=1:n
        sourceAudit{k}=readOne(files(k));
        send(queue,k);
    end
else
    for k=1:n
        sourceAudit{k}=readOne(files(k));
        if mod(k,250)==0 || k==n, fprintf('PWS read %d/%d\n',k,n); end
    end
end
frequency=sourceAudit{1}.Frequency_Hz;
means=NaN(numel(day),numel(frequency)); count=zeros(size(means));
hasFile=false(numel(day),1);
for k=1:n
    s=sourceAudit{k};
    assert(isequal(s.Frequency_Hz,frequency),'PWS frequency grid changes: review source.');
    [found,row]=ismember(s.DayUTC,day);
    assert(found && ~hasFile(row),'Duplicate or out-of-window PWS source day.');
    hasFile(row)=true;
    count(row,:)=s.ValidCount;
    good=s.ValidCount>0;
    means(row,good)=s.Sum_Vm(good)./s.ValidCount(good);
end
wave=table(day+hours(12),means,count,hasFile, ...
    'VariableNames',{'EpochUTC','EMean_Vm','ValidSampleCount','OfficialCDFPresent'});
wave.Properties.UserData=struct('Frequency_Hz',frequency, ...
    'TimeBin','UTC [00:00, next 00:00); unweighted finite sample mean in each channel.', ...
    'Units','V/m','NoInterpolation',true,'NoFrequencyIntegration',true);
    function reportProgress(~)
        completed=completed+1;
        if mod(completed,250)==0 || completed==n
            fprintf('PWS read %d/%d\n',completed,n);
        end
    end
end

function s=readOne(file)
obj=dataobj(char(file));
t=getv(obj,'epoch'); e=getv(obj,'electric_field'); f=getv(obj,'frequency');
assert(~isempty(t) && ~isempty(e) && ~isempty(f),'Missing required PWS variable.');
assert(strcmpi(strtrim(e.UNITS),'V/m'),'Unexpected PWS electric-field units.');
assert(strcmp(e.DEPEND_0,'epoch') && strcmp(e.DEPEND_1,'frequency'));
assert(contains(lower(string(t.type)),'tt2000'),'Unexpected PWS epoch type.');
tt=int64(t.data(:));
assert(numel(unique(tt))==numel(tt),'Duplicate PWS Epoch requires review.');
epoch=datetime(EpochTT(tt).epochUnix,'ConvertFrom','posixtime','TimeZone','UTC');
frequency=double(f.data(:)).';
assert(numel(frequency)==16 && all(diff(frequency)>0));
values=double(e.data);
if isempty(epoch), values=zeros(0,16); end
assert(isequal(size(values),[numel(epoch) 16]),'Unexpected PWS array dimensions.');
dateToken=regexp(char(file),'_(\d{8})_v','tokens','once');
assert(~isempty(dateToken));
day=datetime(dateToken{1},'InputFormat','yyyyMMdd','TimeZone','UTC');
assert(all(epoch>=day & epoch<day+caldays(1)), ...
    'PWS source contains records outside its nominal day; stop for review.');
% 仅采用实际 CDF 给出的质量范围；不增加人为阈值。
badFill=false(size(values)); badLow=badFill; badHigh=badFill;
if isfield(e,'FILLVAL'), badFill=values==double(e.FILLVAL); end
if isfield(e,'VALIDMIN'), badLow=values<double(e.VALIDMIN); end
if isfield(e,'VALIDMAX'), badHigh=values>double(e.VALIDMAX); end
bad=badFill|badLow|badHigh|~isfinite(values);
values(bad)=NaN;
count=sum(isfinite(values),1);
summed=sum(values,1,'omitnan');
% 独立核对已筛选原记录的直接均值与 sum/count。
average=summed./count; average(count==0)=NaN;
direct=mean(values,1,'omitnan');
assert(all(abs(average-direct)<=max(abs(direct),1)*1e-12 | (isnan(average)&isnan(direct))));
metadata=e; metadata=rmfield(metadata,'data');
s=struct('File',file,'SHA256',string(Case1_File_SHA256(file)), ...
    'DayUTC',day,'Records',numel(epoch),'FirstEpochUTC',NaT('TimeZone','UTC'), ...
    'LastEpochUTC',NaT('TimeZone','UTC'),'Frequency_Hz',frequency, ...
    'ValidCount',count,'Sum_Vm',summed,'RejectedCount',sum(bad,1), ...
    'FillCount',sum(badFill,1),'BelowValidMinCount',sum(badLow&~badFill,1), ...
    'AboveValidMaxCount',sum(badHigh&~badFill,1), ...
    'VariableMetadata',metadata,'GlobalAttributes',obj.GlobalAttributes);
if ~isempty(epoch)
    s.FirstEpochUTC=min(epoch); s.LastEpochUTC=max(epoch);
end
end
