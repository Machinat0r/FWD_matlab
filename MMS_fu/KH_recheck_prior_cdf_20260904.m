function KH_recheck_prior_cdf_20260904
% 复核单变量CDF接口的返回类型，修复第一次抽查的时间定位错误。
root='Z:\SPART-WORK\Data\MMS\derived\KH\catalog_audit_20260903';
addpath('C:\Users\Administrator\Documents\irfu-matlab-master\contrib\nasa_cdf_patch');
cd(root);
fp=fullfile(root,'cdf_presence.json');
S=jsondecode(fileread(fp));
T=readtable('C:\Users\Administrator\Documents\KH\MMS_KH_published_event_catalog.csv','TextType','string');
for k=1:numel(S.entries)
    r=S.entries(k);
    if isempty(r.Errors),continue;end
    row=find(T.EventID==string(r.EventID),1);
    t1=datenum(char(T.StartUTC(row)),'yyyy-mm-dd HH:MM:SS');
    t2=datenum(char(T.EndUTC(row)),'yyyy-mm-dd HH:MM:SS');
    oldErrors=r.Errors;r.Errors={};
    for j=1:numel(oldErrors)
        parts=strsplit(oldErrors{j},' : ');f=parts{1};
        try
            tt=spdfcdfread(f,'Variables',{r.EpochVariable},'CombineRecords',true);
            if iscell(tt),tt=tt{1};end
            hit=find(tt>=t1 & tt<t2,32,'first');
            if isempty(hit),r.Status='候选CDF时刻不在事件内';continue;end
            a=spdfcdfread(f,'Variables',{r.EpochVariable,r.Variable},'Records',hit(:)'-1,'CombineRecords',false,'ConvertEpochToDatenum',true);
            for jr=1:size(a,1)
                val=double(a{jr,2});
                if all(isfinite(val(:))) && all(abs(val(:))<1e29)
                    r.Status='有有效记录';r.EvidenceFile=f;r.EvidenceUTC=datestr(double(a{jr,1}),'yyyy-mm-dd HH:MM:SS.FFF');break;
                end
            end
        catch ME
            r.Errors{end+1}=[f ' : ' ME.message];
        end
    end
    S.entries(k)=r;
end
S.checkedUTC=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss'));
fid=fopen(fp,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(S));fclose(fid);
fprintf('Recheck complete; remaining read errors: %d\n',sum(arrayfun(@(r)~isempty(r.Errors),S.entries)));
end
