function KH_check_cdf_presence_20260903(catalogFile)
% 直接读取原CDF确认事件时段内B/Vi/Ve/E记录存在；抽查不等于连续/共同覆盖。
%% 1. 目录、输入
root='Z:\SPART-WORK\Data\MMS';
out=fullfile(root,'derived','KH','catalog_audit_20260903');
addpath('C:\Users\Administrator\Documents\irfu-matlab-master\contrib\nasa_cdf_patch');
cd(out);
if nargin<1, catalogFile='C:\Users\Administrator\Documents\KH\MMS_KH_published_event_catalog.csv'; end
T=readtable(catalogFile,'TextType','string');
products=["B","Vi","Ve","E"]; entries=struct([]);
%% 2. 各产品逐星读取真实记录
firstEvent=1;
if isfile(fullfile(out,'cdf_presence.json'))
    old=jsondecode(fileread(fullfile(out,'cdf_presence.json')));entries=old.entries(:).';firstEvent=old.processed+1;
end
for ie=firstEvent:height(T)
    t1=datenum(char(T.StartUTC(ie)),'yyyy-mm-dd HH:MM:SS');
    t2=datenum(char(T.EndUTC(ie)),'yyyy-mm-dd HH:MM:SS');
    for ip=1:4
        for ic=1:4
            sc=sprintf('mms%d',ic);ep='Epoch';
            switch products(ip)
                case "B", parts={sc,'fgm','brst','l2'};v=[sc '_fgm_b_gse_brst_l2'];
                case "Vi", parts={sc,'fpi','brst','l2','dis-moms'};v=[sc '_dis_bulkv_gse_brst'];
                case "Ve", parts={sc,'fpi','brst','l2','des-moms'};v=[sc '_des_bulkv_gse_brst'];
                case "E", parts={sc,'edp','brst','l2','dce'};v=[sc '_edp_dce_gse_brst_l2'];ep=[sc '_edp_epoch_brst_l2'];
            end
            files={};times=[];
            for day=floor(t1)-1:floor(t2)
                folder=fullfile(root,parts{:},datestr(day,'yyyy'),datestr(day,'mm'),datestr(day,'dd'));
                if ~isfolder(folder),continue;end
                d=dir(fullfile(folder,'*.cdf'));
                for j=1:numel(d)
                    token=regexp(d(j).name,'_(\d{14})_v','tokens','once');
                    if isempty(token),continue;end
                    files{end+1}=fullfile(d(j).folder,d(j).name); %#ok<AGROW>
                    times(end+1)=datenum(token{1},'yyyymmddHHMMSS'); %#ok<AGROW>
                end
            end
            [~,idx]=sort(string(files));files=files(idx);times=times(idx);
            [~,idx]=unique(times,'last');files=files(idx);times=times(idx);
            [times,idx]=sort(times);files=files(idx);
            inside=find(times>=t1 & times<t2);prior=find(times<t1,1,'last');
            candidates=[inside(:);prior(:)].';
            rec=struct('EventID',char(T.EventID(ie)),'Spacecraft',ic,'Product',char(products(ip)),...
                'Status','未确认','FilesStartingInside',numel(inside),'CheckedFiles',0,...
                'EvidenceFile','','EvidenceUTC','','Variable',v,'EpochVariable',ep,...
                'Errors',{{}},'Note','抽查CDF有效记录；不代表完整连续或共同覆盖');
            for jf=candidates
                f=files{jf};rec.CheckedFiles=rec.CheckedFiles+1;
                try
                    inds=0:31;
                    if times(jf)<t1
                        tt=spdfcdfread(f,'Variables',{ep},'CombineRecords',true);if iscell(tt),tt=tt{1};end
                        hit=find(tt>=t1 & tt<t2,32,'first');
                        if isempty(hit),continue;end
                        inds=hit(:)'-1;
                    end
                    a=spdfcdfread(f,'Variables',{ep,v},'Records',inds,'CombineRecords',false,'ConvertEpochToDatenum',true);
                    for jr=1:size(a,1)
                        t=double(a{jr,1});val=double(a{jr,2});
                        if isfinite(t) && t>=t1 && t<t2 && all(isfinite(val(:))) && all(abs(val(:))<1e29)
                            rec.Status='有有效记录';rec.EvidenceFile=f;rec.EvidenceUTC=datestr(t,'yyyy-mm-dd HH:MM:SS.FFF');break;
                        end
                    end
                    if strcmp(rec.Status,'有有效记录'),break;end
                catch ME
                    rec.Errors{end+1}=[f ' : ' ME.message];
                end
            end
            if isempty(candidates),rec.Status='本地未找到候选CDF';end
            if isempty(entries), entries=rec; else, entries(end+1)=rec; end %#ok<AGROW>
        end
    end
    fid=fopen(fullfile(out,'cdf_presence.json'),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(struct('entries',entries,'processed',ie,'total',height(T),...
        'complete',ie==height(T),'checkedUTC',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')))));
    fclose(fid);
    fprintf('%s: valid %d/16\n',T.EventID(ie),sum(strcmp({entries(end-15:end).Status},'有有效记录')));
end
end
