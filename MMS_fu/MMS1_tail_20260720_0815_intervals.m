%% MMS1磁尾四小时分段：2026-07-20至2026-08-15，UTC / GSM
% 只读取官方MEC轨道；X=0穿越用相邻轨道点作线性插值。
% 每次进入X<0后每4h分段，末尾不足4h保留；与请求日期范围取交集。
clear;clc
%% 路径及IRFU初始化
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
ParentDir='Z:\SPART-WORK\Data\MMS\';
RecordDir=[ParentDir 'derived\MMS1_tail_20260720_0815'];
OutputDir='C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\MMS1_tail_20260720_0815';
addpath(CodeDir,IRFDir);irf('check_path');
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
mms.db_init('local_file_db',ParentDir);
if ~isfolder(OutputDir),mkdir(OutputDir);end
%% load R [GSM, km]
tint=irf.tint('2026-07-19T00:00:00Z/2026-08-17T23:59:59Z');
R1_ts=mms.db_get_ts('mms1_mec_srvy_l2_epht89d','mms1_mec_r_gsm',tint);
R1=irf.ts2mat(R1_ts);
assert(all(isfinite(R1),'all'),'MEC轨道存在非有限数据，需检查');
assert(all(diff(R1(:,1))>0)&max(diff(R1(:,1)))<=60,'MEC轨道有缺口，需先补齐');
Range=irf.tint('2026-07-20T00:00:00Z/2026-08-16T00:00:00Z');
RangeStart=Range.start.epochUnix;RangeStop=Range.stop.epochUnix;
%% 找出进入及离开X<0的轨道段
Negative=R1(:,2)<0;
first=find(diff([false;Negative])==1);
last=find(diff([Negative;false])==-1);
Windows=struct([]);Passages=struct([]);iw=0;ip=0;
for j=1:numel(first)
    a=first(j);b=last(j);
    Entry=R1(a,1);Exit=R1(b,1);
    if a>1,Entry=interp1(R1([a-1 a],2),R1([a-1 a],1),0,'linear');end
    if b<size(R1,1),Exit=interp1(R1([b b+1],2),R1([b b+1],1),0,'linear');end
    if Exit<=RangeStart||Entry>=RangeStop,continue;end
    assert(a>1,'首个磁尾段进入时间超出已下载轨道范围，需前向补齐MEC');
    assert(b<size(R1,1),'最后磁尾段离开时间超出已下载轨道范围，需后向补齐MEC');
    ip=ip+1;
    Passage=struct('passage',ip,'entryEpoch',Entry,'exitEpoch',Exit,...
        'entryUTC',char(EpochTT(EpochUnix(Entry)).utc),'exitUTC',char(EpochTT(EpochUnix(Exit)).utc),...
        'entryBracketUTC',{cellstr(EpochTT(EpochUnix(R1([a-1 a],1))).utc)},...
        'exitBracketUTC',{cellstr(EpochTT(EpochUnix(R1([b b+1],1))).utc)});
    if ip==1,Passages=Passage;else,Passages(ip)=Passage;end
    for t0=Entry:14400:Exit-1e-6
        tStart=max(t0,RangeStart);tStop=min([t0+14400,Exit,RangeStop]);
        if tStop<=tStart,continue;end
        iw=iw+1;tMid=(tStart+tStop)/2;
        Rmid=irf_resamp(R1,tMid,'linear');
        assert(Rmid(2)<0,'分段中点不在磁尾');
        StartUTC=char(EpochTT(EpochUnix(tStart)).utc);StopUTC=char(EpochTT(EpochUnix(tStop)).utc);
        Tag=regexprep(StartUTC(1:19),'[-:T]','');
        Window=struct('number',iw,'id',sprintf('W%03d_%s',iw,Tag),'passage',ip,...
            'startUTC',StartUTC,'endUTC',StopUTC,'startEpoch',tStart,'endEpoch',tStop,...
            'durationHours',(tStop-tStart)/3600,'midpointXYZ_RE',Rmid(2:4)/6372,...
            'entryEpoch',Entry,'exitEpoch',Exit);
        if iw==1,Windows=Window;else,Windows(iw)=Window;end
    end
end
%% 保存分段清单；只保存必要元数据，不复制原始轨道
fid=fopen(fullfile(RecordDir,'windows.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(Windows));fclose(fid);
fid=fopen(fullfile(RecordDir,'passages.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(Passages));fclose(fid);
T=struct2table(Windows);writetable(T,fullfile(OutputDir,'MMS1_tail_windows.csv'));
fprintf('PASSAGES=%d WINDOWS=%d FIGURES=%d CADENCE=%gs\n',ip,iw,2*iw,median(diff(R1(:,1))));
disp(T(:,{'number','passage','startUTC','endUTC','durationHours'}));

