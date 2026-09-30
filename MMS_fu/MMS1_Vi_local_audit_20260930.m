%% MMS1 磁尾图件离子流速本地完整性审计（2026-09-30）
% 只读原始 CDF 和原图记录；使用原绘图的 IRFU 接口及坐标转换。
% 时间 UTC，速度 km/s；不平滑、不补缺、不新增质量筛选。
% 结果为文件/窗口统计与时间缺口，写入 Z 盘 derived；不保存全量速度副本。
clearvars -except AuditPilot
if ~exist('AuditPilot','var'), AuditPilot=false; end
%% 路径及现成 IRFU 接口
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
ParentDir='Z:\SPART-WORK\Data\MMS';
RecordDir=fullfile(ParentDir,'derived','MMS1_tail_20260720_0815');
AuditDir=fullfile(RecordDir,'Vi_completeness_20260930','local');
if AuditPilot, AuditDir=fullfile(AuditDir,'pilot'); end
if ~isfolder(AuditDir), mkdir(AuditDir); end
addpath(IRFDir); irf('check_path');
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
clear global MMS_DB
mms.db_init('local_file_db',[ParentDir filesep]);
mms.db_init('db_cache_size_max',2048); mms.db_init('db_cache_enabled',true);
Windows=jsondecode(fileread(fullfile(RecordDir,'windows.json')));
Files=[dir(fullfile(ParentDir,'mms1','fpi','fast','l2','dis-moms','2026','07','*.cdf'));...
    dir(fullfile(ParentDir,'mms1','fpi','fast','l2','dis-moms','2026','08','*.cdf'))];
[~,FileOrder]=sort({Files.name}); Files=Files(FileOrder);
Files=Files(startsWith({Files.name},'mms1_fpi_fast_l2_dis-moms_'));
if AuditPilot, Files=Files(1); Windows=Windows(2); end
%% 逐文件直接读 CDF：明确采样时间、有效点与原始积分时间宽度
FileRows=cell(numel(Files),1); AllGSE=zeros(0,4); AllHalfWidth=zeros(0,1);
for jf=1:numel(Files)
    FilePath=fullfile(Files(jf).folder,Files(jf).name);
    F=struct('name',Files(jf).name,'path',FilePath,'bytes',Files(jf).bytes,...
        'error','','rows',0,'anyFinite',0,'allFinite',0,'partialFinite',0,...
        'firstUTC','','lastUTC','','firstValidUTC','','lastValidUTC','',...
        'cadenceMedianSeconds',NaN,'cadenceMinSeconds',NaN,'cadenceMaxSeconds',NaN,...
        'integrationWidthMinSeconds',NaN,'integrationWidthMaxSeconds',NaN,...
        'integrationMinusUnits','','integrationPlusUnits','','sampleSpanGapsOverCadence',0);
    try
        dobj=dataobj(FilePath);
        Raw=get_variable(dobj,'mms1_dis_bulkv_gse_fast');
        assert(~isempty(Raw),'CDF 中缺少 mms1_dis_bulkv_gse_fast');
        Vi1_ts=mms.variable2ts(Raw); V1=irf.ts2mat(Vi1_ts);
        ValidAny=any(isfinite(V1(:,2:4)),2); ValidAll=all(isfinite(V1(:,2:4)),2);
        F.rows=size(V1,1); F.anyFinite=sum(ValidAny); F.allFinite=sum(ValidAll);
        F.partialFinite=sum(ValidAny & ~ValidAll);
        F.firstUTC=char(datetime(V1(1,1),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
        F.lastUTC=char(datetime(V1(end,1),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
        if any(ValidAny)
            tv=V1(ValidAny,1);
            F.firstValidUTC=char(datetime(tv(1),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
            F.lastValidUTC=char(datetime(tv(end),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
        end
        dt=diff(V1(:,1)); F.cadenceMedianSeconds=median(dt);
        F.cadenceMinSeconds=min(dt); F.cadenceMaxSeconds=max(dt);
        F.sampleSpanGapsOverCadence=sum(dt>F.cadenceMedianSeconds+0.001);
        Minus=Raw.DEPEND_0.DELTA_MINUS_VAR; Plus=Raw.DEPEND_0.DELTA_PLUS_VAR;
        F.integrationMinusUnits=Minus.UNITS; F.integrationPlusUnits=Plus.UNITS;
        dm=double(Minus.data); dp=double(Plus.data);
        if strcmp(Minus.UNITS,'ms'), dm=dm/1000; else, assert(strcmp(Minus.UNITS,'s')); end
        if strcmp(Plus.UNITS,'ms'), dp=dp/1000; else, assert(strcmp(Plus.UNITS,'s')); end
        width=dm(:)+dp(:);
        if isscalar(width), width=repmat(width,size(V1,1),1); end
        assert(numel(width)==size(V1,1));
        F.integrationWidthMinSeconds=min(width); F.integrationWidthMaxSeconds=max(width);
        % IRFU 将粒子时间移到积分中心；覆盖范围采用原 CDF 时宽的一半。
        % 此处只测量覆盖与空档，不删除任何原始有效样本。
        AllGSE=[AllGSE;V1]; AllHalfWidth=[AllHalfWidth;width/2]; %#ok<AGROW>
    catch ME
        F.error=getReport(ME,'extended','hyperlinks','off');
    end
    FileRows{jf}=F;
    fprintf('CDF %d/%d %s rows=%d valid=%d error=%d\n',jf,numel(Files),F.name,F.rows,F.anyFinite,~isempty(F.error));
end
FileRows=vertcat(FileRows{:});
fid=fopen(fullfile(AuditDir,'cdf_files.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(FileRows)); fclose(fid);
writetable(struct2table(FileRows,'AsArray',true),fullfile(AuditDir,'cdf_files.csv'));
[~,ix]=sort(AllGSE(:,1)); AllGSE=AllGSE(ix,:); AllHalfWidth=AllHalfWidth(ix);
%% 原生有效采样的积分覆盖并集：记录所有大于 1 ms 的空档
% 1 ms 仅消除双精度时标舍入，不作为科学数据质量门槛。
Valid=any(isfinite(AllGSE(:,2:4)),2); tv=AllGSE(Valid,1); hw=AllHalfWidth(Valid);
lo=tv-hw; hi=tv+hw; Coverage=zeros(0,2);
if ~isempty(lo)
    begin=lo(1); finish=hi(1);
    for j=2:numel(lo)
        if lo(j)>finish+0.001
            Coverage(end+1,:)=[begin finish]; begin=lo(j); finish=hi(j); %#ok<SAGROW>
        else
            finish=max(finish,hi(j));
        end
    end
    Coverage(end+1,:)=[begin finish];
end
%% 逐窗口重复原 IRFU 读取及 GSM 转换，并核对原图计数和逐 CDF 样本
WindowRows=cell(numel(Windows),1); GapRows=cell(0,1);
for iw=1:numel(Windows)
    Event=Windows(iw); tint=irf.tint([Event.startUTC '/' Event.endUTC]);
    tStart=tint.start.epochUnix; tStop=tint.stop.epochUnix;
    Previous=jsondecode(fileread(fullfile(RecordDir,[Event.id '_overview.json'])));
    W=struct('id',Event.id,'startUTC',Event.startUTC,'endUTC',Event.endUTC,...
        'error','','oldBurstAnyFinite',Previous.nativeCounts(2,1),'oldFastAnyFinite',Previous.nativeCounts(2,2),...
        'burstRows',0,'burstGseAnyFinite',0,'burstGsmAnyFinite',0,...
        'fastRows',0,'fastGseAnyFinite',0,'fastGsmAnyFinite',0,'fastGseAllFinite',0,'fastGsmAllFinite',0,...
        'directRows',0,'directAnyFinite',0,'databaseMatchesDirectTimes',false,'databaseMatchesDirectValues',false,...
        'matchesOriginalCounts',false,'coordinateValidCountLoss',0,...
        'coverageSeconds',0,'uncoveredSeconds',0,'gapCount',0,'maxGapSeconds',0,...
        'gapCategory','','firstValidUTC','','lastValidUTC','');
    try
        for im=1:2
            mode='brst'; if im==2, mode='fast'; end
            Vi1_ts=mms.get_data(['Vi_gse_fpi_' mode '_l2'],tint,1);
            GSE=zeros(0,4); GSM=zeros(0,4);
            if ~isempty(Vi1_ts)
                Vi1_ts=Vi1_ts.tlim(tint); GSE=irf.ts2mat(Vi1_ts);
                gsmVi1_ts=irf_gse2gsm(Vi1_ts); GSM=irf.ts2mat(gsmVi1_ts);
            end
            anyGse=sum(any(isfinite(GSE(:,2:4)),2)); anyGsm=sum(any(isfinite(GSM(:,2:4)),2));
            if im==1
                W.burstRows=size(GSE,1); W.burstGseAnyFinite=anyGse; W.burstGsmAnyFinite=anyGsm;
            else
                W.fastRows=size(GSE,1); W.fastGseAnyFinite=anyGse; W.fastGsmAnyFinite=anyGsm;
                W.fastGseAllFinite=sum(all(isfinite(GSE(:,2:4)),2)); W.fastGsmAllFinite=sum(all(isfinite(GSM(:,2:4)),2));
                Direct=AllGSE(AllGSE(:,1)>=tStart & AllGSE(:,1)<=tStop,:);
                W.directRows=size(Direct,1); W.directAnyFinite=sum(any(isfinite(Direct(:,2:4)),2));
                W.databaseMatchesDirectTimes=isequal(GSE(:,1),Direct(:,1));
                W.databaseMatchesDirectValues=isequaln(GSE,Direct);
                if anyGse>0
                    vv=GSE(any(isfinite(GSE(:,2:4)),2),1);
                    W.firstValidUTC=char(datetime(vv(1),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
                    W.lastValidUTC=char(datetime(vv(end),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
                end
            end
        end
        W.matchesOriginalCounts=isequal([W.burstGsmAnyFinite W.fastGsmAnyFinite],Previous.nativeCounts(2,:));
        W.coordinateValidCountLoss=W.burstGseAnyFinite+W.fastGseAnyFinite-W.burstGsmAnyFinite-W.fastGsmAnyFinite;
        LocalCoverage=Coverage(Coverage(:,2)>tStart & Coverage(:,1)<tStop,:);
        LocalCoverage(:,1)=max(LocalCoverage(:,1),tStart); LocalCoverage(:,2)=min(LocalCoverage(:,2),tStop);
        W.coverageSeconds=sum(LocalCoverage(:,2)-LocalCoverage(:,1));
        W.uncoveredSeconds=tStop-tStart-W.coverageSeconds;
        cursor=tStart; Gaps=zeros(0,2);
        for j=1:size(LocalCoverage,1)
            if LocalCoverage(j,1)>cursor+0.001, Gaps(end+1,:)=[cursor LocalCoverage(j,1)]; end %#ok<SAGROW>
            cursor=LocalCoverage(j,2);
        end
        if tStop>cursor+0.001, Gaps(end+1,:)=[cursor tStop]; end
        W.gapCount=size(Gaps,1);
        if ~isempty(Gaps), W.maxGapSeconds=max(Gaps(:,2)-Gaps(:,1)); end
        if W.fastGsmAnyFinite+W.burstGsmAnyFinite==0
            W.gapCategory='entire_window_blank';
        elseif W.gapCount>0
            W.gapCategory='partial_gap';
        else
            W.gapCategory='covered';
        end
        for j=1:size(Gaps,1)
            GapRows{end+1,1}=struct('id',Event.id,...
                'startUTC',char(datetime(Gaps(j,1),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")),...
                'endUTC',char(datetime(Gaps(j,2),'ConvertFrom','posixtime','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")),...
                'durationSeconds',Gaps(j,2)-Gaps(j,1)); %#ok<SAGROW>
        end
    catch ME
        W.error=getReport(ME,'extended','hyperlinks','off');
    end
    WindowRows{iw}=W;
    fid=fopen(fullfile(AuditDir,[Event.id '.json']),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(W)); fclose(fid);
    fprintf('WINDOW %d/%d %s old=%d new=%d direct=%d match=%d gaps=%d error=%d\n',iw,numel(Windows),W.id,W.oldFastAnyFinite,W.fastGsmAnyFinite,W.directAnyFinite,W.matchesOriginalCounts,W.gapCount,~isempty(W.error));
end
WindowRows=vertcat(WindowRows{:}); GapRows=vertcat(GapRows{:});
fid=fopen(fullfile(AuditDir,'windows_audit.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(WindowRows)); fclose(fid);
writetable(struct2table(WindowRows,'AsArray',true),fullfile(AuditDir,'windows_audit.csv'));
fid=fopen(fullfile(AuditDir,'window_gaps.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(GapRows)); fclose(fid);
writetable(struct2table(GapRows,'AsArray',true),fullfile(AuditDir,'window_gaps.csv'));
%% 汇总及审计依据
Summary=struct('auditUTC',char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss'Z'")),...
    'fileCount',numel(FileRows),'fileReadErrors',sum(~cellfun('isempty',{FileRows.error})),...
    'fileRows',sum([FileRows.rows]),'fileAnyFinite',sum([FileRows.anyFinite]),'filePartialFinite',sum([FileRows.partialFinite]),...
    'windows',numel(WindowRows),'windowReadErrors',sum(~cellfun('isempty',{WindowRows.error})),...
    'originalCountMismatches',sum(~[WindowRows.matchesOriginalCounts]),...
    'databaseDirectTimeMismatches',sum(~[WindowRows.databaseMatchesDirectTimes]),...
    'databaseDirectValueMismatches',sum(~[WindowRows.databaseMatchesDirectValues]),...
    'coordinateValidCountLoss',sum([WindowRows.coordinateValidCountLoss]),...
    'burstValidSamples',sum([WindowRows.burstGsmAnyFinite]),'fastValidSamples',sum([WindowRows.fastGsmAnyFinite]),...
    'entirelyBlankWindows',sum(strcmp({WindowRows.gapCategory},'entire_window_blank')),...
    'partialGapWindows',sum(strcmp({WindowRows.gapCategory},'partial_gap')),...
    'coveredWindows',sum(strcmp({WindowRows.gapCategory},'covered')),...
    'coverageHours',sum([WindowRows.coverageSeconds])/3600,'uncoveredHours',sum([WindowRows.uncoveredSeconds])/3600,...
    'cadenceSecondsRange',[min([FileRows.cadenceMedianSeconds]) max([FileRows.cadenceMedianSeconds])],...
    'notes','Finite counts follow original any-component rule; direct CDF and DB values compared exactly. Coverage is union of finite sample integration supports using CDF delta variables; 1 ms tolerance only for numeric time representation. No scientific quality filter or interpolation.');
fid=fopen(fullfile(AuditDir,'summary.json'),'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(Summary)); fclose(fid);
disp(Summary); fprintf('LOCAL_VI_AUDIT_FINISHED\n');
