%% Case_for_SMILE.pptx中5个事件的AE指数
% 原PPT第2、3、4、5、8页明确标出的时间，UTC，不额外延长。
% 沿用原overview的直接脚本及IRFU接口；没有新增function。
% 输入：NASA OMNI 1分钟原始CDF，AE_INDEX，单位nT（Kyoto quicklook）。
% 输出：5张PNG、每个事件一页的PDF；原始CDF和运行记录保存在Z盘。
close all
clear;clc

%% 路径和事件
CodeDir='C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
IRFDir='C:\Users\Administrator\Documents\irfu-matlab-master';
DataDir='Z:\SPART-WORK\Data\MMS\ancillary\omni\hro_1min\2026';
RecordDir='Z:\SPART-WORK\Data\MMS\derived\Case_SMILE_AE_20260929';
OutputDir='C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\Case_SMILE_AE_20260929';
addpath(CodeDir,IRFDir);
irf('check_path');
setenv('CDF_LEAPSECONDSTABLE',fullfile(IRFDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));
for p={DataDir,RecordDir,OutputDir}
    if ~isfolder(p{1}),mkdir(p{1});end
end
TimeStart={'2026-07-21T00:00:00Z','2026-07-24T08:00:00Z','2026-07-28T05:10:00Z','2026-08-04T01:50:00Z','2026-08-11T06:50:00Z'};
TimeStop={'2026-07-21T04:00:00Z','2026-07-25T05:00:00Z','2026-07-28T06:00:00Z','2026-08-04T03:00:00Z','2026-08-11T08:30:00Z'};
SlideNumbers=[2 3 4 5 8];
MonthIndex=[1 1 1 2 2];
Files={'omni_hro_1min_20260701_v01.cdf','omni_hro_1min_20260801_v01.cdf'};
SourceRoot='https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/2026/';
PDFFile=fullfile(OutputDir,'Case_for_SMILE_AE.pdf');
AE_month=cell(1,2);
FileInfo=cell(1,2);

%% 下载和读取原始CDF：复用已下载文件
for im=1:2
    FilePath=fullfile(DataDir,Files{im});URL=[SourceRoot Files{im}];
    if ~isfile(FilePath)
        fprintf('DOWNLOAD %s\n',URL);
        PartFile=[FilePath '.part'];
        websave(PartFile,URL,weboptions('Timeout',120));
        movefile(PartFile,FilePath);
    end
    dobj=dataobj(FilePath);
    AE_ts=get_ts(dobj,'AE_INDEX');
    AE_month{im}=irf.ts2mat(AE_ts);
    Raw=get_variable(dobj,'AE_INDEX');Fill=getfillval(dobj,'AE_INDEX');
    Expected=double(Raw.data);Expected(Expected==double(Fill))=NaN;
    % AE为整数CDF变量：以原CDF的FILLVAL掩码在double矩阵中保留缺测。
    AE_month{im}(double(Raw.data)==double(Fill),2)=NaN;
    assert(isequaln(AE_month{im}(:,2),Expected(:)),'AE读取结果与CDF原变量不一致');
    assert(all(abs(diff(AE_month{im}(:,1))-60)<1e-3),'原CDF时间轴存在非1分钟间隔');
    Info=dir(FilePath);
    FileInfo{im}=struct('file',FilePath,'url',URL,'bytes',Info.bytes,'variable','AE_INDEX',...
        'units',Raw.UNITS,'fillValue',double(Fill),'cadenceSeconds',60,...
        'version','v01','indexStatus','Kyoto quicklook');
    fprintf('READ %s valid=%d missing=%d\n',Files{im},sum(isfinite(Expected)),sum(~isfinite(Expected)));
end

%% 各事件分别画AE
Results=cell(1,5);
for ie=1:5
    tint=irf.tint([TimeStart{ie} '/' TimeStop{ie}]);
    AE=irf_tlim(AE_month{MonthIndex(ie)},tint.epochUnix);
    ExpectedCount=round(diff(tint.epochUnix)/60); % irf_tlim采用[start,end)
    assert(size(AE,1)==ExpectedCount,'事件时间范围内的分钟记录数不完整');
    assert(any(isfinite(AE(:,2))),'事件没有有效AE值');
    fprintf('PLOT Case %d rows=%d missing=%d\n',ie,size(AE,1),sum(~isfinite(AE(:,2))));

    %% Init figure / AE plot：沿用原overview绘图格式
    fn=figure('Visible','off','Color','w','Position',[10 10 1280 430]);
    h=irf_subplot(1,1,1);
    set(h,'Position',[0.085 0.25 0.89 0.57]);
    irf_plot(h,AE,'color','k','LineWidth',1.1);
    irf_zoom(h,'x',tint);
    ylabel(h,'AE [nT]','FontSize',13);
    set(h,'FontName','Arial','FontSize',12,'TickDir','out','Box','on');
    grid(h,'off');
    axes(h);set(gca,"XTickLabelRotation",0);
    if ~strcmp(TimeStart{ie}(1:10),TimeStop{ie}(1:10))
        xlabel(h,[TimeStart{ie}(1:10) ' / ' TimeStop{ie}(1:10) ' UTC']);
    end
    MaxAE=max(AE(:,2),[],'omitnan');
    ylim(h,[0 max(100,ceil(1.05*MaxAE/100)*100)]);
    StartText=strrep(TimeStart{ie}(1:16),'T',' ');
    StopText=strrep(TimeStop{ie}(1:16),'T',' ');
    title(h,sprintf('Case %d: %s - %s UTC',ie,StartText,StopText),...
        'FontSize',14,'FontWeight','normal');
    annotation(fn,'textbox',[0.085 0.018 0.9 0.06],...
        'String','AE: WDC Kyoto quicklook via NASA OMNI | 1 min',...
        'FontName','Arial','FontSize',10,'Color',[0.35 0.35 0.35],...
        'EdgeColor','none','FitBoxToText','off');
    drawnow;

    %% 出图保存部分
    Base=sprintf('Case%d_AE_%s',ie,strrep(TimeStart{ie}(1:10),'-',''));
    PNGFile=fullfile(OutputDir,[Base '.png']);
    exportgraphics(fn,PNGFile,'Resolution',200,'BackgroundColor','white');
    exportgraphics(fn,PDFFile,'ContentType','vector','BackgroundColor','white','Append',ie>1);
    Results{ie}=struct('caseNumber',ie,'pptSlide',SlideNumbers(ie),...
        'startUTC',TimeStart{ie},'endUTC',TimeStop{ie},'units','nT',...
        'cadenceSeconds',60,'rowCount',size(AE,1),'validCount',sum(isfinite(AE(:,2))),...
        'missingCount',sum(~isfinite(AE(:,2))),'minAE',min(AE(:,2),[],'omitnan'),...
        'maxAE',MaxAE,'sourceFile',Files{MonthIndex(ie)},'png',PNGFile);
    close(fn);
end

%% 保存来源和运行记录；不转换或复制全量CDF为MAT/CSV
Record=struct('sourcePresentation','C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\Case_for_SMILE.pptx',...
    'sourceDescription','NASA OMNI_HRO_1MIN; AE from WDC Kyoto, quicklook',...
    'sourceDocumentation','https://omniweb.gsfc.nasa.gov/html/omni_min_data.html',...
    'processing','No smoothing or interpolation; CDF FILLVAL becomes NaN; irf_tlim uses start <= t < end; axes match PPT intervals',...
    'files',{FileInfo},'events',{Results},'pdf',PDFFile,'complete',true);
fid=fopen(fullfile(RecordDir,'AE_run_record.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(Record));fclose(fid);
fprintf('AE_CASES_FINISHED\n');
