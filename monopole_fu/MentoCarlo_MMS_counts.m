clear; clc;
global ParentDir OutputDir
%% 
ParentDir = '/Volumes/SPART-WORK/Data/MMS/'; 
DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
TempDir     = [DownloadDir, 'temp/']; 
if ~isfolder(TempDir)
    mkdir(TempDir);
end

%% 
Date = '2016-01-01T00:00:00.000Z/2024-01-01T00:00:00.000Z';
ic  = 1:4; 

load([pwd,'/NameTags.mat']);

splitDate = regexp(Date,'/','split');
startDate = splitDate{1};
endDate   = splitDate{2};
startDT = datetime(startDate, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
endDT   = datetime(endDate,   'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
NameTagsDT = datetime(NameTags, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
mask = (NameTagsDT >= startDT) & (NameTagsDT <= endDT);
NameTags = NameTags(mask);

%% 
splitDate{1} = erase(splitDate{1},{':','.'});splitDate{2} = erase(splitDate{2},{':','.'});
OutputDir = [DownloadDir,'TetraVolume/',splitDate{1},'To',splitDate{2},'/'];
if ~isfolder(OutputDir)
    mkdir(OutputDir);
end

% summaryFile = [OutputDir, 'TetraVolumeSummary.csv'];   % 结果汇总
% detailDir   = [OutputDir, 'Detail/'];                  % 存放每个时间段的详细数据
% if ~isfolder(detailDir)
%     mkdir(detailDir);
% end
% 
% % 若已有旧文件，可以根据需要选择是否删除
% if exist(summaryFile, 'file')
%     delete(summaryFile);
% end
% 
% % 结果表头
% fidSummary = fopen(summaryFile, 'w', 'n', 'UTF-8');
% fprintf(fidSummary, ['TimeTag, StartUTC, EndUTC, ',...
%     'Npts_SC1, Npts_SC2, Npts_SC3, Npts_SC4, ', ...
%     'Vmean_km3, Vmin_km3, Vmax_km3, Vstd_km3\n']);
%%
mms.db_init('local_file_db', ParentDir);
formatDate = @(s) [s(2:5), '-', s(6:7), '-', s(8:9), 'T', ...
                   s(10:11), ':', s(12:13), ':', s(14:15), '.000Z'];

%% 
%1
for TDT = 1:length(NameTags)-1
    clc;clear B1 B2 B3 B4 R1 R2 R3 R4
    fprintf(['当前处理时间段: ', NameTags{TDT}(1:end-1), ' ~ ', ...
            NameTags{TDT+1}(1:end-1), '\n']);
    
    tempDate = [char(NameTags(TDT)) '/' char(NameTags(TDT+1))];
    temptint     = irf.tint(tempDate);
    
    try
        c_eval("B?_ts = mms.get_data('B_gsm_brst', temptint, ?);", 1);
        tint = B1_ts.time([1,end]);

    
    c_eval("B?_ts = mms.get_data('B_gsm_brst', tint, ?);", ic);
    c_eval("R?_ts = mms.get_data('R_gsm', tint, ?);", ic);
    c_eval('B? = irf.ts2mat(B?_ts);',ic)
    c_eval('R? = irf.ts2mat(R?_ts);',ic)
    c_eval('B? = irf_resamp(B?, B1);',2:4)
    c_eval('R? = irf_resamp(R?, B1);',ic)
    RR = [R1(1,2:4);R2(1,2:4);R3(1,2:4);R4(1,2:4);...
        R1(end,2:4);R2(end,2:4);R3(end,2:4);R4(end,2:4);];

    N = numel(B1_ts.time.epoch);

    catch
        writematrix([NameTags{TDT}(1:end-1), ' 的 B 数据导入出现问题'], ...
            [OutputDir, 'errorlog.txt'], ...
            'WriteMode','append', 'Encoding','UTF-8');
        continue;
    end

    %2
    try
    tri = delaunayTriangulation(RR);
    [~, V] = convexHull(tri);
    catch
    writematrix([NameTags{TDT}(1:end-1), ' 的三角剖分出现问题'], ...
            [OutputDir, 'errorlog.txt'], ...
            'WriteMode','append', 'Encoding','UTF-8');
        continue
    end
    
    writematrix([NameTags{TDT}(1:end-1),' 计数为',num2str(N),'， 体积为', num2str(V), ' km^3']...
        ,[OutputDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
end
