%------written by Wending Fu, May.2025 in Beijing------------
clear;clc;
global ParentDir 
% global OutputDir
ParentDir = '/gs/home/by2430102/MMS/'; 
DownloadDir = ParentDir;
TempDir = [DownloadDir,'temp/']; mkdir(TempDir);
 
ic = 1:4; iic = 1:4;
load(['/gs/home/by2430102/fwd_matlab_patch/monopole_fu/NameTags.mat']);

% Date = '2019-09-24/2019-09-25';
Date = '2017-01-01T05:16:00.000Z/2017-01-01T06:28:00.000Z';
splitDate = regexp(Date,'/','split');
startDate = splitDate{1};
endDate   = splitDate{2};
startDT = datetime(startDate, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
endDT   = datetime(endDate,   'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
NameTagsDT = datetime(NameTags, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSS'Z'", 'TimeZone', 'UTC');
mask = (NameTagsDT >= startDT) & (NameTagsDT <= endDT);
NameTags = NameTags(mask);

%% output directory & DB init
splitDate{1} = erase(splitDate{1},{':','.'});
splitDate{2} = erase(splitDate{2},{':','.'});
OutputDir = [ParentDir,'Monopole_Search/',splitDate{1},'To',splitDate{2},'/'];
clear mask NameTagsDT
mkdir(OutputDir)
mms.db_init('local_file_db',ParentDir);

units = irf_units;
if isempty(gcp('nocreate'))
parpool('Processes', 8); 
end

for TDT = 1:length(NameTags)-1  % 主循环
    clc; clear B1 B2 B3 B4 R1 R2 R3 R4 Pos;
    disp(NameTags(TDT))
    PlotFlag = 0;
    
    tempDate = [char(NameTags(TDT)) '/' char(NameTags(TDT+1))];
    tempTint = irf.tint(tempDate);

    %% 尝试导入位置数据
    tryTimes = 0;
    while tryTimes <= 10
        try
            mms.db_init('local_file_db',ParentDir);
            Pos = mms.get_data('R_gsm',tempTint,1);
            Pos = irf.ts2mat(Pos);
            tryTimes = 666;
        catch
            tryTimes = tryTimes + 1;
        end
    end

    if tryTimes ~= 666
        % 记录数据导入1错误
        msg = [NameTags{TDT}(2:end-2),'的数据导入1出现问题'];
        fid = fopen([OutputDir,'errorlog.txt'], 'a', 'n', 'UTF-8');
        if fid == -1
            warning('无法打开 %s 进行写入', [OutputDir,'errorlog.txt']);
        else
            fprintf(fid, '%s\n', msg);
            fclose(fid);
        end
        continue
    end

    %% 数据处理与判据
    try
        tint = irf.tint(irf_time(Pos(1,1),'epoch>epochTT'), irf_time(Pos(end,1),'epoch>epochTT'));
        B1_ts = mms.get_data('B_gsm_brst',tint,1);
        B2_ts = mms.get_data('B_gsm_brst',tint,2);
        B3_ts = mms.get_data('B_gsm_brst',tint,3);
        B4_ts = mms.get_data('B_gsm_brst',tint,4);
        B1 = irf.ts2mat(B1_ts);
        B2 = irf.ts2mat(B2_ts);
        B3 = irf.ts2mat(B3_ts);
        B4 = irf.ts2mat(B4_ts);
        c_eval('B? = irf_resamp(B?,B1);');
        
        Pos = mms.get_data('R_gsm',tint);
        R1 = Pos.gsmR1; R2 = Pos.gsmR2; R3 = Pos.gsmR3; R4 = Pos.gsmR4;
        R1 = [Pos.time.epochUnix, R1(:,1:3)];
        R2 = [Pos.time.epochUnix, R2(:,1:3)];
        R3 = [Pos.time.epochUnix, R3(:,1:3)];
        R4 = [Pos.time.epochUnix, R4(:,1:3)];
        c_eval('R? = irf_resamp(R?,B1);');

        CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
        c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');
        
        N = size(B1,1);
        LocPoint = nan(N,3); LocRes = cell(N,1);
        Q = nan(N,1); resQ = cell(N,1);
        Qerror = 1000*ones(N,1);
        dLoc = 1000*ones(N,15);

        gradB = c_4_grad(R1,R2,R3,R4,B1,B2,B3,B4,'grad');
        divB = [gradB(:,1), sum([gradB(:,2),gradB(:,6),gradB(:,10)],2)];
        
        RR12 = irf_abs(R1-R2); RR13 = irf_abs(R1-R3); RR14 = irf_abs(R1-R4);
        RR23 = irf_abs(R2-R3); RR24 = irf_abs(R2-R4); RR34 = irf_abs(R3-R4);
        RR_mean = (RR12(:,5)+RR13(:,5)+RR14(:,5)+RR23(:,5)+RR24(:,5)+RR34(:,5))/6;
        if ceil(max(log10(RR_mean))) > 3, continue; end

        id = nchoosek(1:6,2);
        warning('off')
        tic
        % parfor_progress(N-1);
        parfor i = 1:200
            % parfor_progress;
            [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError(R1,R2,R3,R4,B1,B2,B3,B4,i,i*sign(divB(i,2)),RR_mean(i),1);
            tempd = irf_abs(LocRes{i}(id(:,1),:)-LocRes{i}(id(:,2),:));
            dLoc(i,:) = (tempd(:,4)/RR_mean(i))'*100;
            if ~isnan(resQ{i})
                Qerror(i) = abs(100*std(resQ{i})/Q(i));
            end
        end
        % parfor_progress(0);
        toc

        meand = mean(dLoc,2);
        LocPoint = irf_abs(LocPoint);

        %% 写入 DateList
        msg = [irf_time(B1(1),'epoch>utc'),' ', num2str(min(meand)), ' ', num2str(min(Qerror)), ' ', num2str(min(LocPoint(:,4)))];
        fid = fopen([OutputDir,'DateList.txt'], 'a', 'n', 'UTF-8');
        if fid == -1
            warning('无法打开 %s 进行写入', [OutputDir,'DateList.txt']);
        else
            fprintf(fid, '%s\n', msg);
            fclose(fid);
        end
        
        TimeUTC = erase(irf_time(B1(1),'epoch>utc'),{':','.'});
        save([OutputDir, TimeUTC, '.mat'], 'Q','resQ','LocPoint','LocRes');

        %% 符合判据则写入 CaseList
        try
            tag = find(meand <= 66.1 & Qerror <= 100 & LocPoint(:,4) <= 3*mean(RR_mean));
            if ~isempty(tag)
                PlotFlag = 1;
                flagTimeUTC = erase(irf_time(B1(tag(1)),'epoch>utc'),{':','.'});
                msg = [flagTimeUTC,' ', num2str(meand(tag)), ' ', num2str(Qerror(tag)), ' ', num2str(LocPoint(tag,4))];
                fid = fopen([OutputDir,'CaseList.txt'], 'a', 'n', 'UTF-8');
                if fid == -1
                    warning('无法打开 %s 进行写入', [OutputDir,'CaseList.txt']);
                else
                    fprintf(fid, '%s\n', msg);
                    fclose(fid);
                end
            end
        catch
            % 记录 caselist 错误
            msg = [NameTags{TDT}(2:end-2),'的caselist导入出现问题'];
            fid = fopen([OutputDir,'errorlog.txt'], 'a', 'n', 'UTF-8');
            if fid == -1
                warning('无法打开 %s 进行写入', [OutputDir,'errorlog.txt']);
            else
                fprintf(fid, '%s\n', msg);
                fclose(fid);
            end
        end
    catch
        % 记录数据导入2错误
        msg = [NameTags{TDT}(2:end-2),'的数据导入2出现问题'];
        fid = fopen([OutputDir,'errorlog.txt'], 'a', 'n', 'UTF-8');
        if fid == -1
            warning('无法打开 %s 进行写入', [OutputDir,'errorlog.txt']);
        else
            fprintf(fid, '%s\n', msg);
            fclose(fid);
        end
        continue
    end
end
