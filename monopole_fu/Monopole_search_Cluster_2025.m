%------written by Wending Fu, May.2025 in Beijing------------
clear;clc;
cd /Volumes/SPART-NAS/Data/Cluster/
ParentDir = '/Volumes/SPART-NAS/Data/Cluster/';
ic=1:4;

% TT = '2002-03-17\2003-01-01';
TT = '2004-01-12\2004-01-13';
Datelist = regexp(TT,'\d+-\d+-\d+','match');
TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'\']; mkdir(TaskDir)
Datelist = datenum(Datelist,'yyyy-mm-dd');
Datelist = datestr(Datelist(1):Datelist(2),'yyyy-mm-dd');

mms.db_init('local_file_db',ParentDir);
%%
units = irf_units;
% parfor_progress(length(NameTags)-1);
for tempDate = 1:size(Datelist,1)-1 %This is a distinctive temp  (๑ˉ∀ˉ๑)
clc;clear B1 B2 B3 B4 R1 R2 R3 R4 RR;

Tsta = [Datelist(tempDate,:) 'T00:00:000Z'];
Tend = [Datelist(tempDate+1,:) 'T00:00:000Z'];
tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH
disp(['当前日期:' Datelist(tempDate,:)])

try 
    caa_load_changed_by_fwd('CL_SP_AUX',Tsta,Tend);
catch
    caa_download(tint,'CL_SP_AUX');
    caa_load_changed_by_fwd('CL_SP_AUX',Tsta,Tend);
end

R = c_caa_var_get('sc_r_xyz_gse__CL_SP_AUX','mat');
tint = [R(1,1) R(end,1)];
Tsta = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(1))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];
Tend = [datestr(datenum(1970,1,1,0,0,0)+mean(tint(2))/86400,'yyyy-mm-ddTHH:MM:SS.FFF') 'Z'];

% if min(sqrt(R(:,2).^2+R(:,3).^2+R(:,4).^2)) < 20000+6371
%     continue
% end
errorflag = 1;
try
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
catch
    try
    %    Magnetic fields
    c_eval("caa_download(tint,'C*_CP_FGM_FULL')",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    catch
        errorflag = 666;
    end
end

try
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
catch 
    try
    c_eval("caa_download(tint,'C*_CP_AUX_POSGSE_1M')",ic);  % position & velocity for each sc
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
    catch
    errorflag = 666;
    end
end
  
if errorflag == 666
    writematrix(['Data incompleted at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],[TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
continue
end

try
% % % dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
% % % varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
c_eval('B? = irf_gse2gsm(B?_gse);',ic);
% c_eval('B?=irf_abs(B?_gsm);',ic);
% c_eval('B? = irf_resamp(B?,B1);',2:4);
c_eval('B? = irf_resamp(B?,B1);',2:4);
catch
        writematrix(['Data incompleted at: ',datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
        continue
end
%% 判据
try
    c_eval('R?_gse = c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');',ic);
    c_eval('R? = irf_gse2gsm(R?_gse);',ic);
    c_eval('R? = irf_resamp(R?,B1);')

    c_eval('R? = irf_resamp(R?,B1);')
    CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
    c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

    idx = checkLineIntersect(R1,R2,R3,R4, B1,B2,B3,B4, 3, 1);
    c_eval('R?_copy = R?;');
    c_eval('B?_copy = B?;');
    c_eval('R? = R?(idx,:);');
    c_eval('B? = B?(idx,:);');
    
    LocPoint = zeros(length(B1),3)*nan;
    LocRes = cell(length(B1),1);
    Q = zeros(length(B1),1)*nan;
    resQ = cell(length(B1),1);
    
    Qerror = ones(length(B1),1)*1000;
    dLoc = ones(length(B1),15)*1000;
    
    % div
    gradB=c_4_grad(R1, R2, R3, R4, B1, B2, B3, B4, 'grad');
    divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度

    flag_m = 0;
    time_flagm = 0;    
    
    RR12 = irf_abs(R1-R2); RR13 = irf_abs(R1-R3); RR14 = irf_abs(R1-R4); 
    RR23 = irf_abs(R2-R3); RR24 = irf_abs(R2-R4); RR34 = irf_abs(R3-R4); 
    RR_mean = (RR12(:,5) + RR13(:,5) + RR14(:,5) + RR23(:,5) + RR24(:,5) + RR34(:,5))/6;
    MultiPower = ceil(max([log10(RR_mean)]));
    if MultiPower > 3, continue; end
    
    id = nchoosek(1:6,2);
    warning('off')

    % solve
    tic
    parfor_progress(size(B1, 1)-1);
    parfor i =1:size(B1, 1)
    % parfor i =1:200
    parfor_progress;
    [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError(R1(i,:), R2(i,:), R3(i,:), R4(i,:),...
        B1(i,:), B2(i,:), B3(i,:), B4(i,:), i,i*sign(divB(i,2)),RR_mean(i),1);

    tempd = irf_abs(LocRes{i}(id(:,1),:)-LocRes{i}(id(:,2),:));
    tempd = tempd(:,4)/RR_mean(i);
    tempd = tempd';
    dLoc(i,:) = tempd * 100;

    if ~isnan(resQ{i})
        Qerror(i) = abs(100*std(resQ{i})/Q(i));
    else
        Qerror(i) = 1000;
    end
    end
    parfor_progress(0);
    toc

    meand = mean(dLoc,2);
    LocPoint = irf_abs(LocPoint);

    writematrix([datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF'),' ',num2str(min(meand)),' ', num2str(min(Qerror))...
        ,' ', num2str(min(LocPoint(:,4)))],[TaskDir,'DateList.txt'],'WriteMode','append','Encoding','UTF-8')
    TimeUTC = datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF');TimeUTC = erase(TimeUTC,{':','.'});
    save([TaskDir, TimeUTC, '.mat'], 'Q', 'resQ', 'LocPoint', 'LocRes')
    
    try
    tag = find(meand <= 66.1 & Qerror <= 100 & LocPoint(:,4) <= 3*mean(RR_mean));
    if ~isempty(tag)
    PlotFlag = 1;
    flagTime = B1(tag, 1);
    flagTimeUTC = datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF');flagTimeUTC = erase(flagTimeUTC,{':','.'});
    writematrix([flagTimeUTC,' ',num2str(meand(tag)),' ', num2str(Qerror(tag))...
    ,' ', num2str(LocPoint(tag,4))],[TaskDir,'CaseList.txt'],'WriteMode','append','Encoding','UTF-8')
    end
    catch
    writematrix([datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF'),'的caselist导入出现问题'],[TaskDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    end
catch
    writematrix([datestr(datenum(1970,1,1,0,0,0)+R(1,1)/86400,'yyyymmdd HH:MM:SS.FFF'),'的数据导入2出现问题'],[TaskDir,'errorlog.txt'],...
        'WriteMode','append','Encoding','UTF-8')
    continue
end

%% 符合判据的继续下载并出图
% if PlotFlag == 1
% try
%     SDCFilesDownload(FileGroups{TDT},TempDir) 
%     SDCDataMove(TempDir,ParentDir)
%     mms.db_init('local_file_db',ParentDir);
%     PlotTint = irf_time([flagTime(end)-20,flagTime(end)+20],'epoch>epochTT');
%     id_flagTime = SDCPlot(PlotTint,desmoms,ic,NameTags{TDT},flagTime(end));
% catch
%     writematrix([irf_time(flagTime(end),'epoch>utc'),'的画图出现问题'],[OutputDir,'errorlog.txt'],...
%     'WriteMode','append','Encoding','UTF-8')
%     continue
% end
% end
%% 删除文件夹并生成记录文件
% % % try
% % %     cd(OutputDir)
% % %     rmdir(TempDir,'s');    
% % % catch
% % %     writematrix(['删除文件夹',NameTags{TDT}(2:end-2),'失败'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
% % % end
end