clear;clc
close all
% created by fwd on 2021.12.27
%%
global ParentDir 
ParentDir = 'E:\THEMIS\'; 
TempDir = 'E:\THEMIS\temp\'; 
if ~isfolder(TempDir)
    mkdir(TempDir);
end
%% Time Option
TT = '2011-02-01\2011-04-30';
% 2010 0426
% TT = '2009-03-01T06:00:00Z\2009-03-01T08:00:00Z';
Datelist = regexp(TT,'\d+-\d+-\d+','match');
TaskDir = [ParentDir,Datelist{1},'T',Datelist{2},'\']; mkdir(TaskDir)
Datelist = datenum(Datelist,'yyyy-mm-dd');
Datelist = str2num(datestr(Datelist(1):Datelist(2),'yyyymmdd'));
%% Files Download
ic = {'a','e'};
for tempDate = 1:length(Datelist)
    tmpDay = num2str(Datelist(tempDate));
    fprintf(['Begin ',tmpDay,':\n'])
    for i = 1:length(ic)
        THEMISDownload(tmpDay,['th',ic{i}],'fgm',TempDir)
        THEMISDownload(tmpDay,['th',ic{i}],'ssc',TempDir)
    end
    THEMISDataMove(TempDir,ParentDir)
%% Case Selection
    flag = [];
    tint = iso2epoch([tmpDay(1:4),'-',tmpDay(5:6),'-',tmpDay(7:8),'T','00:00:00Z']) + [0 86400];
    try
    c_eval("B_? = th_read_l2_change_by_fwd('th?_fgl_gsm',tint);",ic);
%     c_eval("Vi_? = th_read_l2_change_by_fwd('th?_peir_velocity_gsm',tint);",ic);
    c_eval("R_? = THEMISLocation('?',str2num(TT(1:4)),ParentDir);",ic)
    catch
        writematrix(['Load data Failed at: ',datestr(datenum(1970,1,1,0,0,0)+mean(tint)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
continue
    end
    
    if length(B_a) ~= length(B_e)
        try
        B_e = irf_resamp(B_e,B_a);
        catch
        B_a = irf_resamp(B_a,B_e);
        end
    end
    
    if length(R_a) ~= length(R_e)
        try
        R_e = irf_resamp(R_e,R_a);
        catch
        R_a = irf_resamp(R_a,R_e);
        end
    end

    for tmpB = 1:length(B_a)-2800
        clc
        flagB = [];
        
        disp(['Date:',tmpDay])
        disp(['B:',num2str(tmpB),'/',num2str(length(B_a)-2800)])
        disp([repmat('■',1,round(10*tmpB/(length(B_a)-2800))),repmat('□',1,10-round(10*tmpB/(length(B_a)-2800)))])
        
        if mean(abs(B_a(tmpB:tmpB+2400,4))) >=  mean(abs(B_a(tmpB:tmpB+2400,3))) && ...
                mean(abs(B_a(tmpB:tmpB+2400,4))) >=  mean(abs(B_a(tmpB:tmpB+2400,2)))
            flagB = [flagB B_a(tmpB,1)];
        elseif mean(abs(B_e(tmpB:tmpB+2400,4))) >=  mean(abs(B_e(tmpB:tmpB+2400,3))) && ...
                mean(abs(B_e(tmpB:tmpB+2400,4))) >=  mean(abs(B_e(tmpB:tmpB+2400,2)))
            flagB = [flagB B_e(tmpB,1)];
        end
    end
    
    for tmp = 1:length(flagB)
        c_eval('flagR? = find(R_?(:,1)==round(flagB(tmp)/60)*60);',ic)
        c_eval('Loc? = R_?(flagR?,:);',ic);
    
        if ~isempty(flag) && abs(Loca(tmp,1)-flag(end)) <= 600
            continue
        else
            if abs(Loca(tmp,4)-Loce(tmp,4)) > abs(Loca(tmp,3)-Loce(tmp,3)) && abs(Loca(tmp,2)-Loce(tmp,2)) <= 1
                flag = [flag,Loca(tmp,1)];
            end
        end
    end

%% Plot
    if ~isempty(flag)
            for i = 1:length(ic)
                THEMISDownload(tmpDay,['th',ic{i}],'efi',TempDir)
                THEMISDownload(tmpDay,['th',ic{i}],'esa',TempDir)
            end
            THEMISDataMove(TempDir,ParentDir)
        for len = 1:length(flag)
            tt = flag(len)+[-3600 3600];
            try
                THEMISPlot('a',tt,TaskDir);
                THEMISPlot('e',tt,TaskDir);
            catch
                writematrix(['Plot Failed at: ',datestr(datenum(1970,1,1,0,0,0)+flag(len)/86400,'yyyymmdd HH:MM:SS.FFF')],...
                    [TaskDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
            end
            disp("Successful Find Eligible Case ~~~///(^v^)\\\~~~")    
        end
    else
        clc
        disp(['No case in: ',tmpDay])
    end
    try
        cd(TempDir)
        rmdir(TempDir,'s');
    catch
        fprintf('Failed to distory the temp folder\n')
    end
end