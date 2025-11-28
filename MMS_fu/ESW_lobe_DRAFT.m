%------written by Wending Fu, Sept.2025 in Beijing------------
clear;clc;
ParentDir = 'Z:\SPART-WORK\Data\MMS/'; 
OutputDir = 'C:\MMS\ESWSearch\2019-01-01To2025-01-01\';
ic = 1:4;
%%
CaseListDir = [OutputDir, 'caselist.txt'];
fid = fopen(CaseListDir,'r');
lines = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
lines = lines{1};
rawToken = regexp(lines, '(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+)', 'match');
rawTime = cellfun(@(x) x{1}, rawToken, 'UniformOutput', false);
timeMs = regexp(rawTime, '^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}', 'match', 'once');
CaseTimeList = strcat(timeMs, 'Z');
%%
mms.db_init('local_file_db',ParentDir);
for iCase = 1:length(CaseTimeList)
clc;
s1 = [num2str(iCase),'/',num2str(length(CaseTimeList))]; s2 = ['正在处理:',CaseTimeList{iCase}];
disp([s1,repmat('■',1,round(10*iCase/length(CaseTimeList))),repmat('□',1,10-round(10*iCase/length(CaseTimeList))),s2])
try
    CaseTime = irf_time(CaseTimeList{iCase}, 'utc>epoch');
    tint = irf_time([CaseTime-0.1,CaseTime+0.1],'epoch>epochTT');
    
    % load B
    c_eval(['B?_ts=mms.get_data(''B_gsm_brst'',tint,?);'],ic);
    % c_eval(['Bt?_ts=B?_ts.abs;'],ic); 
    c_eval(['B?=irf.ts2mat(B?_ts);'],ic);
    %  c_eval(['B?_gsm=irf_gse2gsm(B?,-1);'],ic);


    % load E
    c_eval(['E?_ts=mms.get_data(''E_gse_edp_brst_l2'',tint,?);'],ic);
    %%%%%c_eval(['E?_ts=mms.get_data(''E_gse_edp_fast_l2'',tint,?);'],ic);
    c_eval(['Et?_ts=E?_ts.abs;'],ic); 
    c_eval(['E?_gsm=irf_gse2gsm(E?_ts);'],ic);
    c_eval(['E?=irf.ts2mat(E?_gsm);'],ic);
    % c_eval('Et? = irf_abs(E?);',ic);
    c_eval(['Efac?=irf_convert_fac(E?,B?,[1,0,0]);'],ic);
    c_eval('Efac?_ts = irf.ts_vec_xyz(E?_ts.time,Efac?(:,2:4));',ic);

    % load N
    Ne1_ts=mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint);
    Ne1=irf.ts2mat(Ne1_ts);

    % load R
    tintlong = tint+[-60 60];
    R  = mms.get_data('R_gsm',tintlong);
    c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gsmR?(:,1:3));',ic);

catch 
    writematrix([irf_time(CaseTime,'epoch>utc'),'的数据导入出现问题'],[OutputDir,'DRAFTerrorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
    continue
end
%% wave analysis
try
Pfrange = [2 1000];
[Vp,kmag, waveL, waveE, waveThe, Fre, W1, W2, W3, W4,kx,ky,kz] = WaveAna_4SC_fast('Efac?_ts.z','Rxyz?','B?_ts',tint,'numf',400,...
    'wwidth',1,'frange',Pfrange,'cav',1);  

P=waveThe(:,2:end);
P(P>90)=180-P(P>90);
waveThe(:,2:end)=P;
specE=struct('t',waveE(:,1));
specE.f=Fre;
specE.p=waveE(:,2:end);
specE.f_label='';
specE.p_label={'log_{10} B_{z}^2 (nT^2 Hz^-1)'};


specV=struct('t',Vp(:,1));
specV.f=Fre;
specV.p=Vp(:,2:end);
specV.f_label='';
specV.p_label={'V_p (km/s)'};

specThe=struct('t',waveThe(:,1));
specThe.f=Fre;
specThe.p=waveThe(:,2:end);
specThe.f_label='';
specThe.p_label={'\Theta (deg)'};


speckmag=struct('t',waveThe(:,1));
speckmag.f=Fre;
speckmag.p=kmag(:,1:end).*1000;  % 换单位
speckmag.f_label='';
speckmag.p_label={'k_{} (km^{-1})'};
%% Init figure
        close all;
        h=irf_plot(4,'newfigure');
        xSize=800; ySize=800;
        set(0,'DefaultAxesFontSize',8);
        set(0,'DefaultLineLineWidth',2);
        set(gcf,'Position',[10 10 xSize ySize]);
        i_subplot=4;
        
        for ii=1:i_subplot-1
            pospanel(ii,:)=get(h(ii),'pos');
        end
        for ii=1:i_subplot-2
            pospanel(ii,2)=pospanel(ii,2)+(i_subplot-1-ii)*0.005;
            set(h(ii),'pos',pospanel(ii,:));
        end              

        h(1)=irf_panel('E');
        irf_spectrogram(h(1),specE,'log');
        hold(h(1),'on');
        hold(h(1),'off');
        % caxis(h(1),[-4 -1]);
        set(h(1),'yscale','log');
        set(h(1),'Ylim',Pfrange);
        set(h(1),'ytick', [ 0.02 0.1 1 4 ]);
        grid(h(1),'off');
        ylabel(h(1),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
        irf_legend(h(1),'(a)',[0.02 0.95],'color','k','fontsize',10)
        set(h(1),'FontSize',10); 
        colormap(h(1), jet);
        
        h(2)=irf_panel('Vp');
        irf_spectrogram(h(2),specV,'lin');
        hold(h(2),'on');
        hold(h(2),'off');
        caxis(h(2),[0 300]);
        set(h(2),'yscale','log');
        set(h(2),'Ylim',Pfrange)
        grid(h(2),'off');
        ylabel(h(2),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
        irf_legend(h(2),'(b)',[0.02 0.95],'color','k','fontsize',10)
        % set(h(2),'ytick', [ 0.02 0.1 1 4 ]);
        set(h(2),'FontSize',10);        
        colormap(h(2), jet);
        
        h(3)=irf_panel('V');
        [h(3), hcb3]=irf_spectrogram(h(3),specThe,'lin');
        hold(h(3),'on');
        hold(h(3),'off');
        % caxis(h(3),[0 90]);
        set(hcb3,'xtick', [0 45 90]);
        set(h(3),'yscale','log');
        set(h(3),'Ylim',Pfrange)
        grid(h(3),'off');
        ylabel(h(3),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
        irf_legend(h(3),'(c)',[0.02 0.95],'color','k','fontsize',10) 
        % set(h(3),'ytick', [ 0.02 0.1 1 4]);  
        set(h(3),'FontSize',10);        
        colormap(h(3), jet);      
        
        
        h(4)=irf_panel('kmag');
        irf_spectrogram(h(4),speckmag,'lin');
        hold(h(4),'on');
        hold(h(4),'off');
        % caxis(h(4),[0 0.05]);
        set(h(4),'yscale','log');
        set(h(4),'Ylim',Pfrange)
        grid(h(4),'off');
        ylabel(h(4),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
        % set(h(4),'ytick', [ 0.02 0.1 1 4 ]);  
        irf_legend(h(4),'(d)',[0.02 0.95],'color','k','fontsize',10) 
        set(h(4),'FontSize',10);  
        colormap(h(4), jet);
%%
irf_zoom(tint,'x',h(1:i_subplot));
irf_plot_axis_align(i_subplot)
% for ii = 1:n  
%     tempTime = flagTime(1); %2s之内相近的符合判据的时间只画一条线
%     irf_pl_mark(h(ii),tempTime,'k')
%     id_flagTime = [1];
%     for jj = 2:length(flagTime)
%         if flagTime(jj) - tempTime >= 2
%             tempTime = flagTime(jj);
%             irf_pl_mark(h(ii),tempTime,'k');
%             id_flagTime(end+1) = jj;
%         end
%     end
% end
%%  出图保存部分
set(gca,"XTickLabelRotation",0)
% set(gcf,'render','painters');
% set(gcf,'paperpositionmode','auto')
% figname = [OutputDir,'OverviewFig\',Name(2:end-2)];   
mkdir([OutputDir,'DRAFT\']);
tempName = strrep(CaseTimeList{iCase},':','');tempName = strrep(tempName,'.','');
figname = [OutputDir,'DRAFT\', tempName];  
print(gcf, '-dpng', [figname '.png']);    
%     pause(1)
close all

save([OutputDir,'DRAFT\', tempName, '.mat'], 'specE','specV')
writematrix([irf_time(CaseTime,'epoch>utc'),'的最大速度为',num2str(max(Vp(:,2:end),[],'all'))],...
    [OutputDir,'DRAFTCaseList.txt'],'WriteMode','append','Encoding','UTF-8')
catch
    writematrix([irf_time(CaseTime,'epoch>utc'),'的画图出现问题'],[OutputDir,'DRAFTerrorlog.txt'],...
    'WriteMode','append','Encoding','UTF-8')
    continue
end
end
