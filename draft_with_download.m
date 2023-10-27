%% 初始化
clc;
clear;close all
global ParentDir 
ParentDir = '/Volumes/FWD-T7Disk/MMS/'; 
mms.db_init('local_file_db','/Volumes/FWD-T7Disk/MMS/')
TempDir = '/Volumes/FWD-T7Disk/MMS/temp/';mkdir(TempDir);
TT = '2019-08-05T16:24:32.00Z/2019-08-05T16:24:37.00Z';
ic=1:4;
%% 下载数据
Tint = irf.tint(TT);
tint=Tint;
Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce,scpot');
filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy');
% filenames = [filenames1,filenames2,filenames3,filenames4];
filenames = [filenames1,filenames4];

expr1 = '_\d+\_v';
NameTags = regexp(filenames,expr1,'match');
NameTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NameTags)),'UniformOutput',false);

TTlist = strjoin(regexp(TT,'\d+','match'),'');
i = 1;flag = 0;%若flag=0，说明整段时间都在第i-1个Tag里
while str2double(TTlist(17:30)) > NameTags{i}  % 如果时间段刚好仅在某天的最后一个文件里会出bug，可以把时间往前调1ms
    if str2double(TTlist(1:14)) < NameTags{i}
        flag=1; break  %若flag=1，说明时间段的开始在第i-1个Tag里，结束在第i个里
    else
        i=i+1;
    end
    if i > length(NameTags)
        break
    end
end

% '/' for MacOs, '\' for Windows
if flag == 0
    tempTag = num2str(NameTags{i-1});
    filenames = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag)));
    desmoms =  [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag(1:4),'/',tempTag(5:6),'/',...
            tempTag(7:8),'/',filenames{cellfun(@(x)(~isempty(x)),strfind(filenames,'des-moms'))}];
    desmoms1 = desmoms; desmoms2 = desmoms1;
else
    if i == 1
        errordlg('时间起始处无brst数据，请检查时间范围,或使用Overview_srvydownload程序')
    end
    tempTag1 = num2str(NameTags{i-1});
    tempTag2 = num2str(NameTags{i});
    filenames1 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag1)));
    filenames2 = filenames(cellfun(@(x)(~isempty(x)),strfind(filenames,tempTag2)));
    filenames = [filenames1,filenames2];
    desmoms1 = [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag1(1:4),'/',tempTag1(5:6),'/',...
            tempTag1(7:8),'/',filenames1{cellfun(@(x)(~isempty(x)),strfind(filenames1,'des-moms'))}];
    desmoms2 = [ParentDir,'mms',num2str(ic),'/fpi/brst/l2/des-moms/',tempTag2(1:4),'/',tempTag2(5:6),'/',...
            tempTag2(7:8),'/',filenames2{cellfun(@(x)(~isempty(x)),strfind(filenames2,'des-moms'))}];
end
SDCFilesDownload(filenames,TempDir)
SDCFilesDownload(filenames_srvy,TempDir)
SDCDataMove(TempDir,ParentDir)
%% 读取数据
c_eval('Bxyz?=mms.get_data(''B_gse_fgm_brst_l2'',Tint,?);',ic);
c_eval('Exyz?=mms.get_data(''E_gse_edp_brst_l2'',Tint,?);',ic);
Tintlong = Tint+[-60 60];
R  = mms.get_data('R_gse',Tintlong);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

%% 特性分析
frange = [100 1000];
[Vp,kmag, waveL, waveE, waveThe, Fre, W1, W2, W3, W4,kx,ky,kz] = WaveAna_4SC('Exyz?.abs','Rxyz?','Bxyz?',Tint,'numf',400,...
    'wwidth',2,'frange',frange,'sn',10,'cav',8);     

%% 检验
[Rcor, Ratio_L_sc, Ratio_M] = Checkwave_4SC('Exyz?.abs','Rxyz?', waveL, Fre, W1, W2, W3, W4);

%% 角度修正

P=waveThe(:,2:end);
P(P>90)=180-P(P>90);
waveThe(:,2:end)=P;
%% 截断无波区域  
% waveD=waveE(:,2:end);
% waveD(waveD<0.05)=nan;
% waveD(:,280:end)=nan;
% waveD(waveD>0)=1;
% Vp(:,2:end)=Vp(:,2:end).*waveD;
% waveThe(:,2:end)=waveThe(:,2:end).*waveD;
% waveE(:,2:end)=waveE(:,2:end).*waveD;
% kmag=kmag.*waveD;


%% 结构体
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

specL=struct('t',Ratio_L_sc(:,1));
specL.f=Fre;
specL.p=Ratio_L_sc(:,2:end);
specL.f_label='';
specL.p_label={'\lambda/2L'};


specR=struct('t',Rcor(:,1));
specR.f=Fre;
specR.p=Rcor(:,2:end);
specR.f_label='';
specR.p_label={'cc'};

specM=struct('t',Ratio_M(:,1));
specM.f=Fre;
specM.p=Ratio_M(:,2:end);
specM.f_label='';
specM.p_label={'ma'};



%% Plot figure all quantities

close all;
h=irf_plot(7,'newfigure');
xSize=800; ySize=800;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth',2);
set(gcf,'Position',[10 10 xSize ySize]);
i_subplot=7;

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
caxis(h(1),[-4 -1]);
set(h(1),'yscale','log');
set(h(1),'Ylim',frange);
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
% caxis(h(2),[0 300]);
set(h(2),'yscale','log');
set(h(2),'Ylim',frange)
grid(h(2),'off');
ylabel(h(2),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
irf_legend(h(2),'(b)',[0.02 0.95],'color','k','fontsize',10)
set(h(2),'ytick', [ 0.02 0.1 1 4 ]);
set(h(2),'FontSize',10);        
colormap(h(2), jet);

h(3)=irf_panel('V');
[h(3), hcb3]=irf_spectrogram(h(3),specThe,'lin');
hold(h(3),'on');
hold(h(3),'off');
caxis(h(3),[0 90]);
set(hcb3,'xtick', [0 45 90]);
set(h(3),'yscale','log');
set(h(3),'Ylim',frange)
grid(h(3),'off');
ylabel(h(3),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
irf_legend(h(3),'(c)',[0.02 0.95],'color','k','fontsize',10) 
set(h(3),'ytick', [ 0.02 0.1 1 4]);  
set(h(3),'FontSize',10);        
colormap(h(3), jet);      


h(4)=irf_panel('kmag');
irf_spectrogram(h(4),speckmag,'lin');
hold(h(4),'on');
hold(h(4),'off');
caxis(h(4),[0 0.05]);
set(h(4),'yscale','log');
set(h(4),'Ylim',frange)
grid(h(4),'off');
ylabel(h(4),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
set(h(4),'ytick', [ 0.02 0.1 1 4 ]);  
irf_legend(h(4),'(d)',[0.02 0.95],'color','k','fontsize',10) 
set(h(4),'FontSize',10);  
colormap(h(4), jet);

h(5)=irf_panel('L');
irf_spectrogram(h(5),specL,'lin');
hold(h(5),'on');        
hold(h(5),'off');
caxis(h(5),[0 2]);
set(h(5),'yscale','log');
set(h(5),'Ylim',frange)
grid(h(5),'off');
ylabel(h(5),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
set(h(5),'ytick', [ 0.02 0.1 1 4 ]);
irf_legend(h(5),'(f)',[0.02 0.95],'color','k','fontsize',10) 
set(h(5),'FontSize',10);             
colormap(h(5), jet);


h(6)=irf_panel('R');
irf_spectrogram(h(6),specR,'lin');
hold(h(6),'on');        
hold(h(6),'off');
caxis(h(6),[0 1]);
set(h(6),'yscale','log');
set(h(6),'Ylim',frange)
grid(h(6),'off');
ylabel(h(6),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
set(h(6),'ytick', [ 0.02 0.1 1 4 ]);
irf_legend(h(6),'(g)',[0.02 0.95],'color','k','fontsize',10) 
set(h(6),'FontSize',10);  
colormap(h(6), jet);

h(7)=irf_panel('M');
irf_spectrogram(h(7),specM,'lin');
hold(h(7),'on');      
hold(h(7),'off');
caxis(h(7),[0 1]);
set(h(7),'yscale','log');
set(h(7),'Ylim',frange)
grid(h(7),'off');
ylabel(h(7),{'f (Hz)'},'fontsize',9,'Interpreter','tex');
set(h(7),'ytick', [0.02 0.1 1 4 ]); 
irf_legend(h(7),'(h)',[0.02 0.95],'color','k','fontsize',10) 
set(h(7),'FontSize',10); 
colormap(h(7), jet);


 irf_plot_axis_align(h);          
 % tint = irf.tint('2018-04-16T10:20:00.00Z/2018-04-16T11:00:00.00Z');  % T3
 irf_zoom(h,'x',tint);
 
 
 
 
 %% 
 
 
 
set(gcf,'paperpositionmode','auto')

set(gcf,'render','zbuffer');  
figname=[''];
% print(gcf, '-dpdf','-r400',figname);

        
        