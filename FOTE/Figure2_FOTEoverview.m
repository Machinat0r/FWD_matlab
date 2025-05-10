clear;close all
clc
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);

ic = 1:4;
% tint =irf.tint('2021-08-15T03:34:42Z/2021-08-15T03:44:15Z');

Tsta = '2017-07-12T11:54:33.600Z';
Tend = '2017-07-12T11:54:35.400Z';
tint = irf.tint(Tsta,Tend);

TT = '2017-07-12T11:54:00.000Z/2017-07-12T11:55:00.000Z';
% % % TT = '2017-07-06T17:46:30.000Z/2017-07-06T17:48:00.000Z';
% TT = '2018-07-06T12:36:30.000Z/2018-07-06T12:38:00.000Z';
% TT = '2017-05-22T08:23:30.000Z/2017-05-22T08:24:30.000Z';
% TT = '2017-05-25T04:40:00.000Z/2017-05-25T04:41:00.000Z';

tint2=irf.tint(TT);

Datelist = regexp(TT,'\d+-\d+-\d+','match');
Datelist{2} = datestr(datenum(Datelist{2},'yyyy-mm-dd')+1,'yyyy-mm-dd');
Date = [Datelist{1},'/',Datelist{2}];

filenames1 = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
filenames4 = SDCFilenames(Date,ic,'inst','edp','drm','brst','dpt','dce');
filenames = [filenames1, filenames2,filenames4];
% % % 
[filenames,desmoms1,desmoms2] = findFilenames(TT,filenames,'brst',ic);

SDCFilesDownload_NAS(filenames,TempDir, 'Threads', 48, 'CheckSize', 0)
SDCDataMove(TempDir,ParentDir)
mms.db_init('local_file_db',ParentDir)

thresold = 0.1;
hh = 20;
for ic = 1:4
    
    %% load Fields data
    c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
    % c_eval('Bxyz?=mms.get_data(''B_gse_srvy_l2'',tint,?);',ic);
    c_eval('B?=irf.ts2mat(Bxyz?);',ic);
    c_eval('Bt?=Bxyz?.abs;',ic);

    %% smooth
    c_eval('BS?(:,2)=smooth(B?(:,2),hh);',ic);
    c_eval('BS?(:,3)=smooth(B?(:,3),hh);',ic);
    c_eval('BS?(:,4)=smooth(B?(:,4),hh);',ic);
    c_eval('BS?(:,1)=B?(:,1);',ic);
%     c_eval('BS?(:,1)=smooth(B?(:,1),hh);',ic);
    c_eval('B?=BS?;',ic);
    c_eval('B?=irf_abs(B?);',ic);
    
    c_eval('Ne? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',tint);',ic);
    c_eval('Ni? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_brst'',tint);',ic);
    c_eval('Vi?= mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkv_gse_brst'',tint);',ic);
    c_eval('Ve?= mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkv_gse_brst'',tint);',ic);
    
    c_eval('Exyz?=mms.get_data(''E_gse_edp_fast_l2'',tint,?);',ic);
    c_eval('E?=irf.ts2mat(Exyz?);',ic);
    c_eval('Rgse?=mms.get_data(''R_gse'',tint2,?);',ic);
    c_eval('R?=irf.ts2mat(Rgse?);',ic);
end

%% FOTE method
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
tint=[iso2epoch(Tsta) iso2epoch(Tend)];
for ic=1:4
  c_eval(['R?=irf_resamp(R?,B?);'],ic);
  c_eval(['B?=irf_tlim(B?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
  c_eval(['BmagC?=irf_abs(B?);'],ic);
end
gradB=c_4_grad('R?','B?','grad');
d_dot_B=c_4_grad('R?','B?','div');
d_cros_B=c_4_grad('R?','B?','curl');


%error of curolmeter
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

temp=irf_abs(d_cros_B);
d_cros_B_mag=temp(:,[1 5]);
err_curlmeter=irf_multiply(1,d_dot_B,1,d_cros_B_mag,-1);
err_curlmeter(:,2)=abs(err_curlmeter(:,2))*100;


%Null type identification
for ii=1:length(d_dot_B(:,1))
    mksizSim(ii)=4;
    
    deltB_null=reshape(gradB(ii,2:end),3,3);
    [V,D] = eig(deltB_null);
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='w';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            type(ii)='X'; clr(ii)='k'; faceclr(ii)='w';
        end
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='r';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            type(ii)='o'; clr(ii)='k'; faceclr(ii)='w';
        end
    end
    %=========================================================
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='w';
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if min(abs([D(1,1) D(2,2) D(3,3)]))/max(abs([D(1,1) D(2,2) D(3,3)]))<thresold
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
                typeSim(ii)='X'; clrSim(ii)='b'; faceclrSim(ii)='w'; mksizSim(ii)=7;
            end
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='X'; clrSim(ii)='r'; faceclrSim(ii)='w'; mksizSim(ii)=7;
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            typeSim(ii)='X'; clrSim(ii)='k'; faceclrSim(ii)='w'; mksizSim(ii)=7;
        end
        %------------------------------------------------------------------
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='r';
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))/max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) < (thresold*2)
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
                typeSim(ii)='o'; clrSim(ii)='b'; faceclrSim(ii)='w';
            end
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='o'; clrSim(ii)='r'; faceclrSim(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            typeSim(ii)='o'; clrSim(ii)='k'; faceclrSim(ii)='w';
        end
        %------------------------------------------------------------------
    end
    %=========================================================
    
    eigVal_err(ii,2)=abs(real(D(1,1)+D(2,2)+D(3,3)))/max(abs([real(D(1,1)), real(D(2,2)), real(D(3,3))])) * 100;
    sumeigVal(ii,2)=D(1,1)+D(2,2)+D(3,3);
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;
end
eigVal_err(:,1)=d_dot_B(:,1);
sumeigVal(:,1)=d_dot_B(:,1);
eigVal_err_v2(:,1)=d_dot_B(:,1);


%find null position
for ii=1:length(B1(:,1))
    dBeach=reshape(gradB(ii,2:end),3,3);
    dR1(ii,2:4)=B1(ii,2:4)*inv(dBeach');
    dR2(ii,2:4)=B2(ii,2:4)*inv(dBeach');
    dR3(ii,2:4)=B3(ii,2:4)*inv(dBeach');
    dR4(ii,2:4)=B4(ii,2:4)*inv(dBeach');
end
dR1(:,1)=B1(:,1);
dR2(:,1)=B1(:,1);
dR3(:,1)=B1(:,1);
dR4(:,1)=B1(:,1);

dRmag1=irf_abs(dR1);
dRmag2=irf_abs(dR2);
dRmag3=irf_abs(dR3);
dRmag4=irf_abs(dR4);

dRmin(:,2)=min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1)=dRmag1(:,1);

pause(1)

%% figure1
h1=irf_plot(4,'newfigure');
lnwid = 0.5; fonts = 8;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth',lnwid);
xSize=600; ySize=900; panal_ad = 0.95;
set(gcf,'Position',[100 100 xSize ySize]);


%% magnetic field Bz
% h1(i_plot)=irf_panel('B');
i_plot = 1;
hold(h1(i_plot),'on');
c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,5)],''Linewidth'',lnwid,''color'',''k'')',ic);
c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,2)],''Linewidth'',lnwid,''color'',''b'')',ic);
c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,3)],''Linewidth'',lnwid,''color'',''g'')',ic); %[0,0.3906,0]
c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,4)],''Linewidth'',lnwid,''color'',''r'')',ic);

c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,2)*0],''k--'', ''Linewidth'',0.5);',ic);
hold(h1(i_plot),'off');
ylabel(h1(i_plot),{'B [nT]'},'Interpreter','tex','fontsize',fonts);
set(h1(i_plot),'ylim',[-5 10],'fontsize',fonts);
irf_legend(h1(i_plot),'(c)',[0.99 0.98],'color','k','fontsize',fonts)
grid(h1(i_plot),'off');box(h1(i_plot),'on');
Pos_c = get(h1(i_plot),'position'); Pos_s = get(h1(1),'position');
set(h1(i_plot),'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);

%% ion velocity Vi
% h1(i_plot)=irf_panel('B');
i_plot = i_plot+1;
hold(h1(i_plot),'on');
% c_eval('irf_plot(h1(i_plot),Vi?.abs,''Linewidth'',lnwid,''color'',''k'')',ic);
c_eval('irf_plot(h1(i_plot),Vi?.x,''Linewidth'',lnwid,''color'',''b'')',ic);
c_eval('irf_plot(h1(i_plot),Vi?.y,''Linewidth'',lnwid,''color'',''g'')',ic);
c_eval('irf_plot(h1(i_plot),Vi?.z,''Linewidth'',lnwid,''color'',''r'')',ic);

c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,2)*0+80],''k--'', ''Linewidth'',0.5);',ic);
hold(h1(i_plot),'off');
ylabel(h1(i_plot),{'V_{ion} [km/s]'},'Interpreter','tex','fontsize',fonts);
set(h1(i_plot),'ylim',[-200 350],'ytick',[-200:100:300],'fontsize',fonts);
irf_legend(h1(i_plot),'(d)',[0.99 0.98],'color','k','fontsize',fonts)
grid(h1(i_plot),'off');box(h1(i_plot),'on');
Pos_c = get(h1(i_plot),'position'); Pos_s = get(h1(1),'position');
set(h1(i_plot),'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);


% %% electron velocity Ve
% % h1(i_plot)=irf_panel('B');
% i_plot = i_plot+1;
% hold(h1(i_plot),'on');
% c_eval('irf_plot(h1(i_plot),Ve?.x,''Linewidth'',lnwid,''color'',''b'')',ic);
% c_eval('irf_plot(h1(i_plot),Ve?.y,''Linewidth'',lnwid,''color'',''g'')',ic);
% c_eval('irf_plot(h1(i_plot),Ve?.z,''Linewidth'',lnwid,''color'',''r'')',ic);
% 
% c_eval('irf_plot(h1(i_plot),[B?(:,1) B?(:,2)*0],''k--'', ''Linewidth'',0.5);',ic);
% 
% % c_eval('irf_plot(h1(i_plot),[E?(:,1) E?(:,2)*0],''k--'', ''Linewidth'',0.5);',ic);
% hold(h1(i_plot),'off');
% ylabel(h1(i_plot),{'V_{electron} [km/s]'},'Interpreter','tex','fontsize',fonts);
% % set(h1(i_plot),'ylim',[-1150 1250],'ytick',[-1000:500:1000],'fontsize',fonts);
% irf_legend(h1(i_plot),'(e)',[0.99 0.98],'color','k','fontsize',fonts)
% grid(h1(i_plot),'off');box(h1(i_plot),'on');
% Pos_c = get(h1(i_plot),'position'); Pos_s = get(h1(1),'position');
% set(h1(i_plot),'position',[Pos_c(1) Pos_c(2) Pos_s(3) Pos_c(4)]);

%% distance plot
i_plot = i_plot+1;
hold(h1(i_plot),'on');
% irf_plot(h1(i_plot),dRmag1(:,[1 5]), 'color','k', 'Linewidth',lnwid); hold on;
% irf_plot(h1(i_plot),dRmag2(:,[1 5]), 'color','r', 'Linewidth',lnwid); hold on;
% irf_plot(h1(i_plot),dRmag3(:,[1 5]), 'color','g', 'Linewidth',lnwid); hold on;
% irf_plot(h1(i_plot),dRmag4(:,[1 5]), 'color','b', 'Linewidth',lnwid); hold on;
dRmin(dRmin(:,2)>=2000,2)=3000;
irf_plot(h1(i_plot),dRmin(:,[1 2]), 'color','b', 'Linewidth',lnwid); hold on;
for ii=1:2:length(dRmin(:,1))
    irf_plot(h1(i_plot), dRmin(ii,:), [typeSim(ii) clrSim(ii)], 'MarkerSize',mksizSim(ii),'MarkerFaceColor',faceclrSim(ii), 'Linewidth',lnwid); hold on;
end
set(h1(i_plot),'Ylim',[0 2e3]);
% set(h1(i_plot), 'Ylim',[0 35],'Ytick',[0:10:30]);
ylabel(h1(i_plot),{'|r|','(km)'},'Interpreter','tex');
irf_legend(h1(i_plot),'(h)',[0.99 0.94],'color','k','fontsize',12)
grid(h1(i_plot),'off');

%% error plot
i_plot = i_plot+1;
hold(h1(i_plot),'on');
a1=irf_plot(h1(i_plot),[err_4C(:,1) err_4C(:,2)],'k','Linewidth',lnwid);
a2=irf_plot(h1(i_plot), eigVal_err_v2, 'c','Linewidth',lnwid); 

% fill(h1(i_plot),[a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
fill(h1(i_plot),[a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
fill(h1(i_plot), [a1.XData,fliplr(a1.XData)],[zeros(1,length(a1.YData)),fliplr(a1.YData)],[0.7 0.7 0.7]);


% irf_plot(h1(i_plot), [err_4C(:,1) err_4C(:,2)*0+40], 'k--','Linewidth',lnwid); 
irf_plot(h1(i_plot),[eigVal_err(:,1) eigVal_err(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;
ylim([0, 100]);

grid off;
% set(h1(i_plot),'Ylim',[0 160]);
%ylabel('\fontcolor{k}|\nabla\cdot\bf{B}|/|\nabla\times\bf{B}|,\fontcolor{b}|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}| [%]');
ylabel('\eta,\xi [%]','Interpreter','tex');
irf_legend(h1(i_plot),'(i)',[0.99 0.94],'color','k','fontsize',12)
irf_legend(h1(i_plot),{'|\nabla\cdot\bf{B}|/|\nabla\times\bf{B}|'},[0.75 0.9],'color','k','fontsize',10)
irf_legend(h1(i_plot),{'|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}| [%]'},[0.75 0.75],'color','c','fontsize',10)

irf_plot_ylabels_align(h1)
% title(h(1),'the overview of accident','fontname','times new roman','fontweight','normal')

irf_zoom(h1,'x',tint)
% irf_zoom(h1,'x',irf.tint('2017-01-27T12:05:42.20Z/2017-01-27T12:05:44.20Z'))

