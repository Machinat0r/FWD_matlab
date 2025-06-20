clear
clc

global ParentDir 
ParentDir = 'D:/MMS/'; 
DownloadDir = 'C:/MMS/';
TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
%%
ic=2;
% Tsta='2018-08-22T15:34:32.60Z';   
% Tend='2018-08-22T15:34:33.80Z';
% Tint=irf.tint(Tsta,Tend);
% 
% Tsta1='2018-08-22T15:34:34.30Z';   
% Tend1='2018-08-22T15:34:35.20Z';
% Tint1=irf.tint(Tsta1,Tend1);
Tint=irf.tint('2020-08-03T01:45:27.500Z/2020-08-03T01:45:28.800Z');
% Tint=irf.tint('2017-07-26T07:21:13.00Z/2017-07-26T07:33:42.00Z');
% Tint=irf.tint('2017-07-26T07:21:13.00Z/2017-07-26T07:28:00.00Z');
% Tint=irf.tint('2017-07-26T07:30:00.00Z/2017-07-26T07:32:00.00Z');

%%
c_eval('Bxyz=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tint);',ic);

c_eval(['Bts=Bxyz.abs;'],ic); 
c_eval(['Bt=irf.ts2mat(Bts);'],ic);

% c_eval('Exyz_gse=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);
% Exyz=irf_gse2gsm(Exyz_gse);
c_eval('Bscm_gse=mms.get_data(''B_gse_scm_brst_l2'',Tint,?);',ic);

c_eval(['Btscm=Bscm_gse.abs;'],ic); 
c_eval(['Btscm1=irf.ts2mat(Btscm);'],ic);

% Bscm=irf_gse2gsm(Bscm_gse);

Bscm=irf_gse2gsm(Bscm_gse);


%%

c_eval('ne = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_brst'',Tint);',ic);   %数据版本不一样
% magBh=Bh.abs;
magB=Bxyz.abs;
magBscm=Bscm.abs;

Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);


meanB=mean(magB.data);
meanne=mean(ne.data);
%% Tint
Bt_for_psd=[Bts.time.epochUnix,double(Bt(:,2))];
[psd_bt,f_t] = irf_psd(Bt_for_psd,3000,1/(Bts.time(2)-Bts.time(1)),'mean');

Bt_for_psd1=[magBscm.time.epochUnix,double(magBscm.data(:,1))];
[psd_bt1,f_t1] = irf_psd(Bt_for_psd1,2048,1/(magBscm.time(2)-magBscm.time(1)),'none');

 
% Bx_for_psd=[Bscm.time.epochUnix,double(Bscm.data(:,1))];
% [psd_bx,f_bx] = irf_psd(Bx_for_psd,2048,1/(Bscm.time(2)-Bscm.time(1)),'mean');
% 
% Bx_for_psd1=[Bxyz.time.epochUnix,double(Bxyz.data(:,1))];
% [psd_bx1,f_bx1] = irf_psd(Bx_for_psd1,4000,1/(Bxyz.time(2)-Bxyz.time(1)),'mean');
% 
% By_for_psd=[Bscm.time.epochUnix,double(Bscm.data(:,2))];
% [psd_by,f_by] = irf_psd(By_for_psd,2048,1/(Bscm.time(2)-Bscm.time(1)));
% 
% By_for_psd1=[Bxyz.time.epochUnix,double(Bxyz.data(:,2))];
% [psd_by1,f_by1] = irf_psd(By_for_psd1,2048,1/(Bxyz.time(2)-Bxyz.time(1)));
% 
% Bz_for_psd=[Bscm.time.epochUnix,double(Bscm.data(:,3))];
% [psd_bz,f_bz] = irf_psd(Bz_for_psd,2048,1/(Bscm.time(2)-Bscm.time(1)));
% 
% Bz_for_psd1=[Bxyz.time.epochUnix,double(Bxyz.data(:,3))];
% [psd_bz1,f_bz1] = irf_psd(Bz_for_psd1,2048,1/(Bxyz.time(2)-Bxyz.time(1)));

%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;

B_SI=meanB*1e-9;
Wpe = sqrt(meanne*1e6*e^2/Me/epso);
Wce = e*B_SI/Me;
Wci = e*B_SI/Mp;
Wpp = sqrt(meanne*1e6*e^2/Mp/epso);
Fce = Wce/2/pi;
Fci = Wci/2/pi;
Fce01=Fce*0.1;
Fce05=Fce*0.5;
Fpe = Wpe/2/pi;
Fcp = Fce/Mp_Me;
Fpp = Wpp/2/pi;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);

%%
fn=figure;
set(fn,'Position',[100 100 400 400])
    h(1)=axes('position',[0.30 0.30 0.4 0.60]); 
%     h(2)=axes('position',[0.41 0.40 0.15 0.30]);
    
    
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',0.75); 
    

%% 
% f_t(1,:)=[];
% psd_bt(1,:)=[];
x1=log10(f_t);
x1(1,:)=[];
y1=log10(psd_bt);
y1(1,:)=[];
%%

x2=log10(f_t1);
x2(1,:)=[];
y2=log10(psd_bt1);
y2(1,:)=[];

% f1=x1(1:3,:);
% p1=y1(1:3,:);
% f2=x1(4:151,:);
% p2=y1(4:151,:);
% f3=x1(152:end,:);
% p3=y1(152:end,:);

f1=x1(1:4,:);
p1=y1(1:4,:);
f2=x1(4:151,:);
p2=y1(4:151,:);
f3=x1(152:end,:);
p3=y1(152:end,:);

f4=x2(1:300,:);
p4=y2(1:300,:);
k1=polyfit(f1,p1,1);
k2=polyfit(f2,p2,1);
k3=polyfit(f3,p3,1);
k4=polyfit(f4,p4,1);
yy1=polyval(k1,f1);
yy2=polyval(k2,f2);
yy3=polyval(k3,f3);
yy4=polyval(k4,f4);
f4=10.^f4;yy4=10.^yy4;
f4=f4(8:end,:);yy4=yy4(8:end,:);


% plot(h(1),f_t,psd_bt,'color','[1.00,0.73,0.08]','linewidth',2);hold on;
% plot(h(1),f_t1(8:end,:),psd_bt1(8:end,:),'color','[1.00,0.73,0.08]','linewidth',2);hold on;

%%
plot(h(1),f_t,psd_bt,'color','k','linewidth',2);hold on;
plot(h(1),f_t1(8:end,:),psd_bt1(8:end,:),'color','y','linewidth',2);hold on;
% plot(h(1),10.^f1,10.^yy1,'color','g','linewidth',0.5);hold on;
% plot(h(1),10.^f2,10.^yy2,'color','g','linewidth',0.5);hold on;
% plot(h(1),10.^f3,10.^yy3,'color','r','linewidth',0.5);hold on;
% plot(h(1),f4,yy4,'color','r','linestyle','--','linewidth',0.5);hold on;
%%

% plot(h(1),[Fce Fce],yrange,'color','r','linestyle','--');
% plot(h(1),[Fce05 Fce05],yrange,'color','k','linestyle','--');
% plot(h(1),[Fce01 Fce01],yrange,'color','b','linestyle','--');

plot(h(1),[Flh Flh],[10^(-10) 10^(4)],'color','k','linestyle','--');
plot(h(1),[Fci Fci],[10^(-10) 10^(4)],'color','k','linestyle','--');
hold (h(1),'off');
set(h(1),'yscale','log');
set(h(1),'xscale','log');

xlabel(h(1),'frequency (Hz)','fontsize',10,'fontname','times new roman');
ylabel(h(1),{'B_{t PSD} [(nT)^{2}/Hz]'},'Interpreter','tex','fontsize',10);
% irf_zoom(gca,'x',[10 3.2e4]);
set(h(1),'ylim',[10^(-9.5) 10^(3)],'ytick',[10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4)]);
set(h(1),'xlim',[10^(-2) 10^4],'xtick',[10^(-2) 10^(-1) 10^(0) 10^(1) 10^2  10^3 10^4]);

% irf_legend(gca,'0.1f_{ce}',[0.05 0.05],'color','b','fontsize',8)
% irf_legend(gca,'f_{ce}',[0.05 0.15],'color','r','fontsize',8)
% irf_legend(gca,'0.5f_{ce}',[0.05 0.25],'color','k','fontsize',8)
% irf_legend(gca,'f_{lh}',[0.05 0.35],'color','g','fontsize',8)
% irf_zoom(gca,'y',yrange1);
%%

% plot(h(2),f_bx,psd_bx,'color','k');
% hold (h(2),'on');
% plot(h(2),f_bx1,psd_bx1,'color','b');
% 
% 
% % plot(h(1),[Fce Fce],yrange,'color','r','linestyle','--');
% % plot(h(1),[Fce05 Fce05],yrange,'color','k','linestyle','--');
% % plot(h(1),[Fce01 Fce01],yrange,'color','b','linestyle','--');
% plot(h(2),[8 8],[10^(-10) 10^(4)],'color','g','linestyle','--');
% plot(h(2),[0.1 0.1],[10^(-10) 10^(4)],'color','g','linestyle','--');
% hold (h(2),'off');
% set(h(2),'yscale','log');
% set(h(2),'xscale','log');
% xlabel(h(2),'frequency (Hz)','fontsize',10,'fontname','times new roman');
% ylabel(h(2),{'B_{x PSD} [(nT)^{2}/Hz]'},'Interpreter','tex','fontsize',10);
% % irf_zoom(gca,'x',[10 3.2e4]);
% set(h(2),'ylim',[10^(-10) 10^(4)],'ytick',[10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2)]);
% set(h(2),'xlim',[10^(-2) 10^4],'xtick',[10^(-2) 10^(0) 10^2 10^4]);

%%
% plot(h(3),f_bz,psd_bz,'color','k');
% hold (h(3),'on');
% plot(h(3),f_bz1,psd_bz1,'color','b');
% 
% % plot(h(3),[Fce Fce],yrange,'color','r','linestyle','--');
% % plot(h(3),[Fce05 Fce05],yrange,'color','k','linestyle','--');
% % plot(h(3),[Fce01 Fce01],yrange,'color','b','linestyle','--');
% % plot(h(3),[Flh Flh],yrange,'color','g','linestyle','--');
% plot(h(3),[8 8],[10^(-10) 10^(4)],'color','g','linestyle','--');
% plot(h(3),[0.1 0.1],[10^(-10) 10^(4)],'color','g','linestyle','--');
% hold (h(3),'off');
% set(h(3),'yscale','log');
% set(h(3),'xscale','log');
% xlabel(h(3),'frequency (Hz)','fontsize',10,'fontname','times new roman');
% ylabel(h(3),{'B_{z PSD} [(nT)^{2}/Hz]'},'Interpreter','tex','fontsize',10);
% % irf_zoom(gca,'x',[10 3.2e4]);
% set(h(3),'ylim',[10^(-10) 10^(4)],'ytick',[10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2)]);
% set(h(3),'xlim',[10^(-2) 10^4],'xtick',[10^(-2) 10^(0) 10^2 10^4]);
% irf_legend(gca,'(a)',[0.99 0.98],'color','k','fontsize',10)
% irf_legend(gca,'0.1f_{ce}',[0.05 0.05],'color','b','fontsize',8)
% irf_legend(gca,'f_{ce}',[0.05 0.15],'color','r','fontsize',8)
% irf_legend(gca,'0.5f_{ce}',[0.05 0.25],'color','k','fontsize',8)
% irf_legend(gca,'f_{lh}',[0.05 0.35],'color','g','fontsize',8)


% irf_legend(gca,'(a)',[0.99 0.98],'color','k','fontsize',10)
% irf_legend(gca,'0.1f_{ce}',[0.05 0.05],'color','b','fontsize',8)
% irf_legend(gca,'f_{ce}',[0.05 0.15],'color','r','fontsize',8)
% irf_legend(gca,'0.5f_{ce}',[0.05 0.25],'color','k','fontsize',8)
% irf_legend(gca,'f_{lh}',[0.05 0.35],'color','g','fontsize',8)

% set(h(1:3),'ylim',yrange,'ytick',[10^(-10) 10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2)]);
% set(h(1:3),'xlim',[10^(0) 10^4],'xtick',[10^(0) 10^2 10^4]);
% set(h(4:6),'ylim',yrange,'ytick',[10^(-8) 10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2)]);
% set(h(4:6),'xlim',[10^(0) 10^4],'xtick',[10^(0) 10^2 10^4]);

% irf_zoom(gca,'y',yrange3);
% title(h(2),strcat(Tsta(12:22),' to',32,Tend(12:22),' UT'),'fontname','times new roman','fontweight','normal');  %32是空格

% set(gcf,'render','painters');
% 
% figname=['wave_psd' '_' Tsta(1:4) Tsta(6:7) Tsta(9:10) '-' ...
%     Tsta(12:13) Tsta(15:16) Tsta(18:19) Tsta(21:23) '_' Tend(12:13) Tend(15:16) Tend(18:19) Tend(21:23)];
% print(gcf, '-dpng', [figname '.png']);


