%% 找涡旋之速度线性
%% 1/3程序
clc; close all
clear
mms.db_init('local_file_db','C:/MMS/')
% load err_4CV2
Tsta = '2017-07-06T17:31:57.000Z';
Tend = '2017-07-06T17:31:59.900Z';

% Tsta = '2017-07-06T01:44:55.000Z'; 
% Tend = '2017-07-06T01:45:04.000Z';
% tintStr = '06-July-2017';
tint=irf.tint('2017-07-06T17:31:30.000Z/2017-07-06T17:32:30.000Z');

Time = '2017-07-06T17:31:58.700Z';%别删
T_resta = '2017-07-06T17:31:58.000Z';
T_reend = '2017-07-06T17:31:59.500Z';

% Time = '2017-01-27T12:05:30.99Z';%别删
% T_resta = '2017-01-27T12:05:30.94Z';
% T_reend = '2017-01-27T12:05:31.04Z';

%--------------------------------------------
time = irf_time(Time,'utc>epochtt');
nh = 15;
hh1=15;
hh2=15;
for ic=1:4
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);

c_eval('Vegse?=mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('Ve?=irf.ts2mat(Vegse?);',ic);
c_eval('V?=irf.ts2mat(Vegse?);',ic);

c_eval('Rgse?=mms.get_data(''R_gse'',tint,?);',ic);
c_eval('R?=irf.ts2mat(Rgse?);',ic);

c_eval('e_r? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energy_brst'',tint);',ic);
%% smooth
% c_eval('BS?(:,2)=smooth(B?(:,2),nh);',ic);
% c_eval('BS?(:,3)=smooth(B?(:,3),nh);',ic);
% c_eval('BS?(:,4)=smooth(B?(:,4),nh);',ic);
% c_eval('BS?(:,1)=B?(:,1);',ic);
% c_eval('B?=BS?;',ic);
end

B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
tint=[iso2epoch(Tsta) iso2epoch(Tend)];

lth = 4;

kk = length (e_r1.data);
if mod(kk,2) == 1
   le = kk-2; lo = kk-3;
else
   le = kk-3; lo = kk-2;
end
for ic =1:4
    c_eval('e_r = e_r?;',ic)
    if e_r.data(1)>e_r.data(2)
        for ii = 1:2:le
%             c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('V?(ii+1,2:4)=(V?(ii+2,2:4)+V?(ii,2:4))/2;',ic)
        end
    else
        for ii = 2:2:lo
%             c_eval('ne?(ii+1,2)=(ne?(ii+2,2)+ne?(ii,2))/2;',ic)
            c_eval('V?(ii+1,2:4)=(V?(ii+2,2:4)+V?(ii,2:4))/2;',ic)
        end
    end
end

for ic=1:4
  c_eval(['R?=irf_resamp(R?,B?);'],ic);
  c_eval(['V?=irf_resamp(V?,B?);'],ic);
end
for ic =1:4
    c_eval('VeS?(:,1:4)=V?(:,1:4);',ic)
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh1);',ic);
    c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh1);',ic);
    c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh1);',ic);
    
%     c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh2);',ic);
%     c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh2);',ic);
%     c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh2);',ic);
%     c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh2);',ic);
%     c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh2);',ic);
%     c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh2);',ic);
%     c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh2);',ic);
%     c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh2);',ic);
%     c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh2);',ic);
%     c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh2);',ic);
%     c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh2);',ic);
%     c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh2);',ic);
%     c_eval('VeS?(:,2)=smooth(VeS?(:,2),hh2);',ic);
%     c_eval('VeS?(:,3)=smooth(VeS?(:,3),hh2);',ic);
%     c_eval('VeS?(:,4)=smooth(VeS?(:,4),hh2);',ic);
    c_eval('V?(:,1:4)=VeS?(:,1:4);',ic);
end

for ic=1:4
  c_eval(['B?=irf_tlim(B?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
  c_eval(['V?=irf_tlim(V?,tint);'],ic);
  
  c_eval(['Ro?=irf_tlim(R?,tint);'],ic);
  c_eval(['Rts_sc?=ones(lth,3);'],ic);
end

a=mean(V1,1);
b=mean(V2,1);
c=mean(V3,1);
d=mean(V4,1);
Vbg=(a(2:4)+b(2:4)+c(2:4)+d(2:4))/4;
M1=zeros(length(V1(:,1)),4);
M1(:,2)=Vbg(1);
M1(:,3)=Vbg(2);
M1(:,4)=Vbg(3);
for n=1:4
    c_eval('V?=V?-M1;',n);
end

gradB = c_4_grad('R?','B?','grad');
gradV = c_4_grad('R?','V?','grad');

idxnull=find(gradB(:,1,1)>=iso2epoch(Time)); idxnull=idxnull(1);
idxsta=find(gradB(:,1,1)>=iso2epoch(T_resta)); idxsta=idxsta(1)-1;
idxend=find(gradB(:,1,1)>=iso2epoch(T_reend)); idxend=idxend(1)+1;
% idxnull = 330;
dV_null=reshape(gradV(idxnull,2:end),3,3);

% for ii=idxsta:idxend
%     dBeach=reshape(gradB(ii,2:end),3,3);
%     dR1(ii,2:4)=(B1(ii+1,2:4)-B1(ii,2:4))*inv(dBeach');
%     dR2(ii,2:4)=(B2(ii+1,2:4)-B2(ii,2:4))*inv(dBeach');
%     dR3(ii,2:4)=(B3(ii+1,2:4)-B3(ii,2:4))*inv(dBeach');
%     dR4(ii,2:4)=(B4(ii+1,2:4)-B4(ii,2:4))*inv(dBeach');
%     dR1(ii,1)=B1(ii,1);
%     dR2(ii,1)=B1(ii,1);
%     dR3(ii,1)=B1(ii,1);
%     dR4(ii,1)=B1(ii,1);
% end

if 0
    load GV
for k = 1:34
    GVt(4*k-2:4*k,1:3) = reshape(GradV(k,2:end),3,3);
    [GVt(4*k-2:4*k,5:7),GVt(4*k-2:4*k,9:11)] = eig(GVt(4*k-2:4*k,1:3));   
end
    GVt(find(GVt==0))=0;
end

for ii=idxsta:idxend
    dVeach=reshape(gradV(ii,2:end),3,3);
    dR1(ii,2:4)=(V1(ii+1,2:4)-V1(ii,2:4))*inv(dVeach');
    dR2(ii,2:4)=(V2(ii+1,2:4)-V2(ii,2:4))*inv(dVeach');
    dR3(ii,2:4)=(V3(ii+1,2:4)-V3(ii,2:4))*inv(dVeach');
    dR4(ii,2:4)=(V4(ii+1,2:4)-V4(ii,2:4))*inv(dVeach');
    dR1(ii,1)=V1(ii,1);
    dR2(ii,1)=V1(ii,1);
    dR3(ii,1)=V1(ii,1);
    dR4(ii,1)=V1(ii,1);
    
    Vr1(ii,1)=V1(ii,1);
    Vr2(ii,1)=V1(ii,1);
    Vr3(ii,1)=V1(ii,1);
    Vr4(ii,1)=V1(ii,1);
    Verror(ii,1)=V1(ii,1);
end


for ii=idxsta:idxnull-1
    dR1(ii,6:8)=-sum(dR1(ii:idxnull,2:4));
    dR2(ii,6:8)=-sum(dR2(ii:idxnull,2:4));
    dR3(ii,6:8)=-sum(dR3(ii:idxnull,2:4));
    dR4(ii,6:8)=-sum(dR4(ii:idxnull,2:4));
end

for ii=idxnull+1:idxend
    dR1(ii,6:8)=sum(dR1(idxnull:ii,2:4));
    dR2(ii,6:8)=sum(dR2(idxnull:ii,2:4));
    dR3(ii,6:8)=sum(dR3(idxnull:ii,2:4));
    dR4(ii,6:8)=sum(dR4(idxnull:ii,2:4));
end

for ii=idxsta:idxend
    dR1(ii,10:12) = V1(idxnull,2:4)+(dV_null * (dR1(ii,6:8))')'; 
    dR1(ii,13) = sqrt(sum(dR1(ii,10:12).^2));
    dR1(ii,14) = sqrt(sum(dR1(ii,6:8).^2));
    V1(ii,5) = sqrt(sum(V1(ii,2:4).^2));
    Vr1(ii,2:4) = dR1(ii,10:12)-V1(ii,2:4);
    Vr1(ii,5) = sqrt(sum(Vr1(ii,2:4).^2));
    
    dR2(ii,10:12) = V2(idxnull,2:4)+(dV_null * (dR2(ii,6:8))')'; 
    dR2(ii,13) = sqrt(sum(dR2(ii,10:12).^2));
    dR2(ii,14) = sqrt(sum(dR2(ii,6:8).^2));
    V2(ii,5) = sqrt(sum(V2(ii,2:4).^2));    
    Vr2(ii,2:4) = dR2(ii,10:12)-V2(ii,2:4);
    Vr2(ii,5) = sqrt(sum(Vr2(ii,2:4).^2));
    
    dR3(ii,10:12) = V3(idxnull,2:4)+(dV_null * (dR3(ii,6:8))')'; 
    dR3(ii,13) = sqrt(sum(dR3(ii,10:12).^2));
    dR3(ii,14) = sqrt(sum(dR3(ii,6:8).^2));
    V3(ii,5) = sqrt(sum(V3(ii,2:4).^2));
    Vr3(ii,2:4) = dR3(ii,10:12)-V3(ii,2:4);
    Vr3(ii,5) = sqrt(sum(Vr3(ii,2:4).^2));

    dR4(ii,10:12) = V4(idxnull,2:4)+(dV_null * (dR4(ii,6:8))')'; 
    dR4(ii,13) = sqrt(sum(dR4(ii,10:12).^2));
    dR4(ii,14) = sqrt(sum(dR4(ii,6:8).^2));
    V4(ii,5) = sqrt(sum(V4(ii,2:4).^2));
    Vr4(ii,2:4) = dR4(ii,10:12)-V4(ii,2:4);
    Vr4(ii,5) = sqrt(sum(Vr4(ii,2:4).^2));   
    
    Verror(ii,2) = 100*(Vr1(ii,5)/V1(ii,5)+Vr2(ii,5)/V2(ii,5)+Vr3(ii,5)/V3(ii,5)+Vr4(ii,5)/V4(ii,5))/4;
end

h=irf_plot(5,'newfigure');
xSize=750; ySize=640;
set(gcf,'Position',[100 100 xSize ySize]);
mmscolors=[1 0 0; 0 1 0; 0 0 1]; %红绿蓝
set(h,'ColorOrder',mmscolors)
lnwid =1.5;

%% SC plot
ii =1;
% h(ii)=irf_panel('Sc_location');
% hold(h(ii),'on');
% for ic = 1:4
%     c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,14]),''Linewidth'',lnwid)',ic);
% end 
% hold(h(ii),'off');
% ylabel(h(ii),{'dR','(km)'},'Interpreter','tex');
% % set(h(ii),'ylim',[0 50],'ytick',[0:20:40]);
% % c_eval('set(h(ii),''ylim'',[V?(idxnull,2)-150 V?(idxnull,2)+150]);',ic);
% irf_legend(h(ii),{'MMS1','MMS2','MMS3','MMS4'},[0.98 0.8],'fontsize',12)
% irf_legend(h(ii),'(k)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(ii),'off');
% ii=ii+1;
%% MMS1
ic =1;
size=22;
h(ii)=irf_panel('Ve1');
hold(h(ii),'on');
% c_eval('p1=irf_plot(h(ii),V?(idxsta:idxend,[1,5]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,2]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,3]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,4]),''Linewidth'',lnwid)',ic);
hold(h(ii),'off');
ylabel(h(ii),{'MMS1','Ve','(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-5 40],'ytick',[0:20:40]);
% c_eval('set(h(ii),''ylim'',[V?(idxnull,2)-150 V?(idxnull,2)+150]);',ic);
irf_legend(h(ii),{'Vx','Vy','Vz'},[0.98 0.03],'fontsize',size)
irf_legend(h(ii),'(c)',[0.99 0.98],'color','k','fontsize',size)
grid(h(ii),'off');
ii=ii+1;

h(ii)=irf_panel('Ve1');
hold(h(ii),'on');
% c_eval('p2=irf_plot(h(ii),dR?(idxsta:idxend,[1,13]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,10]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,11]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,12]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
% irf_plot(h(2),[B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.5);
hold(h(ii),'off');
% ylabel(h(2),{'MMS1','Ve','(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-150 110],'ytick',[-150:50:100]);%31712
% set(h(ii),'ylim',[-300 500],'ytick',[-200:200:400]);%33485
% c_eval('set(h(2),''ylim'',[V?(idxnull,3)-150 V?(idxnull,3)+150]);',ic);
% irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',12)
% legend(h(ii),[p1,p2],'Observation','FOTE','location','northwest')
grid(h(ii),'off');
ii=ii+1;

%% MMS2
ic =2;
h(ii)=irf_panel('Ve2');
hold(h(ii),'on');
% c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,5]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,2]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,3]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,4]),''Linewidth'',lnwid)',ic);
hold(h(ii),'off');
ylabel(h(ii),{'MMS2','Ve','(km/s)'},'Interpreter','tex');
% c_eval('set(h(ii),''ylim'',[V?(idxnull,2)-150 V?(idxnull,2)+150]);',ic);
irf_legend(h(ii),{'Vx','Vy','Vz'},[0.98 0.03],'fontsize',size)
irf_legend(h(ii),'(d)',[0.99 0.98],'color','k','fontsize',size)
grid(h(ii),'off');
ii=ii+1;

h(ii)=irf_panel('Ve2');
hold(h(ii),'on');
% c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,13]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,10]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,11]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,12]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
% irf_plot(h(2),[B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.5);
hold(h(ii),'off');
% ylabel(h(4),{'MMS2','Ve','(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-210 40],'ytick',[-200:50:0]);%31712
% set(h(ii),'ylim',[0 480],'ytick',[100:100:400]);%33485
% c_eval('set(h(2),''ylim'',[V?(idxnull,3)-150 V?(idxnull,3)+150]);',ic);
% irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',12)
grid(h(ii),'off');
ii=ii+1;

%% MMS3
ic =3;
h(ii)=irf_panel('Ve3');
hold(h(ii),'on');
% c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,5]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,2]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,3]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,4]),''Linewidth'',lnwid)',ic);
hold(h(ii),'off');
ylabel(h(ii),{'MMS3','Ve','(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-5 40],'ytick',[0:20:40]);
% c_eval('set(h(ii),''ylim'',[V?(idxnull,2)-150 V?(idxnull,2)+150]);',ic);
irf_legend(h(ii),{'Vx','Vy','Vz'},[0.98 0.03],'fontsize',size)
irf_legend(h(ii),'(e)',[0.99 0.98],'color','k','fontsize',size)
grid(h(ii),'off');
ii+ii+1;

h(ii)=irf_panel('Ve3');
hold(h(ii),'on');
% c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,13]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,10]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,11]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,12]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
% irf_plot(h(2),[B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.5);
hold(h(ii),'off');
% ylabel(h(ii),{'MMS3_{Ve}','(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-180 130],'ytick',[-150:50:100]);%31712
% set(h(ii),'ylim',[100 350],'ytick',[150:100:350]);%33485
% c_eval('set(h(2),''ylim'',[V?(idxnull,3)-150 V?(idxnull,3)+150]);',ic);
% irf_legend(h(ii),'(f)',[0.99 0.98],'color','k','fontsize',12)
grid(h(ii),'off');
ii=ii+1;

%% MMS4
ic = 4;
h(ii)=irf_panel('Ve4');
hold(h(ii),'on');
% c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,5]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,2]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,3]),''Linewidth'',lnwid)',ic);
c_eval('irf_plot(h(ii),V?(idxsta:idxend,[1,4]),''Linewidth'',lnwid)',ic);
hold(h(ii),'off');
ylabel(h(ii),{'MMS4';'Ve';'(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-5 40],'ytick',[0:20:40]);
% c_eval('set(h(ii),''ylim'',[V?(idxnull,2)-150 V?(idxnull,2)+150]);',ic);
irf_legend(h(ii),{'Vx','Vy','Vz'},[0.98 0.03],'fontsize',size)
irf_legend(h(ii),'(f)',[0.99 0.98],'color','k','fontsize',size)
grid(h(ii),'off');
ii=ii+1;

h(ii)=irf_panel('Ve4');
hold(h(ii),'on');
% c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,13]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,10]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,11]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
c_eval('irf_plot(h(ii),dR?(idxsta:idxend,[1,12]),''Linewidth'',lnwid,''Linestyle'',''--'')',ic);
% irf_plot(h(2),[B1(:,1) B1(:,2)*0],'k--', 'Linewidth',0.5);
hold(h(ii),'off');
% ylabel(h(ii),{'MMS4_{Ve}','(km/s)'},'Interpreter','tex');
% set(h(ii),'ylim',[-130 90],'ytick',[-100:50:50]);%31712
% set(h(ii),'ylim',[-200 600],'ytick',[-300:200:500]);%33485
% c_eval('set(h(2),''ylim'',[V?(idxnull,3)-150 V?(idxnull,3)+150]);',ic);
% irf_legend(h(ii),'(h)',[0.99 0.98],'color','k','fontsize',12)
grid(h(ii),'off');
ii=ii+1;
% %% Error
h(ii)=irf_panel('error');
hold(h(ii),'on');
a2=irf_plot(h(ii),Verror(idxsta:idxend,[1,2]),'Linewidth',lnwid);
f2=fill([a2.XData,fliplr(a2.XData)],[zeros(1,length(a2.YData)),fliplr(a2.YData)],[0.8706 0.9216 0.9804]);
set(f2,'parent',h(ii))
hold(h(ii),'off');
ylabel(h(ii),{'\beta','[%]'},'Interpreter','tex');
set(h(ii),'ylim',[0 100],'ytick',[-20:20:40]);
% c_eval('set(h(ii),''ylim'',[V?(idxnull,2)-150 V?(idxnull,2)+150]);',ic);
irf_legend(h(ii),'(p)',[0.99 0.98],'color','k','fontsize',12)
grid(h(ii),'off');
ii=ii+1;

% %% error eigen value
% h(ii)=irf_panel('error2');
% hold(h(ii),'on');
% %irf_plot([eigVal_err(:,1) eigVal_err(:,2)], 'color','b', 'Linewidth',lnwid); hold on;
% % irf_plot(h(6),eigVal_err_v2, 'b', 'Linewidth',0.5);
% a3=irf_plot(h(ii),[err_4CV(:,1) err_4CV(:,2)], 'k', 'Linewidth',lnwid);
% % irf_plot(h(6),[err_4C(:,1) err_4C(:,2)], 'r', 'Linewidth',lnwid);
% f3=fill([a3.XData,fliplr(a3.XData)],[zeros(1,length(a3.YData)),fliplr(a3.YData)],[0.8706 0.9216 0.9804]);
% set(f3,'parent',h(ii))
% % irf_plot(h(6),eigVal_err_v2B, 'c', 'Linewidth',0.5);
% % irf_plot(h(6),eigVal_errB, 'c', 'Linewidth',0.5);
% irf_plot(h(ii),[err_4CV(:,1) err_4CV(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;
% hold(h(ii),'off');
% set(h(ii),'Ylim',[0 75],'ytick',[-20:20:60]);
% %ylabel('|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}| [%]');
% ylabel({'|\nabla\cdot\bf{nV}|','|\nabla\times\bf{nV}|', '[%]'});
% % ylabel('\xi [%]','fontname','times new roman','fontweight','normal');
% irf_legend(h(ii),'(q)',[0.99 0.98],'color','k','fontsize',12)
% grid(h(ii),'off');

%% an
irf_plot_ylabels_align(h)
tint=[iso2epoch(T_resta) iso2epoch(T_reend)];
T0=dR1(idxnull);
for ii = 1:7
    c_eval('irf_pl_mark(h(?),T0,''k'',''Linestyle'',''--'')',ii)
end
% tint2=irf.tint('2015-10-16T13:06:31.07Z/2015-10-16T13:06:31.22Z');
% for ii = 1:8
%     c_eval('irf_pl_mark(h(?),tint2,[0.7,0.7,0.7])',ii)
% end
irf_zoom(h,'x',tint)
xSize = 200; ySize=400;
set(gcf,'position',[10,300,xSize,ySize])
set(gcf,'render','painters');
Paper_X = 20; Paper_Y = 40; 
coef=floor(max(xSize/Paper_X,ySize/Paper_Y));
FigSize_X = xSize/coef; FigSize_Y = ySize/coef;
xLeft2 = (Paper_X- FigSize_X)/2;  yTop2 = (Paper_Y- FigSize_Y)/2; 
set(gcf,'PaperSize', [Paper_X Paper_Y]); 
set(gcf,'PaperPosition',[xLeft2 yTop2 FigSize_X FigSize_Y])

irf_plot_ylabels_align(h)
irf_zoom(h,'x',tint);
irf_zoom('y',ylim)
% set(gcf,'paperpositionmode','auto')
set(h,'fontsize',size);
figname = ['Vxyz_31712_try'];
% print(gcf, '-dpdf','-r600',[figname '.pdf']);