clc
clear 
Tcomputsta=clock; 
%--------------------------------------------
mms.db_init('local_file_db','E:\MMS\')
Tsta='2018-08-27T12:15:30.000Z'; 
Tend='2018-08-27T12:15:50.000Z'; 
tint=irf.tint('2018-08-27T12:15:30.000Z/2018-08-27T12:15:50.000Z');
Tnull='2018-08-27T12:15:43.045Z';
% Tnull='2015-09-19T07:43:30.555Z';
% Tnull='2015-09-19T07:43:30.945Z';
% Tnull='2015-09-19T07:43:30.76Z';
Time = Tnull;
time = irf_time(Time,'utc>epochtt');

Tnull2='2018-08-27T12:15:43.045Z';

Tssta = '2018-08-27T12:15:30.000Z';
Tsend = '2018-08-27T12:15:50.000Z';



ic=1:4;
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
% c_eval('Rgse?=mms.get_data(''R_gse'',tint,?);',ic);
Pos = mms.get_data('R_gse',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gseR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
% c_eval('R?=irf.ts2mat(Rgse?);',ic);


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


B2=irf_resamp(B2,B1);%将B2采集到B1的时间点上
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
tint=[iso2epoch(Tsta) iso2epoch(Tend)];

lth = 4;

for ic=1:4
  c_eval(['R?=irf_resamp(R?,B?);'],ic);
  c_eval(['B?=irf_tlim(B?,tint);'],ic);%返回指定时间时刻的数据
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
  c_eval(['Ro?=R?;'],ic);
  c_eval(['Bo?=B?;'],ic);
  c_eval(['Rts_sc?=ones(lth,3);'],ic);
end

gradB=c_4_grad('R?','B?','grad');
gradB2 = gradB;
idxnull=find(gradB(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
idxnull2=find(gradB(:,1,1)>=iso2epoch(Tnull2)); idxnull2=idxnull2(1);
idxsta=find(gradB(:,1,1)>=iso2epoch(Tssta)); idxsta=idxsta(1);
% idxend=find(gradB(:,1,1)>=iso2epoch(Tsend)); idxend=idxend(1);
dB_null=reshape(gradB(idxnull,2:end),3,3);

[V,D] = eig(dB_null);

    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
         [m,n]=find(D==min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))])));
         [mm,nn]=find(D==-(min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]))));
         if isempty(m)
             m=mm;
         else m=m;
         end
    else [m,n]=find(imag(D)==0&D~=0);
    end  
    e1=V(:,m)';
    e1=irf_norm(e1);
    etmp=[0 -1 0];%
    e3=cross(e1,etmp);          e3=irf_norm(e3);
    e2=cross(e3,e1);
    tm1 = [e1' e2' e3'];
    tm2 = inv(tm1);%求矩阵的逆

for ic=1:4
    c_eval(['B?=irf_newxyz(B?, e1,e2,e3);'],ic);
    c_eval(['R?=irf_newxyz(R?, e1,e2,e3);'],ic);
end


gradB=c_4_grad('R?','B?','grad');
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');

%construct B around null
idxnull=find(gradB(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
dB_null=reshape(gradB(idxnull,2:end),3,3);
[V,D] = eig(dB_null);
a1 = V(:,2)';          a2 = V(:,3)';%特征值较大的两个特征向量
as1=[0,a1(2),a1(3)];   as2=[0,a2(2),a2(3)];%为何不取第一个特征向量，
b1=as1/norm(as1);%e2在y-z面的投影反向
b2=as2/norm(as2);%
o1 = (a1+a2)/norm(a1+a2);
i1 = (a1-a2)/norm(a1-a2);
m1 = cross(o1,i1);
b3=(b1+b2)/norm(b1+b2);
c=dot(b3,[0,1,0]);
% c2=dot(a1,a2/(norm(a1)*norm(a2)));
angle=acosd(c);

dR1=inv(dB_null) * (B1(idxnull,2:4))';   
R_null=R1(idxnull,2:4)-dR1';

Rsc1=R1(idxnull,2:end);
Rsc2=R2(idxnull,2:end);
Rsc3=R3(idxnull,2:end);
Rsc4=R4(idxnull,2:end);

% Rbarycenter=(Rsc1+Rsc2+Rsc3+Rsc4)./[4.0 4.0 4.0]    %重心位置

Rsc1=Rsc1-R_null;
Rsc2=Rsc2-R_null;
Rsc3=Rsc3-R_null;
Rsc4=Rsc4-R_null;


j_null=j(idxnull,2:end) .* 1e9;



%% construct B around null
set(0,'defaultLineLineWidth', 0.5);
set(0,'defaultAxesFontSize', 12);
set(0,'defaultTextFontSize', 12);
set(0,'defaultAxesFontUnits', 'pixels');
fig1=figure( ...
          'Name','Dataset coverage', ...
          'Tag','XYXgsm');clf;
set(fig1,'PaperUnits','centimeters')

xSize = 800; ySize=800;
set(gcf,'position',[10,300,xSize,ySize])
set(gcf,'render','painters');
Paper_X = 20; Paper_Y = 20; 
coef=floor(max(xSize/Paper_X,ySize/Paper_Y));
FigSize_X = xSize/coef; FigSize_Y = ySize/coef;
xLeft2 = (Paper_X- FigSize_X)/2;  yTop2 = (Paper_Y- FigSize_Y)/2; 
set(gcf,'PaperSize', [Paper_X Paper_Y]); 
set(gcf,'PaperPosition',[xLeft2 yTop2 FigSize_X FigSize_Y])

% set(fig1,'Position',[10 10 xSize*coef ySize*coef])
% xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
% xLeft = 0; yTop = -1;
% set(fig1,'PaperPosition',[xLeft yTop xSize ySize])
% set(fig1,'Position',[10 10 xSize*coef ySize*coef])


h=[];

h(1)=axes('position',[0.1 0.1 0.8 0.8]); % [x y dx dy]
aaa=1;

BoxWid=200; 
Flimx=200;
Flimy=200;
Flimz=200;
Flim_ax=200;
Flim_ay=200;
Flim_az=200;
Flim_tx=50;
Flim_ty=50;
Flim_tz=50;


% BoxWid=250; 
% Flimx=200;
% Flimy=100;
% Flimz=75;
% Flim_ax=200;
% Flim_ay=100;
% Flim_az=75;
% Flim_tx=50;
% Flim_ty=50;
% Flim_tz=50;

% BoxWid=200;
Gray = 144;
clmax = 5;
%% inflow region
for t_mg=0
%     for l_mg=[-45 -40 -35 -29 -22 -14 -5 -2 2 5 14 22 29 35 40 45]   %30.945
%     for l_mg=[-50 -45 -40 -35 -29 -22 -14 -5 -2 2 5 14 22 29 35 40 45 50]   % 30.555(1)
    for l_mg=[-62:8:-14 -5 5 14:8:62]   %31.3 
        Xgrid=l_mg*i1(1)+t_mg*m1(1);
        Ygrid=l_mg*i1(2)+t_mg*m1(2); 
        Zgrid=l_mg*i1(3)+t_mg*m1(3);
        %-----inverse trace-----
        ii=1;
        Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;
        while abs(Xprev)<BoxWid & abs(Yprev)<BoxWid & abs(Zprev)<BoxWid
            Xcurt=Xprev;
            Ycurt=Yprev;
            Zcurt=Zprev;
            
            Br = dB_null * ([Xcurt Ycurt Zcurt])';%根据磁零点计算零点周围的磁场
%             Br = Br+([Bxbc Bybc Bzbc])';   %加上重心处磁场
            Bxcurt=Br(1);
            Bycurt=Br(2);
            Bzcurt=Br(3);
            Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
%             step=Bmcurt;
            step=1;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
            Xprev=Xcurt-stepvec(1);
            Yprev=Ycurt-stepvec(2);
            Zprev=Zcurt-stepvec(3);

            Xline(ii)=Xcurt;
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
%               plot3(gca, Xline, Yline, Zline,'color',[Gray,Gray,Gray]/256); hold on;
%          plot3(gca, Xline, Yline, Zline,'r'); hold on;
             Nlin=length(Xline);
             cline(Xline, Yline, Zline, Bmline, 0, clmax, jet); view(3); hold on;
%              arrP=fix(Nlin*0.7);
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP-1) Yline(arrP-1) Zline(arrP-1)],'b',6,12); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
        
        
        %-----positive trace-----
        ii=1;
        Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
        while abs(Xnext)<BoxWid & abs(Ynext)<BoxWid & abs(Znext)<BoxWid
            Xcurt=Xnext;
            Ycurt=Ynext;
            Zcurt=Znext;
            
            Br = dB_null * ([Xcurt Ycurt Zcurt])';
%             Br = Br+([Bxbc Bybc Bzbc])';  %加上重心处磁场
            Bxcurt=Br(1);
            Bycurt=Br(2);
            Bzcurt=Br(3);
            Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
%             step=Bmcurt;
              step=1;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
            Xnext=Xcurt+stepvec(1);
            Ynext=Ycurt+stepvec(2);
            Znext=Zcurt+stepvec(3);

            Xline(ii)=Xcurt;
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
%                plot3(gca, Xline, Yline, Zline,'color',[Gray,Gray,Gray]/256); hold on;
%         plot3(gca, Xline, Yline, Zline,'r'); hold on;
             Nlin=length(Xline);
             cline(Xline, Yline, Zline, Bmline, 0, clmax, jet); view(3); hold on;
%              arrP=fix(Nlin*0.7);
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+1) Yline(arrP+1) Zline(arrP+1)],'b',6,12); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
     
    end
end
%% outflow region
% for l_mg=[-62:8:-14 -5 5 14:8:62]   %30.555(1)
for l_mg=[-40 -35 -29 -22 -14 -5 5 14 22 29 35 40]%31.3 
        Xgrid=l_mg*o1(1);
        Ygrid=l_mg*o1(2); 
        Zgrid=l_mg*o1(3);
        %-----inverse trace-----
        ii=1;
        Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;
        while abs(Xprev)<BoxWid & abs(Yprev)<BoxWid & abs(Zprev)<BoxWid
            Xcurt=Xprev;
            Ycurt=Yprev;
            Zcurt=Zprev;
            
            Br = dB_null * ([Xcurt Ycurt Zcurt])';
%             Br = Br+([Bxbc Bybc Bzbc])';   %加上重心处磁场
            Bxcurt=Br(1);
            Bycurt=Br(2);
            Bzcurt=Br(3);
            Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
%             step=Bmcurt;
            step=1;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
            Xprev=Xcurt-stepvec(1);
            Yprev=Ycurt-stepvec(2);
            Zprev=Zcurt-stepvec(3);

            Xline(ii)=Xcurt;
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
%               plot3(gca, Xline, Yline, Zline,'color',[Gray,Gray,Gray]/256); hold on;
%              plot3(gca, Xline, Yline, Zline,'b'); hold on;
%             Nlin=length(Xline);
             cline(Xline, Yline, Zline, Bmline, 0, clmax, jet); view(3); hold on;
%             arrP=fix(Nlin*0.7);
%             arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP-1) Yline(arrP-1) Zline(arrP-1)],'b',6,12); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
        
        
        %-----positive trace-----
        ii=1;
        Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
        while abs(Xnext)<BoxWid & abs(Ynext)<BoxWid & abs(Znext)<BoxWid
            Xcurt=Xnext;
            Ycurt=Ynext;
            Zcurt=Znext;
            
            Br = dB_null * ([Xcurt Ycurt Zcurt])';
%             Br = Br+([Bxbc Bybc Bzbc])';  %加上重心处磁场
            Bxcurt=Br(1);
            Bycurt=Br(2);
            Bzcurt=Br(3);
            Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
%             step=Bmcurt;
              step=1;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
            Xnext=Xcurt+stepvec(1);
            Ynext=Ycurt+stepvec(2);
            Znext=Zcurt+stepvec(3);

            Xline(ii)=Xcurt;
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
%                plot3(gca, Xline, Yline, Zline,'color',[Gray,Gray,Gray]/256); hold on;
%              plot3(gca, Xline, Yline, Zline,'b'); hold on; 
%             Nlin=length(Xline);
             cline(Xline, Yline, Zline, Bmline, 0, clmax, jet); view(3); hold on;
%             arrP=fix(Nlin*0.7);
%             arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+1) Yline(arrP+1) Zline(arrP+1)],'b',6,12); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
     
end
% Jlinelngth=BoxWid;
% Jendp1=Jlinelngth*irf_norm(j_null);
% Jendp2=-Jendp1;
% arrow3(Jendp2, Jendp1,'s5',3,5); hold on;
% light('position',Jendp1); lighting gouraud;

%--------------------------------------------
% calculate the location series
lsn = 1;
% for ti = [idxsta:5:idxend] 
   
% idxnull2 = 4;
plot3(gca, [Rsc1(1) Rsc2(1) Rsc3(1) Rsc4(1)], [Rsc1(2) Rsc2(2) Rsc3(2) Rsc4(2)], ...
           [Rsc1(3) Rsc2(3) Rsc3(3) Rsc4(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc2(1) Rsc4(1) Rsc1(1) Rsc3(1)], [Rsc2(2) Rsc4(2) Rsc1(2) Rsc3(2)], ...
           [Rsc2(3) Rsc4(3) Rsc1(3) Rsc3(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc1(1)], [Rsc1(2)],[Rsc1(3)], 'ks', 'Linewidth',1, ...
           'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12); hold on;
plot3(gca, [Rsc2(1)], [Rsc2(2)],[Rsc2(3)], 'rs', 'Linewidth',1, ...
           'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12); hold on;
plot3(gca, [Rsc3(1)], [Rsc3(2)],[Rsc3(3)], 'gs', 'Linewidth',1, ...
           'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12); hold on;
plot3(gca, [Rsc4(1)], [Rsc4(2)],[Rsc4(3)], 'bs', 'Linewidth',1, ...
           'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12); hold on;
       
       
% for ii = 1:length(Rts_sc4)       
%     text(Rts_sc4(ii,1)+10, Rts_sc4(ii,2)+10,Rts_sc4(ii,3)+10,num2str(ii));
% end
%        
% quiver3(0,0,0,100*V(1,1),100*V(1,2),100*V(1,3)); hold on
% quiver3(0,0,0,100*V(2,1),100*V(2,2),100*V(2,3)); hold on
% quiver3(0,0,0,100*V(3,1),100*V(3,2),100*V(3,3)); hold off
% quiver3(0,0,0,0,100*b3(2),100*b3(3)); hold on
% quiver3(0,0,0,100*a1(1),100*a1(2),100*a1(3)); hold on
% quiver3(0,0,0,100*a2(1),100*a2(2),100*a2(3)); hold on
% quiver3(0,0,0,100*o1(1),100*o1(2),100*o1(3)); hold on
% quiver3(0,0,0,100*i1(1),100*i1(2),100*i1(3)); hold on
maxB=max(Bmax);
minB=min(Bmin);
       
% hcb=colorbar('peer',h(1),'East', 'TickDirection','out');
% posFig=get(h(1),'Position'); 
% left=posFig(1)+posFig(3)*1; low=posFig(2)+posFig(4)*0.87; width=posFig(3)/40; height=0.1;
% set(hcb,'Position',[left low width height]);
% ylabel(hcb,'|B| (nT)');
% caxis([0, 5]);   %here is derived from minB & maxB. it should be consistent with cline range

% colormap('jet')
% hcb=colorbar('peer',h(1),'East', 'TickDirection','out');
% posFig=get(h(1),'Position'); 
% left=posFig(1)+posFig(3)*1; low=posFig(2)+posFig(4)*0.87; width=posFig(3)/40; height=0.1;
% set(hcb,'Position',[left low width height]);
% % ylabel(hcb,'log_{10} B^2','nT^2 Hz^{-1}','fontsize',10);
% ylabel(hcb,{'log_{10} B^2','nT^2 Hz^{-1}'},'fontsize',7.5);
% caxis([-8, -2]); 




set(gca,'Xlim',[-Flimx Flimx], 'Ylim',[-Flimy Flimy], 'Zlim',[-Flimz Flimz]);
set(gca,'xtick',[-Flim_ax:Flim_tx:Flim_ax], 'ytick',[-Flim_ay:Flim_ty:Flim_ay], 'ztick',[-Flim_az:Flim_tz:Flim_az],'fontsize',13);
% set(gca,'Xlim',[-150 150], 'Ylim',[-150 150], 'Zlim',[-150 150]);
% set(gca,'xtick',[-150:50:150], 'ytick',[-150:50:150], 'ztick',[-150:50:150],'fontsize',12);
% set(gca,'Xlim',[-200 200], 'Ylim',[-200 200], 'Zlim',[-200 200]);
% set(gca,'xtick',[-200:50:200], 'ytick',[-200:50:200], 'ztick',[-200:50:200],'fontsize',12);
set(gca,'DataAspectRatio',[1 1.0 1]);
xlabel(gca,'e_{1} [km]');
ylabel(gca,'e_{2} [km]');
zlabel(gca,'e_{3} [km]');
% grid off;


box on
angles=get(gca,'view');
set(gca,'view',angles)


%% save figure
set(fig1,'render','painters');  %矢量
% figname=['-e1touying'];
% print(fig1, '-dpdf', [figname '.pdf']);
% set(fig1,'renderer','opengl');   %位图
% % axis off;
% figname=['Btopology_3D'];
% print(fig1, '-dpng','-r400',[figname '.png']);

% set(gca,'view',[0 90])
% figname=['Btopology_xy_plane'];
% print(fig1, '-dpng','-r400',[figname '.png']);

% set(gca,'view',[0 0])
% figname=['Btopology_xz_plane'];
% print(fig1, '-dpng','-r400',[figname '.png']);

  set(gca,'view',[90 0])
 figname=['Btopology_2D_2'];
% print(fig1, '-dpng','-r400',[figname '.png']);
% print(fig1, '-dpdf','-r400',[figname '.pdf']);
%  set(gca,'view',[135 26])

Tcomputend=clock;
Telaps=etime(Tcomputend, Tcomputsta)/60  %minute