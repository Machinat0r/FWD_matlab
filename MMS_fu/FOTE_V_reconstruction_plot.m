%% 找涡旋之重构彩图
clear ;close all;clc
Tcomputsta=clock; 
%--------------------------------------------
mms.db_init('local_file_db','/Volumes/SPART-NAS/Data/MMS/')
%Tsta='2015-10-16T13:07:00Z'; 
%Tend='2015-10-16T13:07:03Z';
%Tnull='2015-10-16T13:07:02.250Z';
% e1=[-0.552, 0.834, -0.008];
hh1=10;
% % % Tsta = '2017-07-06T17:31:58.000Z';
% % % Tend = '2017-07-06T17:31:59.500Z';
Tsta = '2019-08-05T16:24:00.000Z';
Tend = '2019-08-05T16:25:00.000Z';
% Tnull='2017-01-27T12:05:31.712Z';
Tnull='2019-08-05T16:24:31.137Z';

% e1=[0.928580204157578 -0.145037102078109 -0.341618271565582];%31.712，别删
% e1=[0.384889128652657 -0.890735064677563 -0.241767250054376];%30.99,别删
% e1=[0.440732154171349 -0.213297770361358 -0.871928454311680];
e1 = [-0.5079 0.5042 0.6984];
% e1=[1,0,0];
% e1 = [-0.4172 0.6567 0.3237];
%--------------------------------------------

tint=irf.tint('2019-08-05T16:24:00.000Z/2019-08-05T16:25:00.000Z');

%--------------------------------------------
%background magnetic field
ic=1:4;
c_eval('Vegse?=mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
% c_eval('Vegse? = irf_gse2gsm(Vegse?);',ic);
c_eval('VeS?=irf.ts2mat(Vegse?);',ic);
c_eval('Vet?=Vegse?.abs;',ic);
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
c_eval('B?=irf.ts2mat(Bxyz?);',ic);
c_eval('B? = irf_gse2gsm(B?);',ic);
c_eval('Bt?=Bxyz?.abs;',ic);
c_eval('Rgse?=mms.get_data(''R_gse'',tint,?);',ic);
c_eval('R?=irf.ts2mat(Rgse?);',ic);
c_eval('R? = R?(:,1:4);',ic);
c_eval('e_r? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_energy_brst'',tint);',ic);

%% smooth_step1
kk = size(VeS1,1);
if mod(kk,2) == 1
   le = kk-2; lo = kk-3;%even偶数（这里当它是奇数），odd奇数
else
   le = kk-3; lo = kk-2;%even偶数，odd奇数
end
for ic =1:4
    c_eval('e_r = e_r1;',ic)
    if e_r.data(1)>e_r.data(2)
        for ii = 1:2:le
            c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
        end
    else
        for ii = 2:2:lo
            c_eval('VeS?(ii+1,2:4)=(VeS?(ii+2,2:4)+VeS?(ii,2:4))/2;',ic)
        end
    end
end


VeS1=irf_resamp(VeS1,B1);
VeS2=irf_resamp(VeS2,VeS1);
VeS3=irf_resamp(VeS3,VeS1);
VeS4=irf_resamp(VeS4,VeS1);
%% smooth_step2
for ic =1:4
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
    c_eval('Ve?(:,1:4)=VeS?(:,1:4);',ic);
end

tint=[iso2epoch(Tsta) iso2epoch(Tend)];
for ic=1:4
  c_eval(['R?=irf_resamp(R?,Ve?);'],ic);
  c_eval(['Ve?=irf_tlim(Ve?,tint);'],ic);
  c_eval(['R?=irf_tlim(R?,tint);'],ic);
end


%--------------------------------------------
% change the coordinates to eigen vector corrdinate
e1=irf_norm(e1);
etmp=[0 1 0];
e3=cross(e1,etmp);          e3=irf_norm(e3);
e2=cross(e3,e1);

%L=[0.79 0.58 0.19]; L=irf_norm(L);
%etmp=[0 0 1];
%M=cross(etmp,L);          M=irf_norm(M);
%N=cross(L,M);

for ic=1:4
   c_eval(['Ve?=irf_newxyz(Ve?, e1,e2,e3);'],ic);
   c_eval(['R?=irf_newxyz(R?, e1,e2,e3);'],ic);
end
%--------------------------------------------
%for ic=1:4
%   c_eval('R?(:,2:end)=R?(:,2:end)-R3(:,2:end);',ic)
%end
gradVe=c_4_grad('R?','Ve?','grad');
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','Ve?');

%construct B around null
idxnull=find(gradVe(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
dVe_null=reshape(gradVe(idxnull,2:end),3,3);%选定零点时刻的Ve梯度矩阵

dR1=inv(dVe_null) * (Ve1(idxnull,2:4))';% '为转置，inv为逆矩阵，dR1为卫星距离磁零点的距离【null速度为0＋一阶泰勒展开假设】
R_null=R1(idxnull,2:4)-dR1';%R1=R_null+dR1

Rsc1=R1(idxnull,2:end);%磁零点时刻卫星的位置——绝对
Rsc2=R2(idxnull,2:end);
Rsc3=R3(idxnull,2:end);
Rsc4=R4(idxnull,2:end);

Rsc1=Rsc1-R_null;%磁零点时刻卫星的位置——相对磁零点
Rsc2=Rsc2-R_null;
Rsc3=Rsc3-R_null;
Rsc4=Rsc4-R_null;

c_eval('RR? = R?(:,2:end) - R_null;')
RR_mean = 0.25*(RR1 + RR2 + RR3 + RR4);

tint3 = '2019-08-05T16:24:31.000Z/2019-08-05T16:24:31.100Z';
% c_eval('Vegse? = irf.ts2mat(Vegse?);');
gradVe_3=c_4_grad('R?','Ve?','grad');
gradVe_3 = irf_tlim(gradVe_3,tint3);
Ve1_3 = irf_tlim(Ve1,tint3);
c_eval('R? = irf_tlim(R?,tint3);')
R_sc_all = zeros(size(gradVe_3,1),3);
for i_grad= 1:size(gradVe_3,1)
    dVe_null_i=reshape(gradVe_3(i_grad,2:end),3,3);
    dR1_i=inv(dVe_null_i) * (Ve1_3(i_grad,2:4))';
    R_sc_all(i_grad,:) = dR1_i';
end
%% construct B around null

set(0,'defaultLineLineWidth', 0.5);
set(0,'defaultAxesFontSize', 12);
set(0,'defaultTextFontSize', 12);
set(0,'defaultAxesFontUnits', 'pixels');
fig1=figure( ...
          'Name','Dataset coverage', ...
          'Tag','XYXgsm');clf;
set(fig1,'PaperUnits','centimeters')
xSize = 1000; ySize = 5000; coef=floor(min(1200/xSize,1200/ySize));
% xLeft = 0; yTop = -1;
% set(fig1,'PaperPosition',[xLeft yTop xSize ySize]);
set(fig1,'Position',[100 100 xSize ySize]);

% h=[];

% h(1)=axes('position',[0.1 0.1 0.8 0.8]); % [x y dx dy]
aaa=1;
BoxWid=3e4; %图坐标上下限[-100 100]
Vup=500;

%for Xgrid=[-40 -38 -35 -31 -26 -20 -13 -5 5 13 20 26 31 35 38 40]
% for Xgrid=[-2 -1.3 -0.5 0.5 1.3 2 ]%三个循环，Xgrid变，theta变，X/Y/Z prev/curt变
i_plot = 1;
for Xgrid=[-37.5,  -30,-25]%31712
% for Xgrid=[-25]%31712
    for theta=[0]*pi/180

    fig1=figure(i_plot);clf;
    set(fig1,'PaperUnits','centimeters')
    xSize = 1000; ySize = 5000; coef=floor(min(1200/xSize,1200/ySize));
    % xLeft = 0; yTop = -1;
    % set(fig1,'PaperPosition',[xLeft yTop xSize ySize]);
    set(fig1,'Position',[100 100 xSize ySize]);

%     for theta=[0:120:240]*pi/180
    %for theta=[0:90:270]*pi/180
% for theta=[0:90:350]*pi/180
        Ygrid=Xgrid*cos(theta);%暂时不理解他在图上是哪一段，theta让人特别迷惑，先往下看，别纠结
        Zgrid=Xgrid*sin(theta);
        
        %-----inverse trace-----
        ii=1;
        Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;%在此之前为Xgrid和theta变,控制XYZgrid;此后为
        Bmcurt=0;
        while abs(Xprev)<BoxWid & abs(Yprev)<BoxWid & abs(Zprev)<BoxWid & Bmcurt<Vup%速度的绝对值小于阈值
            Xcurt=Xprev;
            Ycurt=Yprev;
            Zcurt=Zprev;
            Ver =dVe_null * ([Xcurt Ycurt Zcurt])';%一阶泰勒+零点速度为0+零点坐标位置为0，算速度分量
            Bxcurt=Ver(1);
            Bycurt=Ver(2);
            Bzcurt=Ver(3);
%             Bxcurt=Br(1)+B3(idxnull,2);
%             Bycurt=Br(2)+B3(idxnull,3);
%             Bzcurt=Br(3)+B3(idxnull,4);
            Bmcurt=sqrt(Ver(1).^2+Ver(2).^2+Ver(3).^2);%速度绝对值
            %step=Bmcurt;
            step=1;
            stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;%步长向量，速度矢量方向，“逆”着磁力线画
            Xprev=Xcurt-stepvec(1);%下一次while循环的XYZcurt为逆着磁力线往前推1个单位长度
            Yprev=Ycurt-stepvec(2);
            Zprev=Zcurt-stepvec(3);

            Xline(ii)=Xcurt;%XYZ点的|B|对应存储
            Yline(ii)=Ycurt;
            Zline(ii)=Zcurt;
            Bmline(ii)=Bmcurt;
            ii=ii+1;
        end
        if exist('Xline')
            % plot3(gca, Xline, Yline, Zline,'b'); hold on;
            Nlin=length(Xline);
%           cline(Xline, Yline, Zline, Bmline, 0, 1500, cool); view(3); hold on;%就是画个图，字面意思
            % cline(Xline, Yline, Zline, Bmline, 0, Vup, jet); view(3); hold on;%就是画个图，字面意思
            arrP=fix(Nlin*0.65);%315%舍入到最接近的整数,向0取整%完全没懂，但似乎影响不大
%             daspect([1,1,1]); 
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP-1) Yline(arrP-1) Zline(arrP-1)],'r',2.5,3.5); hold on;
            arrP=fix(Nlin*0.5);%115
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP-1) Yline(arrP-1) Zline(arrP-1)],'r',2.5,3.5); hold on;
            arrP=fix(Nlin*0.7);%250
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP-1) Yline(arrP-1) Zline(arrP-1)],'r',2.5,3.5); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;%计算这条磁力线上磁场强度最大值，最小值
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
        
        
        %-----positive trace-----
        ii=1;
        Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
        Bmcurt=0;
        while abs(Xnext)<BoxWid & abs(Ynext)<BoxWid & abs(Znext)<BoxWid & Bmcurt<Vup
            Xcurt=Xnext;
            Ycurt=Ynext;
            Zcurt=Znext;
            Ver = dVe_null * ([Xcurt Ycurt Zcurt])';
            Bxcurt=Ver(1);
            Bycurt=Ver(2);
            Bzcurt=Ver(3);
%             Bxcurt=Br(1)+B3(idxnull,2);
%             Bycurt=Br(2)+B3(idxnull,3);
%             Bzcurt=Br(3)+B3(idxnull,4);
            Bmcurt=sqrt(Ver(1).^2+Ver(2).^2+Ver(3).^2);
            %step=Bmcurt;
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
            % plot3(gca, Xline, Yline, Zline,'r'); hold on;
            Nlin=length(Xline);
            if Xline~=0
                % cmap = 'jet';
                cmap = othercolor('RdYlGn11');
                cmap = flip(cmap);

            n = 256;  % 想要的渐变级数
            
            % 先定义那 6 个关键颜色（已从图中取近似 RGB 并归一化到 [0,1]）
            cols = [
                18, 128, 144;   % 深青  #128090
                39, 167, 176;   % 青    #27A7B0
                252,182, 41;    % 金黄  #FCB629
                217,120,105;    % 暗粉  #D97869
                215,142,165;    % 浅紫  #D78EA5
                172, 59,125;    % 洋红  #AC3B7D
            ] / 255;
            
            % 在位置 [0,1] 均匀布置这 6 点
            pos = linspace(0,1,size(cols,1))';
            
            % 对 [0,1] 做线性插值，生成 n×3 的 cmap
            t = linspace(0,1,n)';
            cmap = interp1(pos, cols, t, 'linear');

            % n = 256;  % 渐变级数
            % pos = [0; 1/3; 2/3; 1];
            % cols = [
            %     0.6, 0.0, 0.0;   % 深红  (RGB [153, 0,   0  ]/255)
            %     1.0, 0.6, 0.6;   % 浅红  (RGB [255,153,153]/255)
            %     0.6, 0.8, 1.0;   % 浅蓝  (RGB [153,204,255]/255)
            %     0.0, 0.0, 0.6;   % 深蓝  (RGB [  0,  0,153]/255)
            % ];
            % 
            % t = linspace(0,1,n)';
            % cmap = interp1(pos, cols, t, 'linear');
            % cmap = flip(cmap);

          cline_arrow(Xline, Yline, Zline, Bmline, 0, Vup, cmap);hold on;
          i_plot=i_plot+1;
            end

          % cline(Xline, Yline, Zline, Bmline, 0, Vup, jet); view(3); hold on;
            % % % % % cline(Xline, Yline, Zline, Bmline, 0, Vup, jet); view(3); hold on;
            arrP=fix(Nlin*0.55);%315
%             arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+1) Yline(arrP+1) Zline(arrP+1)],'r',2.5,3.5); hold on;
            arrP=fix(Nlin*0.7);%115
%             arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+1) Yline(arrP+1) Zline(arrP+1)],'r',2.5,3.5); hold on;
            arrP=fix(Nlin*0.55);%250
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+1) Yline(arrP+1) Zline(arrP+1)],'r',2.5,3.5); hold on;
%              arrow3([Xline(arrP) Yline(arrP) Zline(arrP)], [Xline(arrP+4) Yline(arrP+4) Zline(arrP+4)],'r',2.5,3.5); hold on;
            Bmax(aaa)=max(Bmline); Bmin(aaa)=min(Bmline); aaa=aaa+1;
        end
        clear Xline
        clear Yline
        clear Zline
        clear Bmline
        clear Nlin
        
        box_bd1 = 1000;box_bd2 = 250;box_bd3 = 300;
        caxis([0, Vup]);   %here is derived from minB & maxB. 【it should be consistent with cline range】
        set(gca, 'Ylim',[-box_bd1 box_bd1], 'Ylim',[-box_bd2 box_bd2], 'Zlim',[-box_bd3 box_bd3]);
        set(gca,'xtick',[-box_bd1:box_bd1:box_bd1], 'ytick',[-box_bd2:box_bd2:box_bd2], 'ztick',[-box_bd3:box_bd3:box_bd3],'fontsize',23);
        set(gca,'DataAspectRatio',[5 1.0 1]);
        set(gca,'view',[90 0])
        set(gcf,'render','painters');
        set(gcf,'paperpositionmode','auto')
        print(gcf, '-dpdf', ['/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/new/打死不改了版/Ve_resconstruction/Ver2/3', ...
            num2str(i_plot),'.pdf']);    

     
    end
end


% Jlinelngth=BoxWid;
% Jendp1=Jlinelngth*irf_norm(j_null);
% Jendp2=-Jendp1;
% arrow3(Jendp2, Jendp1,'s5',3,5); hold on;
% light('position',Jendp1); lighting gouraud;


plot3(gca, [Rsc1(1) Rsc2(1) Rsc3(1) Rsc4(1)], [Rsc1(2) Rsc2(2) Rsc3(2) Rsc4(2)], ...%画线，四点出三线
           [Rsc1(3) Rsc2(3) Rsc3(3) Rsc4(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc2(1) Rsc4(1) Rsc1(1) Rsc3(1)], [Rsc2(2) Rsc4(2) Rsc1(2) Rsc3(2)], ...%同画线，四点出另外三线
           [Rsc2(3) Rsc4(3) Rsc1(3) Rsc3(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc1(1)], [Rsc1(2)],[Rsc1(3)], 'ks', 'Linewidth',1, ...%画点
           'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12); hold on;
plot3(gca, [Rsc2(1)], [Rsc2(2)],[Rsc2(3)], 'rs', 'Linewidth',1, ...
           'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12); hold on;
plot3(gca, [Rsc3(1)], [Rsc3(2)],[Rsc3(3)], 'gs', 'Linewidth',1, ...
           'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12); hold on;
plot3(gca, [Rsc4(1)], [Rsc4(2)],[Rsc4(3)], 'bs', 'Linewidth',1, ...
           'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12); hold on;   
RSC(1,1)=(Rsc1(1,1)+Rsc2(1,1)+Rsc3(1,1)+Rsc4(1,1))/4;
RSC(2,1)=(Rsc1(1,2)+Rsc2(1,2)+Rsc3(1,2)+Rsc4(1,2))/4;
RSC(3,1)=(Rsc1(1,3)+Rsc2(1,3)+Rsc3(1,3)+Rsc4(1,3))/4;
plot3(R_sc_all(:,1),R_sc_all(:,2),R_sc_all(:,3));hold off;

maxB=max(Bmax);
minB=min(Bmin);
       
% hcb=colorbar('peer',h(1),'North', 'XDir','reverse', 'TickDirection','out', 'XAxisLocation','top');
% posFig=get(h(1),'Position'); 
% left=posFig(1)+posFig(3)*0.7; low=posFig(2)+posFig(4)*4/5; width=posFig(3)/4; height=0.015;

% hcb=colorbar('peer',h(1),'East','TickDirection','out');
% posFig=get(h(1),'Position'); 
% left=posFig(1)+posFig(3)*1; low=posFig(2)+posFig(4)*0.87; width=posFig(3)/40; height=0.1;
% set(hcb,'Position',[left low width height]);
% ylabel(hcb,'|Ve|');


box_bd1 = 1000;box_bd2 = 250;box_bd3 = 300;
caxis([0, Vup]);   %here is derived from minB & maxB. 【it should be consistent with cline range】
set(gca, 'Ylim',[-box_bd1 box_bd1], 'Ylim',[-box_bd2 box_bd2], 'Zlim',[-box_bd3 box_bd3]);
set(gca,'xtick',[-box_bd1:box_bd1:box_bd1], 'ytick',[-box_bd2:box_bd2:box_bd2], 'ztick',[-box_bd3:box_bd3:box_bd3],'fontsize',23);
% set(gca,'DataAspectRatio',[1 1.0 1]);
% xlabel(gca,'e_{1} [km]','fontsize',20);
% ylabel(gca,'e_{2} [km]','fontsize',20);
% zlabel(gca,'e_{3} [km]','fontsize',20);
xlabel(gca,'x [km]','fontsize',20);
ylabel(gca,'y [km]','fontsize',20);
zlabel(gca,'z [km]','fontsize',20);
%irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.98 0.26],'color','cluster')
grid off;

set(gca,'YTickLabel',[]);set(gca,'ZTickLabel',[])

%angles=get(gca,'view');
%set(gca,'view',angles)



%% save figure
% set(fig1,'render','painters');
% figname=['B_around_null_cycle10000_spiral_cluster'];
% print(fig1, '-dpdf', [figname '.pdf']);

% set(fig1,'renderer','opengl');
%axis off;
% figname=['Vetopology'];
% set(gca,'view',[-46.5949303652596,20.7818750000000])%31712角度
% set(gca,'view',[11.3094596268930	28.0706250000000])%3099角度
%  print(fig1, '-dpng','-r400',[figname '.png']);
% set(fig1, 'Position', [100, 100, 600, 600]);
% colorbar;
% viewInfo=get(gca,'view');

% set(gca,'view',[0 90])
% figname=['Vetopology_xy_33485'];
% print(fig1, '-dpng','-r400',[figname '.png']);
% 
% set(gca,'view',[0 0])
% figname=['Vetopology_xz_33485'];
% print(fig1, '-dpng','-r400',[figname '.png']);

% % % set(gca,'view',[90 0])
% % % figname=['Vetopology'];
% % % set(gcf,'render','painters');
% % % set(gcf,'paperpositionmode','auto')
% % % print(gcf, '-dpdf', ['/Users/fwd/Documents/Ti~mor~/M/Sandglass/Nat/submission/Figures/new/打死不改了版/Ve_res', '.pdf']);    
