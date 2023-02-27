%% 1.判断通量曲线大于100eV处是否存在斜率=0
%% 2.大于10/20/30eV 存在连续五个斜率绝对值小于0的点




clc;
clear;
Tsta='2017-11-15T00:00:00.000Z';
Tend='2017-11-15T04:00:00.000Z';


% Tsta='2017-11-15T18:20:00.000Z';
% Tend='2017-11-15T18:30:00.000Z';

% Tsta='2017-11-15T05:55:00.000Z';
% Tend='2017-11-15T06:00:00.000Z';

% Tsta='2017-11-15T23:50:00.000Z';
% Tend='2017-11-15T23:59:00.000Z';

% Tsta='2017-11-15T22:20:00.000Z';
% Tend='2017-11-15T22:30:00.000Z';

% Tsta='2017-11-15T20:00:00.000Z';
% Tend='2017-11-15T20:10:00.000Z';

Tmia='/';
tintStr=[Tsta,Tmia,Tend];

Xsvyspec_NAME ='C:\Matlab\bin\新建文件夹\fwd\zcl\zcl\mvn_swe_l2_svyspec_20171115_v04_r04.cdf';
Xsvyspec_file= dataobj(Xsvyspec_NAME);%读文件
espectrum=get_variable(Xsvyspec_file,'diff_en_fluxes');%单位为eV/eV cm2 s sr，?

Energy=espectrum.DEPEND_1.data;%能道
T_espec=espectrum.DEPEND_0.data;%时间
T_espec=irf_time(T_espec,'ttns>epoch');

espectrum=[T_espec double(espectrum.data)];
espectrum=irf_tlim(espectrum,tintStr);

Energy=double(fliplr(Energy')); %矩阵转置、能道
espectrum(:,2:end)=fliplr(espectrum(:,2:end));%数组左右翻转（大到小）变（小到大）

%% 相空间密度的计算
ePSD=espectrum;
% Energy2=Energy.^2; %能道的平方 单位为eV
% Energy2=repmat(Energy2,size(espectrum,1),1);
% ePSD(:,2:end)=(ePSD(:,2:end)./Energy2).*0.1593;% s^3 km(-6)


%% plot
% TIME1='2017-11-15T21:28:30.07Z'
% TIME1='2017-11-15T18:25:30.07Z'
TIME1=Tsta;
TIME2=[Tend(1:18),num2str(1),'.07Z'];

AA_t1=iso2epoch(TIME1);
t_dif1=abs(ePSD(:,1)-AA_t1);
At1=find(min(min(t_dif1))==t_dif1);
TIME_sure1=ePSD(At1,1);
cc1=find(ePSD(:,1)==TIME_sure1);
TIME_sure1=epoch2iso(TIME_sure1);
AA_t2=iso2epoch(TIME2);
t_dif2=abs(ePSD(:,1)-AA_t2);
At2=find(min(min(t_dif2))==t_dif2);
TIME_sure2=ePSD(At2,1);
cc2=find(ePSD(:,1)==TIME_sure2);



%YY_PSD1=ePSD(cc1,2:end);
ePSD=log10(ePSD(cc1:cc2,2:end));%
Energy_lg=log10(Energy); 
%  通量关于能量的对数的斜率
gradient=zeros(cc2,64);
for i =2:63
gradient(:,i)=(ePSD(:,i-1)-ePSD(:,i+1))./(Energy_lg(:,i-1)-Energy_lg(:,i+1));
end   %求斜率方法
gradient(:,1)=gradient(:,2);
gradient(:,64)=gradient(:,63)
abs_gra=abs(gradient);


%% 1.判断80eV后是否存在连续几个点斜率绝对值小于1
   
%        judge=zeros(cc2,1);
% for d=1:cc2
%        num=0;     %斜率小于1连续点个数
% 
%  for j=31:51
%   k=find(abs_gra(d,j)<=1);
% 
%   if  length(k)==1
%      num=num+1
%   else
%      num=0;
%   end
%   
%   
%   if num==3
%      judge(d)=1;    %满足连续d个点斜率小于1
%      break
%   else
%       judge(d)=NaN;
%   end
%   
%  end 
%  
% end
%% 2.1.100eV后曲线是否有极大值（导数为0）
judge=zeros(cc2,1)
for d =1:cc2
    for j=32:51
        if gradient(d,j-1)>=0  &&  gradient(d,j-2)>=0 && gradient(d,j+1)<=0  &&  gradient(d,j+2)<=0 && gradient(d,j-3)>=0&&gradient(d,j+3)<=0  ...
%                      abs_gra(d,j)<=2   && abs_gra(d,j+1)<=2 && abs_gra(d,j-1)<=2 ...
%              %尽可能使导数严格单调递减

            judge(d)=1;
            break
        else
            judge(d)=NaN;
        end
    end
end

%% 2.2 100eV后曲线有极大值（某能道通量大于两侧）
% judge=zeros(cc2,1)
% for d =1:cc2
%     for j=31:51
%         if  ePSD(d,j)>=ePSD(d,j-1) && ePSD(d,j-1)>=ePSD(d,j-2)&& ePSD(d,j)>=ePSD(d,j+1) && ePSD(d,j+1)>=ePSD(d,j+2)
%             judge(d)=1;
%             break
%         else
%             judge(d)=NaN;
%         end
%     end
% end
%%  2.3差分算法  不是很行
% ePSD=ePSD';
% dY=diff(ePSD(:,1:cc2));
% for d =1:cc2
%     for j=30:50
%         if dY(j,d)>=0  && dY(j+1,d)<=0 &&dY(j-1,d)>=0 &&dY(j+2)<=0
%              judge(d)=1;
%             break
%         else
%             judge(d)=NaN;
%         end
%     end
% end
%% 最大值在100e之后取到，并左右两边值都小于他
% judge=zeros(cc2,1);
% for d =1:cc2
%    [m,index]=max(ePSD(d,31:51));
%     
%     if ePSD(d,index+29)<m &&ePSD(d,index+31)<m &&ePSD(d,index+28)<ePSD(d,index+29) && ePSD(d,index+31)>ePSD(d,index+32)
%         judge(d)=1;
%         
%     else
%         judge(d,1)=NaN;
%     end
% end
%  [m,index]=max(ePSD(555,31:51));
%% 肉眼观察的两带时间录入  
%  judge(1:cc2,1)=1

%% 把judge=1的时间对应的空间位置存入矩阵     
%1.用mvn数据
% Xk_NAME='D:\maven\OVERVIEW\mvn_insitu_kp-4sec_20171115_v15_r05.cdf';
% Xk_File=dataobj(Xk_NAME);
% tint1=irf.tint(Tsta,Tend);
% tint = irf_time(tint1,'epochtt>epoch');
% % tintS=tint1;
% 
% sp_MSO= get_variable(Xk_File,'SPICE_spacecraft_MSO');
% R_mso= getmat(Xk_File,'SPICE_spacecraft_MSO');
% R_mso=irf_tlim(R_mso,tint);
% R=R_mso;



%2.用MAG的1秒数据
dataBR=importdata('C:\Matlab\bin\新建文件夹\fwd\zcl\zcl\mvn_mag_l2_2017319ss1s_20171115_v01_r01.txt');
Tsmag='2017-11-15T00:00:01.307';
Temag='2017-11-16T00:00:00.390';
T1=iso2epoch(Tsmag);
T2=iso2epoch(Temag);
tint1=irf.tint(Tsta,Tend);
tint = irf_time(tint1,'epochtt>epoch');
TT = linspace(T1,T2,length(dataBR(:,1)))';
R_mso=[TT,dataBR(:,12),dataBR(:,13),dataBR(:,14)];% km
R_mso=irf_tlim(R_mso,tint);
R=R_mso;
loca=espectrum(:,1).*judge    %判据判断发生两带现象的时间

location=zeros(length(loca),3);
for z=1:length(loca)
   
    if isnan(loca(z,1))==1
        location(z,:)=NaN;
       
   else 
   t_dif3=abs(loca(z,1)-R(:,1));
   At3=find(min(min(t_dif3))==t_dif3);
   location(z,:)=[R(At3,2),R(At3,3),R(At3,4)];
   end
end     %发生两带的空间位置
location=location./3396
[xt,yt,zt]=sphere(25);
surf(xt,yt,zt)
hold on 
Xm=location(:,1)
Ym=location(:,2)
Zm=location(:,3)
%  plot3(Xm,Ym,Zm,'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);
 plot3(Xm,Ym,Zm,'LineWidth',1.5,'color',[0.9290 0.6940 0.1250]);  %画出判断出来的现象空间位置
ylabel({'Y_{MSO}'},'fontsize',8);
xlabel({'X_{MSO}'},'fontsize',8);
zlabel({'Z_{MSO}'},'fontsize',8)

% %% 模型的计算
% xshock0=[1.6285:-0.1:-4];
% % zshock0=[1.6285:-0.1:-4];
% a0=1.026;b=0.6;L=2.081;
% y0=((a0^2)-1).*((xshock0-b).^2)-2.*a0.*L.*(xshock0-b)+L^2;%方程
% y0=y0.^(0.5);
% y0=real(y0);
% 
% %磁层边界层
% x1=[0:0.0001:1.253];
% L1=1.08;a1=0.77;b1=0.64;
% y1=((a1^2)-1)*((x1-b1).^2)-2*a1*L1*(x1-b1)+L1^2;
% y1=real(y1.^(0.5));
% 
% x2=[0:-0.0001:-4];
% L2=0.528;a2=1.009;b2=1.6;
% y2=((a2^2)-1)*((x2-b2).^2)-2*a2*L2*(x2-b2)+L2^2;
% y2=real(y2.^(0.5));
% 
% x12=[fliplr(x1) x2];
% y12=[fliplr(y1) y2];
% 
% 
% %% plor 3D
% X_Y12=[x12' y12'];%磁层边界层
% 
% X_Y=[xshock0' y0'];%弓激波
% 
% angle=[0:20:360]./180.*pi;
%   
% for i_n=1:length(angle)
%     X_Y_Z{1,i_n}=[X_Y(:,1) X_Y(:,2).*cos(angle(i_n)) X_Y(:,2).*sin(angle(i_n))];
%     X_Y_Z12{1,i_n}=[X_Y12(:,1) X_Y12(:,2).*cos(angle(i_n)) X_Y12(:,2).*sin(angle(i_n))];
% end
% 
% 
% %% 画每一个圆
% nn_cc=5;
% number_c=floor(length(X_Y_Z{1,1}(:,1))/nn_cc);
% angle_c=([0:1:360]./180).*pi;
% for i_c=1:number_c
%     X_Y_Z_c{1,i_c}=X_Y_Z{1,1}(i_c*nn_cc,:);
%     x_c(1:length(angle_c),1)=X_Y_Z_c{1,i_c}(1,1);
%     for i_angle_c=1:length(angle_c)
%         y_c=(X_Y_Z_c{1,i_c}(1,2).*cos(angle_c))';
%         z_c=(X_Y_Z_c{1,i_c}(1,2).*sin(angle_c))';
%         X_Y_Z_c1{1,i_c}=[x_c y_c z_c];      
%     end 
% end
% for i_c=1:number_c
%     X_Y_Z_c12{1,i_c}=X_Y_Z12{1,1}(i_c*nn_cc,:);
%     x_c(1:length(angle_c),1)=X_Y_Z_c12{1,i_c}(1,1);
%     for i_angle_c=1:length(angle_c)
%         y_c=(X_Y_Z_c12{1,i_c}(1,2).*cos(angle_c))';
%         z_c=(X_Y_Z_c12{1,i_c}(1,2).*sin(angle_c))';
%         X_Y_Z_c112{1,i_c}=[x_c y_c z_c];      
%     end 
% end
% %% 画每一条曲线
% for i_nn=1:length(angle)   
% plot3(X_Y_Z{1,i_nn}(:,1),X_Y_Z{1,i_nn}(:,2),X_Y_Z{1,i_nn}(:,3),'-','LineWidth',1,'color',[127/255 255/255 0/255]);
% hold on
% end
% hold on
% for i_nc=1:number_c
% plot3(X_Y_Z_c1{1,i_nc}(:,1),X_Y_Z_c1{1,i_nc}(:,2),X_Y_Z_c1{1,i_nc}(:,3),'-','LineWidth',1,'color',[127/255 255/255 0/255]);
% hold on
% end
% hold on
% for i_nn=1:length(angle)   
% plot3(X_Y_Z12{1,i_nn}(:,1),X_Y_Z12{1,i_nn}(:,2),X_Y_Z12{1,i_nn}(:,3),'-','LineWidth',1,'color','k');
% hold on
% end
% hold on
% for i_nc=1:number_c
% plot3(X_Y_Z_c112{1,i_nc}(:,1),X_Y_Z_c112{1,i_nc}(:,2),X_Y_Z_c112{1,i_nc}(:,3),'-','LineWidth',1,'color','k');
% hold on
% end
% [xt,yt,zt]=sphere(25);
% mesh(xt,yt,zt);
% alpha(0.5)
% axis equal
% axis on
% 
% 
% set(gca,'XDir','reverse'); 
%% 
% fn=figure;
% set(fn,'Position',[200 200 800 400])
% %     h(1)=axes('position',[0.08 0.12 0.3 0.80]); 
% %     h(1)=axes('position',[0.30 0.10 0.30 0.80]); 
% h(1)=axes('position',[0.10 0.10 0.80 0.80]); 
%     ud=get(fn,'userdata');
%     ud.subplot_handles=h;
%     set(fn,'userdata',ud);
%     set(fn,'defaultLineLineWidth',2); 
% % plot(h(1),log10(Energy),smooth(log10(YY_PSD1),5),'-','LineWidth',0.8,'color','b');
% plot(h(1),(Energy_lg),gradient,'-','LineWidth',0.8,'color','b');
% % plot(h(1),(Energy),YY_PSD1,'-','LineWidth',0.8,'color','b');
% hold(h(1),'on');
% % ylabel(h(1),'log 10(f_e (s^3 km^{-6}))');
% ylabel(h(1),'斜率k');
% %ylabel(h(1),'log(f_e) (s^3 km^{-6})');
% xlabel(h(1),'log(Ee) (eV)')
% set(h(1),'yscale');
% set(h(1),'xscale');
% % set(h(1),'xlim',[1 6000]);
% % set(h(1),'ylim',[0.0001 100000000])
% set(h(1),'xtick',[1 2 3])
% irf_legend(gca,{TIME1(12:19)},[0.1 0.12],'fontsize',15);
% grid on