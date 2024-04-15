clear;clc;close all
%------written by Wending Fu, Jul.8.2022 in Beijing------------
%------modified by Wending Fu, Nov.2023 in Beijing for dipole------------
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       南无电子阿弥陀佛驱散仿生bug
%                                _ooOoo_
%                               o8888888o
%                               88" . "88
%                               (| -_- |)
%                               O\  =  /O
%                            ____/`---'\____
%                          .'  \\|     |//  `.
%                         /  \\|||  :  |||//  \
%                        /  _||||| -:- |||||-  \
%                        |   | \\\  -  /// |   |
%                        | \_|  ''\-/''  |   |
%                        \  .-\__  `-`  ___/-. /
%                      ___`. .'  /-.-\  `. . __
%                   ."" '<  `.___\_<|>_/___.'  >'"".
%                  | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%                  \  \ `-.   \_ __\ /__ _/   .-` /  /
%             ======`-.____`-.___\_____/___.-`____.-'======
% 	                   `=-='
%                 天地玄宗，万气本根。广修亿劫，证吾神通。
%                 三界内外，惟道独尊。体有金光，覆映吾身。
%                 视之不见，听之不闻。包罗天地，养育群生。
%                 受持万遍，身有光明。三界侍卫，五帝司迎。
%                 万神朝礼，役使雷霆。鬼妖丧胆，精怪忘形。
%                 内有霹雳，雷神隐名。洞慧交彻，五炁腾腾。
%                金光速现，覆护真人。急急如律令，bug全去除！
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter
units = irf_units;
k0 = 1/(4*pi);
Q1 = 1e5;Q2 = 0;
r1 = [-1,-1,-1];r2 = [-1,0,5];
% a1 = 10*[0,0,sqrt(3)]; a2 = 10*[1, sqrt(2), 0]; a3 = 10*[0, -sqrt(2), -1]; a4 = 10*[-1, sqrt(2), 0];
a1 = 10*[0,0,3]; a2 = 10*[-2*sqrt(2)-1,2,-1]; a3 = 10*[sqrt(2)+2,sqrt(6)+1,-2]; a4 = 10*[sqrt(2)-1,-sqrt(6)+1,-3];

%% Calculate B
c_eval('d1? = sqrt((a?(:,1)-r1(1)).^2+(a?(:,2)-r1(2)).^2+(a?(:,3)-r1(3)).^2);',1:4);
c_eval('d2? = sqrt((a?(:,1)-r2(1)).^2+(a?(:,2)-r2(2)).^2+(a?(:,3)-r2(3)).^2);',1:4);
c_eval("Bt1? = transpose(k0*Q1./(d1?'.^2));",1:4);
c_eval("Bt2? = transpose(k0*Q2./(d2?'.^2));",1:4);
c_eval('Bt1?(d1?<=0,1)=0;',1:4);c_eval('Bt2?(d2?<=0,1)=0;',1:4);
c_eval('B?(:,1) = Bt1?.*((a?(:,1)-r1(1))./d1?) + Bt2?.*((a?(:,1)-r2(1))./d2?);',1:4);
c_eval('B?(:,2) = Bt1?.*((a?(:,2)-r1(2))./d1?) + Bt2?.*((a?(:,2)-r2(2))./d2?);',1:4);
c_eval('B?(:,3) = Bt1?.*((a?(:,3)-r1(3))./d1?) + Bt2?.*((a?(:,3)-r2(3))./d2?);',1:4);
c_eval('Bt? = irf_abs(B?);',1:4);
c_eval('Bt? = Bt?(:,4);',1:4);
%% Rot B
theta = transpose([-180:1:180]);
phi = 135;
B1 = zeros(length(theta),3);
c_eval('B? = repmat(B?,length(theta),1);',2:4);
c_eval('a? = repmat(a?,length(theta),1);');
% n = a2(1,:)/30;n = [-n(3),0,n(1)];
for i = 1:length(theta)
    B1(i,:) = [Bt1*sind(theta(i))*sind(phi),Bt1*sind(theta(i))*cosd(phi),Bt1*cosd(theta(i))];
    %向量v绕向量u旋转θ，v'=vcosθ + (u×v)sinθ + (u·v)u(1-cosθ) 
    % tempB2 = B2(i,:)/Bt2;
    % B2(i,:) = Bt2*(tempB2*cosd(theta(i))+cross(tempB2,n)*sind(theta(i))+dot(n,tempB2)*n*(1-cosd(theta(i))));
end

%% Poincare Index
indices=c_fgm_poincare_index(B1,B2,B3,B4);
indices(abs(indices)<0.5) = 0;
indices(indices>=0.5) = 1;
indices(indices<=-0.5) = -1;

%% div & FOTE err
gradB=c_4_grad('a?','B?','grad');
% gradB=c_4_grad(a3,a2,a1,a4,B3,B2,B1,B4,'grad');
divB=sum([gradB(:,1) gradB(:,5) gradB(:,9)],2); %nT/km^2

% FOTE err
eigVal_err_v2 = B1(:,1);
monopole_index = zeros(length(theta),1);

for ii=2:length(B1(:,1))  
deltB_null=reshape(gradB(ii,1:end),3,3);
if ~isnan(deltB_null)
    [V,D] = eig(deltB_null);
    % Figure 1o    以最大特征值归一化 |Δ·B|/|ΔxB|
    eigVal_err_v2(ii)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]);  
    if indices(ii)~=0 && isreal(D(1,1)) && isreal(D(2,2)) && isreal(D(3,3)) && D(1,1) > 0 && D(2,2) > 0 && D(3,3) > 0
        monopole_index(ii) = 1;
    elseif indices(ii)~=0 && isreal(D(1,1)) && isreal(D(2,2)) && isreal(D(3,3)) && D(1,1) < 0 && D(2,2) < 0 && D(3,3) < 0
        monopole_index(ii) = -1;
    end
else
    eigVal_err_v2(ii,2)=nan;
end
end

[j,divB,~,~,~,~] = c_4_j('a?','B?');
temp=irf_abs(j);
jmag=temp(:,4);
err_4C=irf_multiply(1,divB,1,jmag,-1);  %% η
err_4C=abs(err_4C);    

%% solve monopole
c_eval('a? = [transpose(1:size(a?,1)),a?];',1:4);
c_eval('B? = [transpose(1:size(B?,1)),B?];',1:4);
[coeff, res] = VSH_Expand('a?', 'B?', 2);

% % % c_eval('a?(:,2:4) = a?(:,2:4)-r1;',1:4);
% % % [coeff, res] = VSH_Expand('a?', 'B?', [0,0,0], 1);

coef00 = coeff(:,[1,5,9]);
coef12 = coeff(:,[2,6,10]);
coef10 = coeff(:,[3,7,11]);
coef11 = coeff(:,[4,8,12]);
%% movie1
% % % figname = ['/Users/fwd/Documents/Ti~mor~/M/Monopole/VSH/VSH_test2_para.mp4'];
% % % v = VideoWriter(figname, 'MPEG-4');
% % % v.FrameRate = 60;
% % % v.Quality = 100;
% % % open(v)
% % % for f = 1:1:361
f = 361;
%% Init figure 1
figure(1)
n=5;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 120; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% B1
h(i)=irf_subplot(n,1,-i);
irf_plot([theta(1:f) B1(1:f,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) B1(1:f,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) B1(1:f,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([theta(1:f) B1(1:f,5)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= 361
irf_plot([theta(f) B1(f,2)], 'ow','MarkerFaceColor','b'); hold on;
irf_plot([theta(f) B1(f,3)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([theta(f) B1(f,4)], 'ow','MarkerFaceColor','r'); hold on;
% irf_plot([theta(f) B1(f,5)], 'ow','MarkerFaceColor','k'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,2)) min(B1(:,3)) min(B1(:,4))])-2 ...
    max([max(B1(:,2)) max(B1(:,3)) max(B1(:,4))])+2]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'B_x','B_y','B_z'},[0.97 0.92]);
irf_legend(gca,{'b'},[0.03 0.92],'color','k');
ylabel('B1 [nT]','fontsize',12);
i=i+1;
%% divB
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([theta(1:f) divB(1:f)], 'color','k', 'Linewidth',0.75); hold on;
% % % if f ~= 361
% % % irf_plot([theta(f) divB(f)], 'ow','MarkerFaceColor','k'); hold on;
% % % end
% % % irf_plot([theta(1:f) 0*divB(1:f)],'k--', 'Linewidth',0.75); hold off;
% % % 
% % % grid off;
% % % set(gca,'Ylim',[min(divB) max(divB)]);
% % % % set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % irf_legend(gca,{'c'},[0.03 0.92],'color','k');
% % % ylabel('\nabla·B [nT/km]','fontsize',12);
% % % i=i+1;
%% PI
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([theta(1:f) indices(1:f)], 'color','k', 'Linewidth',0.75); hold on;
% % % if f ~= 361
% % % irf_plot([theta(f) indices(f)], 'ow','MarkerFaceColor','k'); hold on;
% % % end
% % % 
% % % % c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
% % % set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % irf_legend(gca,{'d'},[0.03 0.92],'color','k');
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('PI ','fontsize',12);
% % % i=i+1;
%% Coeffcients 00
h(i)=irf_subplot(n,1,-i);
irf_plot([theta(1:f) coef00(1:f,1)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef00(1:f,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef00(1:f,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([theta(1:f) B1(1:f,5)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= 361
irf_plot([theta(f) coef00(f,1)], 'ow','MarkerFaceColor','b'); hold on;
irf_plot([theta(f) coef00(f,2)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([theta(f) coef00(f,3)], 'ow','MarkerFaceColor','r'); hold on;
% irf_plot([theta(f) B1(f,5)], 'ow','MarkerFaceColor','k'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% % % set(gca,'Ylim',[min([min(coef00(:,1)) min(coef00(:,2)) min(coef00(:,3))])-2 ...
% % %     max([max(coef00(:,1)) max(coef00(:,2)) max(coef00(:,3))])+2]);
set(gca,'Ylim',[min(coeff,[],'all')-2 , max(coeff,[],'all')+2]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\alpha','\beta','\gamma'},[0.97 0.92]);
irf_legend(gca,{'c'},[0.03 0.92],'color','k');
ylabel('[l,m]=[0,0]','fontsize',12);
i=i+1;
%% Coeffcients 12
h(i)=irf_subplot(n,1,-i);
irf_plot([theta(1:f) coef12(1:f,1)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef12(1:f,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef12(1:f,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([theta(1:f) B1(1:f,5)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= 361
irf_plot([theta(f) coef12(f,1)], 'ow','MarkerFaceColor','b'); hold on;
irf_plot([theta(f) coef12(f,2)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([theta(f) coef12(f,3)], 'ow','MarkerFaceColor','r'); hold on;
% irf_plot([theta(f) B1(f,5)], 'ow','MarkerFaceColor','k'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% % % set(gca,'Ylim',[min([min(coef12(:,1)) min(coef12(:,2)) min(coef12(:,3))])-2 ...
% % %     max([max(coef12(:,1)) max(coef12(:,2)) max(coef12(:,3))])+2]);
set(gca,'Ylim',[min(coeff,[],'all')-2 , max(coeff,[],'all')+2]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\alpha','\beta','\gamma'},[0.97 0.92]);
irf_legend(gca,{'d'},[0.03 0.92],'color','k');
ylabel('[l,m]=[1,-1]','fontsize',12);
i=i+1;
%% Coeffcients 10
h(i)=irf_subplot(n,1,-i);
irf_plot([theta(1:f) coef10(1:f,1)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef10(1:f,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef10(1:f,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([theta(1:f) B1(1:f,5)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= 361
irf_plot([theta(f) coef10(f,1)], 'ow','MarkerFaceColor','b'); hold on;
irf_plot([theta(f) coef10(f,2)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([theta(f) coef10(f,3)], 'ow','MarkerFaceColor','r'); hold on;
% irf_plot([theta(f) B1(f,5)], 'ow','MarkerFaceColor','k'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% % % set(gca,'Ylim',[min([min(coef10(:,1)) min(coef10(:,2)) min(coef10(:,3))])-2 ...
% % %     max([max(coef10(:,1)) max(coef10(:,2)) max(coef10(:,3))])+2]);
set(gca,'Ylim',[min(coeff,[],'all')-2 , max(coeff,[],'all')+2]);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\alpha','\beta','\gamma'},[0.97 0.92]);
irf_legend(gca,{'e'},[0.03 0.92],'color','k');
ylabel('[l,m]=[1,0]','fontsize',12);
i=i+1;
%% Coeffcients 11
h(i)=irf_subplot(n,1,-i);
irf_plot([theta(1:f) coef11(1:f,1)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef11(1:f,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) coef11(1:f,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([theta(1:f) B1(1:f,5)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= 361
irf_plot([theta(f) coef11(f,1)], 'ow','MarkerFaceColor','b'); hold on;
irf_plot([theta(f) coef11(f,2)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([theta(f) coef11(f,3)], 'ow','MarkerFaceColor','r'); hold on;
% irf_plot([theta(f) B1(f,5)], 'ow','MarkerFaceColor','k'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([theta(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
% % % set(gca,'Ylim',[min([min(coef11(:,1)) min(coef11(:,2)) min(coef11(:,3))])-2 ...
% % %     max([max(coef11(:,1)) max(coef11(:,2)) max(coef11(:,3))])+2]);
set(gca,'Ylim',[min(coeff,[],'all')-2 , max(coeff,[],'all')+2]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
set(gca,'ColorOrder',[[0 0 1];[0 1 0];[1 0 0]]);
irf_legend(gca,{'\alpha','\beta','\gamma'},[0.97 0.92]);
irf_legend(gca,{'f'},[0.03 0.92],'color','k');
ylabel('[l,m]=[1,1]','fontsize',12);
i=i+1;
%% Adjust the position
irf_zoom(h(1:end),'x',[theta(1),theta(end)]);
irf_plot_axis_align;
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
irf_pl_mark(h(1:5),[theta(f)],[0.5,0.5,0.5],'linestyle','-');
set(gcf,'color','w')

% figname = [OutputDir,'OverviewFig\',NameTags{TDT}(2:end-2)];    
% print(gcf, '-dpng', [figname '.png']);
% % % frame = getframe(figure(1));
% % % writeVideo(v,frame)
% % % % % % print(gcf, '-dpdf', [figname '.pdf']);
% % % end
% % % close(v)

%% Init figure 2
% % % figname = ['/Users/fwd/Documents/Ti~mor~/M/Monopole/VSH/VSH_test2_configuration.mp4'];
% % % v = VideoWriter(figname, 'MPEG-4');
% % % v.FrameRate = 60; % 默认 30
% % % v.Quality = 100; % 默认 75
% % % open(v)

% % % for i = 1:1:361
i = 181;
cla(figure(2))
figure(2)
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% S/C configuration
cor = {'k','r','g','b'};
c_eval("plot3(a?(i,2),a?(i,3),a?(i,4),'s' ,'color',cor{?},'linewidth',5); hold on;");

RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[a',num2str(ii),'(i,2),a',num2str(ii),'(i,3),a',num2str(ii),'(i,4);',...
    'a?(i,2),a?(i,3),a?(i,4)];'],ii+1:4);
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
c_eval(['plot3(RR',num2str(ii),'?(:,1),RR',num2str(ii),'?(:,2),RR',num2str(ii),'?(:,3),''-- '',''color'',[0.5,0.5,0.5]);hold on;'],ii+1:4)
end
RR_mean = RR_mean(4)/6;
% plot3(RR12(:,1),RR12(:,2),RR12(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  
% plot3(RR14(:,1),RR14(:,2),RR14(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  
% plot3(RR34(:,1),RR34(:,2),RR34(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
irf_legend(gca,{['\theta=',num2str(theta(i)),'°']},[0 0.92],'fontsize',20);hold on;
irf_legend(gca,{'a'},[0.98 1],'fontsize',15);hold on;
xlabel('X [km]','fontsize',15);
ylabel('Y [km]','fontsize',15);
zlabel('Z [km]','fontsize',15);
%% Quiver
Rmean = [0,0,0];
Bt5 = max([Bt1,Bt2,Bt3,Bt4]);

c_eval('Rmean = Rmean+ a?(i,2:4);');Rmean = Rmean/4;
% for i = 1:10:359
c_eval("quiver3(a?(i,2),a?(i,3),a?(i,4),RR_mean*B?(i,2)/Bt5"+...
",RR_mean*B?(i,3)/Bt5,RR_mean*B?(i,4)/Bt5,0.3 ,'color',cor{?},'linewidth',3,'MaxHeadSize',1.5); hold on;")
% end
%% Loc res
% try
% plotPolyhedron(LocRes{i}(:,1),LocRes{i}(:,2),LocRes{i}(:,3),'#FFB8CE',0.3);
% catch
% end
% c_eval("PointRes = plot3(LocRes{i}(?,1),LocRes{i}(?,2),LocRes{i}(?,3),'p','color','#FFB8CE');hold on;",1:6)
% text(LocPoint(1)+5,LocPoint(2),LocPoint(3),['Q = ',num2str(Q,3),'nT·km^2'])  
% c_eval("PointMonopole = plot3(LocPoint(i,1),LocPoint(i,2),LocPoint(i,3),'*' ,'color','m','linewidth',2); hold on;")
PointZero = plot3(0,0,0,'ow','MarkerFaceColor','#7E2F8E','MarkerSize',10);hold on;
% legend([PointRes,PointZero,PointMonopole],{'residual position of the magnetic monopole','modeled magnetic monopole',...
%     'solved magnetic monopole'},'Position',[0.085,0.82,0.4,0.075],'FontSize',8);hold on;
legend([PointZero],{'theoretical magnetic monopole'},'Position',[0.10,0.85,0.4,0.075],'FontSize',12);hold on;
legend('boxoff')
%% output figure 2
RR_max = RR_mean;
axis equal 
axis([-RR_max RR_max -RR_max RR_max -RR_max RR_max]);
view([-30,30])

box on;
set(gca,'LineWidth',0.75);
set(gca,'ZTick',[-50,0,50]);
set(gca,'FontSize',12);
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gcf,'color','w')

% % % frame = getframe(figure(2));
% % % writeVideo(v,frame)
% % % % print(gcf, '-dpdf', [figname '.pdf']);
% % % end
% % % close(v);

% print(gcf, '-dpng', [figname '.png']);