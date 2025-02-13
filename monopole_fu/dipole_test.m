clear;clc;close all
%------written by Wending Fu, Jan.2024 in Beijing------------
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
%% Parameter (two monopoles)
units = irf_units;
k0 = 1/(4*pi);
Q1 = 1e6;Q2 = -1e6;
r1 = [0,0,0.1];r2 = [0,0,-0.1];
x = -50:0.5:50;x=x';

a1 = 10*[0,0,3]; a2 = 10*[-2*sqrt(2)-1,2,-1]; a3 = 10*[sqrt(2)+2,sqrt(6)+1,-2]; a4 = 10*[sqrt(2)-1,-sqrt(6)+1,-3];
c_eval('a? = x.*[1,0,0] + ones(length(x),1).*a?;',1:4);
temp = 1:length(a1); temp = temp';
%% calculate B (two monopoles)
c_eval('d1? = sqrt((a?(:,1)-r1(1)).^2+(a?(:,2)-r1(2)).^2+(a?(:,3)-r1(3)).^2);',1:4);
c_eval('d2? = sqrt((a?(:,1)-r2(1)).^2+(a?(:,2)-r2(2)).^2+(a?(:,3)-r2(3)).^2);',1:4);
c_eval("Bt1? = transpose(k0*Q1./(d1?'.^2));",1:4);
c_eval("Bt2? = transpose(k0*Q2./(d2?'.^2));",1:4);
c_eval('Bt1?(d1?<=0,1)=0;',1:4);c_eval('Bt2?(d2?<=0,1)=0;',1:4);
c_eval('B?(:,1) = Bt1?.*((a?(:,1)-r1(1))./d1?) + Bt2?.*((a?(:,1)-r2(1))./d2?);',1:4);
c_eval('B?(:,2) = Bt1?.*((a?(:,2)-r1(2))./d1?) + Bt2?.*((a?(:,2)-r2(2))./d2?);',1:4);
c_eval('B?(:,3) = Bt1?.*((a?(:,3)-r1(3))./d1?) + Bt2?.*((a?(:,3)-r2(3))./d2?);',1:4);
c_eval('B? = B?*2.5;',1:4);
c_eval('Bt? = irf_abs(B?);',1:4);
c_eval('Bt? = Bt?(:,4);',1:4);
c_eval("a? = [temp,a?];")
c_eval('B? = [temp,B?];',1:4);
c_eval('Bt? = [temp,Bt?];',1:4);
%% Parameter (dipole moment)
% % % units = irf_units;
% % % k0 = 1/(4*pi);
% % % m = 1e5;
% % % r1 = [0,0,0];
% % % x = -50:0.5:50;x=x';
% % % 
% % % a1 = 10*[0,0,3]; a2 = 10*[-2*sqrt(2)-1,2,-1]; a3 = 10*[sqrt(2)+2,sqrt(6)+1,-2]; a4 = 10*[sqrt(2)-1,-sqrt(6)+1,-3];
% % % c_eval('a? = x.*[1,0,0] + ones(length(x),1).*a?;',1:4);
% % % temp = 1:length(a1); temp = temp';
%% calculate B (dipole moment)

%% PI
% PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
% PI=c_4_poincare_index(B1(:,2:4),B4(:,2:4),B3(:,2:4),B2(:,2:4));
% [PI,~,~,~,~]=poincare_index_and_solid_angle(B1_gsm(2:4),B2_gsm(2:4),B3_gsm(2:4),B4_gsm(2:4));
PI(abs(PI)<0.5) = 0;
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
%% divB
gradB=c_4_grad('a?(:,2:4)','B?(:,2:4)','grad');
divB=[gradB(:,1) sum([gradB(:,1) gradB(:,5) gradB(:,9)],2)];      %% 未归一化散度
PI(isnan(divB(:,2)))=0;
% tint = [iso2epoch(Tsta2) iso2epoch(Tend2)];
%% FOTE 误差
LocPoint = zeros(length(PI),3)*nan;
LocRes = cell(length(PI),1);
Q = zeros(length(PI),1)*nan;
resQ = cell(length(PI),1);

Qerror = ones(length(PI),1)*1000;
Locerror = ones(length(PI),1)*200;

c_eval("Btest? = [temp,B?(:,2:4)];")

PI_id = find(PI~=0)';
for i = 1:length(B1)
% for i = PI_id
clc;
disp(['current calculate:',num2str(i),'/',num2str(length(PI))]);

RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[a',num2str(ii),'(i,2),a',num2str(ii),'(i,3),a',num2str(ii),'(i,4);',...
    'a?(i,2),a?(i,3),a?(i,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
% MultiPower = ceil(max([log10(abs(R1(i,2:4))),log10(abs(R2(i,2:4))),log10(abs(R3(i,2:4))),log10(abs(R4(i,2:4)))]));
% solve

[Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('a?','Btest?',i,i*sign(divB(i,2)),RR_mean,1);
% [Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('a?','Btest?',i,i,RR_mean,1);

if ~isnan(LocRes{i})
id = nchoosek(1:6,2);
c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
tempd = [];
c_eval('tempd = [tempd,tempd?(4)/RR_mean];',1:15);
dLoc(i,:) = tempd;
else
    dLoc(i,:) = 10*ones(1,15);
end

%error
if ~isnan(resQ{i})
    Qerror(i) = abs(100*std(resQ{i})/Q(i));
else
    Qerror(i) = 1000;
end
% % % tri_a = delaunayTriangulation([a1(i,2:4);a2(i,2:4);a3(i,2:4);a4(i,2:4)]);
% % % [~,volume_a] = convexHull(tri_a);
% % % tri = delaunayTriangulation(LocRes{i});%%delaunay三角剖分
% % %     if size(tri.Points,1)==1
% % %         volume = 0;
% % %     else
% % %     [~,volume] = convexHull(tri);%%计算多面体体积
% % %     end
% % % Locerror(i) = 100*volume/volume_a;
end
meand = mean(dLoc,2);
meand(meand>10) = 10;
Qerror(meand==10) = 1e3;
%% Coordinate
% c_eval('a?(:,2:4) = a?(:,2:4)-LocPoint;');
% % NameTags = cellfun(@(x)(str2double(x(2:end-2))),unique(cellfun(@cellstr,NameTags)),'UniformOutput',false);
% % cellfun(@(x)(x-LocPoint))
% % LocRes{tempidx_B1} = LocRes{tempidx_B1}-LocPoint(tempidx_B1,:);
% LocPoint = LocPoint - LocPoint;
%% div error
% for i = 1:length(a1)
% c_eval('[tempBx?,tempBy?,tempBz?,tempBt?]=calculateB(a?(i,2),a?(i,3),a?(i,4),Q(i),LocPoint(i,:));')
% c_eval("Bm?(i,:) = [i,tempBx?,tempBy?,tempBz?,tempBt?];")
% end
% gradBm=c_4_grad('a?(:,2:4)','Bm?(:,2:4)','grad');
% divBm=[gradBm(:,1) sum([gradBm(:,1) gradBm(:,5) gradBm(:,9)],2)]; 
% divErr = [divB(:,1) 100*abs(divB(:,2)-divBm(:,2))./abs(divBm(:,2))]; %|deltB-deltBm|/deltB
%% B error
% c_eval('Bcross? = irf_dot(B?(:,2:4),Bm?(:,2:4));')
% c_eval('Btheta? = asin(sqrt(1-(Bcross?./(B?(:,1).*Bm?(:,5))).^2))*180/pi;')
% c_eval('Blength? = abs(B?(:,1)-Bm?(:,5))./Bm?(:,5);')
% Blength = 0.25*(Blength1+Blength2+Blength3+Blength4);
% c_eval('BErr? = 100*(1-(1-Btheta?).*(1-Blength?));')
% BErr = mean([BErr1,BErr2,BErr3,BErr4],2);
%% movie1
% % % figname = ['C:\Users\fwd\Desktop\Ti~mor~\M\magnetic_monopole\supplementary\illustration\SlideModel(S1)\SlideModel-paras.mp4'];
% % % v = VideoWriter(figname, 'MPEG-4');
% % % v.FrameRate = 60;
% % % v.Quality = 100;
% % % open(v)
% % % temp = x;
nn = length(temp);
% % % for f = 1:1:nn
f = length(temp);
%% plot
n=7;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 60; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Btotal
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([temp(1:f) B1(1:f,1)], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([temp(1:f) B2(1:f,1)], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([temp(1:f) B3(1:f,1)], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([temp(1:f) B4(1:f,1)], 'color','b', 'Linewidth',0.75); hold on;
% % % if f ~= nn
% % % irf_plot([temp(f) B1(f,1)], 'ow','MarkerFaceColor','k'); hold on;
% % % irf_plot([temp(f) B2(f,1)], 'ow','MarkerFaceColor','r'); hold on;
% % % irf_plot([temp(f) B3(f,1)], 'ow','MarkerFaceColor','g'); hold on;
% % % irf_plot([temp(f) B4(f,1)], 'ow','MarkerFaceColor','b'); hold on;
% % % end
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([temp(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
% % % grid off;
% % % set(gca,'Ylim',[min([min(B1(:,1)) min(B2(:,1)) min(B3(:,1)) min(B4(:,1))])-2 ...
% % %     max([max(B1(:,1)) max(B2(:,1)) max(B3(:,1)) max(B4(:,1))])+2]);
% % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% % % % set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
% % % set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
% % % irf_legend(gca,{'B_1','B_2','B_3','B_4'},[0.97 0.92]);
% % % irf_legend(gca,{'b'},[0.03 0.92],'color','k');
% % % ylabel('|B| [nT]','fontsize',12);
% % % set(gca,'XtickLabel',[])
% % % i=i+1;
%% Bx
h(i)=irf_subplot(n,1,-i);
irf_plot([temp(1:f) B1(1:f,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B2(1:f,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B3(1:f,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B4(1:f,2)], 'color','b', 'Linewidth',0.75); hold on;
if f ~= nn
irf_plot([temp(f) B1(f,2)], 'ow','MarkerFaceColor','k'); hold on;
irf_plot([temp(f) B2(f,2)], 'ow','MarkerFaceColor','r'); hold on;
irf_plot([temp(f) B3(f,2)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([temp(f) B4(f,2)], 'ow','MarkerFaceColor','b'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,2)) min(B2(:,2)) min(B3(:,2)) min(B4(:,2))])-2 ...
    max([max(B1(:,2)) max(B2(:,2)) max(B3(:,2)) max(B4(:,2))])+2]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'b'},[0.03 0.92],'color','k');
irf_legend(gca,{'B_1','B_2','B_3','B_4'},[0.97 0.92]);
ylabel('B_x [nT]','fontsize',12);
set(gca,'XtickLabel',[])
i=i+1;
%% By
h(i)=irf_subplot(n,1,-i);
irf_plot([temp(1:f) B1(1:f,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B2(1:f,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B3(1:f,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B4(1:f,3)], 'color','b', 'Linewidth',0.75); hold on;
if f ~= nn
irf_plot([temp(f) B1(f,3)], 'ow','MarkerFaceColor','k'); hold on;
irf_plot([temp(f) B2(f,3)], 'ow','MarkerFaceColor','r'); hold on;
irf_plot([temp(f) B3(f,3)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([temp(f) B4(f,3)], 'ow','MarkerFaceColor','b'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,3)) min(B2(:,3)) min(B3(:,3)) min(B4(:,3))])-2 ...
    max([max(B1(:,3)) max(B2(:,3)) max(B3(:,3)) max(B4(:,3))])+2]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'c'},[0.03 0.92],'color','k');
% irf_legend(gca,{'B_1','B_2','B_3','B_4'},[0.97 0.92]);
ylabel('B_y [nT]','fontsize',12);
set(gca,'XtickLabel',[])
i=i+1;
%% Bz
h(i)=irf_subplot(n,1,-i);
irf_plot([temp(1:f) B1(1:f,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B2(1:f,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B3(1:f,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) B4(1:f,4)], 'color','b', 'Linewidth',0.75); hold on;
if f ~= nn
irf_plot([temp(f) B1(f,4)], 'ow','MarkerFaceColor','k'); hold on;
irf_plot([temp(f) B2(f,4)], 'ow','MarkerFaceColor','r'); hold on;
irf_plot([temp(f) B3(f,4)], 'ow','MarkerFaceColor','g'); hold on;
irf_plot([temp(f) B4(f,4)], 'ow','MarkerFaceColor','b'); hold on;
end
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([temp(1:f) 0*B1(1:f,2)],'k--', 'Linewidth',0.75); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,4)) min(B2(:,4)) min(B3(:,4)) min(B4(:,4))])-2 ...
    max([max(B1(:,4)) max(B2(:,4)) max(B3(:,4)) max(B4(:,4))])+2]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
set(gca,'ColorOrder',[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]]);
irf_legend(gca,{'d'},[0.03 0.92],'color','k');
% irf_legend(gca,{'B_1','B_2','B_3','B_4'},[0.97 0.92]);
ylabel('B_z [nT]','fontsize',12);
set(gca,'XtickLabel',[])
i=i+1;
%% divB
h(i)=irf_subplot(n,1,-i);
irf_plot([temp(1:f) divB(1:f,2)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= nn
irf_plot([temp(f) divB(f,2)], 'ow','MarkerFaceColor','k'); hold on;
end
irf_plot([temp(1:f) 0*divB(1:f,2)],'k--', 'Linewidth',0.75); hold off;

grid off;
set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);
% set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
irf_legend(gca,{'e'},[0.03 0.92],'color','k');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
ylabel('\nabla·B [nT/km]','fontsize',12);
set(gca,'XtickLabel',[])
i=i+1;
%% PI
h(i)=irf_subplot(n,1,-i);
irf_plot([temp(1:f) PI(1:f)], 'color','k', 'Linewidth',0.75); hold on;
if f ~= nn
irf_plot([temp(f) PI(f)], 'ow','MarkerFaceColor','k'); hold on;
end

% c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
irf_legend(gca,{'f'},[0.03 0.92],'color','k');
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('PI ','fontsize',12);
set(gca,'XtickLabel',[])
i=i+1;
%% Location Error
h(i)=irf_subplot(n,1,-i);

irf_plot([temp,100*meand], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Locerror)]);
% set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
% set(gca,'Xlim',[1,21])
pos1=get(gca,'pos');
irf_legend(gca,{'g'},[0.03 0.92],'color','k');
ylabel('mean','fontsize',12);
set(gca,'XtickLabel',[])
i=i+1;
%% Q Error
h(i)=irf_subplot(n,1,-i);
irf_plot([temp,Qerror], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Qerror)]);
% set(gca,'Ylim',[0,200], 'ytick',[0 100 200],'fontsize',9);
pos1=get(gca,'pos');
% set(gca,'Xlim',[1,21])
irf_legend(gca,{'h'},[0.03 0.92],'color','k');
ylabel('Q Err [%]','fontsize',12);
xlabel('Slide Distance [km]')
i=i+1;
%% Q
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) Q], 'color','k', 'Linewidth',0.75); hold on;
% % % grid off;
% % % set(gca,'Ylim',[0 max(Q)]);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 500 1000],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Q [nT·km^2]','fontsize',12);
% % % set(gca,'Yscale','log')
% % % i=i+1;
%% plot div err
% % % h(i)=subplot(n,1,i);
% % % plot(temp,divErr(:,2), 'color','k', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
% % % %     max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('div Error [%]','fontsize',12);
% % % set(gca,'Xlim',[1,21])
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% plot B theta err
% % % h(i)=subplot(n,1,i);
% % % plot(temp,Btheta1, 'color','k', 'Linewidth',0.75); hold on;
% % % plot(temp,Btheta2, 'color','r', 'Linewidth',0.75); hold on;
% % % plot(temp,Btheta3, 'color','g', 'Linewidth',0.75); hold on;
% % % plot(temp,Btheta4, 'color','b', 'Linewidth',0.75); hold on;
% % % 
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
% % % %     max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% % % set(gca,'Ylim',[0 90], 'ytick',[0 45 90],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% % % irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
% % % ylabel('B theta Error [°]','fontsize',12);
% % % set(gca,'Xlim',[1,21])
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% plot B length err
% % % h(i)=subplot(n,1,i);
% % % % % % irf_plot([B1(:,1) BErr], 'color','k', 'Linewidth',0.75); hold on;
% % % plot(temp,100*Blength1, 'color','k', 'Linewidth',0.75); hold on;
% % % plot(temp,100*Blength2, 'color','r', 'Linewidth',0.75); hold on;
% % % plot(temp,100*Blength3, 'color','g', 'Linewidth',0.75); hold on;
% % % plot(temp,100*Blength4, 'color','b', 'Linewidth',0.75); hold on;
% % % % irf_plot([B1(:,1) Blength], 'color','k', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
% % % %     max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % set(gca,'Xlim',[1,21])
% % % pos1=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% % % % irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
% % % ylabel('B length Error [%]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% plot B com err
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) BErr1], 'color','k', 'Linewidth',0.75); hold on;
% % % irf_plot([B2(:,1) BErr2], 'color','r', 'Linewidth',0.75); hold on;
% % % irf_plot([B3(:,1) BErr3], 'color','g', 'Linewidth',0.75); hold on;
% % % irf_plot([B4(:,1) BErr4], 'color','b', 'Linewidth',0.75); hold on;
% % % 
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
% % % %     max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% % % irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
% % % ylabel('B Error [%]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% plot mean B com err
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) BErr], 'color','k', 'Linewidth',0.75); hold on;
% % % % irf_plot([B1(:,1) Blength], 'color','k', 'Linewidth',0.75); hold on;
% % % %irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
% % % % c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
% % % %     max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % % set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% % % % irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
% % % ylabel('B Error [%]','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Annotation
irf_zoom(h(1:end),'x',[temp(1),temp(end)]);
irf_plot_axis_align;

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
irf_pl_mark(h(1:5),[temp(f)],[0.5,0.5,0.5],'linestyle','-');
set(gcf,'color','w')

% figname = [OutputDir,'OverviewFig\',NameTags{TDT}(2:end-2)];    
% print(gcf, '-dpng', [figname '.png']);
% % % frame = getframe(figure(1));
% % % writeVideo(v,frame)
% % % % print(gcf, '-dpdf', [figname '.pdf']);
% % % end
% % % close(v)

%% Init figure 2
% % % figname = ['C:\Users\fwd\Desktop\Ti~mor~\M\magnetic_monopole\supplementary\illustration\SlideModel(S1)\SlideModel-configuration.mp4'];
% % % v = VideoWriter(figname, 'MPEG-4');
% % % v.FrameRate = 60; % 默认 30
% % % v.Quality = 100; % 默认 75
% % % open(v)

% % % for i = 1:1:nni = 101;
if ~isnan(LocRes{i})
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
c_eval("line(a?(:,2),a?(:,3),a?(:,4),'LineStyle','--','color','#898989');hold on;")

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
% irf_legend(gca,{['\theta:',num2str(theta(i)),'°']},[0.98 0.92],'fontsize',20);hold on;
xlabel('X [km]','fontsize',15);
ylabel('Y [km]','fontsize',15);
zlabel('Z [km]','fontsize',15);
%% Quiver
Rmean = [0,0,0];
c_eval('Bt? = irf_abs(B?);Bt? = Bt?(:,5);')
% Bt5 = max([Bt1(i),Bt2(i),Bt3(i),Bt4(i)]);
Bt5 = 5;
c_eval('Rmean = Rmean+ a?(i,2:4);');Rmean = Rmean/4;
% for i = 1:10:359
c_eval("quiver3(a?(i,2),a?(i,3),a?(i,4),RR_mean*B?(i,2)/Bt5"+...
",RR_mean*B?(i,3)/Bt5,RR_mean*B?(i,4)/Bt5,2 ,'color',cor{?},'linewidth',3,'MaxHeadSize',1.5); hold on;")
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
if i <=floor(nn/2)
legend([PointZero],{'theoretical magnetic monopole'},'Position',[0.085,0.73,0.4,0.075],'FontSize',12);hold on;
irf_legend(gca,{'a'},[0.98 1.05],'fontsize',15);hold on;
else
legend([PointZero],{'theoretical magnetic monopole'},'Position',[0.53,0.73,0.4,0.075],'FontSize',12);hold on;
irf_legend(gca,{'a'},[0.03 1.05],'fontsize',15);hold on;
end
legend('boxoff')
%% Point in Tetrahedron
% % % LocTetra = [a1(i,2:4);a2(i,2:4);a3(i,2:4);a4(i,2:4)];
% % % if PointInTetrahedron([0,0,0],a1(i,2:4),a2(i,2:4),a3(i,2:4),a4(i,2:4))
% % %     plotPolyhedron(LocTetra(:,1),LocTetra(:,2),LocTetra(:,3),'#FFB8CE',0.2);hold on;
% % %     legend({'Inside'},'Position',[0.085,0.76,0.4,0.075],'FontSize',12);hold on;
% % % else
% % %     plotPolyhedron(LocTetra(:,1),LocTetra(:,2),LocTetra(:,3),'#CCDDFF',0.2);hold on;
% % %     legend({'Outside'},'Position',[0.53,0.76,0.4,0.075],'FontSize',12);hold on;
% % % end
% % % legend('boxoff')
%% output figure 2
RR_max = RR_mean;
axis equal 
axis([-120,120,-40,60,-60,60]);
view([-20+i*0.2,35])

box on;
set(gca,'LineWidth',0.75);
set(gca,'ZTick',[-50,0,50]);
set(gca,'FontSize',12);
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
set(gcf,'color','w')

% % % frame = getframe(figure(2));
% % % writeVideo(v,frame)
% print(gcf, '-dpdf', [figname '.pdf']);
end
% % % end
% % % close(v);
%%
% % % %% streamline
% % % linestep = 20;
% % % for i = ['X','Y','Z']
% % % eval(['grid',i,' = grid',i,'(1:linestep:end,1:linestep:end,1:linestep:end);'])
% % % eval(['gridB',i,' = gridB',i,'(1:linestep:end,1:linestep:end,1:linestep:end);'])
% % % end
% % % lines = streamline(stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,gridX,gridY,gridZ));
% % % set(lines,'color',[0.8,0.58,1]);
% % % for i = 1:length(lines)
% % %     xLine = lines(i).XData;
% % %     yLine = lines(i).YData;
% % %     zLine = lines(i).ZData;
% % %     if length(xLine) >= 3
% % %         LineVec = [xLine(end)-xLine(end-1),yLine(end)-yLine(end-1),zLine(end)-zLine(end-1)];
% % %         NormPara = sqrt(sum(LineVec.^2,2)/sum((Rmax-Rmin).^2,2));
% % % %         LineVec = LineVec*NormPara;
% % %         quiver3(xLine(end),yLine(end),zLine(end),LineVec(1),LineVec(2),LineVec(3),0,'MaxHeadsize',5,'LineWidth',1,'color',[0.8,0.58,1]); hold on;
% % %     end
% % % end
% % % view(3)
%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % idx = find(indices==1);
% % % n=8;
% % % i=1;
% % % set(0,'DefaultAxesFontSize',8);
% % % set(0,'DefaultLineLineWidth', 0.5);
% % % fn=figure(2);clf;
% % % set(gcf,'PaperUnits','centimeters')
% % % xSize = 6; ySize = 9; coef=floor(min(800/xSize,800/ySize));
% % % xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% % % set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% % % set(gcf,'Position',[10 10 xSize*coef ySize*coef])
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,divE(:,2), 'Linewidth',0.75); hold on;
% % % c_eval("plot(x,0*divE(:,2),'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(divE(:,2))-0.1*abs(min(divE(:,2))) max(divE(:,2))+0.1*max(divE(:,2))]);");
% % % pos1=get(gca,'pos');
% % % ylabel('divE [V/m^2]','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,indices, 'Linewidth',0.75); hold on;
% % % grid off;
% % % set(gca, 'ytick',[-1 0 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('PI','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,Ang12, 'Linewidth',0.75); hold on;
% % % plot(x(idx),Ang12(idx), 'color','r', 'Linewidth',0.75)
% % % c_eval("plot(x,0*Ang12,'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(Ang12)-1 max(Ang12)+1]);");
% % % pos1=get(gca,'pos');
% % % ylabel('Ang12 [°]','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,Ang13, 'Linewidth',0.75); hold on;
% % % plot(x(idx),Ang13(idx), 'color','r', 'Linewidth',0.75)
% % % c_eval("plot(x,0*Ang13,'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(Ang13)-1 max(Ang13)+1]);");
% % % pos1=get(gca,'pos');
% % % ylabel('Ang13 [°]','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,Ang14, 'Linewidth',0.75); hold on;
% % % plot(x(idx),Ang14(idx), 'color','r', 'Linewidth',0.75)
% % % c_eval("plot(x,0*Ang14,'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(Ang14)-1 max(Ang14)+1]);");
% % % pos1=get(gca,'pos');
% % % ylabel('Ang14 [°]','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,Ang23, 'Linewidth',0.75); hold on;
% % % plot(x(idx),Ang23(idx), 'color','r', 'Linewidth',0.75)
% % % c_eval("plot(x,0*Ang23,'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(Ang23)-1 max(Ang23)+1]);");
% % % pos1=get(gca,'pos');
% % % ylabel('Ang23 [°]','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,Ang24, 'Linewidth',0.75); hold on;
% % % plot(x(idx),Ang24(idx), 'color','r', 'Linewidth',0.75)
% % % c_eval("plot(x,0*Ang24,'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(Ang24)-1 max(Ang24)+1]);");
% % % pos1=get(gca,'pos');
% % % ylabel('Ang24 [°]','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,Ang34, 'Linewidth',0.75); hold on;
% % % plot(x(idx),Ang34(idx), 'color','r', 'Linewidth',0.75)
% % % c_eval("plot(x,0*Ang34,'k--', 'Linewidth',0.75);"); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(Ang34)-1 max(Ang34)+1]);");
% % % pos1=get(gca,'pos');
% % % ylabel('Ang34 [°]','fontsize',10);
% % % xlabel('x [m]','fontsize',10);
% % % i = i+1;
% % % set(gcf,'render','painters');
% % % set(gcf,'paperpositionmode','auto')
% % % %%%
% % % % % subplot(2,1,1)
% % % % % plot(x,divE(:,2))
% % % % % ylabel('divE','fontsize',12);
% % % % % subplot(2,1,2)
% % % % % plot(x,indices)
% % % % % xlabel('X','fontsize',12);
% % % % % ylabel('PI','fontsize',12);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3
figure(3)
subplot(3,1,1)
irf_plot([x B1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([x B2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([x B3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([x B4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
set(gca,'Ylim',[min([min(B1(:,2)) min(B2(:,2)) min(B3(:,2)) min(B4(:,2))]) ...
    max([max(B1(:,2)) max(B2(:,2)) max(B3(:,2)) max(B4(:,2))])]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bx [V/m]','fontsize',12);

subplot(3,1,2)
irf_plot([x B1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([x B2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([x B3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([x B4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
set(gca,'Ylim',[min([min(B1(:,3)) min(B2(:,3)) min(B3(:,3)) min(B4(:,3))]) ...
    max([max(B1(:,3)) max(B2(:,3)) max(B3(:,3)) max(B4(:,3))])]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('By [V/m]','fontsize',12);

subplot(3,1,3)
irf_plot([x B1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([x B2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([x B3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([x B4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
set(gca,'Ylim',[min([min(B1(:,4)) min(B2(:,4)) min(B3(:,4)) min(B4(:,4))]) ...
    max([max(B1(:,4)) max(B2(:,4)) max(B3(:,4)) max(B4(:,4))])]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bz [V/m]','fontsize',12);