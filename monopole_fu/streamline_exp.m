clear;clc;close all
%------written by Wending Fu, Mar.09.2022 in Beijing------------
%The program simulates a 3d vector field with three hyperbolic zeros
%There is a saddle zero near the space zero
%Four-point detection is used to simulate the measurements of four satellites in the field
%Figure1 shows the running track of the four-point detector in space and the vector direction of the middle position
%Figure2 shows the divergence of the field measured by the four-point detector, poincare index, and the unipolar index under the linearization hypothesis
%Figure3 shows the vector measured by the four-point detector
%Citation: https://www.ljll.math.upmc.fr/~frey/cours/Roscoff/slidesAusoni2014_3.pdf
%see also div_exp
%% Parameter
[x,y,z] = meshgrid(-10:1:10);
u = -y;
v = x;
w = 0*z;
% u = 10*(y-x);
% v = 28*x-y-x.*z;
% w = x.*y-8*z/3;
% [sx,sy,sz] = meshgrid(0,-10:2:10,-10:2:10);
% quiver3(x,y,z,u,v,w)
quiver3(x,y,z,u,0*u,0*u)
% streamline(x,y,z,u,v,w,sx,sy,sz)
% view(3)
xlabel('e1','fontsize',12);
ylabel('e2','fontsize',12);
zlabel('e3','fontsize',12);
%%
a1 = [0,0,3]; a2 = [-2*sqrt(2),0,-1]; a3 = [sqrt(2),sqrt(6),-2]; a4 = [sqrt(2),-sqrt(6),-1];
xx = -10:1:10;xx = xx';
zz = -10:1:10;zz = zz';
c_eval('a? = xx.*[1,0,0] + ones(length(xx),1).*a?;',1:4);
% c_eval('a? = zz.*[0,0,1] + xx.*[1,0,0] + ones(length(xx),1).*a?;',1:4);
c_eval('E? = ones(length(xx),3);');
for i = 1:length(xx)
%     c_eval('[~,tempx?] = sort(abs(xx-a?(i,1)));');c_eval('tempx? = tempx?(1);');
%     c_eval('[~,tempy?] = sort(abs(xx-a?(i,2)));');c_eval('tempy? = tempy?(1);');
%     c_eval('[~,tempz?] = sort(abs(xx-a?(i,3)));');c_eval('tempz? = tempz?(1);');
%     c_eval('E?(i,:) = [u(tempx?,tempy?,tempz?),v(tempx?,tempy?,tempz?),w(tempx?,tempy?,tempz?)];')
c_eval('E?(i,:) = [-a?(i,2),a?(i,1),0];')
end
curlE=c_4_grad('a?','E?','curl');
c_eval('E? = E? - curlE;')
%% div
gradE=c_4_grad('a?','E?','grad');
units = irf_units;

divE=[gradE(:,1) sum([gradE(:,1) gradE(:,5) gradE(:,9)],2)];      %% 未归一化散度
% divE(:,2) = divE(:,2)*units.eps0;

%% Poincare Index
% temp = E1;E1 = E4;E4 = temp;
% E1 = -E1;
% E3 = -E3;
% [indices,sa1,sa2,sa3,sa4]=liner_vector_poincare_index(E1(:,2:4),E2(:,2:4),E3(:,2:4),E4(:,2:4));
[indices,sa1,sa2,sa3,sa4]=poincare_index_and_solid_angle(E1,E2,E3,E4);
% [indices,sa1,sa2,sa3,sa4] = monopole_index(E1(:,2:4),E2(:,2:4),E3(:,2:4),E4(:,2:4));
% indices=c_fgm_poincare_index(E1(:,2:4),E2(:,2:4),E3(:,2:4),E4(:,2:4));
% [indices,sa1,sa2,sa3,sa4]=c_4_poincare_index(E1(:,2:4),E2(:,2:4),E3(:,2:4),E4(:,2:4));
% c_eval('sa? = sa?/(4*pi);')
indices(abs(indices)<0.5) = 0;
indices(indices>=0.5) = 1;
indices(indices<=-0.5) = -1;
% plot(xx,indices)
%% FOTE 误差
eigVal_err_v2=E1(:,1);
monopole_index = zeros(length(x),1);
for ii=2:length(E1(:,1))  
deltB_null=reshape(gradE(ii,1:end),3,3);
if ~isnan(deltB_null)
    [V,D] = eig(deltB_null);
    % Figure 1o    以最大特征值归一化 |Δ·B|/|ΔxB|
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]);  
    if  indices(ii)~=0 && isreal(D(1,1)) && isreal(D(2,2)) && isreal(D(3,3)) && D(1,1) > 0 && D(2,2) > 0 && D(3,3) > 0
        monopole_index(ii) = 1;
    elseif indices(ii)~=0 && isreal(D(1,1)) && isreal(D(2,2)) && isreal(D(3,3)) && D(1,1) < 0 && D(2,2) < 0 && D(3,3) < 0
        monopole_index(ii) = -1;
    end
else
    eigVal_err_v2(ii,2)=nan;
end
end
%% plot
n=3;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(2);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 60; ySize = 40; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%%
h(i)=irf_subplot(n,1,-i);
plot(xx,divE(:,2), 'Linewidth',0.75); hold on;
c_eval("plot(xx,0*divE(:,2),'k--', 'Linewidth',0.75);"); hold off;
grid off;
c_eval("set(gca,'Ylim',[min(divE(:,2))-0.1*abs(min(divE(:,2))) max(divE(:,2))+0.1*max(divE(:,2))]);");
pos1=get(gca,'pos');
ylabel('divE [V/m^2]','fontsize',10);
set(gca,'xtick',[])
i = i+1;
%%
h(i)=irf_subplot(n,1,-i);
plot(xx,indices, 'Linewidth',0.75); hold on;
grid off;
set(gca,'ylim',[-1.1,1.1], 'ytick',[-1 0 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('PI','fontsize',10);
set(gca,'xtick',[])
i = i+1;
%%
h(i)=irf_subplot(n,1,-i);
plot(xx,monopole_index, 'Linewidth',0.75); hold on;
grid off;
set(gca,'ylim',[-1.1,1.1], 'ytick',[-1 0 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('monopole index','fontsize',10);
set(gca,'xtick',[])
i = i+1;
%%
% h(i)=irf_subplot(n,1,-i);
% plot(xx,sa1, 'Linewidth',0.75); hold on;
% grid off;
% set(gca, 'ytick',[-2*pi -pi 0 pi 2*pi],'yticklabels',{'-2\pi' '-\pi' '0' '\pi' '2\pi'},'fontsize',9);
% pos1=get(gca,'pos');
% ylabel('SA 123','fontsize',10);
% set(gca,'xtick',[])
% i = i+1;
%%
% h(i)=irf_subplot(n,1,-i);
% plot(xx,sa2, 'Linewidth',0.75); hold on;
% grid off;
% set(gca, 'ytick',[-2*pi -pi 0 pi 2*pi],'yticklabels',{'-2\pi' '-\pi' '0' '\pi' '2\pi'},'fontsize',9);
% pos1=get(gca,'pos');
% ylabel('SA 142','fontsize',10);
% set(gca,'xtick',[])
% i = i+1;
%%
% h(i)=irf_subplot(n,1,-i);
% plot(xx,sa3, 'Linewidth',0.75); hold on;
% grid off;
% set(gca, 'ytick',[-2*pi -pi 0 pi 2*pi],'yticklabels',{'-2\pi' '-\pi' '0' '\pi' '2\pi'},'fontsize',9);
% pos1=get(gca,'pos');
% ylabel('SA 134','fontsize',10);
% set(gca,'xtick',[])
% i = i+1;
%%
% h(i)=irf_subplot(n,1,-i);
% plot(xx,sa4, 'Linewidth',0.75); hold on;
% grid off;
% set(gca, 'ytick',[-2*pi -pi 0 pi 2*pi],'yticklabels',{'-2\pi' '-\pi' '0' '\pi' '2\pi'},'fontsize',9);
% pos1=get(gca,'pos');
% ylabel('SA 243','fontsize',10);
% % set(gca,'xtick',[])
% xlabel('x [m]','fontsize',10);
% i = i+1;

%%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,err_4C(:,2), 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca, 'ytick',[-1 0 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Solid Angle 243','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
% % % %%
% % % h(i)=irf_subplot(n,1,-i);
% % % plot(x,eigVal_err_v2(:,2), 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca, 'ytick',[-1 0 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Solid Angle 243','fontsize',10);
% % % set(gca,'xtick',[])
% % % i = i+1;
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')


%%
figure(1)
plot3(a1(:,1),a1(:,2),a1(:,3), 'color','k'); hold on;
plot3(a2(:,1),a2(:,2),a2(:,3), 'color','r'); hold on;
plot3(a3(:,1),a3(:,2),a3(:,3), 'color','g'); hold on;
plot3(a4(:,1),a4(:,2),a4(:,3), 'color','b'); hold on;
% [~,tempidx] = min(divE(:,2));
% tempidx = find(indices==1);
% tempidx = median(tempidx);
tempidx = 11;
plot3(a1(tempidx,1),a1(tempidx,2),a1(tempidx,3),'o' ,'color','k','linewidth',5); hold on;
plot3(a2(tempidx,1),a2(tempidx,2),a2(tempidx,3),'o', 'color','r','linewidth',5); hold on;
plot3(a3(tempidx,1),a3(tempidx,2),a3(tempidx,3),'o', 'color','g','linewidth',5); hold on;
plot3(a4(tempidx,1),a4(tempidx,2),a4(tempidx,3),'o', 'color','b','linewidth',5); hold on;


for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[a',num2str(ii),'(tempidx,1),a',num2str(ii),'(tempidx,2),a',num2str(ii),'(tempidx,3);',...
    'a?(tempidx,1),a?(tempidx,2),a?(tempidx,3)];'],ii+1:4);  %% ♥
end
plot3(RR12(:,1),RR12(:,2),RR12(:,3),'--k');hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'--k');hold on;  
plot3(RR14(:,1),RR14(:,2),RR14(:,3),'--k');hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'--k');hold on;  
plot3(RR34(:,1),RR34(:,2),RR34(:,3),'--k');hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'--k');hold on;  

c_eval('E? = irf_abs(E?);');
quiver3(a1(tempidx,1),a1(tempidx,2),a1(tempidx,3),2*E1(tempidx,1)/E1(tempidx,4),...
    2*E1(tempidx,2)/E1(tempidx,4),2*E1(tempidx,3)/E1(tempidx,4) ,'color','k','MaxHeadSize','3'); hold on;
quiver3(a2(tempidx,1),a2(tempidx,2),a2(tempidx,3),2*E2(tempidx,1)/E2(tempidx,4),...
    2*E2(tempidx,2)/E2(tempidx,4),2*E2(tempidx,3)/E2(tempidx,4) ,'color','r','MaxHeadSize','3'); hold on;
quiver3(a3(tempidx,1),a3(tempidx,2),a3(tempidx,3),2*E3(tempidx,1)/E3(tempidx,4),...
    2*E3(tempidx,2)/E3(tempidx,4),2*E3(tempidx,3)/E3(tempidx,4) ,'color','g','MaxHeadSize','3'); hold on;
quiver3(a4(tempidx,1),a4(tempidx,2),a4(tempidx,3),2*E4(tempidx,1)/E4(tempidx,4),...
    2*E4(tempidx,2)/E4(tempidx,4),2*E4(tempidx,3)/E4(tempidx,4) ,'color','b','MaxHeadSize','3'); hold on;
axis equal

%%
figure(3)
subplot(3,1,1)
irf_plot([xx E1(:,1)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([xx E2(:,1)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([xx E3(:,1)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([xx E4(:,1)], 'color','b', 'Linewidth',0.75); hold on;
set(gca,'Ylim',[min([min(E1(:,1)) min(E2(:,1)) min(E3(:,1)) min(E4(:,1))]) ...
    max([max(E1(:,1)) max(E2(:,1)) max(E3(:,1)) max(E4(:,1))])]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ex [V/m]','fontsize',12);

subplot(3,1,2)
irf_plot([xx E1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([xx E2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([xx E3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([xx E4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
set(gca,'Ylim',[min([min(E1(:,2)) min(E2(:,2)) min(E3(:,2)) min(E4(:,2))]) ...
    max([max(E1(:,2)) max(E2(:,2)) max(E3(:,2)) max(E4(:,2))])]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ey [V/m]','fontsize',12);

subplot(3,1,3)
irf_plot([xx E1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([xx E2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([xx E3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([xx E4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
% set(gca,'Ylim',[min([min(E1(:,3)) min(E2(:,3)) min(E3(:,3)) min(E4(:,3))]) ...
%     max([max(E1(:,3)) max(E2(:,3)) max(E3(:,3)) max(E4(:,3))])]);
set(gca,'Ylim',[-1 1]);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Ez [V/m]','fontsize',12);