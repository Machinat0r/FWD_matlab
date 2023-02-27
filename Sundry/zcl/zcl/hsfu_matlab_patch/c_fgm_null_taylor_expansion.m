

function null=c_fgm_null_taylor_expansion(varargin)
%This function plots the null position and null type
%
%*************************************************************************
%**if you find any bugs in this function, pls contact huishanf@gmail.com**
%*************************************************************************
%Input parameter (8 parameters or 9)
%1-4th parameters: B1-B4 (1st colomn is time; 2nd to 4th colomn is Bx, By, Bz)
%5-8th parameters: R1-R4 (1st colomn is time; 2nd to 4th colomn is X, Y, Z)
%9th parameter (optional): threshold for simplify 3D null to 2D null (0<threshold<0.5), default is 0.25
%Notice: the period should not be too long. Otherwise, the script runs slow and the figure is messy
%
%see Fu et al. [2014JGR] for details
%--------------------------------------------------------------------------------------
%Example: 
%c_fgm_null_taylor_expansion(B1,B2,B3,B4, R1,R2,R3,R4);
%c_fgm_null_taylor_expansion(B1,B2,B3,B4, R1,R2,R3,R4, 0.25);
%--------------------------------------------------------------------------------------

%----writen by Huishan Fu at BUAA (2014-05-27)----



[ax,args,nargs] = axescheck(varargin{:});
if nargs==0, % nothing to plot
	return;
end
if isempty(ax), % if no axis input, then plot in current axis
	ax=gca;
end
if nargs==9, 
    thresold=args{9};
else
    thresold=0.25; 
end


B1=args{1}; B2=args{2}; B3=args{3}; B4=args{4};
R1=args{5}; R2=args{6}; R3=args{7}; R4=args{8};


B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);

for ic=1:4
  c_eval(['R?=irf_resamp(R?,B?);'],ic);
end

indPoincare=c_fgm_poincare_index(B1,B2,B3,B4);



for ic=1:4
  c_eval('Bmag?=irf_abs(B?);',ic);
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



%Null type identification
for ii=1:length(d_dot_B(:,1))
    mksizSim(ii)=5;
    
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
end
eigVal_err(:,1)=d_dot_B(:,1);



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

R0_1=irf_add(1,R1,-1,dR1);
R0_2=irf_add(1,R2,-1,dR2);
R0_3=irf_add(1,R3,-1,dR3);
R0_4=irf_add(1,R4,-1,dR4);

minX=min(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
maxX=max(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
minY=min(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
maxY=max(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
minZ=min(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);
maxZ=max(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);



%% Init figure
n_subplots=8;

i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(61);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 12; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = 5; yTop = -1;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

lnwid=0.5;
Re=6372;

%% Poincare index plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(gca, [indPoincare(:,1) indPoincare(:,2)*0], 'k--', 'Linewidth',lnwid); hold on;
irf_plot(gca, indPoincare, 'color','k', 'Linewidth',0.75); hold off;

ylabel('Poincare Index');
set(gca,'Ylim',[-1.5 1.5])
grid off;


%% Xnull plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([R0_4(:,1) R0_4(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([R0_4(:,1) minX], 'b--', 'Linewidth',lnwid); hold on;
irf_plot([R0_4(:,1) maxX], 'b--', 'Linewidth',lnwid); hold off;

set(gca, 'Ylim',[min(minX)-400 max(maxX)+600]);
ylabel('Xnull [km]');
grid off;


%% Ynull plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([R0_4(:,1) R0_4(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([R0_4(:,1) minY], 'b--', 'Linewidth',lnwid); hold on;
irf_plot([R0_4(:,1) maxY], 'b--', 'Linewidth',lnwid); hold off;

set(gca, 'Ylim',[min(minY)-300 max(maxY)+400]);
ylabel('Ynull [km]');
grid off;


%% Znull plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([R0_4(:,1) R0_4(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([R0_4(:,1) minZ], 'b--', 'Linewidth',lnwid); hold on;
irf_plot([R0_4(:,1) maxZ], 'b--', 'Linewidth',lnwid); hold off;

set(gca, 'Ylim',[min(minZ)-100 max(maxZ)+100]);
ylabel('Znull [km]');
grid off;


%% distance plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(dRmag1(:,[1 5]), 'color','k', 'Linewidth',lnwid); hold on;
irf_plot(dRmag2(:,[1 5]), 'color','r', 'Linewidth',lnwid); hold on;
irf_plot(dRmag3(:,[1 5]), 'color','g', 'Linewidth',lnwid); hold on;
irf_plot(dRmag4(:,[1 5]), 'color','b', 'Linewidth',lnwid); hold off;

set(gca, 'Ylim',[0 1000],'Ytick',[0:200:1000]);
ylabel('\deltaR [km]');
grid off;
irf_legend({'C1','C2','C3','C4'},[0.02, 0.95],'color','cluster')


%% minmium distance plot 1
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(gca, dRmin, 'color','k', 'Linewidth',lnwid); hold on;
for ii=1:length(dRmin(:,1))
    irf_plot(gca, dRmin(ii,:), [type(ii) clr(ii)], 'MarkerSize',5,'MarkerFaceColor',faceclr(ii)); hold on;
end
hold off

set(gca, 'Ylim',[0 1000],'Ytick',[0:200:1000]);
ylabel('\deltaR_{min} [km]');
grid off;
set(gca,'ColorOrder',[0 0 0]);
irf_legend(gca,{'3D null'},[0.02, 0.95]);


%% minmium distance plot 2
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot(gca, dRmin, 'color','k', 'Linewidth',lnwid); hold on;
for ii=1:length(dRmin(:,1))
    irf_plot(gca, dRmin(ii,:), [typeSim(ii) clrSim(ii)], 'MarkerSize',mksizSim(ii),'MarkerFaceColor',faceclrSim(ii)); hold on;
end
hold off

set(gca, 'Ylim',[0 1000],'Ytick',[0:200:1000]);
ylabel('\deltaR_{min} [km]');
grid off;
set(gca,'ColorOrder',[0 0 0]);
irf_legend(gca,{'3D -> 2D'},[0.02, 0.95]);


%% error plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([err_4C(:,1) err_4C(:,2)],'k', 'Linewidth',lnwid); hold on;
irf_plot([eigVal_err(:,1) eigVal_err(:,2)], 'b', 'Linewidth',lnwid); hold on;
irf_plot([err_4C(:,1) err_4C(:,2)*0+40], 'k--', 'Linewidth',lnwid); hold off;

grid off;
set(gca,'Ylim',[0 200]);
legend('|\nabla\cdot\bf{B}|/|\nabla\times\bf{B}|','|(\lambda_1+\lambda_2+\lambda_3)/\lambda_{am}|','Location','NorthEast');
ylabel('null type error [%]');



%% Annotation
irf_adjust_panel_position
irf_zoom([B1(1,1) B1(end,1)],'x',h)
set(h(1:end-1),'XTickLabe','');


% set(gcf,'render','painters');
% print(gcf, '-dpdf', ['tmp.pdf']);


