function monopoleSearch_MMS_plot(varargin)
% tint
%% Load data
narginchk(1,2)
OutputDir = varargin{end};
if nargin == 1
tint = evalin('caller','tint');
elseif nargin == 2
    tint = evalin('caller',varargin{1});
else
disp('Incorrect number of input parameters. See usage:')
help CalError
return
end
%%
ic = 1:4;
c_eval("B?_ts=mms.get_data('B_gsm_brst',tint,?);");
c_eval('B?_gsm = irf.ts2mat(B?_ts);'); 
c_eval('B? = irf_abs(B?_gsm);');
c_eval('B? = irf_resamp(B?,B1);',2:4);

Pos = mms.get_data('R_gsm',tint);
c_eval('R? = Pos.gsmR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')
c_eval('R? = irf_resamp(R?,B1);')
CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

%% Parameaters
PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
PI(abs(PI)<0.5) = 0;

LocPoint = zeros(length(PI),3)*nan;
LocRes = cell(length(PI),1);
Q = zeros(length(PI),1)*nan;
resQ = cell(length(PI),1);

Qerror = ones(length(PI),1)*1000;
Locerror = ones(length(PI),1)*200;
dLoc = ones(length(PI),15)*5;

% div
gradB=c_4_grad('R?','B?_gsm','grad');
divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度

PI_id = find(PI~=0)';
for i = PI_id
flag_m = 0;
time_flagm = 0;
clc;
disp(['current calculate:',num2str(i),'/',num2str(length(PI))]);

MultiPower = ceil(max([log10(abs(R1(i,2:4))),log10(abs(R2(i,2:4))),log10(abs(R3(i,2:4))),log10(abs(R4(i,2:4)))]));

if MultiPower > 3
    continue
end

RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
    'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;

% solve
[Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?_gsm',i,i*sign(divB(i,2)),RR_mean,1);

id = nchoosek(1:6,2);
c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
tempd = [];
c_eval('tempd = [tempd,tempd?(4)/RR_mean];',1:15);
dLoc(i,:) = tempd;
end
%% calculate error
Qerror = zeros(length(PI),1);
Locerror = zeros(length(PI),1);

for i = 1:length(PI)
if ~isnan(resQ{i})
    Qerror(i) = abs(100*std(resQ{i})/Q(i));
else
    Qerror(i) = 1001;
end

if ~isnan(LocRes{i})
tri_a = delaunayTriangulation([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)]);
[~,volume_a] = convexHull(tri_a);
tri = delaunayTriangulation(LocRes{i});%%delaunay三角剖分
    if size(tri.Points,1)==1
        volume = 0;
    else
    [~,volume] = convexHull(tri);%%计算多面体体积
    end
Locerror(i) = 100*volume/volume_a;
else
Locerror(i) = 200;
end
end              
%% FOTE err & div
gradB=c_4_grad('R?','B?_gsm','grad');
divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度
%% Init figure 1
n=9;
i=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Btotal
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) B1(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2(:,1) B2(:,5)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3(:,1) B3(:,5)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4(:,1) B4(:,5)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[min([min(B1(:,5)) min(B2(:,5)) min(B3(:,5)) min(B4(:,5))])-5 ...
    max([max(B1(:,5)) max(B2(:,5)) max(B3(:,5)) max(B4(:,5))])+5]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('|B| [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bx
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gsm(:,1) B1_gsm(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gsm(:,1) B2_gsm(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gsm(:,1) B3_gsm(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gsm(:,1) B4_gsm(:,2)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,2)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
set(gca,'Ylim',[min([min(B1_gsm(:,2)) min(B2_gsm(:,2)) min(B3_gsm(:,2)) min(B4_gsm(:,2))])-5 ...
    max([max(B1_gsm(:,2)) max(B2_gsm(:,2)) max(B3_gsm(:,2)) max(B4_gsm(:,2))])+5]);
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,2))])-100 max([max(B?_gsm(:,2))])+100]);",ic);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bx [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% By
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gsm(:,1) B1_gsm(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gsm(:,1) B2_gsm(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gsm(:,1) B3_gsm(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gsm(:,1) B4_gsm(:,3)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,3)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,3))])-100 max([max(B?_gsm(:,3))])+100]);",ic);
set(gca,'Ylim',[min([min(B1_gsm(:,3)) min(B2_gsm(:,3)) min(B3_gsm(:,3)) min(B4_gsm(:,3))])-5 ...
    max([max(B1_gsm(:,3)) max(B2_gsm(:,3)) max(B3_gsm(:,3)) max(B4_gsm(:,3))])+5]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('By [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% Bz
h(i)=irf_subplot(n,1,-i);
irf_plot([B1_gsm(:,1) B1_gsm(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([B2_gsm(:,1) B2_gsm(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([B3_gsm(:,1) B3_gsm(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([B4_gsm(:,1) B4_gsm(:,4)], 'color','b', 'Linewidth',0.75); hold on;
%irf_plot([B1(:,1) B1], 'color','k', 'Linewidth',0.75); hold on;
c_eval("irf_plot([B?_gsm(:,1) 0*B?_gsm(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,4))])-100 max([max(B?_gsm(:,4))])+100]);",ic);
set(gca,'Ylim',[min([min(B1_gsm(:,4)) min(B2_gsm(:,4)) min(B3_gsm(:,4)) min(B4_gsm(:,4))])-5 ...
    max([max(B1_gsm(:,4)) max(B2_gsm(:,4)) max(B3_gsm(:,4)) max(B4_gsm(:,4))])+5]);
% set(gca,'Ylim',[-40 60], 'ytick',[-40:20:60],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
ylabel('Bz [nT]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% η
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([err_4C(:,1) err_4C(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % % % c_eval("set(gca,'Ylim',[min([min(B?_gsm(:,4))])-3 max([max(B?_gsm(:,4))])+3]);",ic);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % c_eval("set(gca,'Ylim',[min(err_4C(:,2))-0.1 max(err_4C(:,2))+0.1]);",ic);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('η','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% ξ
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([eigVal_err_v2(:,1) eigVal_err_v2(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % c_eval("set(gca,'Ylim',[min(eigVal_err_v2(:,2))-0.1 max(eigVal_err_v2(:,2))+0.1]);",ic);
% % % % set(gca,'Ylim',[0 100], 'ytick',[0 50 100],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('ξ','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% divB
h(i)=irf_subplot(n,1,-i);
irf_plot([divB(:,1) divB(:,2)], 'color','k', 'Linewidth',0.75); hold on;

% c_eval("irf_plot([eigVal_err_v2(:,1) 0*eigVal_err_v2(:,4)],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
% set(gca,'Ylim',[0 0.1], 'ytick',[0 0.05 0.1],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('divB [nT/km^2]','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% PI
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) PI], 'color','k', 'Linewidth',0.75); hold on;

c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
grid off;
% c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
pos1=get(gca,'pos');
set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
ylabel('PI ','fontsize',12);
% irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
i=i+1;
%% monopole_index
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) monopole_index], 'color','k', 'Linewidth',0.75); hold on;
% % % 
% % % % c_eval("irf_plot([B1(:,1) PI],'k--', 'Linewidth',0.75);",ic); hold off;
% % % grid off;
% % % % c_eval("set(gca,'Ylim',[min(divB(:,2)) max(divB(:,2))]);",ic);
% % % set(gca,'Ylim',[-1.2 1.2], 'ytick',[-1 0 1],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % set(gca,'ColorOrder',[[0 0 1];[1 0 0]]);
% % % % irf_legend(gca,{'B_x','B_z'},[0.97 0.92]);
% % % ylabel('monopole index ','fontsize',12);
% % % % irf_legend(gca,{'B_N'},[pos2(1)+1.15*pos2(3),pos2(2)]);
% % % % irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',12)
% % % i=i+1;
%% Location Error
% % % h(i)=irf_subplot(n,1,-i);
% % % irf_plot([B1(:,1) Locerror], 'color','k', 'Linewidth',0.75); hold on;
% % % grid off;
% % % % set(gca,'Ylim',[0 max(Locerror)]);
% % % set(gca,'Ylim',[0,200], 'ytick',[0 100 200],'fontsize',9);
% % % pos1=get(gca,'pos');
% % % ylabel('Loc Err [%]','fontsize',12);
% % % i=i+1;
%% dLoc Error
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) max(dLoc,[],2)], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Locerror)]);
set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('Loc Err','fontsize',12);
i=i+1;
%% meandLoc Error
h(i)=irf_subplot(n,1,-i);

irf_plot([B1(:,1) mean(dLoc,2)], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Locerror)]);
set(gca,'Ylim',[0,1], 'ytick',[0 0.5 1],'fontsize',9);
pos1=get(gca,'pos');
ylabel('mean','fontsize',12);
i=i+1;
%% Q Error
h(i)=irf_subplot(n,1,-i);
irf_plot([B1(:,1) Qerror], 'color','k', 'Linewidth',0.75); hold on;
grid off;
% set(gca,'Ylim',[0 max(Qerror)]);
set(gca,'Ylim',[0,1000], 'ytick',[0 500 1000],'fontsize',9);
pos1=get(gca,'pos');
ylabel('Q Err [%]','fontsize',12);
i=i+1;
%% Adjust the position
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align;
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
colormap(jet)
figname = [OutputDir,datestr(datenum(1970,1,1,0,0,0)+B1(1)/86400,'yyyymmddTHHMMSS')];    
print(gcf, '-dpng', [figname '.png']);
clf

for tempidx_B1 = PI_id
%% id
[~,tempidx_R] = sort(abs(R1(:,1)-B1(tempidx_B1,1)));
tempidx_R = tempidx_R(1);
%% Init Figure 2
figure(2)
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% Coordinate
c_eval('tempR?(:,2:4) = R?(:,2:4)-LocPoint(tempidx_B1,:);');
tempLocRes{tempidx_B1} = LocRes{tempidx_B1}-LocPoint(tempidx_B1,:);
tempLocPoint = LocPoint - LocPoint(tempidx_B1,:);
%% Location 
cor = 'krgb';
c_eval("plot3(tempR?(tempidx_R,2),tempR?(tempidx_R,3),tempR?(tempidx_R,4),'s' ,'color',cor(?),'linewidth',5);hold on;")
%% Tetrahedron configuration
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[tempR',num2str(ii),'(tempidx_R,2),tempR',num2str(ii),'(tempidx_R,3),tempR',num2str(ii),'(tempidx_R,4);',...
    'tempR?(tempidx_R,2),tempR?(tempidx_R,3),tempR?(tempidx_R,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
plot3(RR12(:,1),RR12(:,2),RR12(:,3),'--k');hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'--k');hold on;  
plot3(RR14(:,1),RR14(:,2),RR14(:,3),'--k');hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'--k');hold on;  
plot3(RR34(:,1),RR34(:,2),RR34(:,3),'--k');hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'--k');hold on;  
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
xlabel('e_1 [km]','fontsize',12);
ylabel('e_2 [km]','fontsize',12);
zlabel('e_3 [km]','fontsize',12);

%% Quiver
c_eval('B?_gsm = irf_abs(B?_gsm);');
c_eval('B?_gsm = irf_resamp(B?_gsm,B1_gsm);',2:4);
maxB = max([B1_gsm(tempidx_B1,5),B2_gsm(tempidx_B1,5),B3_gsm(tempidx_B1,5),B4_gsm(tempidx_B1,5)]);
c_eval("quiver3(tempR?(tempidx_R,2),tempR?(tempidx_R,3),tempR?(tempidx_R,4),RR_mean*B?_gsm(tempidx_B1,2)/maxB,RR_mean*B?_gsm(tempidx_B1,3)/maxB,RR_mean*B?_gsm(tempidx_B1,4)/maxB,'color',cor(?));hold on;")

%% Loc res
plotPolyhedron(tempLocRes{tempidx_B1}(:,1),tempLocRes{tempidx_B1}(:,2),tempLocRes{tempidx_B1}(:,3),'#FFB8CE',0.3);
% idx = [166:198];
% plot3(LocPoint(idx,1),LocPoint(idx,2),LocPoint(idx,3),'*','color','#FFBAF1');
% line(LocPoint(173:195,1),LocPoint(173:195,2),LocPoint(173:195,3),'color','#FFBAF1');
plot3(tempLocPoint(tempidx_B1,1),tempLocPoint(tempidx_B1,2),tempLocPoint(tempidx_B1,3),'*','color','m');
plot3(tempLocRes{tempidx_B1}(:,1),tempLocRes{tempidx_B1}(:,2),tempLocRes{tempidx_B1}(:,3),'*','color','#FFB8CE');
%%
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname = [OutputDir,datestr(datenum(1970,1,1,0,0,0)+B1(tempidx_B1,1)/86400,'yyyymmddTHHMMSS'),'-Configuration'];
savefig([figname '.fig']);
clf
end
end