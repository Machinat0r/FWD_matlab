%% Init figure 3
h = figure(3);
set(gcf,'PaperUnits','centimeters')
xSize = 50; ySize = 100; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])

%% plot distance
h(1) = subplot(5,1,1);
for i = 1:length(resQ)
    RR_temp = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
    'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
c_eval(['RR_temp=RR_temp+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean(i) = RR_temp(4)/6;
end

id = nchoosek(1:6,2);
for i = 1:length(LocRes)
    c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
    ttd = [];
    c_eval('ttd = [ttd,tempd?(4)/RR_mean(i)];',1:15);
    dLoc(i,:) = ttd;
end
    meand = mean(dLoc,2);
    
% 15 errorbars
resd = dLoc-meand;
resdpos = resd; resdneg = resd;
resdpos(resdpos<0) = 0;resdneg(resdneg>=0) = 0;
c_eval('d?neg = errorbar(B1(:,1),meand,resdneg(:,?),resdpos(:,?));hold on;',1:15);
line(B1(:,1),meand,'color','#BEBEBE','LineWidth',1);
scatter(B1(:,1),meand,'d','filled','MarkerFaceColor','#EDB120')
% c_eval('d?neg = errorbar(B1(:,1),meand,resdneg(:,?),zeros(size(resdneg(:,?))));hold on;',1:6);
% c_eval('d?pos = errorbar(B1(:,1),meand,zeros(size(resdpos(:,?))),resdpos(:,?));hold on;',1:6);
% c_eval("d?pos.Marker = 'o';",1:6)
% c_eval("d?pos.Color = 'k';",1:6)
% c_eval("d?neg.Marker = 'o';",1:6)
c_eval("d?neg.Color = '#696969';",1:15)
c_eval("d?neg.LineWidth = 0.8;",1:15)
c_eval("d?neg.CapSize = 10;",1:15)

% errorbar
% % % resd = dLoc-meand;
% % % resdpos = max(resd,[],2); resdneg = min(resd,[],2);
% % % e1 = errorbar(B1(:,1),meand,resdneg,resdpos);hold on;
% % % e1.LineWidth = 1;
% % % e1.Color = '#696969';
% % % quart = quantile(dLoc,[.25,.75],2);
% % % c_eval("line([B1(?,1),B1(?,1)],[quart(?,1),quart(?,2)],'color','#BEBEBE','linewidth',5);hold on;",1:length(Q))
% % % line(B1(:,1),meand,'linewidth',2,'color','#BEBEBE');hold on;
% % % scatter(B1(:,1),meand,'d','filled','MarkerFaceColor','#EDB120')

ax = gca;
ax.YLim = [0,2];
irf_timeaxis(ax,'date')
ylabel('Loc ','fontsize',12);
%% subplot Q
h(2) = subplot(5,1,2);

% QErrorbar(B1(:,1),Q,resQ)

% 6 errorbars
for i = 1:length(resQ)
c_eval('resQ?(i) = resQ{i}(?)-Q(i);',1:6)
end
c_eval('resQ?pos = resQ?;',1:6);c_eval('resQ?neg = resQ?;',1:6);
c_eval('resQ?pos(resQ?<0) = 0;',1:6);
c_eval('resQ?neg(resQ?>=0) = 0;',1:6);

c_eval('e?pos = errorbar(B1(:,1),Q,resQ?neg,resQ?pos);hold on;',1:6)
% c_eval('e?pos = errorbar(B1(:,1),Q,resQ?neg,zeros(size(resQ?pos)));hold on;',1:6)
% c_eval('e?neg = errorbar(B1(:,1),Q,zeros(size(resQ?neg)),resQ?pos);hold on;',1:6)
% c_eval("e?pos.LineWidth = 1.5;",1:6)
scatter(B1(:,1),Q,'d','filled','MarkerFaceColor','#EDB120')
c_eval("e?pos.Color = 'k';",1:6)
% c_eval("e?neg.Marker = 'o';",1:6)
% c_eval("e?neg.Color = 'k';",1:6)

% errorbar
% % % c_eval('resQmat(?,:) = resQ{?}-Q(?);',1:length(Q))
% % % 
% % % resQpos = max(resQmat,[],2); resQneg = min(resQmat,[],2);
% % % e2 = errorbar(B1(:,1),Q,resQneg,resQpos);hold on;
% % % e2.LineWidth = 1;
% % % e2.Color = '#696969';
% % % e2.CapSize = 0.01;
% % % quart = quantile(transpose(reshape(cell2mat(resQ),6,[])),[.25,.75],2);
% % % c_eval("line([B1(?,1),B1(?,1)],[quart(?,1),quart(?,2)],'color','#BEBEBE','linewidth',5);hold on;",1:length(Q))
% % % line(B1(:,1),Q,'linewidth',2,'color','#BEBEBE');hold on;
% % % scatter(B1(:,1),Q,'d','filled','MarkerFaceColor','#EDB120')
% % % scatter(B1(:,1),resQpos+Q,'filled','MarkerFaceColor','k')

ylabel('Q [nT·km^2] ','fontsize',12);
ax = gca;
ax.YScale = 'log';
ax.YLim = [1e5,1e8];
irf_timeaxis(ax,'date')
%% plot div err
h(3) = subplot(5,1,3);

line(B1(:,1),divErr(:,2),'color','#8470FF');hold on;
s1 = scatter(B1(:,1),divErr(:,2),100,'p','filled','MarkerFaceColor','#7E2F8E');hold on;

ax = gca;
ax.YLim = [0,100];
irf_timeaxis(ax,'date')
ylabel('div Error [%]','fontsize',12);
%% Model direction Error
h(4) = subplot(5,1,4);
cor = {'#000000','#D95319','#77AC30','#0072BD'};
c_eval("line(B?(:,1),Btheta?,'color',cor{?});hold on;")
s1 = scatter(B1(:,1),Btheta1,50,'k','filled');hold on;
s2 = scatter(B2(:,1),Btheta2,50,'filled');hold on;
s2.MarkerFaceColor = '#D95319';
s3 = scatter(B3(:,1),Btheta3,50,'filled');hold on;
s3.MarkerFaceColor = '#77AC30';
s4 = scatter(B4(:,1),Btheta4,50,'filled');hold on;
s4.MarkerFaceColor = '#0072BD';


ax = gca;
ax.YLim = [0,100];
irf_timeaxis(ax,'date')
ylabel('Direction Error [%]','fontsize',12);
%% Model strength Error
h(5) = subplot(5,1,5);
c_eval("line(B?(:,1),100*Blength?,'color',cor{?});hold on;")
% cor = {'#000000','#D95319','#77AC30','#0072BD'};
s1 = scatter(B1(:,1),100*Blength1,50,'k','filled');hold on;
s2 = scatter(B2(:,1),100*Blength2,50,'filled');hold on;
s2.MarkerFaceColor = '#D95319';
s3 = scatter(B3(:,1),100*Blength3,50,'filled');hold on;
s3.MarkerFaceColor = '#77AC30';
s4 = scatter(B4(:,1),100*Blength4,50,'filled');hold on;
s4.MarkerFaceColor = '#0072BD';


ax = gca;
ax.YLim = [0,100];
irf_timeaxis(ax,'date')
ylabel('Strength Error [%]','fontsize',12);
%%
set(gcf,'render','painters');
irf_zoom(tint,'x',h(1:end));
irf_plot_axis_align(h);
% ax.YLimMode = 'auto';
set(gcf,'paperpositionmode','auto')
