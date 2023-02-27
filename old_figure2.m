%% Init Figure 2
h2 = figure(2);
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 80; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Index id
units = irf_units;
tempidx_B1 = find(abs(monopole_index)==0.5);
% tempidx_B1 = find(PI~=0);
tempidx_B1 = tempidx_B1(1);
% tempidx_B1 = round((tempidx_B1(1)+tempidx_B1(end))/2);
 
%     [~,tempidx_B] = max(abs(divB(:,2)));
c_eval('[~,tempidx_B?] = sort(abs(B?_gsm(:,1)-B1_gsm(tempidx_B1,1)));',2:4);
c_eval('tempidx_B? = tempidx_B?(1);')
[~,tempidx_R] = sort(abs(R1_gsm(:,1)-B1(tempidx_B1,1)));
tempidx_R = tempidx_R(1);
 
% [~,curlB]=HelmholtzDecomposition('R?_gsm(tempidx_R,:)','B?_gsm(tempidx_B?,:)');
[Br1,Br2,Br3,Br4,curlB] =  EliminateCurl('R?_gsm','B?_gsm',tempidx_B1,tempidx_R);
% c_eval('Br? = [B?_gsm(:,1) Br?];')
 
PlotFlag = 3;
 
% [Loc,res,LocPoint,gridPoint] = CalError('R?_gsm','B?_gsm',tempidx_R,sign(monopole_index(tempidx_B1))*tempidx_B1,PlotFlag);
[Loc,res,LocPoint,gridPoint] = CalError('R?_gsm','Br?',tempidx_R,sign(monopole_index(tempidx_B1))*1,PlotFlag);

c_eval('Rtemp? = R?_gsm(tempidx_R,:);')
gradB=c_4_grad('Rtemp?','Br?','grad');
% gradB=c_4_grad('R?_gsm','B?_gsm','grad');
deltB_null=reshape(gradB(1,2:end),3,3);
[V,D] = eig(deltB_null);
 
 
%% transform coordinate
TimeZero = datestr(datenum(1970,1,1,0,0,0)+B1(tempidx_B1,1)/86400,'yyyy-mm-ddTHH:MM:SS.FFFZ');
[~, Trans_mat,~] = FOTE('R?_gsm','B?_gsm',TimeZero);
% [~, Trans_mat,~] = FOTE('R?_gsm','Br?',TimeZero);
if abs(monopole_index(tempidx_B1)) == 0.5
    c_eval('e? = Trans_mat.e?;',1:3);
else
    c_eval("e? = V(:,?)';",1:3);
    e3 = cross(e1,e2);
end
 
% Loc = irf_newxyz(Loc,e1,e2,e3);
c_eval('Loc{?} = irf_newxyz(Loc{?},e1,e2,e3);')
c_eval(['B?_gsm(:,2:4)=irf_newxyz(B?_gsm(:,2:4),e1,e2,e3);'],ic);
c_eval(['R?_gsm(:,2:4)=irf_newxyz(R?_gsm(:,2:4),e1,e2,e3);'],ic);
c_eval(['Br?(:,2:4)=irf_newxyz(Br?(:,2:4),e1,e2,e3);'],ic);
% % % 
% % for i = 1:size(R1_gsm,1)
% % CenterPoint(i,:) = (R1_gsm(i,2:4)+R2_gsm(i,2:4)+R3_gsm(i,2:4)+R4_gsm(i,2:4))/4;
% % end
% % c_eval('R?_gsm(:,2:4) = (R?_gsm(:,2:4) - CenterPoint)*units.RE/1000;');
% % c_eval('Loc{?} = (Loc{?} - CenterPoint(tempidx_R,:))*units.RE/1000;');

% [Q,resQ] = solveMonopole('R?_gsm(tempidx_R,:)','Br?',LocPoint,sign(monopole_index(tempidx_B1))*tempidx_B1);

if PlotFlag ~= 1 
LocPoint = irf_newxyz(LocPoint,e1,e2,e3);
gridPoint = irf_newxyz(gridPoint,e1,e2,e3);
CenterPoint = repmat(LocPoint,length(R1_gsm),1);
c_eval('R?_gsm(:,2:4) = (R?_gsm(:,2:4) - CenterPoint)*units.RE/1000;');
c_eval('Loc{?} = (Loc{?} - CenterPoint(tempidx_R,:))*units.RE/1000;');
LocPoint = (LocPoint - CenterPoint(tempidx_R,:))*units.RE/1000;
gridPoint = (gridPoint - CenterPoint(tempidx_R,:))*units.RE/1000;
end

[Q,resQ] = solveMonopole('R?_gsm(tempidx_R,:)','B?_gsm(tempidx_B?,:)',LocPoint,sign(monopole_index(tempidx_B1))*tempidx_B1); 
%% Location 
cor = {'k','r','g','b'};
c_eval("plot3(R?_gsm(tempidx_R,2),R?_gsm(tempidx_R,3),R?_gsm(tempidx_R,4),'o' ,'color',cor{?},'linewidth',5); hold on;");
 
%% Tetrahedron configuration
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'_gsm(tempidx_R,2),R',num2str(ii),'_gsm(tempidx_R,3),R',num2str(ii),'_gsm(tempidx_R,4);',...
    'R?_gsm(tempidx_R,2),R?_gsm(tempidx_R,3),R?_gsm(tempidx_R,4)];'],ii+1:4);  %% ?
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
% plot3(RR12(:,1),RR12(:,2),RR12(:,3),'--k');hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'--k');hold on;  
% plot3(RR14(:,1),RR14(:,2),RR14(:,3),'--k');hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'--k');hold on;  
% plot3(RR34(:,1),RR34(:,2),RR34(:,3),'--k');hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'--k');hold on; 
plot3(RR12(:,1),RR12(:,2),RR12(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  
plot3(RR14(:,1),RR14(:,2),RR14(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  
plot3(RR34(:,1),RR34(:,2),RR34(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'color',[0.5,0.5,0.5],'LineWidth',1.5);hold on;  
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.97 0.92]);
irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
% xlabel('X_G_S_M [R_E]','fontsize',12);
% ylabel('Y_G_S_M [R_E]','fontsize',12);
% zlabel('Z_G_S_M [R_E]','fontsize',12);
xlabel('e1 [km]','fontsize',12);
ylabel('e2 [km]','fontsize',12);
zlabel('e3 [km]','fontsize',12);
 
%% Quiver
c_eval('B?_gsm = irf_abs(B?_gsm);');
Rmean = [0,0,0];
c_eval('Rmean = Rmean+ R?_gsm(tempidx_R,2:4);');Rmean = Rmean/4;
if PlotFlag == 1
c_eval("quiver3(R?_gsm(tempidx_R,2),R?_gsm(tempidx_R,3),R?_gsm(tempidx_R,4),RR_mean*B?_gsm(tempidx_B?,2)/B?_gsm(tempidx_B?,5)"+...
    ",RR_mean*B?_gsm(tempidx_B?,3)/B?_gsm(tempidx_B?,5),RR_mean*B?_gsm(tempidx_B?,4)/B?_gsm(tempidx_B?,5),0.3 ,'color',cor{?},'linewidth',0.5); hold on;")
c_eval('Br? = irf_abs(Br?);');
c_eval("quiver3(R?_gsm(tempidx_R,2),R?_gsm(tempidx_R,3),R?_gsm(tempidx_R,4),RR_mean*Br?(2)/Br?(5),"+...
    "RR_mean*Br?(3)/Br?(5),RR_mean*Br?(4)/Br?(5),0.3 ,'color',cor{?},'linewidth',1.3); hold on;")
elseif PlotFlag == 2
% c_eval('Br? = irf_abs(Br?);');
% c_eval("quiver3(R?_gsm(tempidx_R,2),R?_gsm(tempidx_R,3),R?_gsm(tempidx_R,4),RR_mean*Br?(2)/Br?(5),"+...
%     "RR_mean*Br?(3)/Br?(5),RR_mean*Br?(4)/Br?(5),0.3 ,'color',cor{?},'linewidth',1.3); hold on;")
c_eval("quiver3(R?_gsm(tempidx_R,2),R?_gsm(tempidx_R,3),R?_gsm(tempidx_R,4),RR_mean*B?_gsm(tempidx_B?,2)/B?_gsm(tempidx_B?,5)"+...
",RR_mean*B?_gsm(tempidx_B?,3)/B?_gsm(tempidx_B?,5),RR_mean*B?_gsm(tempidx_B?,4)/B?_gsm(tempidx_B?,5),0.3 ,'color',cor{?},'linewidth',1.3); hold on;")
end
%% axis
% if monopole_index(tempidx_B1) == 0.5
% quiver3(Rmean(1),Rmean(2),Rmean(3),RR_mean*eigV(1),RR_mean*eigV(2),RR_mean*eigV(3),0.5,'color','k'); hold on;
% eigV = -eigV;
% quiver3(Rmean(1),Rmean(2),Rmean(3),RR_mean*eigV(1),RR_mean*eigV(2),RR_mean*eigV(3),0.5,'color','k'); hold on;
% elseif monopole_index(tempidx_B1) == -0.5
% quiver3(Rmean(1)-0.5*RR_mean*eigV(1),Rmean(2)-0.5*RR_mean*eigV(2),Rmean(3)-0.5*RR_mean*eigV(3),RR_mean*eigV(1),RR_mean*eigV(2),RR_mean*eigV(3),0.5,'color','k'); hold on;
% eigV = -eigV;
% quiver3(Rmean(1)-0.5*RR_mean*eigV(1),Rmean(2)-0.5*RR_mean*eigV(2),Rmean(3)-0.5*RR_mean*eigV(3),RR_mean*eigV(1),RR_mean*eigV(2),RR_mean*eigV(3),0.5,'color','k'); hold on;
% end
%% Monopole Solution
if PlotFlag == 1
c_eval("plot3(Loc{?}(1,1),Loc{?}(1,2),Loc{?}(1,3),'o' ,'color','c','linewidth',3); hold on;")
Vertices = [Loc{1};Loc{2};Loc{3};Loc{4}];
Faces = [1 2 3;1 3 4;1 2 4;2 3 4];
a = patch('Faces',Faces,'Vertices',Vertices,'FaceColor','c');
alpha(a,0.2)
elseif PlotFlag == 2
plotPolyhedron(gridPoint(:,1),gridPoint(:,2),gridPoint(:,3));
% % for ii = 1:size(gridPoint,1)
% %     plot3(gridPoint(ii,1),gridPoint(ii,2),gridPoint(ii,3),'ro');
% % end
% % c_eval("plot3(Loc{?}(1,1),Loc{?}(1,2),Loc{?}(1,3),'*' ,'color','k','linewidth',0.5); hold on;")

% scatter3(gridPoint(:,1),gridPoint(:,2),gridPoint(:,3));
% % % cmp = colormap;colorbar;
% % % resQ = [min(resQ),max(resQ)];
% % % caxis(resQ);
% % % cmpid = fix((max(resQ)-Q)*(max(resQ)-min(resQ)));
% % % cmp = cmp(cmpid,:);
c_eval("plot3(LocPoint(1),LocPoint(2),LocPoint(3),'*' ,'color','m','linewidth',2); hold on;")
elseif PlotFlag == 3
    c_eval("plot3(LocPoint(1),LocPoint(2),LocPoint(3),'*' ,'color','m','linewidth',2); hold on;")
end
%% Projection
RR_max = 1.2*max(abs([R1_gsm(tempidx_R,2:4),R2_gsm(tempidx_R,2:4),R3_gsm(tempidx_R,2:4),R4_gsm(tempidx_R,2:4)]));
 
% Projection line
% % % lineProperties = {'linestyle',':','linewidth',0.8};
% % % c_eval('line([R?_gsm(tempidx_R,2) -RR_max],[R?_gsm(tempidx_R,3) R?_gsm(tempidx_R,3)],[R?_gsm(tempidx_R,4) R?_gsm(tempidx_R,4)],lineProperties{:});');
% % % c_eval('line([R?_gsm(tempidx_R,2) R?_gsm(tempidx_R,2)],[R?_gsm(tempidx_R,3) -RR_max],[R?_gsm(tempidx_R,4) R?_gsm(tempidx_R,4)],lineProperties{:});');
% % % c_eval('line([R?_gsm(tempidx_R,2) R?_gsm(tempidx_R,2)],[R?_gsm(tempidx_R,3) R?_gsm(tempidx_R,3)],[R?_gsm(tempidx_R,4) -RR_max],lineProperties{:});');
 
if PlotFlag == 1
% % % %     c_eval("plot3(Loc{?}(1,1),Loc{?}(1,2),Loc{?}(1,3),'o' ,'color','c','linewidth',3); hold on;")
% % %     Vertices1 = [[-RR_max,Loc{1}(2:3)];[-RR_max,Loc{2}(2:3)];[-RR_max,Loc{3}(2:3)];[-RR_max,Loc{4}(2:3)]];
% % %     Vertices2 = [[Loc{1}(1),-RR_max,Loc{1}(3)];[Loc{2}(1),-RR_max,Loc{2}(3)];[Loc{3}(1),-RR_max,Loc{3}(3)];[Loc{4}(1),-RR_max,Loc{4}(3)]];
% % %     Vertices3 = [[Loc{1}(1:2),-RR_max];[Loc{2}(1:2),-RR_max];[Loc{3}(1:2),-RR_max];[Loc{4}(1:2),-RR_max]];
% % %     Faces = [1 2 3;1 3 4;1 2 4;2 3 4];
% % %     patch('Faces',Faces,'Vertices',Vertices1,'FaceColor','c','FaceAlpha',0.3,'EdgeColor','none');
% % %     patch('Faces',Faces,'Vertices',Vertices2,'FaceColor','c','FaceAlpha',0.3,'EdgeColor','none');
% % %     patch('Faces',Faces,'Vertices',Vertices3,'FaceColor','c','FaceAlpha',0.3,'EdgeColor','none');
elseif PlotFlag == 2
text(LocPoint(1)-5,LocPoint(2),LocPoint(3),['Q = ',num2str(Q,3),'nT·km^2'])    
axis equal 
% axis([-5 5 -5 5 -5 5]);
axis([-RR_max RR_max -RR_max RR_max -RR_max RR_max]);

% S/C Projection
R_x = [R1_gsm(tempidx_R,2);R2_gsm(tempidx_R,2);R3_gsm(tempidx_R,2);R4_gsm(tempidx_R,2)];
R_y = [R1_gsm(tempidx_R,3);R2_gsm(tempidx_R,3);R3_gsm(tempidx_R,3);R4_gsm(tempidx_R,3)];
R_z = [R1_gsm(tempidx_R,4);R2_gsm(tempidx_R,4);R3_gsm(tempidx_R,4);R4_gsm(tempidx_R,4)];
R_side = ones(4,1)*(-RR_max);
Vx = [R_side,R_y,R_z];
Vy = [R_x,R_side,R_z];
Vz = [R_x,R_y,R_side];
Faces = [1,2,3,4;1,2,4,3;1,4,2,3;1,4,3,2;1,3,2,4;1,3,4,2];

figure('name','errorbar')
subplot(2,2,1)
plotPolyhedron(-RR_max*ones(size(gridPoint,1),1),gridPoint(:,2),gridPoint(:,3));hold on;
patch('Faces',Faces,'Vertices',Vx,'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1,'EdgeColor','none');hold on;
ylabel('e2 [km]');zlabel('e3 [km]');
axis equal; view(90,0);
ylim([-5,11]);zlim([-5,15]);

subplot(2,2,2)
plotPolyhedron(gridPoint(:,1),-RR_max*ones(size(gridPoint,2),1),gridPoint(:,3));hold on;
patch('Faces',Faces,'Vertices',Vy,'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1,'EdgeColor','none');hold on;
xlabel('e1 [km]');zlabel('e3 [km]');
axis equal; view(0,0);
xlim([-5,11]);zlim([-5,15]);

% Error Region Projection
% % % plotPolyhedron(-RR_max*ones(size(gridPoint,1),1),gridPoint(:,2),gridPoint(:,3));
% % % plotPolyhedron(gridPoint(:,1),-RR_max*ones(size(gridPoint,2),1),gridPoint(:,3));
% % % plotPolyhedron(gridPoint(:,1),gridPoint(:,2),-RR_max*ones(size(gridPoint,3),1));
 
% Monopole Projection
% % % lineProperties = {'linestyle',':','linewidth',0.8,'color','m'};
% % % c_eval('line([LocPoint(1) -RR_max],[LocPoint(2) LocPoint(2)],[LocPoint(3) LocPoint(3)],lineProperties{:});');
% % % c_eval('line([LocPoint(1) LocPoint(1)],[LocPoint(2) -RR_max],[LocPoint(3) LocPoint(3)],lineProperties{:});');
% % % c_eval('line([LocPoint(1) LocPoint(1)],[LocPoint(2) LocPoint(2)],[LocPoint(3) -RR_max],lineProperties{:});');
% % % c_eval("plot3(-RR_max,LocPoint(2),LocPoint(3),'*' ,'color','m','linewidth',2); hold on;")
% % % c_eval("plot3(LocPoint(1),-RR_max,LocPoint(3),'*' ,'color','m','linewidth',2); hold on;")
% % % c_eval("plot3(LocPoint(1),LocPoint(2),-RR_max,'*' ,'color','m','linewidth',2); hold on;")
 
% Q Errorbar
% % % n = [0,1,0];
% % % Center = [LocPoint(1),-RR_max,LocPoint(3)];
% % % R = resQ/max(resQ)*RR_mean/50;
% % % plotCircle3(Center,R(1),n,'w');
% % % plotCircle3(Center,R(2),n,'w');
% % % R0 = log10(Q)/max(resQ)*RR_mean/50;
% % % plotCircle3(Center,R0,n,'m');
 
% C1 = [LocPoint(1),-RR_max,LocPoint(3)];
% C2 = C1 + [0,1,0]*abs(Q)/max(abs(resQ))*RR_mean/4;
% [Cylinder1,EndPlate1,EndPlate2] = cylinder3(C1,C2,RR_mean/50,50,[1,0.58,0.8],1,0);
% alpha(Cylinder1,0.8);
% resC1 = C1 + [0,1,0]*abs(resQ(1))/max(abs(resQ))*RR_mean/4;
% resC2 = C1 + [0,1,0]*abs(resQ(2))/max(abs(resQ))*RR_mean/4;
% [Cylinder2,EndPlate3,EndPlate4] = cylinder3(resC1,resC2,1.2*RR_mean/50,50,[0.58,0.8,1],1,0);
% alpha(Cylinder2,0.6);
end

%% Reconstruction
if PlotFlag ==3  
% % % ss = 2;
% % % % step = (RR_max-0.1)/ss;
% % % % radial = linspace(0,RR_max,ss+1);
radial = [0.25*RR_max,0.75*RR_max];
phi = 0:45:360;
z = linspace(-0.25*RR_max,0.25*RR_max,3);
[~,focus_points] = FOTE_fp_select(radial, phi, z, 'Cylindrical');
% % % % focus_points = [];
% % % % c_eval('focus_points = [focus_points;R?_gsm(tempidx_R,2:4)];');
% % % % focus_points = [focus_points;[focus_points(:,1),focus_points(:,2),-focus_points(:,3)]];
% % % % focus_points = [focus_points;[0,5,0;0,-5,0;5,0,0;-5,0,0]];
% for ii = 1:4
%     c_eval(['focus_points = [focus_points;0.5*(R?_gsm(tempidx_R,2:4)+R',num2str(ii),'_gsm(tempidx_R,2:4))];'],ii:4)
% end
% % % Monopole_res = Monopole_reconstruction(curlB, Q, focus_points, RR_max, [-45,20], 'cline',0,60);
%%
% Q = Q*(units.RE^2/1e6);
% c_eval('Bproj? = R?_gsm(tempidx_R,2:4)*dot(R?_gsm(tempidx_R,2:4),B?_gsm(tempidx_B?,2:4))/(norm(R?_gsm(tempidx_R,2:4))*norm(R?_gsm(tempidx_R,2:4)));')
% c_eval('Bres? = B?_gsm(tempidx_B?,2:4)-Bproj?;');
c_eval('R?_2 = sum(R?_gsm(tempidx_R,2:4).^2,2);');
c_eval('Bres? = Q/(4*pi*R?_2);');
c_eval('Bres? = [Bres?*R?_gsm(tempidx_R,2)/sqrt(R?_2),Bres?*R?_gsm(tempidx_R,3)/sqrt(R?_2),Bres?*R?_gsm(tempidx_R,4)/sqrt(R?_2)];');
c_eval('Bres? = B?_gsm(tempidx_B?,2:4)-Bres?;');
% curlB = c_4_grad('R?_gsm(tempidx_R,2:4)','Bres?','grad');
% curlB=reshape(curlB,3,3);
curlB = curlB/(units.RE/1e3);
% ss = 50;
% step = RR_max/ss;
% [gridX,gridY,gridZ] = meshgrid([-RR_max:step:RR_max],[-RR_max:step:RR_max],[-RR_max:step:RR_max]);
% gridB = cell(size(gridX));gridBt = zeros(size(gridX));gridBX = gridBt;gridBY = gridBX; gridBZ = gridBX;
% for i = 1:numel(gridX)
%     [gridBX(i),gridBY(i),gridBZ(i),gridBt(i)] = CalB(gridX(i),gridY(i),gridZ(i),...
%         Q,curlB,R1_gsm(tempidx_R,2:4),B1_gsm(tempidx_B1,2:4));
%     gridB{i} = [gridBX(i),gridBY(i),gridBZ(i)];
% end
% lenX = length(-RR_max:step:LocPoint(1));
% lenY = length(-RR_max:step:LocPoint(2));
% lenZ = length(-RR_max:step:LocPoint(3));
% strX = '[round(lenY*3/5)+1,lenY+1,1*round((2*ss+1-lenY)/5)+lenY+1]';
% strY = '[round(lenX*3/5)+1,lenX+1,2*round((2*ss+1-lenX)/5)+lenX+1]';
% strZ = '[round(lenZ*2/5)+1,lenZ+1,1*round((2*ss+1-lenZ)/5)+lenZ+1]';
% for i = ['X','Y','Z']
%     eval(['st',i,' = grid',i,'(',strX,',',strY,',',strZ,');'])
% end

% Monopole_reconstruction(curlB, Q,[0,0,0],R1_gsm(tempidx_R,2:4),B1_gsm(tempidx_B1,2:4),1e2, RR_max, [-45,20], 'cline',0,60);
Monopole_res = Monopole_reconstruction(curlB, Q, focus_points,R1_gsm(tempidx_R,2:4),B1_gsm(tempidx_B1,2:4),1e5, [0.5,1.5*RR_max], [-45,20], 'cline',0,80);

% % % [stX,stY,stZ] = meshgrid(linspace(-RR_max*0.6,RR_max*0.6,3),linspace(-RR_max*0.6,RR_max*0.6,3),linspace(-RR_max*0.6,RR_max*0.6,3));
% st = [R1_gsm(tempidx_R,2:4);R2_gsm(tempidx_R,2:4);R3_gsm(tempidx_R,2:4);R4_gsm(tempidx_R,2:4);0,0,0];
% lines = streamline(stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,gridX,gridY,gridZ));
% lines = streamline(stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ));
% % % % lines = stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ);
% % % % % streamline(lines);
% % % % for i = 1:length(lines)
% % % %     lines_len = length(lines{i});
% % % %     Bt = ones(lines_len,1);
% % % %     for j = 1:lines_len
% % % %         [~,~,~,Bt(j)] = CalB(lines{i}(j,1),lines{i}(j,2),lines{i}(j,3),...
% % % %         Q,curlB,R1_gsm(tempidx_R,2:4),B1_gsm(tempidx_B1,2:4));
% % % %     end
% % % %     cline(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),Bt,0,80,jet); hold on;
% % % % end
% caxis([min(gridBt,[],'all'),max(gridBt,[],'all')])
% colorbar;caxis([0,30]);
% lines = streamslice(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ,0.05,'noarrows');
% set(lines,'color',[0.8,0.58,1]);
% for i = 1:length(lines)
%     xLine = lines(i).XData;
%     yLine = lines(i).YData;
%     zLine = lines(i).ZData;
%     if length(xLine) >= 3
%         LineVec = [xLine(end)-xLine(end-1),yLine(end)-yLine(end-1),zLine(end)-zLine(end-1)];
%         NormPara = sqrt(sum(LineVec.^2,2)/RR_max);
% %         LineVec = LineVec*NormPara;
%         quiver3(xLine(end),yLine(end),zLine(end),LineVec(1),LineVec(2),LineVec(3),0,'MaxHeadsize',5,'LineWidth',1,'color',[0.8,0.58,1]); hold on;
%     end
% end
% view([45,20])
end
%%
if PlotFlag == 2
subplot(2,2,3)
QErrorbar(Q,resQ,'m')

Loc_X = [Loc{1}(1),Loc{1}(1)]; Loc_Y = [Loc{1}(2),Loc{1}(2)]; Loc_Z = [Loc{1}(3),Loc{1}(3)];
c_eval('Loc_X = [min(Loc{?}(1),Loc_X(1)),max(Loc{?}(1),Loc_X(2))];',2:4);
c_eval('Loc_Y = [min(Loc{?}(1),Loc_Y(1)),max(Loc{?}(1),Loc_Y(2))];',2:4);
c_eval('Loc_Z = [min(Loc{?}(1),Loc_Z(1)),max(Loc{?}(1),Loc_Z(2))];',2:4);
R_X = [R1_gsm(tempidx_R,2),R1_gsm(tempidx_R,2)]; R_Y = [R1_gsm(tempidx_R,3),R1_gsm(tempidx_R,3)];R_Z = [R1_gsm(tempidx_R,4),R1_gsm(tempidx_R,4)];
c_eval('R_X = [min(R?_gsm(tempidx_R,2),R_X(1)),max(R?_gsm(tempidx_R,2),R_X(2))];',2:4);
c_eval('R_Y = [min(R?_gsm(tempidx_R,3),R_Y(1)),max(R?_gsm(tempidx_R,3),R_Y(2))];',2:4);
c_eval('R_Z = [min(R?_gsm(tempidx_R,4),R_Z(1)),max(R?_gsm(tempidx_R,4),R_Z(2))];',2:4);

% QErrorbar(LocPoint,Loc_Z,Loc_X,R_Z,R_X,'m');
end
%%

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

% figname = [OutputDir,'OverviewFig\',NameTags{TDT}(2:end-2),'-Configuration'];
% print(gcf, '-dpng', [figname '.png']);