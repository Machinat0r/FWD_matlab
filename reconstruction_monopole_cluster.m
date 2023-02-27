clear;clc
close all
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
%% load data
cd C:\Cluster
ic=1:4;
% 
% Tsta='2003-12-30T18:55:00.00Z';
% Tend='2003-12-30T18:56:00.00Z';%15,268
Tsta = '2004-01-17T19:31:00.000Z';
Tend = '2004-01-17T19:32:00.000Z';
% Tsta='2004-01-21T08:56:00.00Z';
% Tend='2004-01-21T08:58:00.00Z';


tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH

% Tsta2='2003-12-30T18:55:06.000Z';
% Tend2='2003-12-30T18:55:09.800Z';
Tsta2 = '2004-01-17T19:31:45.765Z';
Tend2 = '2004-01-17T19:31:45.925Z';
% Tsta2='2004-01-21T08:56:59.275Z';
% Tend2='2004-01-21T08:56:59.435Z';
% tint2=[iso2epoch(Tsta2) iso2epoch(Tend2)];

try
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
%     c_eval("caa_load('C?_CP_FGM_FULL',Tsta,Tend);",ic);
%     c_eval("caa_load('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
catch
    
    %    Magnetic fields
c_eval("caa_download(tint,'C?_CP_FGM_FULL')",ic);
c_eval("caa_download(tint,'C?_CP_AUX_POSGSE_1M')",ic);  % position & velocity for each sc
    c_eval("caa_load_changed_by_fwd('C?_CP_FGM_FULL',Tsta,Tend);",ic);
    c_eval("caa_load_changed_by_fwd('C?_CP_AUX_POSGSE_1M',Tsta,Tend);",ic);
end
% caa_download(tint,'CL_SP_AUX')% position,attitude.. for all sc
% caa_download(tint,'C2_CP_FGM_FULL');
% caa_download(tint,'C4_CP_FGM_FULL');

% caa_load C  %load data from datebase form C

%background magnetic field

% % % dobjname=irf_ssub('C?_CP_FGM_FULL',ic); 
% % % varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic); 
c_eval('B?_gse=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',ic); 
%smooth
% for rec = 2:4
% c_eval('B?_gse(:,rec) = smooth(B?_gse(:,rec),9);')
% end

c_eval('B?_gsm = irf_gse2gsm(B?_gse);',ic);
c_eval('B?=irf_abs(B?_gsm);',ic);
c_eval('B? = irf_resamp(B?,B1);',2:4);
c_eval('B?_gsm = B?;')

c_eval('R?_gse = c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');',ic);
c_eval('R? = irf_gse2gsm(R?_gse);',ic);
c_eval('R? = irf_resamp(R?,B1);')
CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

try
    [~,Tsta_id] = sort(abs(B1(:,1)-iso2epoch(Tsta2)));
    Tsta_id = Tsta_id(1);
    [~,Tend_id] = sort(abs(B1(:,1)-iso2epoch(Tend2)));
    Tend_id = Tend_id(1);
    c_eval('B? = B?(Tsta_id:Tend_id,:);');
    c_eval('B?_gsm = B?;');
    c_eval('R? = R?(Tsta_id:Tend_id,:);');
    tint = [iso2epoch(Tsta2) iso2epoch(Tend2)];
catch
    disp('no tint2,calculate all over tint')
end
%% PI
% PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
% [PI,~,~,~,~]=poincare_index_and_solid_angle(B1_gsm(2:4),B2_gsm(2:4),B3_gsm(2:4),B4_gsm(2:4));
PI(abs(PI)<0.5) = 0;
PI(PI>=0.5) = 1;
PI(PI<=-0.5) = -1;
%% divB
gradB=c_4_grad('R?','B?(:,1:4)','grad');
divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)],2)];      %% 未归一化散度
PI(isnan(divB(:,2)))=0;
% tint = [iso2epoch(Tsta2) iso2epoch(Tend2)];
%% FOTE 误差
LocPoint = zeros(length(PI),3)*nan;
LocRes = cell(length(PI),1);
Q = zeros(length(PI),1)*nan;
resQ = cell(length(PI),1);

Qerror = ones(length(PI),1)*1000;
Locerror = ones(length(PI),1)*200;

PI_id = find(PI~=0)';
% for i = 1:length(B1)
for i = PI_id
clc;
disp(['current calculate:',num2str(i),'/',num2str(length(PI))]);

RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(i,2),R',num2str(ii),'(i,3),R',num2str(ii),'(i,4);',...
    'R?(i,2),R?(i,3),R?(i,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
% MultiPower = ceil(max([log10(abs(R1(i,2:4))),log10(abs(R2(i,2:4))),log10(abs(R3(i,2:4))),log10(abs(R4(i,2:4)))]));
% solve
[Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?',i,i*sign(divB(i,2)),10,1);
id = nchoosek(1:6,2);
c_eval('tempd? = irf_abs(LocRes{i}(id(?,1),:)-LocRes{i}(id(?,2),:));',1:15)
tempd = [];
c_eval('tempd = [tempd,tempd?(4)/RR_mean];',1:15);
dLoc(i,:) = tempd;
%error
Qerror(i) = abs(100*std(resQ{i})/Q(i));
tri_a = delaunayTriangulation([R1(i,2:4);R2(i,2:4);R3(i,2:4);R4(i,2:4)]);
[~,volume_a] = convexHull(tri_a);
tri = delaunayTriangulation(LocRes{i});%%delaunay三角剖分
    if size(tri.Points,1)==1
        volume = 0;
    else
    [~,volume] = convexHull(tri);%%计算多面体体积
    end
Locerror(i) = 100*volume/volume_a;
end
%% Index id
% length(find(Locerror<50 & Qerror<=1000))
tempidx_B1 = find(mean(dLoc,2) == min(mean(dLoc,2)));
% tempidx_B1 = 80;
%     [~,tempidx_B] = max(abs(divB(:,2)));
c_eval('[~,tempidx_B?] = sort(abs(B?_gsm(:,1)-B1_gsm(tempidx_B1,1)));',2:4);
c_eval('tempidx_B? = tempidx_B?(1);')
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
c_eval('R?(:,2:4) = R?(:,2:4)-LocPoint(tempidx_B1,:);');
LocRes{tempidx_B1} = LocRes{tempidx_B1}-LocPoint(tempidx_B1,:);
LocPoint = LocPoint - LocPoint(tempidx_B1,:);
%% Location 
cor = 'krgb';
c_eval("plot3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),'s' ,'color',cor(?),'linewidth',5);hold on;")

%% Tetrahedron configuration
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(tempidx_R,2),R',num2str(ii),'(tempidx_R,3),R',num2str(ii),'(tempidx_R,4);',...
    'R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
c_eval(['plot3(RR',num2str(ii),'?(:,1),RR',num2str(ii),'?(:,2),RR',num2str(ii),'?(:,3),''--k '');hold on;'],ii+1:4)
end
RR_mean = RR_mean(4)/6;
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
% irf_legend(gca,{'C1','C2','C3','C4'},[0.97 0.92]);
% irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
xlabel('e_1 [km]','fontsize',12);
ylabel('e_2 [km]','fontsize',12);
zlabel('e_3 [km]','fontsize',12);

%% Quiver
c_eval('B?_gsm = irf_abs(B?_gsm);');
maxB = max([B1_gsm(tempidx_B1,5),B2_gsm(tempidx_B1,5),B3_gsm(tempidx_B1,5),B4_gsm(tempidx_B1,5)]);
% c_eval("quiver3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),RR_mean*B?_gsm(tempidx_B1,2)/maxB,RR_mean*B?_gsm(tempidx_B1,3)/maxB,RR_mean*B?_gsm(tempidx_B1,4)/maxB,'color',cor(?));hold on;")

axis equal 
%% Reconstruction
% RR_max = 1.2*max(abs([R1(tempidx_R,2:4),R2(tempidx_R,2:4),R3(tempidx_R,2:4),R4(tempidx_R,2:4)]));
RR_max = 250;
ss = 25;
step = 100/ss;
[gridX,gridY,gridZ] = meshgrid([-RR_max:step:RR_max],[-RR_max:step:RR_max],[-RR_max:step:RR_max]);
gridB = cell(size(gridX));gridBt = zeros(size(gridX));gridBX = gridBt;gridBY = gridBX; gridBZ = gridBX;
for i = 1:numel(gridX)
    [gridBX(i),gridBY(i),gridBZ(i),gridBt(i)] = calculateB(gridX(i),gridY(i),gridZ(i),Q(tempidx_B1),LocPoint(tempidx_B1,:));
    gridB{i} = [gridBX(i),gridBY(i),gridBZ(i)];
end


ls = [-1,0,1];
[stX,stY,stZ] = meshgrid(ls,ls,ls);
lines = stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ,[1,80]);
for i = 1:length(lines)
    lines_len = size(lines{i},1);
    Bt = ones(lines_len,1);
    for j = 1:lines_len
        [~,~,~,Bt(j)] = calculateB(lines{i}(j,1),lines{i}(j,2),lines{i}(j,3),Q(tempidx_B1),LocPoint(tempidx_B1,:));
    end
    try
    cline(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),Bt,0,2,'jet',log(Bt)/log(5)); hold on;
    catch
        disp(['error at lines',num2str(i)])
    end
end

ls = [-2:1:2];
[stX,stY,stZ] = meshgrid(ls,ls,ls);
lines = stream3(gridX,gridY,gridZ,gridBX,gridBY,gridBZ,stX,stY,stZ,[1,80]);
for i = 1:length(lines)
    lines_len = size(lines{i},1);
    Bt = ones(lines_len,1);
    for j = 1:lines_len
        [~,~,~,Bt(j)] = calculateB(lines{i}(j,1),lines{i}(j,2),lines{i}(j,3),Q(tempidx_B1),LocPoint(tempidx_B1,:));
    end
    try
    cline(lines{i}(:,1),lines{i}(:,2),lines{i}(:,3),Bt,0,2,'jet',0.5); hold on;
    catch
        disp(['error at lines',num2str(i)])
    end
end

% cmap = othercolor('Blues9');
cmap = colormap('jet');
cmap = cmap(1:240,:);
c = colorbar;
caxis([0,2])
set(c,'YTick',[0,1,2]);
set(c,'YTickLabel',{'10^0','10^1','10^2'}) ;
% set(c,'Color','w')
set(get(c,'ylabel'),'string','|B| [nT]','fontsize',12);
%% adjust the position
axis equal 
box on
axis([-250 250 -250 250 -250 250])
view(40,42)
% set(gca,'Color','k')
% set(gca,'XColor','w');set(gca,'YColor','w');set(gca,'ZColor','w')
% set(gcf,'color','k')
% axis off

set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')
figname = ['C:\Users\fwd\Desktop\Ti~mor~\M\magnetic_monopole\Figure\new\Figure4\' ,'20020218_reconstruction'];
% print(gcf, '-dpdf', [figname '.pdf']);