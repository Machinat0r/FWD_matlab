clear;clc;close all
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
% % % global ParentDir 
% % % ParentDir = '/Volumes/172.17.190.41/Data/MMS/'; 
% % % DownloadDir = '/Users/fwd/Documents/MATLAB/MMS/';
% % % TempDir = [DownloadDir,'temp/'];mkdir(TempDir);
% cd /Volumes/FWD-WorkDisk/Cluster
% cd /Volumes/172.17.190.41/Data/Cluster/
cd /Users/fwd/Documents/MATLAB/Cluster
ic=1:4;
% 
% Tsta='2002-02-18T08:41:00.00Z';
% Tend='2002-02-18T08:42:00.00Z';%1,114
% Tsta='2003-12-30T18:55:00.00Z';
% Tend='2003-12-30T18:56:00.00Z';%7,268
Tsta = '2004-01-17T19:31:00.000Z';
Tend = '2004-01-17T19:32:00.000Z';
% % % Tsta = '2004-01-21T08:56:00.000Z';
% % % Tend = '2004-01-21T08:58:00.000Z';

% % % Tsta = '2000-01-01T00:00:00.000Z';
% % % Tend = '2020-01-01T00:00:00.000Z';


% Tsta='2004-01-01T11:10:00.00Z';
% Tend='2004-01-01T11:11:00.00Z';
% % % Tsta='2004-01-21T08:56:00.00Z';
% % % Tend='2004-01-21T08:58:00.00Z';
% Tsta='2004-05-30T20:22:00.00Z';
% Tend='2004-05-30T20:23:00.00Z';
% Tsta='2003-12-29T13:40:00.00Z';
% Tend='2003-12-29T13:41:00.00Z';

% Tsta='2004-05-03T05:02:00.00Z';
% Tend='2004-05-03T05:03:00.00Z';


tint=[iso2epoch(Tsta) iso2epoch(Tend)]; %ISO time to ISDAT EPOCH
% 
% Tsta2='2002-02-18T08:41:07.000Z';
% Tend2='2002-02-18T08:41:09.000Z';
% Tend2='2002-02-18T08:41:10.000Z';
% Tsta2='2002-06-13T18:41:30.00Z';
% Tend2='2002-06-13T18:41:50.00Z';
% Tsta2='2003-12-17T01:20:05.00Z';
% Tend2='2003-12-17T01:20:09.00Z';
% % % Tsta2='2003-12-30T18:55:08.560Z';
% % % Tend2='2003-12-30T18:55:08.740Z';
% Tsta2='2003-12-30T18:55:06.800Z';
% Tend2='2003-12-30T18:55:09.800Z';
Tsta2 = '2004-01-17T19:31:45.765Z';
Tend2 = '2004-01-17T19:31:45.925Z';
% Tsta2 = '2004-01-17T19:31:44.900Z';
% Tend2 = '2004-01-17T19:31:46.500Z';
% Tsta2='2004-05-03T05:02:28.00Z';
% Tend2='2004-05-03T05:02:30.00Z';

% Tsta2='2004-01-01T11:10:03.500Z';
% Tend2='2004-01-01T11:10:04.500Z';
% % % Tsta2='2004-01-21T08:56:58.200Z';
% % % Tend2='2004-01-21T08:57:00.200Z';
% Tsta2='2004-05-30T20:22:30.00Z';
% Tend2='2004-05-30T20:22:32.00Z';
% Tsta2='2003-12-29T13:40:25.00Z';
% Tend2='2003-12-29T13:40:26.600Z';
% Tsta2='2003-12-29T13:40:26.150Z';
% Tend2='2003-12-29T13:40:26.320Z';
% % % Tsta2='2004-01-21T08:56:59.275Z';
% % % Tend2='2004-01-21T08:56:59.435Z';
tint2=[iso2epoch(Tsta2) iso2epoch(Tend2)];

% % % caa_download(tint,'C?_CP_EFW_L3_E3D_INERT');
% % % caa_download(tint,'C?_CP_EFW_L?_E');
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

% c_eval('B?_gse = irf_gse2gsm(B?_gse);',ic);
c_eval('B?=irf_abs(B?_gse);',ic);
c_eval('B? = irf_resamp(B?,B1);',2:4);

c_eval('R?_gse = c_caa_var_get(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''mat'');',ic);
c_eval('R? = R?_gse;')
% c_eval('R? = irf_gse2gsm(R?_gse);',ic);
c_eval('R? = irf_resamp(R?,B1);')
CenterPoint = (R1(:,2:4)+R2(:,2:4)+R3(:,2:4)+R4(:,2:4))/4;
c_eval('R?(:,2:4) = R?(:,2:4)-CenterPoint;');

try
    [~,Tsta_id] = sort(abs(B1(:,1)-iso2epoch(Tsta2)));
    Tsta_id = Tsta_id(1);
    [~,Tend_id] = sort(abs(B1(:,1)-iso2epoch(Tend2)));
    Tend_id = Tend_id(1);
    c_eval('B? = B?(Tsta_id:Tend_id,:);');
    c_eval('B?_gse = B?;');
    c_eval('R? = R?(Tsta_id:Tend_id,:);');
    tint = [iso2epoch(Tsta2) iso2epoch(Tend2)];
catch
    disp('no tint2,calculate all over tint')
end
%% PI
% PI=c_fgm_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
PI=c_4_poincare_index(B1(:,2:4),B2(:,2:4),B3(:,2:4),B4(:,2:4));
% [PI,~,~,~,~]=poincare_index_and_solid_angle(B1_gse(2:4),B2_gse(2:4),B3_gse(2:4),B4_gse(2:4));
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
for i = 1:size(B1,1)
% for i = PI_id
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
[Q(i),resQ{i},LocPoint(i,:),LocRes{i}] = CalError('R?','B?_gse',i,i*sign(divB(i,2)),RR_mean,10);
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
meand = mean(dLoc,2);
%% Coordinate
% % % c_eval('R?(:,2:4) = R?(:,2:4)-LocPoint;');
maxR = nan*zeros(1,3); minR = nan*zeros(1,3);
c_eval('maxR = max(maxR,max(R?(:,2:4),[],1));',1:4);
c_eval('minR = min(minR,min(R?(:,2:4),[],1));',1:4);
GridX = linspace(minR(1),maxR(1),20);
GridY = linspace(minR(2),maxR(2),20);
GridZ = linspace(minR(3),maxR(3),20);
[MGX, MGY, MGZ] = meshgrid(GridX, GridY, GridZ);
X = MGX(:);Y = MGY(:);Z = MGZ(:);
points = [X, Y, Z];

c_eval('V? = [R1(?,2:4); R2(?,2:4); R3(?,2:4); R4(?,2:4)];',1:size(B1,1));
c_eval('Tri? = delaunayTriangulation(V?);',1:size(B1,1));
c_eval('points? = points(pointLocation(Tri?,points) == 1,:);',1:size(B1,1));

for i = 1:length(PI)
    LocRes{i} = LocRes{i}-LocPoint(i,:);
end
LocPoint = LocPoint - LocPoint;

%% cal B
for i = 1:length(R1)
c_eval('[tempBx?,tempBy?,tempBz?,tempBt?]=calculateB(R?(i,2),R?(i,3),R?(i,4),Q(i),LocPoint(i,:));')
c_eval('Bm?(i,:) = [B1(i,1),tempBx?,tempBy?,tempBz?,tempBt?];')
end
gradBm=c_4_grad('R?','Bm?(:,1:4)','grad');
divBm=[gradBm(:,1) sum([gradBm(:,2) gradBm(:,6) gradBm(:,10)],2)]; 
divErr = [divB(:,1) 100*abs(divB(:,2)-divBm(:,2))./abs(divBm(:,2))]; %|deltB-deltBm|/deltB
%% B error
c_eval('Bcross? = irf_dot(B?(:,2:4),Bm?(:,2:4));')
% % % c_eval('Btheta? = acosd(Bcross?./(B?(:,5).*Bm?(:,5)));')
c_eval('R? = irf_abs(R?);')
c_eval('theta? = acosd(Bcross?./(B?(:,5).*Bm?(:,5)));')
% % Rres = cellfun(@irf_abs,LocRes,'UniformOutput',false);
% % Rres = cellfun(@(x) [max(x([1,4,5],4)),max(x([1,2,6],4)),max(x([2,3,4],4)),max(x([3,5,6],4)),mean(x(:,4))],Rres,'UniformOutput',false);
% % Rres = cell2mat(Rres);
% % c_eval('Btheta? = 100*theta?./RR_mean.*Rres(:,?);')
c_eval('Btheta? = theta?./RR_mean.*R?(:,5);')
% c_eval('Btheta? = asin(sqrt(1-(Bcross?./(B?(:,5).*Bm?(:,5))).^2))*180/pi;')
% % c_eval('Blength? = abs(B?(:,5)-Bm?(:,5))./Bm?(:,5)./RR_mean.*Rres(:,?);')
c_eval('Blength? = abs(B?(:,5)-Bm?(:,5))./Bm?(:,5)./RR_mean.*R?(:,5);')
Blength = 0.25*(Blength1+Blength2+Blength3+Blength4);
c_eval('BErr? = 100*(1-(1-Btheta?).*(1-Blength?));')
BErr = mean([BErr1,BErr2,BErr3,BErr4],2);
%% VSH
c_eval('R?_grid = cell(size(points1,1));', 1:4);
coeff_grid = cell(size(points1,1)); res_grid = cell(size(points1,1));
c_eval('rate1?_grid = zeros(size(B?,1),size(points1,1));', 1:4);
c_eval('rate2?_grid = zeros(size(B?,1),size(points1,1));', 1:4);
for i_grid = 1:size(points1,1)
c_eval('R?_grid{i_grid} = [R?(:,1), R?(:,2:4) - points1(i_grid,:)];',1:4);
[coeff, res] = VSH_Expand('R?_grid{i_grid}', 'B?(:,1:4)', [0,0,0], 1);
coeff_grid{i_grid} = coeff; res_grid{i_grid} = res;


% % % [coeff, res] = VSH_Expand('R?', 'B?(:,1:4)', [0,0,0], 1);
% coeff0 = coeff([1,5,9,13,14,15]);
% coeff1 = coeff([2,3,4,6,7,8,10,11,12]);
alpha_00 = coeff(:,1); beta_00 = coeff(:,5); gamma_00 = coeff(:,9);
alpha_12 = coeff(:,2); beta_12 = coeff(:,6); gamma_12 = coeff(:,10);
alpha_10 = coeff(:,3); beta_10 = coeff(:,7); gamma_10 = coeff(:,11);
alpha_11 = coeff(:,4); beta_11 = coeff(:,8); gamma_11 = coeff(:,12);
c_eval('[r?, theta?, phi?] = Coor_Trans(R?(:,2:4));', 1:4);
[r1, r2, r3, r4, ~] = Normal_R(r1, r2, r3, r4);
for temp_ic = 1:4
temp_ic = num2str(temp_ic);
eval(['Br0',temp_ic,' = Cal_Br0(r',temp_ic,', theta',temp_ic,', alpha_00);']);
eval(['Br1',temp_ic,' = Cal_Br1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', alpha_12, alpha_10, alpha_11);']);
eval(['Bt1',temp_ic,' = Cal_Bt1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', beta_12, beta_10, beta_11, gamma_12, gamma_11);']);
eval(['Bp1',temp_ic,' = Cal_Bp1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', beta_12, beta_11, gamma_12, gamma_10, gamma_11);']);
end
c_eval('[~, B?_sph] = Coor_Trans(R?(:,2:4), B?(:,2:4));', 1:4);
c_eval('B0? = irf_abs([B?(:,1), Br0?, zeros(size(Br0?)), zeros(size(Br0?))]);', 1:4);
c_eval('B1? = irf_abs([B?(:,1), Br1?, Bt1?, Bp1?]);', 1:4);
c_eval('B2? = irf_abs([B?(:,1), B?_sph(:,1) - Br0? - Br1?, B?_sph(:,2) - Bt1?, B?_sph(:,3) - Bp1?]);', 1:4);
c_eval('rate1? = B1?(:,5) ./ B0?(:,5);', 1:4);
c_eval('rate2? = B2?(:,5) ./ B0?(:,5);', 1:4);

c_eval('rate1?_grid(:,i_grid) = rate1?;',1:4);
c_eval('rate2?_grid(:,i_grid) = rate2?;',1:4);
end
%% Init Figure 2
figure(2)
set(gcf,'PaperUnits','centimeters')
xSize = 70; ySize = 60; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 10 xSize*coef ySize*coef])
%% Index id
% tempidx_B1 = find(mean(dLoc,2) == min(mean(dLoc,2)));
tempidx_B1 = 3;

% for tempidx_B1 = tempidx_B'
%     [~,tempidx_B] = max(abs(divB(:,2)));
c_eval('[~,tempidx_B?] = sort(abs(B?_gse(:,1)-B1_gse(tempidx_B1,1)));',2:4);
c_eval('tempidx_B? = tempidx_B?(1);')
[~,tempidx_R] = sort(abs(R1(:,1)-B1(tempidx_B1,1)));
tempidx_R = tempidx_R(1);
%% Location 
cor = 'krgb';
c_eval("plot3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),'s' ,'color',cor(?),'linewidth',5);hold on;")
%% Tetrahedron configuration
RR_mean = zeros(1,4);
for ii = 1:3 
c_eval(['RR',num2str(ii),'?=[R',num2str(ii),'(tempidx_R,2),R',num2str(ii),'(tempidx_R,3),R',num2str(ii),'(tempidx_R,4);',...
    'R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4)];'],ii+1:4);  %% ♥
c_eval(['RR_mean=RR_mean+irf_abs(RR',num2str(ii),'?(2,:)-RR',num2str(ii),'?(1,:));'],ii+1:4);  
end
RR_mean = RR_mean(4)/6;
plot3(RR12(:,1),RR12(:,2),RR12(:,3),'--k');hold on;  plot3(RR13(:,1),RR13(:,2),RR13(:,3),'--k');hold on;  
plot3(RR14(:,1),RR14(:,2),RR14(:,3),'--k');hold on;  plot3(RR23(:,1),RR23(:,2),RR23(:,3),'--k');hold on;  
plot3(RR34(:,1),RR34(:,2),RR34(:,3),'--k');hold on;  plot3(RR24(:,1),RR24(:,2),RR24(:,3),'--k');hold on;  
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'Cluster: 1','2','3','4'},[0.97 0.92]);
irf_legend(gca,{['Separation Distance:',num2str(roundn(RR_mean,-1)),'km']},[0.05 0.92])
xlabel('x [km]','fontsize',12);
ylabel('y [km]','fontsize',12);
zlabel('z [km]','fontsize',12);

%% Quiver
c_eval('B?_gse = irf_abs(B?_gse);');
maxB = max([B1_gse(tempidx_B1,5),B2_gse(tempidx_B1,5),B3_gse(tempidx_B1,5),B4_gse(tempidx_B1,5)]);
% c_eval("quiver3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),RR_mean*B?_gse(tempidx_B1,2)/maxB,RR_mean*B?_gse(tempidx_B1,3)/maxB,RR_mean*B?_gse(tempidx_B1,4)/maxB,'color',cor(?));hold on;")
% c_eval("quiver3(R?(tempidx_R,2),R?(tempidx_R,3),R?(tempidx_R,4),RR_mean*Bm?(tempidx_B1,2)/maxB,RR_mean*Bm?(tempidx_B1,3)/maxB,RR_mean*Bm?(tempidx_B1,4)/maxB,'linewidth',1.5,'color',cor(?));hold on;")
%% Loc res
% % % plotPolyhedron(LocRes{tempidx_B1}(:,1),LocRes{tempidx_B1}(:,2),LocRes{tempidx_B1}(:,3),'#4DBEEE',0.1);
% idx = [166:198];
% plot3(LocPoint(idx,1),LocPoint(idx,2),LocPoint(idx,3),'*','color','#FFBAF1');
% line(LocPoint(173:195,1),LocPoint(173:195,2),LocPoint(173:195,3),'color','#FFBAF1');
% % % plot3(LocPoint(tempidx_B1,1),LocPoint(tempidx_B1,2),LocPoint(tempidx_B1,3),'o','linewidth',5,'color','#0072BD');hold on;
% plot3(LocRes{tempidx_B1}(:,1),LocRes{tempidx_B1}(:,2),LocRes{tempidx_B1}(:,3),'*','color','#FFB8CE');
% plot3(LocRes{tempidx_B1}([1;4;5],1),LocRes{tempidx_B1}([1;4;5],2),LocRes{tempidx_B1}([1;4;5],3),'*','linewidth',1,'color','k')
% text(LocRes{tempidx_B1}([1;4;5],1),LocRes{tempidx_B1}([1;4;5],2),LocRes{tempidx_B1}([1;4;5],3),{'(1,2)','(1,3)','(1,4)'})
% % % iid = [2,3,4];
% % % plot3(LocRes{tempidx_B1}(iid,1),LocRes{tempidx_B1}(iid,2),LocRes{tempidx_B1}(iid,3),'*','linewidth',1,'color','g')
% % % text(LocRes{tempidx_B1}(iid,1),LocRes{tempidx_B1}(iid,2),LocRes{tempidx_B1}(iid,3),{'(2,3)','(3,4)','(1,3)'})
%% Grid rate
rate1_grid = mean([rate11_grid(tempidx_B1,:);rate12_grid(tempidx_B1,:);rate13_grid(tempidx_B1,:);rate14_grid(tempidx_B1,:)]);
rate2_grid = mean([rate21_grid(tempidx_B1,:);rate22_grid(tempidx_B1,:);rate23_grid(tempidx_B1,:);rate24_grid(tempidx_B1,:)]);

cmp = othercolor('BuOr_10'); cmp = flip(cmp);
colormap(cmp);
s = scatter3(points1(:,1), points1(:,2), points1(:,3), 50*ones(size(points1,1),1),rate1_grid*100, 'filled'); hold on;
alpha(s,1./rate1_grid)
clim([67,75])
view(27.662,37.777)
% view(36.537021187953911,-29.347955406451824)
c = colorbar;
c.Label.String = '\Re_2 [%]';
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')


%% two-step component function
function Br0 = Cal_Br0(r, theta, alpha_00)
Br0 = r.^-2 .* alpha_00 .* Cal_SPH(0, 0, theta)';
end

function Bt0 = Cal_Bt0(r, theta, beta_00)
Bt0 = r.^-1 .* beta_00 .* Cal_dPlm(0, 0, theta)';
end

function Bp0 = Cal_Bp0(r, theta, gamma_00)
Bp0 = r.^-1 .* gamma_00 .* Cal_dPlm(0, 0, theta)';
end

function Br1 = Cal_Br1(r, theta, phi, alpha_12, alpha_10, alpha_11)
Br1 = r.^-3 .* (alpha_12 .* cos(-phi) .* Cal_SPH(1, -1, theta)'...
    + alpha_10 .* Cal_SPH(1, 0, theta)' + ...
    alpha_11 .* cos(phi) .* Cal_SPH(1, 1, theta)');
end

function Bt1 = Cal_Bt1(r, theta, phi, beta_12, beta_10, beta_11, gamma_12, gamma_11)
Bt1 = r.^-2 .* (gamma_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta)' + ...
    beta_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta)' + ...
    beta_10 ./ 2 .* Cal_dPlm(1, 0, theta)' + ...
    gamma_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta)' + ...
    beta_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta)');
end

function Bp1 = Cal_Bp1(r, theta, phi, beta_12, beta_11, gamma_12, gamma_10, gamma_11)
Bp1 = r.^-2 .* (gamma_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta)' -...
    beta_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta)' + ...
    gamma_10 ./ 2 .* Cal_dPlm(1, 0, theta)' + ...
    gamma_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta)' -...
    beta_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta)');
end

%% Calculate Spherical Harmonic
function Ylm = Cal_SPH(l, m, theta)
% theta, phi should be rad
% without exp(i*m*phi), it is placed in the magnetic field calculation

if abs(m) > l, Ylm = 0; return; end

% unnormal legendre funciton
Plm = legendre(l, cos(theta));
if m < 0, Plm = Plm(-m + 1,:); else, Plm = Plm(m + 1,:); end

% orthonormalized coefficient 
a = (2*l+1)*factorial(l-m);
b = 4*pi*factorial(l+m);
C = sqrt(a/b);

Ylm = sqrt(2) .* C .* Plm;
end

%% Calculate Difference of Legendre Function
function dPlm = Cal_dPlm(l, m, theta)
Clm1 = 0.5 * sqrt((l + m) * (l - m + 1));
Clm2 = 0.5 * sqrt((l + m + 1) * (l - m));
Plm1 = Cal_SPH(l, m-1, theta);
Plm2 = Cal_SPH(l, m+1, theta);

dPlm = Clm1 .* Plm1 - Clm2 .* Plm2;
end
%% R normalization
function [R1, R2, R3, R4, maxR] = Normal_R(R1, R2, R3, R4)
maxR = 0;
c_eval('maxR = max(maxR,R?(:,1));',1:4);
c_eval('R? = R?./maxR;',1:4);
end