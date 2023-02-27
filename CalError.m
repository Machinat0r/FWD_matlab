function [Q,resQ,LocPoint,resLoc] = CalError(varargin)
% function [LocRes,Loc,LocPoint,Point] = CalError(varargin)
% PlotFlag choose whether to find the intersection of two tetrahedrons, if
% PlotFlag = 2, solute the intersection
%------written by Wending Fu, May.8.2022 in Beijing------------
% lb,ub未启用
%% input data
if length(varargin) == 6 || length(varargin)==12
    LocationSkew = varargin{end};
    varargin = varargin(1:end-1);
    MultiPower = varargin{end};
    varargin = varargin(1:end-1);
    idx = varargin{end};
    if sign(idx) == -1 
        lb = -inf;ub = 0;
        idx = abs(idx);
%         MultiPower = -MultiPower;
    elseif sign(idx) == 1
        lb = 0;ub = inf;  %2022.7.11,修改ModelSolve时调换了该边界；若使用LineIntersection求解，应将+-1边界反过来
    end
    varargin = varargin(1:end-1);
    idx_R = varargin{end};
    varargin = varargin(1:end-1);
end
if length(varargin) == 2
    Rs = varargin{1};
	Bs = varargin{2};
    for id=1:4
        temp = evalin('caller',irf_ssub(Rs,id));
        eval(irf_ssub('R? =temp;',id));
        temp = evalin('caller',irf_ssub(Bs,id));
        eval(irf_ssub('B? =temp;',id)); clear temp
    end
elseif length(varargin) == 8
    c_eval('R? = varargin{?};',1:4);
    c_eval('B? = varargin{?+4};',1:4);
else
    disp('Incorrect number of input parameters. See usage:')
    help CalError
    return
end
%% resample
% units = irf_units;
if ~isempty(idx)
%     c_eval('R?(:,2:4) = units.RE*R?(:,2:4);')
    c_eval('B? = irf_resamp(B?,B1);',2:4)
%     c_eval('R? = irf_resamp(R?,B1);')
    c_eval('B? = B?(idx,:);')
    c_eval('R? = R?(idx_R,:);')
end

%% residual
%% 使用2颗卫星代入模型求解
[t,res] = ModelSolve(R1,R2,R3,R4,B1,B2,B3,B4,lb,ub,MultiPower,idx,LocationSkew);
resQ = [t{1}(1);t{2}(1);t{3}(1);t{4}(1);t{5}(1);t{6}(1)];
% units = irf_units;
% resQ = resQ / units.mu0;
resLoc = [t{1}(2:4);t{2}(2:4);t{3}(2:4);t{4}(2:4);t{5}(2:4);t{6}(2:4)];
Q = mean(resQ);

LocPoint = resLoc(1,:);
try
    LocPoint = centroid3(resLoc);
    if isnan(LocPoint)
        LocPoint = mean(resLoc);
    end
catch
    LocPoint = mean(resLoc);
end
return
%% 无限边界
[t,residual] = LineIntersection(R1,R2,R3,R4,B1,B2,B3,B4,lb,ub);
id1 = [1,2,3,4,1,2];id2 = [2,3,4,1,3,4];
LocRes = zeros(6,3);
for ii = 1:6
eval(['Loc{ii} = (B',num2str(id1(ii)),'(:,2:4).*t(ii,1)+R',num2str(id1(ii)),'(:,2:4)+B',num2str(id2(ii)),...
    '(:,2:4).*t(ii,2)+R',num2str(id2(ii)),'(:,2:4))/2;'])
LocRes(ii,:) = cell2mat(Loc(ii));
end

% 摆！
if length(find(abs(t)<1e-5)) >= 3
    LocPoint = [nan,nan,nan];
    LocRes = nan;
    Point = nan;
    return
else
    Point = [Loc{1};Loc{2};Loc{3};Loc{4};Loc{5};Loc{6}];
    LocPoint = (Loc{1}+Loc{2}+Loc{3}+Loc{4}+Loc{5}+Loc{6})/6;
    return
end

% c_eval('Loc{?} = B?(:,2:4).*t(?)+R?(:,2:4);')
%% 四面体边界
% 如果要仅计算四面体为边界进行计算，可以取消这部分注释（解出单极和其误差位置仅在四面体内）
% % for i = 1:4
% %     ii = 1:4;ii(i) = [];
% %     eval(['[a,b,c,d] = PlaneEquation(R',num2str(ii(1)),'(2:4),R',num2str(ii(2)),'(2:4),R',num2str(ii(3)),'(2:4));'])
% %     eval(['tb(i) = TBoundary(a,b,c,d,B',(num2str(i)),'(2:4),R',(num2str(i)),'(2:4));'])
% % end
% % tb(2) = 0;tb(3) = 0;
% % lb = lb*tb;ub = ub*tb;
% % for i =1:4
% %     if lb(i)>=ub(i) && lb(i) ~=0
% %         lb(i) = -inf;
% %     elseif ub(i)<=lb(i) && ub(i) ~=0
% %         ub(i) = inf;
% %     end
% % end
%% Find the intersection of two tetrahedron
if PlotFlag == 2 || PlotFlag ==3
% 数值求解，不准确
% % % Rmin = Loc{1}; Rmax = Loc{1};
% % % c_eval('Rmin = [min(Rmin(1),Loc{?}(1)),min(Rmin(2),Loc{?}(2)),min(Rmin(3),Loc{?}(3))];')
% % % c_eval('Rmax = [max(Rmax(1),Loc{?}(1)),max(Rmax(2),Loc{?}(2)),max(Rmax(3),Loc{?}(3))];')
% % % % Rmin = Rmin-(Rmax-Rmin);Rmax = Rmax+(Rmax-Rmin);
% % % % Rmin = sign(Rmin).*(abs(Rmin)-abs(Rmax-Rmin));Rmax = sign(Rmax).*(abs(Rmax)+abs(Rmax-Rmin));
% % % step = mean((Rmax-Rmin)/50);
% % % [gridX,gridY,gridZ] = meshgrid(Rmin(1):step:Rmax(1),Rmin(2):step:Rmax(2),Rmin(3):step:Rmax(3));
% % % gridPoint = cell(size(gridX));
% % % for i = 1:numel(gridX(:))
% % %     tempPoint = [gridX(i),gridY(i),gridZ(i)];
% % %     if PointInTetrahedron(tempPoint,R1,R2,R3,R4) && PointInTetrahedron(tempPoint,Loc{1},Loc{2},Loc{3},Loc{4})
% % %         gridPoint{i} = tempPoint;
% % %     end
% % % end
% % % gridPoint(cellfun(@isempty,gridPoint))=[];
% % % LocPoint = [0,0,0];Point = ones(length(gridPoint),3);
% % % for i = 1:numel(gridPoint)
% % %     Point(i,1) = gridPoint{i}(1);
% % %     Point(i,2) = gridPoint{i}(2);
% % %     Point(i,3) = gridPoint{i}(3);
% % %     LocPoint = LocPoint + gridPoint{i};
% % % end
% % % LocPoint =LocPoint/length(gridPoint);

% 解析求交点
Loc_inside = {};Loc_outside = {};id_inside = [];
for i = 1:6
    if PointInTetrahedron(Loc{i},R1,R2,R3,R4)
        Loc_inside{end+1} = Loc{i};
        id_inside(end+1) = i;
    else
        Loc_outside{end+1} = Loc{i};
    end
end
    
if isempty(Loc_outside)
    Point = [Loc{1};Loc{2};Loc{3};Loc{4};Loc{5};Loc{6}];
    LocPoint = (Loc{1}+Loc{2}+Loc{3}+Loc{4}+Loc{5}+Loc{6})/6;
else
    id_choose = nchoosek(1:6,2);
    if length(Loc_inside) >= 2
        id_nchoose = nchoosek(id_inside,2);
        c_eval("id_choose(pdist2(id_nchoose(?,:),id_choose,'hamming')==0,:)=[];",1:size(id_nchoose,1));
    end

    Loc_inter = [];
    for ii = 1:size(id_choose,1)
        for j = 1:4
            jj = 1:4;jj(j) = [];
            eval(['[a,b,c,d] = PlaneEquation(R',num2str(jj(1)),...
                '(2:4),R',num2str(jj(2)),'(2:4),R',num2str(jj(3)),'(2:4));'])
            eval(['Loc_inter(end+1,:) = TBoundary(a,b,c,d,Loc{',...
                (num2str(id_choose(ii,1))),'},Loc{',(num2str(id_choose(ii,2))),'});'])
        end
    end

    Point = unique(Loc_inter,'rows');
    Point(pdist2([0,0,0],Point,'hamming')==0,:)=[];
    Point = [Point;reshape(cell2mat(Loc_inside),3,length(Loc_inside))'];
    id_point = [];
    for lenPoint = 1:size(Point,1)
        if ~PointInTetrahedron(Point(lenPoint,:),R1,R2,R3,R4)
            id_point(end+1) = lenPoint;
        end
    end
    Point(id_point,:)=[];
    LocPoint = sum(Point)/size(Point,1);
end
else
    LocPoint = sum(LocRes{:})/6;Point = nan;
end
%% random find
% TV = 0.001;iter = 0;
% while ~PointInTetrahedron(Loc{1}(1,:),R1,R2,R3,R4)
%     rate = 2*(rand(4,3)-0.5)*5/100+1;
%     c_eval('Brand?(:,2:4) = B?(:,2:4).*rate(?,:);')
%     [t,residual] = LineIntersection(R1,R2,R3,R4,Brand1,Brand2,Brand3,Brand4);
%     iter = iter+1
%     c_eval('Loc{?}(1,:) = B?(:,2:4)*t{?}(1)+R?(:,2:4);')
%     c_eval('Loc{?}(2,:) = B?(:,2:4)*t{?}(2)+R?(:,2:4);')
% end
%% gradient descent to find an acceptable residual & error
% % % TV = 0.001;iter = 0;
% % % step = 0.99;
% % % rate{1} = ones(4,3);rate{2} = 2*(rand(4,3)-0.5)*0.1/100+1;
% % % while max(residual) > TV
% % %     tempres = residual;
% % %     c_eval('B?(:,2:4) = B?(:,2:4).*rate{iter+2}(?,:);')
% % % %     [t,residual] = LineIntersection(R1,R2,R3,R4,Brand1,Brand2,Brand3,Brand4);
% % %     [t,residual] = LineIntersection(R1,R2,R3,R4,B1,B2,B3,B4)
% % %     resgrad = (tempres-residual)/(max(abs(tempres-residual)));
% % %     rate{iter+3} = (rate{iter+1}-rate{iter+2});
% % %     c_eval('rate{iter+3}(?,:) = -step*rate{iter+3}(?,:)*resgrad(?)+1;')
% % %     if ~isempty(residual(residual<=TV))
% % %         i_row = find(residual<=0.1);
% % %         for ii = 1:length(i_row)
% % %             rate{iter+2}(i_row(ii),:) = ones(1,3);
% % %             rate{iter+3}(i_row(ii),:) = ones(1,3);
% % %         end
% % %     end
% % %     iter = iter+1
% % % end
end
%% point of intersection function
function [t,residual]=LineIntersection(R1,R2,R3,R4,B1,B2,B3,B4,lb,ub)
% % % x0 = rand(1,4);
% % % lb = repmat(lb,1,4); ub = repmat(ub,1,4);
% % % [t,~,residual] = lsqnonlin(@(t)myfunc(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub);

x0 = rand(1,2);lb = repmat(lb,1,2); ub = repmat(ub,1,2);
[t1,~,res1] = lsqnonlin(@(t)fun1(t,B1,B2,R1,R2),x0,lb,ub);
[t2,~,res2] = lsqnonlin(@(t)fun2(t,B2,B3,R2,R3),x0,lb,ub);
[t3,~,res3] = lsqnonlin(@(t)fun3(t,B3,B4,R3,R4),x0,lb,ub);
[t4,~,res4] = lsqnonlin(@(t)fun4(t,B4,B1,R4,R1),x0,lb,ub);
[t5,~,res5] = lsqnonlin(@(t)fun5(t,B1,B3,R1,R3),x0,lb,ub);
[t6,~,res6] = lsqnonlin(@(t)fun6(t,B2,B4,R2,R4),x0,lb,ub);

t = [t1;t2;t3;t4;t5;t6];
residual = [res1;res2;res3;res4;res5;res6];
end

function fun1 = fun1(t,B1,B2,R1,R2)
  fun1 = B1(:,2:4).*t(1)+R1(:,2:4) - B2(:,2:4).*t(2)-R2(:,2:4);
end
function fun2 = fun2(t,B2,B3,R2,R3)
  fun2 = B2(:,2:4).*t(1)+R2(:,2:4) - B3(:,2:4).*t(2)-R3(:,2:4);
end
function fun3 = fun3(t,B3,B4,R3,R4)
  fun3 = B3(:,2:4).*t(1)+R3(:,2:4) - B4(:,2:4).*t(2)-R4(:,2:4);
end
function fun4 = fun4(t,B4,B1,R4,R1)
  fun4 = B4(:,2:4).*t(1)+R4(:,2:4) - B1(:,2:4).*t(2)-R1(:,2:4);
end
function fun5 = fun5(t,B1,B3,R1,R3)
  fun5 = B1(:,2:4).*t(1)+R1(:,2:4) - B3(:,2:4).*t(2)-R3(:,2:4);
end
function fun6 = fun6(t,B2,B4,R2,R4)
  fun6 = B2(:,2:4).*t(1)+R2(:,2:4) - B4(:,2:4).*t(2)-R4(:,2:4);
end
%% 2 S/C data into the monopole model to solve the Q & location
function [t,residual]=ModelSolve(R1,R2,R3,R4,B1,B2,B3,B4,lb,ub,MultiPower,idx,LocationSkew)
% x0 = [sign(idx)*1e4*MultiPower+15,0,0,0];
% x0 = [sign(idx)*1e4*MultiPower+15,LocationSkew*ones(1,3)];
x0 = [sign(idx)*1e4*MultiPower,3,-1,2];
% lb = [lb,-inf,-inf,-inf]; ub = [ub,inf,inf,inf];
lb = [];ub = [];
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective');
[t1,~,residual1] = lsqnonlin(@(t)myfunc1(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t2,~,residual2] = lsqnonlin(@(t)myfunc2(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t3,~,residual3] = lsqnonlin(@(t)myfunc3(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t4,~,residual4] = lsqnonlin(@(t)myfunc4(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t5,~,residual5] = lsqnonlin(@(t)myfunc5(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t6,~,residual6] = lsqnonlin(@(t)myfunc6(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);

t = {t1,t2,t3,t4,t5,t6};
residual = {residual1,residual2,residual3,residual4,residual5,residual6};
end

function myfun1 = myfunc1(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];
    c_eval('d? = R?(:,2:4)-loc;');
    c_eval('d? = irf_abs(d?);');
    c_eval('dnorm? = d?(:,1:3)/d?(:,4);');
     
    myfun1 = [t(1)/(4*pi*(d1(4))^2)*dnorm1' - B1(:,2:4)';...
             t(1)/(4*pi*(d2(4))^2)*dnorm2' - B2(:,2:4)';];
end
function myfun2 = myfunc2(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];
    c_eval('d? = R?(:,2:4)-loc;');
    c_eval('d? = irf_abs(d?);');
    c_eval('dnorm? = d?(:,1:3)/d?(:,4);');
     
    myfun2 = [t(1)/(4*pi*(d2(4))^2)*dnorm2' - B2(:,2:4)';...
             t(1)/(4*pi*(d3(4))^2)*dnorm3' - B3(:,2:4)';];
end
function myfun3 = myfunc3(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];
    c_eval('d? = R?(:,2:4)-loc;');
    c_eval('d? = irf_abs(d?);');
    c_eval('dnorm? = d?(:,1:3)/d?(:,4);');
     
    myfun3 = [t(1)/(4*pi*(d3(4))^2)*dnorm3' - B3(:,2:4)';...
             t(1)/(4*pi*(d4(4))^2)*dnorm4' - B4(:,2:4)';];
end
function myfun4 = myfunc4(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];
    c_eval('d? = R?(:,2:4)-loc;');
    c_eval('d? = irf_abs(d?);');
    c_eval('dnorm? = d?(:,1:3)/d?(:,4);');
     
    myfun4 = [t(1)/(4*pi*(d1(4))^2)*dnorm1' - B1(:,2:4)';...
             t(1)/(4*pi*(d3(4))^2)*dnorm3' - B3(:,2:4)';];
end
function myfun5 = myfunc5(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];
    c_eval('d? = R?(:,2:4)-loc;');
    c_eval('d? = irf_abs(d?);');
    c_eval('dnorm? = d?(:,1:3)/d?(:,4);');
     
    myfun5 = [t(1)/(4*pi*(d1(4))^2)*dnorm1' - B1(:,2:4)';...
             t(1)/(4*pi*(d4(4))^2)*dnorm4' - B4(:,2:4)';];
end
function myfun6 = myfunc6(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];
    c_eval('d? = R?(:,2:4)-loc;');
    c_eval('d? = irf_abs(d?);');
    c_eval('dnorm? = d?(:,1:3)/d?(:,4);');
     
    myfun6 = [t(1)/(4*pi*(d2(4))^2)*dnorm2' - B2(:,2:4)';...
             t(1)/(4*pi*(d4(4))^2)*dnorm4' - B4(:,2:4)';];
end
%% plane equation
% % % function [a,b,c,d] = PlaneEquation(R1,R2,R3)
% % % % ax+by+cz+d=0
% % % syms x y z
% % % D=[ones(4,1),[[x,y,z];R1;R2;R3]];%D的行列式等于零即为平面方程
% % % detd=det(D);
% % % coeff = coeffs(detd);
% % % d = double(coeff(1));
% % % c = double(coeff(2));
% % % b = double(coeff(3));
% % % a = double(coeff(4));
% % % end
%% t boundary
% % % function Loc_inter = TBoundary(a,b,c,d,r1,r2)
% % % syms x y z real
% % % f1 = a*x+b*y+c*z+d == 0;
% % % f2 = (x-r1(1))/(r2(1)-r1(1))==(y-r1(2))/(r2(2)-r1(2));
% % % f3 = (y-r1(2))/(r2(2)-r1(2))==(z-r1(3))/(r2(3)-r1(3));
% % % [x,y,z] = solve([f1,f2,f3],x,y,z);
% % % x = double(x);
% % % y = double(y);
% % % z = double(z);
% % % % Loc_inter = [x,y,z];
% % % if (min(r1(1),r2(1)) <= x) &&  (x <= max(r1(1),r2(1)))
% % % Loc_inter = [x,y,z];
% % % else
% % %     Loc_inter = [0,0,0];
% % % end
% % % end