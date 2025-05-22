function [Q,resQ,LocPoint,resLoc] = CalError_Remake(varargin)
% function [LocRes,Loc,LocPoint,Point] = CalError_Remake(varargin)
% PlotFlag choose whether to find the intersection of two tetrahedrons, if
% PlotFlag = 2, solute the intersection
% R, B, idx_flag, MultiPower, LocationSkew
%------written by Wending Fu, May.8.2022 in Beijing------------
%------remaked by Wending Fu, May.19.2025 in Beijing------------
% lb,ub未启用
%% input data
if length(varargin) == 5 || length(varargin)==11
    LocationSkew = varargin{end};
    varargin = varargin(1:end-1);
    MultiPower = varargin{end};
    varargin = varargin(1:end-1);
    idx = varargin{end};        idx_flag = sign(idx);
    if sign(idx) == -1 
        lb = -inf;ub = 0;
        idx = abs(idx);
%         MultiPower = -MultiPower;
    elseif sign(idx) == 1
        lb = 0;ub = inf;  %2022.7.11,修改ModelSolve时调换了该边界；若使用LineIntersection求解，应将+-1边界反过来
    end
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
% if ~isempty(idx)
% %     c_eval('R?(:,2:4) = units.RE*R?(:,2:4);')
%     c_eval('B? = irf_resamp(B?,B1);',2:4)
%     c_eval('R? = irf_resamp(R?,B1);')
%     % c_eval('B? = B?(idx,:);')
%     % c_eval('R? = R?(idx_R,:);')
% end
B2 = irf_resamp(B2, B1); B3 = irf_resamp(B3, B1); B4 = irf_resamp(B4, B1);
R1 = irf_resamp(R1, B1); R2 = irf_resamp(R2, B1); R3 = irf_resamp(R3, B1); R4 = irf_resamp(R4, B1);
%% residual
%% 使用2颗卫星代入模型求解
Q = nan*ones(size(B1,1),1); resQ = cell(size(B1,1),1);
resLoc = cell(size(B1,1),1); LocPoint = nan*ones(size(B1,1),3);
parfor ii = 1:size(B1,1)
tmpR1 = R1(ii,:); tmpR2 = R2(ii,:); tmpR3 = R3(ii,:); tmpR4 = R4(ii,:); 
tmpB1 = B1(ii,:); tmpB2 = B2(ii,:); tmpB3 = B3(ii,:); tmpB4 = B4(ii,:); 
[t,~,flag] = ModelSolve(tmpR1,tmpR2,tmpR3,tmpR4,tmpB1,tmpB2,tmpB3,tmpB4,...
    lb,ub,MultiPower,idx_flag,LocationSkew);
% % % if isempty(find(cell2mat(flag)>=-1, 1))
% % %     Q(ii) = nan;resQ = nan*ones(6,1);
% % %     resLoc = nan*ones(6,3);LocPoint = nan*ones(1,3);
% % % else
resQ{ii} = [t{1}(1);t{2}(1);t{3}(1);t{4}(1);t{5}(1);t{6}(1)];
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
% end
end
return

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
function [t,residual,flag]=ModelSolve(R1,R2,R3,R4,B1,B2,B3,B4,lb,ub,MultiPower,idx_flag,LocationSkew)
% x0 = [idx_flag*1e4*MultiPower+15,10,10,10];
% x0 = [idx_flag*1e4*MultiPower+15,LocationSkew*ones(1,3)];
x0 = [idx_flag*1e4*MultiPower,LocationSkew*ones(1,3)];
% x0 = [idx_flag*1e4*MultiPower,3,-1,2];
% x0 = [idx_flag*1e4*MultiPower,ones(1,3)*MultiPower];
% x0 = [idx_flag*1e4*MultiPower,0,0,0];
% lb = [0.1*idx_flag*1e4*MultiPower,1000*ones(1,3)*MultiPower]; ub = [ub,1000*ones(1,3)*MultiPower];
lb = [];ub = [];
options = optimoptions('lsqnonlin');
% options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
    % 'OptimalityTolerance',1e-10,'FunctionTolerance',1e-10,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);
[t1,~,residual1] = lsqnonlin(@(t)myfunc1(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t2,~,residual2] = lsqnonlin(@(t)myfunc2(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t3,~,residual3] = lsqnonlin(@(t)myfunc3(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t4,~,residual4] = lsqnonlin(@(t)myfunc4(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t5,~,residual5] = lsqnonlin(@(t)myfunc5(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
[t6,~,residual6] = lsqnonlin(@(t)myfunc6(t,B1,B2,B3,B4,R1,R2,R3,R4),x0,lb,ub,options);
% % % [t1,residual1,flag1] = fsolve(@(t)myfunc1(t,B1,B2,B3,B4,R1,R2,R3,R4),x0);
% % % [t2,residual2,flag2] = fsolve(@(t)myfunc2(t,B1,B2,B3,B4,R1,R2,R3,R4),x0);
% % % [t3,residual3,flag3] = fsolve(@(t)myfunc3(t,B1,B2,B3,B4,R1,R2,R3,R4),x0);
% % % [t4,residual4,flag4] = fsolve(@(t)myfunc4(t,B1,B2,B3,B4,R1,R2,R3,R4),x0);
% % % [t5,residual5,flag5] = fsolve(@(t)myfunc5(t,B1,B2,B3,B4,R1,R2,R3,R4),x0);
% % % [t6,residual6,flag6] = fsolve(@(t)myfunc6(t,B1,B2,B3,B4,R1,R2,R3,R4),x0);

t = {t1,t2,t3,t4,t5,t6};
residual = {residual1,residual2,residual3,residual4,residual5,residual6};
flag1 = 2; flag2 = 2; flag3 = 2; flag4 = 2; flag5 = 2; flag6 = 2;
flag = {flag1,flag2,flag3,flag4,flag5,flag6};
end

function myfun1 = myfunc1(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];

    d1 = irf_abs(R1(:,2:4)-loc); d2 = irf_abs(R2(:,2:4)-loc); 
    d3 = irf_abs(R3(:,2:4)-loc); d4 = irf_abs(R4(:,2:4)-loc);

    dnorm1 = d1(:,1:3)/d1(:,4); dnorm2 = d2(:,1:3)/d2(:,4);
    dnorm3 = d3(:,1:3)/d3(:,4); dnorm4 = d4(:,1:3)/d4(:,4);
     
    myfun1 = [t(1)/(4*pi*(d1(4))^2)*dnorm1' - B1(:,2:4)';...
             t(1)/(4*pi*(d2(4))^2)*dnorm2' - B2(:,2:4)';];
end
function myfun2 = myfunc2(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];

    d1 = irf_abs(R1(:,2:4)-loc); d2 = irf_abs(R2(:,2:4)-loc); 
    d3 = irf_abs(R3(:,2:4)-loc); d4 = irf_abs(R4(:,2:4)-loc);

    dnorm1 = d1(:,1:3)/d1(:,4); dnorm2 = d2(:,1:3)/d2(:,4);
    dnorm3 = d3(:,1:3)/d3(:,4); dnorm4 = d4(:,1:3)/d4(:,4);
     
    myfun2 = [t(1)/(4*pi*(d2(4))^2)*dnorm2' - B2(:,2:4)';...
             t(1)/(4*pi*(d3(4))^2)*dnorm3' - B3(:,2:4)';];
end
function myfun3 = myfunc3(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];

    d1 = irf_abs(R1(:,2:4)-loc); d2 = irf_abs(R2(:,2:4)-loc); 
    d3 = irf_abs(R3(:,2:4)-loc); d4 = irf_abs(R4(:,2:4)-loc);

    dnorm1 = d1(:,1:3)/d1(:,4); dnorm2 = d2(:,1:3)/d2(:,4);
    dnorm3 = d3(:,1:3)/d3(:,4); dnorm4 = d4(:,1:3)/d4(:,4);
     
    myfun3 = [t(1)/(4*pi*(d3(4))^2)*dnorm3' - B3(:,2:4)';...
             t(1)/(4*pi*(d4(4))^2)*dnorm4' - B4(:,2:4)';];
end
function myfun4 = myfunc4(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];

    d1 = irf_abs(R1(:,2:4)-loc); d2 = irf_abs(R2(:,2:4)-loc); 
    d3 = irf_abs(R3(:,2:4)-loc); d4 = irf_abs(R4(:,2:4)-loc);

    dnorm1 = d1(:,1:3)/d1(:,4); dnorm2 = d2(:,1:3)/d2(:,4);
    dnorm3 = d3(:,1:3)/d3(:,4); dnorm4 = d4(:,1:3)/d4(:,4);
     
    myfun4 = [t(1)/(4*pi*(d1(4))^2)*dnorm1' - B1(:,2:4)';...
             t(1)/(4*pi*(d3(4))^2)*dnorm3' - B3(:,2:4)';];
end
function myfun5 = myfunc5(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];

    d1 = irf_abs(R1(:,2:4)-loc); d2 = irf_abs(R2(:,2:4)-loc); 
    d3 = irf_abs(R3(:,2:4)-loc); d4 = irf_abs(R4(:,2:4)-loc);

    dnorm1 = d1(:,1:3)/d1(:,4); dnorm2 = d2(:,1:3)/d2(:,4);
    dnorm3 = d3(:,1:3)/d3(:,4); dnorm4 = d4(:,1:3)/d4(:,4);
     
    myfun5 = [t(1)/(4*pi*(d1(4))^2)*dnorm1' - B1(:,2:4)';...
             t(1)/(4*pi*(d4(4))^2)*dnorm4' - B4(:,2:4)';];
end
function myfun6 = myfunc6(t,B1,B2,B3,B4,R1,R2,R3,R4)
% t = [Q,x,y,z]
    loc = [t(2),t(3),t(4)];

    d1 = irf_abs(R1(:,2:4)-loc); d2 = irf_abs(R2(:,2:4)-loc); 
    d3 = irf_abs(R3(:,2:4)-loc); d4 = irf_abs(R4(:,2:4)-loc);

    dnorm1 = d1(:,1:3)/d1(:,4); dnorm2 = d2(:,1:3)/d2(:,4);
    dnorm3 = d3(:,1:3)/d3(:,4); dnorm4 = d4(:,1:3)/d4(:,4);
     
    myfun6 = [t(1)/(4*pi*(d2(4))^2)*dnorm2' - B2(:,2:4)';...
             t(1)/(4*pi*(d4(4))^2)*dnorm4' - B4(:,2:4)';];
end