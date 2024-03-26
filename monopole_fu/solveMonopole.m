function [Q,resQ] =  solveMonopole(varargin)
% 
%------written by Wending Fu, Apr.19.2022 in Beijing------------
%% input data
if length(varargin) == 4 || length(varargin)==10
    idx = varargin{end};
    if sign(idx) == -1
        flag = -1;
        idx = abs(idx);
    elseif sign(idx) == 1
        flag = 1;
    else
        disp('Incorrect index. See usage:')
        help solveMonopole
    return
    end
    Loc = varargin{end-1};
    varargin = varargin(1:end-2);
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
    help solveMonopole
    return
end
%% resamp
if ~isempty(idx)
    c_eval('B? = irf_resamp(B?,B1);',2:4)
    c_eval('R? = irf_resamp(R?,B1);')
%     c_eval('B? = B?(idx,:);')
%     c_eval('R? = R?(idx,:);')
end
%%
if size(R1,2) >= 4, c_eval('tR? = R?(:,1);');c_eval('R? = R?(:,2:4);');end
if size(B1,2) >= 4, c_eval('tB? = B?(:,1);');c_eval('B? = B?(:,2:4);');end
%% B projection
c_eval('DirectionVector? = R? - Loc;')
% c_eval('Bproj? = DirectionVector?*dot(DirectionVector?,B?)/(norm(DirectionVector?)*norm(DirectionVector?));')
c_eval('Bproj? = B?;')
%%
units = irf_units;
% c_eval('Bproj? = Bproj?*1e-2;')
% c_eval('R? = R?*units.RE;')
% Loc = Loc*units.RE;
c_eval('Bt? = irf_abs(Bproj?);')



%solve
% syms Q real
c_eval('deltR? = (R?(1)-Loc(1)).^2 + (R?(2)-Loc(2)).^2 + (R?(3)-Loc(3)).^2;')
c_eval('f?  = @(Q)Q - Bt?(4)*(4*pi*deltR?);')
x0 = 10;
% lb = 0; ub = inf;
[Q] = lsqnonlin(@(Q)myfunc(Q,deltR1,deltR2,deltR3,deltR4,Bt1,Bt2,Bt3,Bt4),x0);%nT*km^2
% Q = Q/((units.RE/1e3)^2); %nT*RE^2


c_eval('resQ? = Bt?(4)*4*pi*deltR?;')
% c_eval('resQ? = log10(resQ?);')
% resQ = [min([resQ1,resQ2,resQ3,resQ4]), max([resQ1,resQ2,resQ3,resQ4])];
resQ = [resQ1,resQ2,resQ3,resQ4];
% c_eval('resQ? = residual(?)+Bt?(4)*4*pi*deltR?;')
% c_eval('resQ? = resQ?*1e-9*units.RE^2;')
% c_eval('logQ? = log10(abs(resQ?));')
% logQ0 = log10(Q);
% c_eval('resQ(?) = resQ?;');c_eval('logQ(?) = logQ?;')
% for i =1:4
%     if resQ(i) >0
%         logQ(i) = logQ(i) + logQ0;
%     elseif resQ(i) <0
%         logQ(i) = logQ0 - logQ(i);
%     end
% end
% resQ = [min([logQ(1),logQ(2),logQ(3),logQ(4)]) max([logQ(1),logQ(2),logQ(3),logQ(4)])];

Q = flag*Q;
resQ = flag*resQ;
end
function myfun = myfunc(Q,deltR1,deltR2,deltR3,deltR4,Bt1,Bt2,Bt3,Bt4)
    myfun = [Q-(4*pi*deltR1) * Bt1(4);...
                Q-(4*pi*deltR2) * Bt2(4);...
                Q-(4*pi*deltR3) * Bt3(4);...
                Q-(4*pi*deltR4) * Bt4(4);];
end