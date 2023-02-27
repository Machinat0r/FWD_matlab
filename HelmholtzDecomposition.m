function [Br,Bp] = HelmholtzDecomposition(varargin)
% wrong function
% cannot be used yet
%------written by Wending Fu, Apr.19.2022 in Beijing------------
if length(varargin) == 2 
    Rs = varargin{1};
	Bs = varargin{2};
    for id=1:4
        temp = evalin('caller',irf_ssub(Rs,id));
        eval(irf_ssub('R? =temp;',id));
        temp = evalin('caller',irf_ssub(Bs,id));
        eval(irf_ssub('B? =temp;',id)); clear ttt
    end
elseif length(varargin) == 8
    c_eval('R? = varargin{?};',1:4);
    c_eval('B? = varargin{?+4};',1:4);
else
    disp('Incorrect number of input parameters. See usage:')
    help HelmholtzDecomposition
    return
end
c_eval('irf_resamp(R?,B?);')
% if length(R1) == 4, c_eval('R? = R?(:,2:4);') ,end
% if length(B1) == 4, c_eval('B? = B?(:,2:4);'),end
c_eval('B? = 1e-9*B?(:,2:4);')
units = irf_units;
c_eval('R? = R?(:,2:4)*units.RE;')

% curl
curlB = 1e9*c_4_grad('R?','B?','curl');
Bp = curlB;
% curlB = curl_Jac(gradB);

Br = 0;
end