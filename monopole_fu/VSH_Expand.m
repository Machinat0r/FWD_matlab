function [coeff, res] = VSH_Expand(varargin)
%------written by Wending Fu, Mar.2024 in Beijing------------
%% Input Data
flag = varargin{end};
if flag == 1
SphereCenter = varargin{end-1}; %#ok
varargin = varargin(1:end-2);
else
varargin = varargin(1:end-1);
end
if length(varargin) == 2
    Rs = varargin{1};
	Bs = varargin{2};
    for id=1:4
        temp = evalin('caller',irf_ssub(Rs,id)); %#ok
        eval(irf_ssub('R? =temp;',id));
        temp = evalin('caller',irf_ssub(Bs,id)); %#ok
        eval(irf_ssub('B? =temp;',id)); clear temp
    end
elseif length(varargin) == 8
    c_eval('R? = varargin{?};',1:4);
    c_eval('B? = varargin{?+4};',1:4);
else
    disp('Incorrect number of input parameters. See usage:')
    help VSH_Expand
    return
end
%% Resample
c_eval('B? = irf_resamp(B?, B1);', 2:4);
c_eval('R? = irf_resamp(R?, B1);', 1:4);
c_eval('B?_sph = B?;', 1:4);
c_eval('R?_sph = R?;', 1:4);
%% Coordinate Transformation
if flag == 1
% % % SphereCenter = zeros(size(R1,1),3);
% % % for temp_t = 1:length(R1)
% % % TR = delaunayTriangulation([R1(temp_t,2:4); R2(temp_t,2:4);R3(temp_t,2:4); R4(temp_t,2:4)]);
% % % SphereCenter(temp_t,:) = circumcenter(TR);
% % % end

c_eval('R?(:,2:4) = R?(:,2:4) - SphereCenter;', 1:4);
[B1_sph, R1_sph] = Coor_Trans(B1(:,2:4), R1(:,2:4)); %#ok
[B2_sph, R2_sph] = Coor_Trans(B2(:,2:4), R2(:,2:4)); %#ok
[B3_sph, R3_sph] = Coor_Trans(B3(:,2:4), R3(:,2:4)); %#ok
[B4_sph, R4_sph] = Coor_Trans(B4(:,2:4), R4(:,2:4)); %#ok
c_eval('B?_sph = [B?(:,1), B?_sph];', 1:4);
c_eval('R?_sph = [R?(:,1), R?_sph];', 1:4);
end
%% Calculate Coefficients
if flag ~= 3
coeff = zeros(size(B1_sph,1),15); %12
res = zeros(size(B1_sph,1),12); %12
else
coeff = zeros(size(B1_sph,1),30); %12
res = zeros(size(B1_sph,1),15); %12
end

parfor temp_t = 1:size(B1_sph,1)
[coeff(temp_t,:), res(temp_t,:)] = SPH_coef(B1_sph(temp_t,2:4), B2_sph(temp_t,2:4), B3_sph(temp_t,2:4),...
    B4_sph(temp_t,2:4), R1_sph(temp_t,2:4), R2_sph(temp_t,2:4), R3_sph(temp_t,2:4),...
    R4_sph(temp_t,2:4), flag); %#ok
end

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

%% Solve Spherical Harmonic Coeffecients
function [coeff, res] = SPH_coef(B1, B2, B3, B4, R1, R2, R3, R4, flag)
if flag == 1
% one-step
c_eval('r? = R?(:,1);',1:4);
[r1, r2, r3, r4, ~] = Normal_R(r1, r2, r3, r4); %#ok

c_eval('theta? = R?(:,2);',1:4);
c_eval('phi? = R?(:,3);',1:4);
x0 = zeros(12,1);

c_eval('theta? = theta? + eps;', 1:4);
c_eval('phi? = phi? + eps;', 1:4);
[coeff,~,res,~] = lsqnonlin(@(coeff) Cal_func(B1, B2, B3, B4, r1, r2, r3, r4,...
    theta1, theta2, theta3, theta4, phi1, phi2, phi3, phi4, coeff), x0);
coeff = [coeff; 0;0;0];


elseif flag == 2
% two-step
x0 = eps*ones(6,1);
x1 = eps*ones(9,1);
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1e4,...
    'MaxIterations',1e4);
% % % lb = [-inf, -inf, -inf, -5, -5, -5]; ub = [inf, inf, inf, 5,5,5];
[coeff0,~,res,~,~] = lsqnonlin(@(coeff) Cal_func0(B1, B2, B3, B4, R1, R2, R3, R4, coeff),...
    x0, [], [], options);
% [coeff0,~,res,~,~] = lsqnonlin(@(coeff) Cal_func0(B1, B2, B3, B4, R1, R2, R3, R4, coeff), x0);
R0 = [coeff0(4),coeff0(5),coeff0(6)];
% % % R0 = [1,0.5,-0.5];

% % % B1 = [res(1), res(5), res(9)]; B2 = [res(2), res(6), res(10)];
% % % B3 = [res(3), res(7), res(11)]; B4 = [res(4), res(8), res(12)];
% % % c_eval('R? = R? - R0;', 1:4);
% % % [~, R1_temp] = Coor_Trans(B1, R1); [~, R2_temp] = Coor_Trans(B2, R2); %#ok
% % % [~, R3_temp] = Coor_Trans(B3, R3); [~, R4_temp] = Coor_Trans(B4, R4); %#ok
% % % c_eval('B?x = B?(1)*sin(R?_temp(2))*cos(R?_temp(3)) + B?(2)*cos(R?_temp(2))*cos(R?_temp(3)) - B?(3)*sin(R?_temp(3));', 1:4);
% % % c_eval('B?y = B?(1)*sin(R?_temp(2))*sin(R?_temp(3)) + B?(2)*cos(R?_temp(2))*sin(R?_temp(3)) + B?(3)*cos(R?_temp(3));', 1:4);
% % % c_eval('B?z = B?(1)*cos(R?_temp(2)) - B?(2)*sin(R?_temp(2));', 1:4);
% % % c_eval('B? = -[B?x, B?y, B?z];', 1:4);

B3_res = [res(1), res(3), res(5)]; B1_res = [res(2), res(4), res(6)]; 
c_eval('R? = R? - R0;', 1:4);
[B1, R1_temp] = Coor_Trans(B1, R1); [B2, R2_temp] = Coor_Trans(B2, R2); %#ok
[B3, R3_temp] = Coor_Trans(B3, R3); [B4, R4_temp] = Coor_Trans(B4, R4); %#ok
[R1_temp_norm, R2_temp_norm, R3_temp_norm, R4_temp_norm, ~] = Normal_R(R1_temp(1), R2_temp(1), R3_temp(1), R4_temp(1)); %#ok
alpha_00 = coeff0(1); beta_00 = coeff0(2); gamma_00 = coeff0(3);
for ic = 1:4
temp_ic = num2str(ic);
eval(['theta',temp_ic,' = R',temp_ic,'_temp(2) + eps;']);
eval(['Br0',temp_ic,' = Cal_Br0(R',temp_ic,'_temp_norm,theta',temp_ic,', alpha_00);']);
eval(['Bt0',temp_ic,' = Cal_Bt0(R',temp_ic,'_temp_norm,theta',temp_ic,', beta_00);']);
eval(['Bp0',temp_ic,' = Cal_Bp0(R',temp_ic,'_temp_norm,theta',temp_ic,', gamma_00);']);
end
B1_res = [Br01, Bt01, Bp01] - B1; B2_res = [Br02, Bt02, Bp02] - B2;
B3_res = [Br03, Bt03, Bp03] - B3; B4_res = [Br04, Bt04, Bp04] - B4;
c_eval('B? = B?_res;',1:4);
c_eval('B?x = B?(1)*sin(R?_temp(2))*cos(R?_temp(3)) + B?(2)*cos(R?_temp(2))*cos(R?_temp(3)) - B?(3)*sin(R?_temp(3));', 1:4);
c_eval('B?y = B?(1)*sin(R?_temp(2))*sin(R?_temp(3)) + B?(2)*cos(R?_temp(2))*sin(R?_temp(3)) + B?(3)*cos(R?_temp(3));', 1:4);
c_eval('B?z = B?(1)*cos(R?_temp(2)) - B?(2)*sin(R?_temp(2));', 1:4);
c_eval('B? = -[B?x, B?y, B?z];', 1:4);

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1e4,...
    'MaxIterations',1e4);
[coeff1,~,res,~] = lsqnonlin(@(coeff) Cal_func1(B1, B2, B3, B4, R1, R2, R3, R4, coeff), x1,[],[],options);
coeff = [coeff0(1); coeff1(1:3); coeff0(2); coeff1(4:6); coeff0(3); coeff1(7:9);coeff0(4:6)];
res = [res;0;0;0];
end
end
%% one-step calculation
% Calculate Br, Btheta, Bphi
% coeff = {alpha00, alpha1-1, alpha10, alpha11, beta00, beta1-1, beta10, 
% beta11, gamma00, gamma1-1, gamma10, gamma11} ----[alpha,beta,gamma][l,m]
function func = Cal_func(B1, B2, B3, B4, r1, r2, r3, r4, theta1, theta2, theta3, theta4, phi1, phi2, phi3, phi4, coeff) %#ok
alpha_00 = coeff(1); alpha_12 = coeff(2); alpha_10 = coeff(3); alpha_11 = coeff(4); %#ok, alpha_12 = alpha1-1 
beta_00 = coeff(5); beta_12 = coeff(6); beta_10 = coeff(7); beta_11 = coeff(8); %#ok
gamma_00 = coeff(9); gamma_12 = coeff(10); gamma_10 = coeff(11); gamma_11 = coeff(12); %#ok
 
% one-step
for ic = 1:4
temp_ic = num2str(ic);
eval(['Br',temp_ic,' = Cal_Br(r',temp_ic,',theta',temp_ic,', phi',temp_ic,', alpha_00, alpha_12, alpha_10, alpha_11);']);
eval(['Bt',temp_ic,' = Cal_Bt(r',temp_ic,',theta',temp_ic,', phi',temp_ic,', beta_00, beta_12, beta_10, beta_11, gamma_12, gamma_11);']);
eval(['Bp',temp_ic,' = Cal_Bp(r',temp_ic,',theta',temp_ic,', phi',temp_ic,', beta_12, beta_11, gamma_00, gamma_12, gamma_10, gamma_11);']);
end

func = [Br1 - B1(1); Br2 - B2(1); Br3 - B3(1); Br4 - B4(1);
        Bt1 - B1(2); Bt2 - B2(2); Bt3 - B3(2); Bt4 - B4(2);
        Bp1 - B1(3); Bp2 - B2(3); Bp3 - B3(3); Bp4 - B4(3);];
end
%% two-step calculation
% l = 0
function func = Cal_func0(B1, B2, B3, B4, R1, R2, R3, R4, coeff) %#ok
alpha_00 = coeff(1); beta_00 = coeff(2); gamma_00 = coeff(3); %#ok
Rx0 = coeff(4); Ry0 = coeff(5); Rz0 = coeff(6);
R0 = [Rx0, Ry0, Rz0]; %#ok
% % % R0 = [1,0.5,-0.5];

for ic = 1:4
temp_ic = num2str(ic);
eval(['R',temp_ic,' = R',temp_ic,' - R0;']);
eval(['[B',temp_ic,', R',temp_ic,'] = Coor_Trans(B',temp_ic,', R',temp_ic,');']);
eval(['r',temp_ic,' = R',temp_ic,'(1);']);
eval(['theta',temp_ic,' = R',temp_ic,'(2) + eps;']); % avoid divide sin(0)
end
[r1, r2, r3, r4, ~] = Normal_R(r1, r2, r3, r4); %#ok
% [B1, B2, B3, B4, ~] = Normal_B(B1, B2, B3, B4);
for ic = 1:4
temp_ic = num2str(ic);
eval(['Br0',temp_ic,' = Cal_Br0(r',temp_ic,',theta',temp_ic,', alpha_00);']);
eval(['Bt0',temp_ic,' = Cal_Bt0(r',temp_ic,',theta',temp_ic,', beta_00);']);
eval(['Bp0',temp_ic,' = Cal_Bp0(r',temp_ic,',theta',temp_ic,', gamma_00);']);
end
% % % 
% % % func = [Br01 - B1(1); Br02 - B2(1); Br03 - B3(1); Br04 - B4(1);
% % %         Bt01 - B1(2); Bt02 - B2(2); Bt03 - B3(2); Bt04 - B4(2);
% % %         Bp01 - B1(3); Bp02 - B2(3); Bp03 - B3(3); Bp04 - B4(3);];
func = [Br03 - B3(1); Br01 - B1(1); 
        Bt03 - B3(2); Bt01 - B1(2);
        Bp03 - B3(3); Bp01 - B1(3);]; 
end


% l = 1
function func = Cal_func1(B1, B2, B3, B4, R1, R2, R3, R4, coeff) %#ok
% alpha_12 = alpha1-1
alpha_12 = coeff(1); alpha_10 = coeff(2); alpha_11 = coeff(3); %#ok
beta_12 = coeff(4); beta_10 = coeff(5); beta_11 = coeff(6); %#ok
gamma_12 = coeff(7); gamma_10 = coeff(8); gamma_11 = coeff(9); %#ok
% Rx0 = coeff(10); Ry0 = coeff(11); Rz0 = coeff(12);
% R0 = [Rx0, Ry0, Rz0]; %#ok
R0 = [0, 0, 0];

for ic = 1:4
temp_ic = num2str(ic);
eval(['R',temp_ic,' = R',temp_ic,' - R0;']);
eval(['[B',temp_ic,', R',temp_ic,'] = Coor_Trans(B',temp_ic,', R',temp_ic,');']);
eval(['r',temp_ic,' = R',temp_ic,'(1);']);
eval(['theta',temp_ic,' = R',temp_ic,'(2) + eps;']);
eval(['phi',temp_ic,' = R',temp_ic,'(3) + eps;']);
end
[r1, r2, r3, r4, ~] = Normal_R(r1, r2, r3, r4); %#ok
for ic = 1:4
temp_ic = num2str(ic);
eval(['Br1',temp_ic,' = Cal_Br1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', alpha_12, alpha_10, alpha_11);']);
eval(['Bt1',temp_ic,' = Cal_Bt1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', beta_12, beta_10, beta_11, gamma_12, gamma_11);']);
eval(['Bp1',temp_ic,' = Cal_Bp1(r',temp_ic,', theta',temp_ic,', phi',temp_ic,', beta_12, beta_11, gamma_12, gamma_10, gamma_11);']);
end

func = [Br11 - B1(1); Br12 - B2(1); Br13 - B3(1); 
        Bt11 - B1(2); Bt12 - B2(2); Bt13 - B3(2); 
        Bp11 - B1(3); Bp12 - B2(3); Bp13 - B3(3);];
% % % func = [Br11 - B1(1); Br12 - B2(1); Br13 - B3(1); Br14 - B4(1); 
% % %         Bt11 - B1(2); Bt12 - B2(2); Bt13 - B3(2); Bt14 - B4(2);
% % %         Bp11 - B1(3); Bp12 - B2(3); Bp13 - B3(3); Bp14 - B4(3);];
end
%% component function
function Br = Cal_Br(r, theta, phi, alpha_00, alpha_12, alpha_10, alpha_11)
Br = r.^-2 .* alpha_00 .* Cal_SPH(0, 0, theta) + ...
    r.^-3 .* (alpha_12 .* cos(-phi) .* Cal_SPH(1, -1, theta)...
    + alpha_10 .* Cal_SPH(1, 0, theta) + ...
    alpha_11 .* cos(phi) .* Cal_SPH(1, 1, theta));
end

function Bt = Cal_Bt(r, theta, phi, beta_00, beta_12, beta_10, beta_11, gamma_12, gamma_11)
Bt =r.^-1 .* beta_00 .* Cal_dPlm(0, 0, theta) + ...
    r.^-2 .* (gamma_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta) + ...
    beta_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta) + ...
    beta_10 ./ 2 .* Cal_dPlm(1, 0, theta) + ...
    gamma_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta) + ...
    beta_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta));
end

function Bp = Cal_Bp(r, theta, phi, beta_12, beta_11, gamma_00, gamma_12, gamma_10, gamma_11)
Bp =r.^-1 .* gamma_00 .* Cal_dPlm(0, 0, theta) + ...
    r.^-2 .* (gamma_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta) -...
    beta_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta) + ...
    gamma_10 ./ 2 .* Cal_dPlm(1, 0, theta) + ...
    gamma_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta) -...
    beta_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta));
end
%% two-step component function
function Br0 = Cal_Br0(r, theta, alpha_00)
Br0 = r.^-2 .* alpha_00 .* Cal_SPH(0, 0, theta);
end

function Bt0 = Cal_Bt0(r, theta, beta_00)
Bt0 = r.^-1 .* beta_00 .* Cal_dPlm(0, 0, theta);
end

function Bp0 = Cal_Bp0(r, theta, gamma_00)
Bp0 = r.^-1 .* gamma_00 .* Cal_dPlm(0, 0, theta);
end

function Br1 = Cal_Br1(r, theta, phi, alpha_12, alpha_10, alpha_11)
Br1 = r.^-3 .* (alpha_12 .* cos(-phi) .* Cal_SPH(1, -1, theta)...
    + alpha_10 .* Cal_SPH(1, 0, theta) + ...
    alpha_11 .* cos(phi) .* Cal_SPH(1, 1, theta));
end

function Bt1 = Cal_Bt1(r, theta, phi, beta_12, beta_10, beta_11, gamma_12, gamma_11)
Bt1 = r.^-2 .* (gamma_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta) + ...
    beta_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta) + ...
    beta_10 ./ 2 .* Cal_dPlm(1, 0, theta) + ...
    gamma_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta) + ...
    beta_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta));
end

function Bp1 = Cal_Bp1(r, theta, phi, beta_12, beta_11, gamma_12, gamma_10, gamma_11)
Bp1 = r.^-2 .* (gamma_12 .* sin(phi) ./ 2 .* Cal_dPlm(1, -1, theta) -...
    beta_12 .* cos(phi) ./ (2*sin(theta)) .* Cal_SPH(1, -1, theta) + ...
    gamma_10 ./ 2 .* Cal_dPlm(1, 0, theta) + ...
    gamma_11 .* cos(phi) ./ 2 .* Cal_dPlm(1, 1, theta) -...
    beta_11 .* -sin(phi) ./ (2*sin(theta)) .* Cal_SPH(1, 1, theta));
end

%% R normalization
function [R1, R2, R3, R4, maxR] = Normal_R(R1, R2, R3, R4)
maxR = 0;
c_eval('maxR = max(maxR,R?(:,1));',1:4);
c_eval('R? = R?./maxR;',1:4);
end
%% B normalization
function [B1, B2, B3, B4, maxB_norm] = Normal_B(B1, B2, B3, B4)
maxB_norm = 0;
c_eval('maxB_norm = max(maxB_norm,norm(B?(:,1)));',1:4);
c_eval('B? = B?./maxB_norm;',1:4);
end
%% coordinate transformation
function [B, R] = Coor_Trans(B, R)
[phi, theta, r] = cart2sph(R(:,1), R(:,2), R(:,3));
theta = pi/2 - theta;
R = [r, theta, phi];

Br = B(:,1).* sin(theta).* cos(phi) + B(:,2).* sin(theta).* sin(phi) + B(:,3).* cos(theta);
Bt = B(:,1).* cos(theta).* cos(phi) + B(:,2).* cos(theta).* sin(phi) - B(:,3).* sin(theta);
Bp = -B(:,1).* sin(phi) + B(:,2).* cos(phi);
B = [Br, Bt, Bp];
end
%% constraint function i.e. boundary point in the tetradron
function [c, ceq] = constraintFunction0(R1, R2, R3, R4, coeff)
R0 = coeff(4:6)';
TR = delaunayTriangulation([R1;R2;R3;R4]);
ID = pointLocation(TR,R0);
c = ID - 1;
ceq = [];
end