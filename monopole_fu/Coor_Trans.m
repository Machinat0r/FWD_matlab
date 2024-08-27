%% coordinate transformation
function varargout = Coor_Trans(varargin)
if length(varargin) == 1
flag = 'R';
R = varargin{1};
elseif length(varargin) == 2
flag = 'BR';
R = varargin{1};
B = varargin{2};
end
%%
[phi, theta, r] = cart2sph(R(:,1), R(:,2), R(:,3));
theta = pi/2 - theta;
R = [r, theta, phi];

if flag == 'R'
varargout{1} = r;
varargout{2} = theta;
varargout{3} = phi;
elseif flag == 'BR'
Br = B(:,1).* sin(theta).* cos(phi) + B(:,2).* sin(theta).* sin(phi) + B(:,3).* cos(theta);
Bt = B(:,1).* cos(theta).* cos(phi) + B(:,2).* cos(theta).* sin(phi) - B(:,3).* sin(theta);
Bp = -B(:,1).* sin(phi) + B(:,2).* cos(phi);
B = [Br, Bt, Bp];
varargout{1} = R;
varargout{2} = B;
end
end