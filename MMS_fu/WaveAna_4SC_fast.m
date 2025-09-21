function [Vp,kmag, waveL, waveE, waveThe, Fre, W1, W2, W3, W4,kx,ky,kz] = WaveAna_4SC_fast(varargin)
%
% [Vp,kpara, kperp, waveL, waveE, waveThe, Fre, W1, W2, W3, W4] = WaveAna_4SC('E?','R?','B?',Tints,options...)
%
%
% Input: (All data must be in TSeries format)
%       E? -        Fields to apply 4SC cross-spectral analysis to. E.g., E
%                   or B fields (if multiple components only the first is used). 
%       R? -        Positions of the four spacecraft
%       B? -        Background magnetic field in the same coordinates as R?. 
%                   Used to determine the parallel and perpendicular wave numebers using 4SC average.
%       Tints -     Time interval over which the power spectrum is calculated. 
%                   To avoid boundary effects use a longer time interval
%                   for E? and B?. 
%
% Options:
%       linear -    Linearly spaced frequencies. Set number to df (default is logarithmic spacing).
%       numf -      Set number of frequencies used in spectrogram.
%       wwidth -    Multiplier for Morlet wavelet width. Default is 1.
%       linear -    Use linear spacing between frequencies. df argument is
%                   required.
%       sn  -    set number for FGM data smoothing，尽量单数
%       cav  -   降低采样率倍数，加快计算速度，FGM数据慎用
% Output:
%       Vp    -     phase velocity km/s
%       k     -     wavenumber   1/km
%       waveL -     wave length, km
%       waveE -     power spectrum
%       waveThe -   propagation angle
%       Fre      -    frequence
%% 初始化数据
if (nargin < 4)
	help mms.fk_powerspec4SC;
	powerxy = NaN;
	xvariable = NaN;
	yvariable = NaN;
	return;
end

ic = 1:4;

c_eval('E?=evalin(''base'',irf_ssub(varargin{1},?));',ic);
c_eval('R?=evalin(''base'',irf_ssub(varargin{2},?));',ic);
c_eval('B?=evalin(''base'',irf_ssub(varargin{3},?));',ic);
Tint = varargin{4};
c_eval('B? = B?.resample(B1);',2:4);
c_eval('E? = E?.resample(E1);',2:4);
c_eval('R? = R?.resample(E1);',1:4);

time = E1.time;

cav = 8;      % 估计相位数据点
numf = 400;   % 频率 分划
uselinear = 0; %  线性化
wwidth = 1;   % 默认小波变换参数
frange = 0;   % 频率范围

args=varargin(5:end);
if numel(args)>0
	haveoptions=1;
  irf.log('notice','Options were passed.');
else
	haveoptions=0;
end

while haveoptions
	l = 2;
	switch(lower(args{1}))
    case 'cav'
      if numel(args)>1 && isnumeric(args{2})
        cav = args{2};
      end
    case 'sn'
      if numel(args)>1 && isnumeric(args{2})
        sn = args{2};      
      end     
    case 'numf'
      if numel(args)>1 && isnumeric(args{2})
        numf = floor(args{2});
      end
    case 'linear'
      if numel(args)>1 && isnumeric(args{2})
        df = args{2}; 
        uselinear = 1;
        irf.log('notice','Using linearly spaced frequencies');
      end
    case 'wwidth'
      if numel(args)>1 && isnumeric(args{2})
        wwidth = args{2};
      end
    case 'frange' 
       if numel(args)>1 && isnumeric(args{2})
        frange = args{2};      
       end   
       
    otherwise
      irf.log('warning',['Unknown flag: ' args{1}]);
      l=1;
      break;
  end
    args = args(l+1:end);
	if isempty(args), haveoptions=0; end
end

%% 平滑磁场
Bav = irf.ts_vec_xyz(B1.time,(B1.data+B2.data+B3.data+B4.data)/4); % 直流 磁场平均
% Bav=Bav.filt(0,0.05,128,5);% 滤波 低通 (0,0.05)hz 0.05hz以上都算波 0.05hz以下当做背景磁场
Bav=Bav.filt(0,0.1,128,5); %20170619

idx = tlim(E1.time,Tint);  %idx 为T2的时间序列个数

% If odd, remove last data point (as is done in irf_wavelet)
if mod(length(idx),2)
    idx(end)=[];
end
%% 小波变换
if ~uselinear
	c_eval('W? = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''nf'',numf,''wavelet_width'',5.36*wwidth,''f'',[frange(1) frange(2)]);',ic);
else
  c_eval('W? = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''linear'',df,''wavelet_width'',5.36*wwidth);',ic);
end
%% 
Fre=W1.f;  % 频率分划
numf = length(W1.f); 

L = length(idx);  % 截止后的时间长度
times=time(idx);
c_eval('W?.p = {W?.p{1,1}(idx,:)};',ic);
c_eval('W?.t = times;',ic);%将W中的t改为秒

%% 主体计算

fkPower = 0.25*(cell2mat(W1.p).*conj(cell2mat(W1.p)) + cell2mat(W2.p).*conj(cell2mat(W2.p)) ...
   + cell2mat(W3.p).*conj(cell2mat(W3.p)) + cell2mat(W4.p).*conj(cell2mat(W4.p))); %将w中的复数 点乘自己的共轭 得到power能量


N = floor(L/cav)-1;
posav = round(cav/2) + (0:1:N)*cav;
avtimes = times(posav);  %% cav个点计算一次，重新定义时间序列 （数据降频）

Bav = Bav.resample(avtimes); 
c_eval('R? = R?.resample(avtimes);',ic);
%% 相位差，能量计算
cx12 = zeros(N+1,numf);
cx13 = zeros(N+1,numf);
cx14 = zeros(N+1,numf);
cx23 = zeros(N+1,numf);
cx24 = zeros(N+1,numf);
cx34 = zeros(N+1,numf);
waveE = zeros(N+1,numf);

% cx12 = double(W1.p{1,1}.*conj(W2.p{1,1}));  % 相位差
% cx13 = double(W1.p{1,1}.*conj(W3.p{1,1}));
% cx14 = double(W1.p{1,1}.*conj(W4.p{1,1}));
% cx23 = double(W2.p{1,1}.*conj(W3.p{1,1}));
% cx24 = double(W2.p{1,1}.*conj(W4.p{1,1}));
% cx34 = double(W3.p{1,1}.*conj(W4.p{1,1}));


% % % for m = 1:N
% % %     for n=1:numf
% % %         cx12(m,n) = irf.nanmean(W1.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n).*conj(W2.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n))); % 相位差
% % %         cx13(m,n) = irf.nanmean(W1.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n).*conj(W3.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n)));
% % %         cx14(m,n) = irf.nanmean(W1.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n).*conj(W4.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n)));
% % %         cx23(m,n) = irf.nanmean(W2.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n).*conj(W3.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n)));
% % %         cx24(m,n) = irf.nanmean(W2.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n).*conj(W4.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n)));
% % %         cx34(m,n) = irf.nanmean(W3.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n).*conj(W4.p{1,1}([posav(m)-floor(cav/2)+1:posav(m)+floor(cav/2)],n)));
% % %     end
% % %         waveE(m,:) = irf.nanmean(fkPower([posav(m)-round(cav/2)+1:posav(m)+round(cav/2)],:)); % 能量时间平均、
% % %         clc
% % %         disp([repmat('■',1,round(10*m/N)),repmat('□',1,10-round(10*m/N)),' 啾咪 ⊂(˃̶͈̀ε ˂̶͈́ ⊂ )'])
% % % end  

P1 = full(double(cell2mat(W1.p)));   % [L x numf]
P2 = full(double(cell2mat(W2.p)));
P3 = full(double(cell2mat(W3.p)));
P4 = full(double(cell2mat(W4.p)));

C12 = P1.*conj(P2);
C13 = P1.*conj(P3);
C14 = P1.*conj(P4);
C23 = P2.*conj(P3);
C24 = P2.*conj(P4);
C34 = P3.*conj(P4);

Cx12_sm = movmean(C12, cav, 1, 'omitnan');    % [L x numf]
Cx13_sm = movmean(C13, cav, 1, 'omitnan');
Cx14_sm = movmean(C14, cav, 1, 'omitnan');
Cx23_sm = movmean(C23, cav, 1, 'omitnan');
Cx24_sm = movmean(C24, cav, 1, 'omitnan');
Cx34_sm = movmean(C34, cav, 1, 'omitnan');

E_sm    = movmean(fkPower, cav, 1, 'omitnan'); % [L x numf]

cx12 = Cx12_sm(posav, :);     % [N+1 x numf]
cx13 = Cx13_sm(posav, :);
cx14 = Cx14_sm(posav, :);
cx23 = Cx23_sm(posav, :);
cx24 = Cx24_sm(posav, :);
cx34 = Cx34_sm(posav, :);

waveE = double(E_sm(posav, :));       % [N+1 x numf]
%% Compute phase differences between each spacecraft pair
th12 = atan2(imag(cx12),real(cx12));
th13 = atan2(imag(cx13),real(cx13));
th14 = atan2(imag(cx14),real(cx14));
th23 = atan2(imag(cx23),real(cx23));
th24 = atan2(imag(cx24),real(cx24));
th34 = atan2(imag(cx34),real(cx34));

wmat = 2*pi*ones(N+1,1)*(W1.f)';

% 换算时间延后
dt12 = th12./wmat;
dt13 = th13./wmat;
dt14 = th14./wmat;
dt23 = th23./wmat;
dt24 = th24./wmat;
dt34 = th34./wmat;

% Weighted averaged time delay using all spacecraft pairs
dt2 = 0.5*dt12 + 0.2*(dt13 - dt23) + 0.2*(dt14 - dt24) + 0.1*(dt14 - dt34 - dt23);
dt3 = 0.5*dt13 + 0.2*(dt12 + dt23) + 0.2*(dt14 - dt34) + 0.1*(dt12 + dt24 - dt34);
dt4 = 0.5*dt14 + 0.2*(dt12 + dt24) + 0.2*(dt13 + dt34) + 0.1*(dt12 + dt23 + dt34);

% Compute phase speeds

R1 = R1.data;
R2 = R2.data;
R3 = R3.data;
R4 = R4.data;

kx = zeros(N+1,numf);
ky = zeros(N+1,numf);
kz = zeros(N+1,numf);
Vp = zeros(N+1, numf);
% ** k的方向和Vp的方向是绝对一致的
for ii = 1:N+1
	dR = [R2(ii,:);R3(ii,:);R4(ii,:)]-[R1(ii,:);R1(ii,:);R1(ii,:)];
  for jj = 1:numf
    % Volumetric tensor with SC1 as center.
    m = dR\[dt2(ii,jj);dt3(ii,jj);dt4(ii,jj)]; %相速度
    Vp(ii,jj)=1/norm(m); %相速度的模分之一 无方向
    kx(ii,jj) = 2*pi*W1.f(jj)*m(1);
    ky(ii,jj) = 2*pi*W1.f(jj)*m(2);
    kz(ii,jj) = 2*pi*W1.f(jj)*m(3);  % v=f/k,k=2pi*w/v
  end
end

kx = kx/1e3; % 标准单位
ky = ky/1e3;
kz = kz/1e3;
kmag = sqrt(kx.*kx + ky.*ky + kz.*kz);  % 波数数值

Bavxmat = Bav.x.data*ones(1,numf);
Bavymat = Bav.y.data*ones(1,numf);
Bavzmat = Bav.z.data*ones(1,numf);
Bavabsmat = Bav.abs.data*ones(1,numf);

kpara = (kx.*Bavxmat + ky.*Bavymat + kz.*Bavzmat)./Bavabsmat;  % 平行 k 
kperp = sqrt(kmag.^2 - kpara.^2);   % 垂直 k

kmax = max(max(kmag))*1.1;
kmin = -kmax;
kvec = linspace(-kmax,kmax,500);
kmagvec = linspace(0,kmax,500);
dkmag = kmax/500;
dk = 2*kmax/500;

%% 计算速度
Vpx = zeros(N+1, numf);
Vpy = zeros(N+1, numf);
Vpz = zeros(N+1, numf);
for ii = 1:N+1
  for jj = 1:numf
    Vpx(ii,jj)=Vp(ii,jj)*kx(ii,jj)/kmag(ii,jj);
     Vpy(ii,jj)=Vp(ii,jj)*ky(ii,jj)/kmag(ii,jj);
      Vpz(ii,jj)=Vp(ii,jj)*kz(ii,jj)/kmag(ii,jj); % 相速度大小*k的单位向量 及相速度的向量
  end
end



%% 计算波长
waveL = zeros(N+1,numf);
for ii=1:N+1
    for jj=1:numf
        waveL(ii,jj)=2*pi/kmag(ii,jj);
    end
end


%% 计算传播角度

waveThe = zeros(N+1,numf);
for ii=1:N+1
    for jj=1:numf
        
%         V_1=Vp1{ii,jj};
        V_1=[Vpx(ii,jj) Vpy(ii,jj) Vpz(ii,jj)];
        V_B=[Bavxmat(ii,jj) Bavymat(ii,jj) Bavzmat(ii,jj)];
        
        V_1=V_1/norm(V_1); %相速度单位方向 
        V_B=V_B/norm(V_B); %背景磁场单位方向 
        waveThe(ii,jj)=acos(dot(V_1, V_B))/pi*180;
    end
%     ii
end


avtimes=irf_time(avtimes,'EpochTT>epoch');


Vp=[avtimes, Vp];
kpara=[avtimes,kpara];
kperp=[avtimes,kperp];
waveL=[avtimes,waveL];
waveE=[avtimes, waveE];
waveThe=[avtimes, waveThe];

W1 = [irf_time(W1.t,'EpochTT>epoch')  double(cell2mat(W1.p))];
W2 = [irf_time(W2.t,'EpochTT>epoch')  double(cell2mat(W2.p))];
W3 = [irf_time(W3.t,'EpochTT>epoch')  double(cell2mat(W3.p))];
W4 = [irf_time(W4.t,'EpochTT>epoch')  double(cell2mat(W4.p))];

W1=W1(posav, :);
W2=W2(posav, :);
W3=W3(posav, :);
W4=W4(posav, :);

return
end
%% 
function Y0 = replaceNaNWithZero(Yin)
    Y0 = Yin;
    Y0(isnan(Y0)) = 0;
end
%%
function Z = nanwinmean(S, Y)
    % S: N x T 稀疏矩阵（double）
    % Y: T x numf，可能是 single/double，实数或复数
    cls = class(Y);               % 记录原始精度 single/double

    % 转 double 并做 NaN->0
    Yd = double(Y);
    mask = ~isnan(Yd);            % logical，NaN 位置为 0
    Yd(~mask) = 0;

    % 窗口求和（分子）与非NaN计数（分母），均为 double
    num = S * Yd;                 % N x numf
    den = S * double(mask);       % N x numf
    den(den == 0) = 1;            % 避免除零（对应全NaN窗口）

    Z = num ./ den;               % double 结果

    % 转回原始精度（保留复数）
    if strcmp(cls,'single')
        Z = single(Z);
    end
end
