%% ================= MDD / MGA Fig.9-like plot =================
% 需要已有变量：
% R1 R2 R3 R4: 四颗卫星位置，TSeries 或 [time x y z]，单位 km
% B1 B2 B3 B4: 四颗卫星磁场，TSeries 或 [time Bx By Bz]，单位 nT
%
% 如果你还没有读取数据，可参考：
% tint = irf.tint('2015-10-16T13:04:28.0Z/2015-10-16T13:04:31.0Z');
% c_eval('B? = mms.get_data(''B_gsm_brst'',tint,?);',1:4);
% c_eval('R? = mms.get_data(''R_gsm'',tint,?);',1:4);

method = 'MGA';        % 'MGA' 更接近文章 Fig.9；'MDD' 为真正 MDD
coordName = 'GSM';     % 改成 'GSE' 也可以
planeStr = 'xy';       % 轨迹投影平面：'xy', 'xz', 'yz'
deltaB = 0.05;         % nT，用于画 deltaB/lmax 参考线，可按仪器误差修改

%% ---------- convert TSeries to matrix ----------
if isa(B1,'TSeries'); B1m = irf.ts2mat(B1); else; B1m = B1; end
if isa(B2,'TSeries'); B2m = irf.ts2mat(B2); else; B2m = B2; end
if isa(B3,'TSeries'); B3m = irf.ts2mat(B3); else; B3m = B3; end
if isa(B4,'TSeries'); B4m = irf.ts2mat(B4); else; B4m = B4; end

if isa(R1,'TSeries'); R1m = irf.ts2mat(R1); else; R1m = R1; end
if isa(R2,'TSeries'); R2m = irf.ts2mat(R2); else; R2m = R2; end
if isa(R3,'TSeries'); R3m = irf.ts2mat(R3); else; R3m = R3; end
if isa(R4,'TSeries'); R4m = irf.ts2mat(R4); else; R4m = R4; end

B1m = B1m(:,1:4);
B2m = B2m(:,1:4);
B3m = B3m(:,1:4);
B4m = B4m(:,1:4);

R1m = R1m(:,1:4);
R2m = R2m(:,1:4);
R3m = R3m(:,1:4);
R4m = R4m(:,1:4);

%% ---------- resample everything to MMS1 B time ----------
B2m = irf_resamp(B2m,B1m);
B3m = irf_resamp(B3m,B1m);
B4m = irf_resamp(B4m,B1m);

R1m = irf_resamp(R1m,B1m);
R2m = irf_resamp(R2m,B1m);
R3m = irf_resamp(R3m,B1m);
R4m = irf_resamp(R4m,B1m);

ok = all(isfinite([B1m(:,2:4),B2m(:,2:4),B3m(:,2:4),B4m(:,2:4), ...
                   R1m(:,2:4),R2m(:,2:4),R3m(:,2:4),R4m(:,2:4)]),2);

B1m = B1m(ok,:);
B2m = B2m(ok,:);
B3m = B3m(ok,:);
B4m = B4m(ok,:);

R1m = R1m(ok,:);
R2m = R2m(ok,:);
R3m = R3m(ok,:);
R4m = R4m(ok,:);

%% ---------- calculate gradient tensor with c_4_grad ----------
% c_4_grad 输出顺序：
% dxBx dxBy dxBz dyBx dyBy dyBz dzBx dzBy dzBz

try
    gradB = c_4_grad(R1m,R2m,R3m,R4m,B1m,B2m,B3m,B4m,'grad');
catch
    gradB = c_4_grad('R?m','B?m','grad');
end

if isa(gradB,'TSeries')
    gradBm = irf.ts2mat(gradB);
else
    gradBm = gradB;
end

t = gradBm(:,1);
Gdata = gradBm(:,2:10);
nt = size(Gdata,1);

B1g = irf_resamp(B1m,gradBm);
R1g = irf_resamp(R1m,gradBm);
R2g = irf_resamp(R2m,gradBm);
R3g = irf_resamp(R3m,gradBm);
R4g = irf_resamp(R4m,gradBm);

%% ---------- MDD or MGA eigen analysis ----------
lambda = nan(nt,3);
Nmax = nan(nt,3);
Nmid = nan(nt,3);
Nmin = nan(nt,3);

divB = nan(nt,1);
curlB = nan(nt,3);
Qcurl = nan(nt,1);
Qgrad = nan(nt,1);

for ii = 1:nt

    G = [Gdata(ii,1) Gdata(ii,2) Gdata(ii,3); ...
         Gdata(ii,4) Gdata(ii,5) Gdata(ii,6); ...
         Gdata(ii,7) Gdata(ii,8) Gdata(ii,9)];

    if any(~isfinite(G),'all')
        continue
    end

    switch upper(method)
        case 'MDD'
            Lmat = G*G.';
        case 'MGA'
            Lmat = G.'*G;
        otherwise
            error('method must be MDD or MGA');
    end

    Lmat = 0.5*(Lmat + Lmat.');

    [V,D] = eig(Lmat);
    lam = real(diag(D));

    [lamSort,idx] = sort(lam,'descend');
    Vsort = real(V(:,idx));

    lambda(ii,:) = lamSort(:).';

    Nmax(ii,:) = Vsort(:,1).';
    Nmid(ii,:) = Vsort(:,2).';
    Nmin(ii,:) = Vsort(:,3).';

    divB(ii) = G(1,1) + G(2,2) + G(3,3);

    curlB(ii,:) = [G(2,3) - G(3,2), ...
                   G(3,1) - G(1,3), ...
                   G(1,2) - G(2,1)];

    curlNorm = sqrt(sum(curlB(ii,:).^2));
    gradMax = max(abs(G),[],'all');

    Qcurl(ii) = abs(divB(ii))./curlNorm;
    Qgrad(ii) = abs(divB(ii))./gradMax;
end

%% ---------- remove arbitrary sign flips of eigenvectors ----------
for ii = 2:nt
    if all(isfinite(Nmax(ii,:))) && all(isfinite(Nmax(ii-1,:)))
        if dot(Nmax(ii,:),Nmax(ii-1,:)) < 0
            Nmax(ii,:) = -Nmax(ii,:);
        end
    end

    if all(isfinite(Nmid(ii,:))) && all(isfinite(Nmid(ii-1,:)))
        if dot(Nmid(ii,:),Nmid(ii-1,:)) < 0
            Nmid(ii,:) = -Nmid(ii,:);
        end
    end

    if all(isfinite(Nmin(ii,:))) && all(isfinite(Nmin(ii-1,:)))
        if dot(Nmin(ii,:),Nmin(ii-1,:)) < 0
            Nmin(ii,:) = -Nmin(ii,:);
        end
    end
end

sqrtLam = sqrt(max(lambda,0));

D1 = (sqrtLam(:,1) - sqrtLam(:,2))./sqrtLam(:,1);
D2 = (sqrtLam(:,2) - sqrtLam(:,3))./sqrtLam(:,1);
D3 = sqrtLam(:,3)./sqrtLam(:,1);

D1(~isfinite(D1)) = nan;
D2(~isfinite(D2)) = nan;
D3(~isfinite(D3)) = nan;

%% ---------- spacecraft separation and error reference line ----------
r12 = sqrt(sum((R1g(:,2:4) - R2g(:,2:4)).^2,2));
r13 = sqrt(sum((R1g(:,2:4) - R3g(:,2:4)).^2,2));
r14 = sqrt(sum((R1g(:,2:4) - R4g(:,2:4)).^2,2));
r23 = sqrt(sum((R2g(:,2:4) - R3g(:,2:4)).^2,2));
r24 = sqrt(sum((R2g(:,2:4) - R4g(:,2:4)).^2,2));
r34 = sqrt(sum((R3g(:,2:4) - R4g(:,2:4)).^2,2));

lmax = max([r12 r13 r14 r23 r24 r34],[],2);
errLine = deltaB ./ lmax;

%% ---------- Fig.9-like plot ----------
h = irf_plot(7,'newfigure');
set(gcf,'Color','w','Position',[80 40 980 1050]);

% 颜色顺序：蓝、绿、红、黑
cBlue  = [0 0 1];
cGreen = [0 0.6 0];
cRed   = [1 0 0];
cBlack = [0 0 0];

%% ---------- Panel a: B ----------
irf_plot(h(1),[B1g(:,1) B1g(:,2:4)]);
hold(h(1),'on');

Babs = sqrt(sum(B1g(:,2:4).^2,2));
irf_plot(h(1),[B1g(:,1) Babs]);

lines = findobj(h(1),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);
set(lines(4),'Color',cBlack);

ylabel(h(1),['B ' coordName ' [nT]']);

leg = irf_legend(h(1),{'B_x','B_y','B_z','|B|'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);
set(leg(4),'Color',cBlack);

%% ---------- Panel b: sqrt eigenvalues ----------
irf_plot(h(2),[t sqrtLam]);
hold(h(2),'on');

% 如果你希望第4条黑线为误差参考线，保留下面这一行
irf_plot(h(2),[t errLine]);

set(h(2),'YScale','lin');
ylabel(h(2),'\surd\lambda [nT/km]');
ylim(h(2),[1e-5 1.5e-1]);

lines = findobj(h(2),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);
set(lines(4),'Color',cBlack);
set(lines(4),'LineStyle','--');

leg = irf_legend(h(2),{'\surd\lambda_{max}','\surd\lambda_{mid}','\surd\lambda_{min}','\deltaB/l_{max}'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);
set(leg(4),'Color',cBlack);

%% ---------- Panel c: dimensionality indices ----------
irf_plot(h(3),[t D1 D2 D3]);
ylabel(h(3),'D index');
ylim(h(3),[0 1]);

lines = findobj(h(3),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);

leg = irf_legend(h(3),{'D_1','D_2','D_3'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);

%% ---------- Panel d: maximum direction ----------
irf_plot(h(4),[t Nmax]);
ylabel(h(4),'N_{max}');
ylim(h(4),[-1 1]);

lines = findobj(h(4),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);

leg = irf_legend(h(4),{'x','y','z'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);

%% ---------- Panel e: intermediate direction ----------
irf_plot(h(5),[t Nmid]);
ylabel(h(5),'N_{mid}');
ylim(h(5),[-1 1]);

lines = findobj(h(5),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);

leg = irf_legend(h(5),{'x','y','z'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);

%% ---------- Panel f: minimum direction ----------
irf_plot(h(6),[t Nmin]);
ylabel(h(6),'N_{min}');
ylim(h(6),[-1 1]);

lines = findobj(h(6),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);

leg = irf_legend(h(6),{'x','y','z'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);

%% ---------- Panel g: quality indicators ----------
irf_plot(h(7),[t Qcurl Qgrad]);
hold(h(7),'on');

irf_plot(h(7),[t 0.4*ones(nt,1)]);
irf_plot(h(7),[t 0.6*ones(nt,1)]);

ylabel(h(7),'quality');
ylim(h(7),[0 1]);

lines = findobj(h(7),'Type','line');
lines = flipud(lines);

set(lines(1),'Color',cBlue);
set(lines(2),'Color',cGreen);
set(lines(3),'Color',cRed);
set(lines(4),'Color',cBlack);

set(lines(3),'LineStyle','--');
set(lines(4),'LineStyle',':');

leg = irf_legend(h(7),{'|\nabla\cdotB|/|\nabla\timesB|','|\nabla\cdotB|/max|dB_i/dj|','0.4','0.6'},[0.98 0.90]);
set(leg(1),'Color',cBlue);
set(leg(2),'Color',cGreen);
set(leg(3),'Color',cRed);
set(leg(4),'Color',cBlack);

%% ---------- figure formatting ----------
for ii = 1:7
    grid(h(ii),'off');
    set(h(ii),'FontSize',10);
end

irf_zoom(h(1:7),'x',[t(1) t(end)]);
irf_plot_axis_align(h(1:7));
irf_timeaxis(h(7),'usefig');

title(h(1),[upper(method) ' analysis using c\_4\_grad, ' coordName]);
irf_pl_number_subplots(h,[0.02 0.90],'fontsize',12);
% 
% % Panel h: spacecraft trajectory projection
% axes(h(8));
% cla(h(8));
% hold(h(8),'on');
% 
% R0 = mean([R1g(1,2:4); R2g(1,2:4); R3g(1,2:4); R4g(1,2:4)],1);
% 
% R1rel = R1g(:,2:4) - R0;
% R2rel = R2g(:,2:4) - R0;
% R3rel = R3g(:,2:4) - R0;
% R4rel = R4g(:,2:4) - R0;
% 
% switch lower(planeStr)
%     case 'xy'
%         ix = 1; iy = 2; xlab = '\DeltaX [km]'; ylab = '\DeltaY [km]';
%     case 'xz'
%         ix = 1; iy = 3; xlab = '\DeltaX [km]'; ylab = '\DeltaZ [km]';
%     case 'yz'
%         ix = 2; iy = 3; xlab = '\DeltaY [km]'; ylab = '\DeltaZ [km]';
%     otherwise
%         error('planeStr must be xy, xz, or yz');
% end
% 
% plot(R1rel(:,ix),R1rel(:,iy),'k','LineWidth',1.2);
% plot(R2rel(:,ix),R2rel(:,iy),'r','LineWidth',1.2);
% plot(R3rel(:,ix),R3rel(:,iy),'g','LineWidth',1.2);
% plot(R4rel(:,ix),R4rel(:,iy),'b','LineWidth',1.2);
% 
% plot(R1rel(1,ix),R1rel(1,iy),'ko','MarkerFaceColor','k');
% plot(R2rel(1,ix),R2rel(1,iy),'ro','MarkerFaceColor','r');
% plot(R3rel(1,ix),R3rel(1,iy),'go','MarkerFaceColor','g');
% plot(R4rel(1,ix),R4rel(1,iy),'bo','MarkerFaceColor','b');
% 
% plot(R1rel(end,ix),R1rel(end,iy),'ks','MarkerFaceColor','k');
% plot(R2rel(end,ix),R2rel(end,iy),'rs','MarkerFaceColor','r');
% plot(R3rel(end,ix),R3rel(end,iy),'gs','MarkerFaceColor','g');
% plot(R4rel(end,ix),R4rel(end,iy),'bs','MarkerFaceColor','b');
% 
% axis equal;
% box on;
% xlabel(xlab);
% ylabel(ylab);
% legend({'MMS1','MMS2','MMS3','MMS4'},'Location','best');
% title(['SC trajectory in ' upper(planeStr) ' plane']);

%% ---------- figure formatting ----------
for ii = 1:7
    grid(h(ii),'off');
    set(h(ii),'FontSize',10);
end

irf_zoom(h(1:7),'x',[t(1) t(end)]);
irf_plot_axis_align(h(1:7));
irf_timeaxis(h(7),'usefig');

title(h(1),[upper(method) ' analysis using c\_4\_grad, ' coordName]);
irf_pl_number_subplots(h,[0.02 0.90],'fontsize',12);

colormap(jet)
set(gca,"XTickLabelRotation",0)
set(gcf,'render','painters');
set(gcf,'paperpositionmode','auto')

%% ---------- output useful values ----------
MDD_MGA_out.method = method;
MDD_MGA_out.time = t;
MDD_MGA_out.sqrtLambda = sqrtLam;
MDD_MGA_out.lambda = lambda;
MDD_MGA_out.Nmax = Nmax;
MDD_MGA_out.Nmid = Nmid;
MDD_MGA_out.Nmin = Nmin;
MDD_MGA_out.D1 = D1;
MDD_MGA_out.D2 = D2;
MDD_MGA_out.D3 = D3;
MDD_MGA_out.Qcurl = Qcurl;
MDD_MGA_out.Qgrad = Qgrad;
MDD_MGA_out.errLine = errLine;
MDD_MGA_out.gradB = gradBm;