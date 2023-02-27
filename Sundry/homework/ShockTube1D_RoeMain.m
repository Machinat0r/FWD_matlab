clc;clear;close all

%初始化变量，dt
n = 1001;%1001个结点，1000个网格
dx = 1;%空间步长
dt = 0.15;%时间步长
t = 100;
global r mu
mu = 4*pi;%磁导率
r = 5/3;
Bx = 0.75*sqrt(4*pi)*ones(1, n);
By = sqrt(4*pi)*ones(1, n);
for i = 501:n
    By(i) = -By(i);
end
Bz = zeros(1, n);
rho = ones(1, n);
rho(501:n) = 0.125;
u = zeros(1, n);
v = zeros(1, n);
w = zeros(1, n);
p = ones(1, n);
p(501:n) = 0.1;
et = p./(r-1)./rho + 0.5*(u.^2 + v.^2 + w.^2)  + (Bx.^2 + By.^2 + Bz.^2)/(2*mu)./rho; 
H=r*p/(r-1)./rho+0.5*(u.^2+v.^2+w.^2)+(Bx.^2+By.^2+Bz.^2)/mu./rho;
U = zeros(7, n);
F = zeros(7, n);
U(1, :) = rho;
U(2, :) = rho.*u;
U(3, :) = rho.*v;
U(4, :) = rho.*w;
U(5, :) = By;
U(6, :) = Bz;
U(7, :) = rho.*et;
F(1, :) = rho.*u;
F(2, :) = rho.*u.*u + p + (-Bx.^2 + By.^2 + Bz.^2)/(2*mu);
F(3, :) = rho.*u.*v - Bx.*By/mu;
F(4, :) = rho.*u.*w - Bx.*Bz/mu;
F(5, :) = u.*By - v.*Bx;
F(6, :) = u.*Bz - w.*Bx;
F(7, :) = (rho.*et + p + (Bx.^2 + By.^2 + Bz.^2)/(2*mu)).*u - Bx.*(u.*Bx + v.*By + w.*Bz)/mu;
f = zeros(7, n);%数值通量

%计算


x = 0:dx:1000;
T1 = 100;
[U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_Roe(U,F,T1,dt,dx,Bx,n);

%%
figure(1)
n=8;
cor = {'#611E1E'};
set(gcf,'Position',[0,-200,600,1000])
i=1;

h(i)=irf_subplot(n,1,-i);
plot(x,Rho1,'color',cor{1});hold on;
xlabel('x')
ylabel('\rho');
ylim([0 1.1])
grid off
i = i+1;

h(i)=irf_subplot(n,1,-i);
plot(x,Ux1,'color',cor{1});hold on;
xlabel('x')
ylabel('U_x');
ylim([-0.8 0.8])
grid off
i = i+1;

h(i)=irf_subplot(n,1,-i);
i = i+1;
plot(x,Uy1,'color',cor{1});hold on;
xlabel('x')
ylabel('U_y');
ylim([-2.2 0.2])
grid off

h(i)=irf_subplot(n,1,-i);
i = i+1;
plot(x,Uz1,'color',cor{1});hold on;
xlabel('x')
ylabel('U_z');
ylim([-1 1])
grid off

h(i)=irf_subplot(n,1,-i);
i = i+1;
plot(x,Bx1,'color',cor{1});hold on;
xlabel('x')
ylabel('B_x');
grid off

h(i)=irf_subplot(n,1,-i);
i = i+1;
plot(x,By1,'color',cor{1});hold on;
xlabel('x')
ylabel('B_y');
grid off

h(i)=irf_subplot(n,1,-i);
i = i+1;
plot(x,Bz1,'color',cor{1});hold on;
xlabel('x')
ylabel('B_z');
ylim([-1 1])
grid off

h(i)=irf_subplot(n,1,-i);
i = i+1;
plot(x,P1,'color',cor{1});hold on;
xlabel('x')
ylabel('P');
ylim([0 1.1])
grid off

%%
% T1 = 50;T2 = 150;T3 = 250;
% [U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_Roe(U,F,T1,dt,dx,Bx,n);
% [U2,F2,Rho2,Ux2,Uy2,Uz2,Bx2,By2,Bz2,P2] = ShockTube1D_Roe(U,F,T2,dt,dx,Bx,n);
% [U3,F3,Rho3,Ux3,Uy3,Uz3,Bx3,By3,Bz3,P3] = ShockTube1D_Roe(U,F,T3,dt,dx,Bx,n);
% 
% figure(1)
% cor = {'#F3D9BE','#EF9163','#CD4432','#611E1E'};
% set(gcf,'Position',[0,-200,800,1000])
% subplot(4,2,1);
% plot(x,rho,'color',cor{1});hold on;
% plot(x,Rho1,'color',cor{2});hold on;
% plot(x,Rho2,'color',cor{3});hold on;
% plot(x,Rho3,'color',cor{4});hold on;
% legend(['Time:',num2str(0)],['Time:',num2str(T1)],['Time:',num2str(T2)],['Time:',num2str(T3)],'Location','southwest')
% xlabel('x')
% ylabel('\rho');
% ylim([0 1.1])
% grid off
% 
% subplot(4,2,2);
% plot(x,u,'color',cor{1});hold on;
% plot(x,Ux1,'color',cor{2});hold on;
% plot(x,Ux2,'color',cor{3});hold on;
% plot(x,Ux3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('U_x');
% ylim([-0.8 0.8])
% grid off
% 
% subplot(4,2,3);
% plot(x,v,'color',cor{1});hold on;
% plot(x,Uy1,'color',cor{2});hold on;
% plot(x,Uy2,'color',cor{3});hold on;
% plot(x,Uy3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('U_y');
% ylim([-2.2 0.2])
% grid off
% 
% subplot(4,2,4);
% plot(x,w,'color',cor{1});hold on;
% plot(x,Uz1,'color',cor{2});hold on;
% plot(x,Uz2,'color',cor{3});hold on;
% plot(x,Uz3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('U_z');
% % ylim([-1e-17 1e-17])
% grid off
% 
% subplot(4,2,5);
% plot(x,Bx,'color',cor{1});hold on;
% plot(x,Bx1,'color',cor{2});hold on;
% plot(x,Bx2,'color',cor{3});hold on;
% plot(x,Bx3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('B_x');
% grid off
% 
% subplot(4,2,6);
% plot(x,By,'color',cor{1});hold on;
% plot(x,By1,'color',cor{2});hold on;
% plot(x,By2,'color',cor{3});hold on;
% plot(x,By3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('B_y');
% grid off
% 
% subplot(4,2,7);
% plot(x,Bz,'color',cor{1});hold on;
% plot(x,Bz1,'color',cor{2});hold on;
% plot(x,Bz2,'color',cor{3});hold on;
% plot(x,Bz3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('B_z');
% % ylim([-1e-17 1e-17])
% grid off
% 
% subplot(4,2,8);
% plot(x,p,'color',cor{1});hold on;
% plot(x,P1,'color',cor{2});hold on;
% plot(x,P2,'color',cor{3});hold on;
% plot(x,P3,'color',cor{4});hold on;
% xlabel('x')
% ylabel('P');
% ylim([0 1.1])
% grid off
