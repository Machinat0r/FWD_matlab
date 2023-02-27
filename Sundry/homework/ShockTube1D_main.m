clc; clear;close all
%% initial condition
% mesh generation
N = 1000;
dx = 1;
x = 0:dx:dx*N;

% initial condition
global gamma mu0
gamma = 5/3;
mu0 = 4*pi;
Xleft = x<=0.5*max(x);Xright = x>0.5*max(x);
Bx0 = 0.75*sqrt(4*pi)*ones(size(x));
Rho0 = 1*Xleft + 0.125*Xright;
Ux0   = 0*Xleft + 0*Xright;
Uy0   = 0*Xleft + 0*Xright;
Uz0   = 0*Xleft + 0*Xright;
P0   = 1*Xleft + 0.1*Xright;
By0 = sqrt(4*pi)*Xleft - sqrt(4*pi)*Xright;
Bz0 = 0*Xleft + 0*Xright;
E = P0./((gamma-1).*Rho0) + 0.5*(Ux0.^2+Uy0.^2+Uz0.^2) + 0.5*(Bx0.^2+By0.^2+...
    Bz0.^2)./(mu0.*Rho0);

% CFL
T1 = 50; T2 = 150; T3 = 250;
CFL = 0.15;
dt = dx*CFL;
Times1 = ceil(T1/dt);Times2 = ceil(T2/dt);Times3 = ceil(T3/dt);
%% evaluation U & F
U(1,:) = Rho0;
U(2,:) = Rho0.*Ux0;
U(3,:) = Rho0.*Uy0;
U(4,:) = Rho0.*Uz0;
U(5,:) = Bx0;
U(6,:) = By0;
U(7,:) = Bz0;
U(8,:) = Rho0.*E;

F(1,:) = Rho0.*Ux0;
F(2,:) = Rho0.*Ux0.^2 + P0 + 0.5*(-Bx0.^2+By0.^2+Bz0.^2)/mu0;
F(3,:) = Rho0.*Ux0.*Uy0 - Bx0.*By0/mu0;
F(4,:) = Rho0.*Ux0.*Uz0 - Bx0.*Bz0/mu0;
F(5,:) = 0;
F(6,:) = Ux0.*By0 - Uy0.*Bx0;
F(7,:) = Ux0.*Bz0 - Uz0.*Bx0;
F(8,:) = (Rho0.*E + P0 + 0.5*(Bx0.^2+By0.^2+Bz0.^2)/mu0).*Ux0 - ...
    Bx0.*(Ux0.*Bx0+Uy0.*By0+Uz0.*Bz0)./mu0;

%% plot
[U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_numsol(U,F,Times1,'Lax-Friedrichs',dt,dx);
[U2,F2,Rho2,Ux2,Uy2,Uz2,Bx2,By2,Bz2,P2] = ShockTube1D_numsol(U,F,Times2,'Lax-Friedrichs',dt,dx);
[U3,F3,Rho3,Ux3,Uy3,Uz3,Bx3,By3,Bz3,P3] = ShockTube1D_numsol(U,F,Times3,'Lax-Friedrichs',dt,dx);
% [U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_numsol(U,F,Times1,'MacCormack',dt,dx,0.3);
% [U2,F2,Rho2,Ux2,Uy2,Uz2,Bx2,By2,Bz2,P2] = ShockTube1D_numsol(U,F,Times2,'MacCormack',dt,dx,0.3);
% [U3,F3,Rho3,Ux3,Uy3,Uz3,Bx3,By3,Bz3,P3] = ShockTube1D_numsol(U,F,Times3,'MacCormack',dt,dx,0.3);
% [U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_numsol(U,F,Times1,'Roe',dt,dx);
% [U2,F2,Rho2,Ux2,Uy2,Uz2,Bx2,By2,Bz2,P2] = ShockTube1D_numsol(U,F,Times2,'Roe',dt,dx);
% [U3,F3,Rho3,Ux3,Uy3,Uz3,Bx3,By3,Bz3,P3] = ShockTube1D_numsol(U,F,Times3,'Roe',dt,dx);
% [U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_numsol(U,F,Times1,'TVD',dt,dx,0.1);
% [U2,F2,Rho2,Ux2,Uy2,Uz2,Bx2,By2,Bz2,P2] = ShockTube1D_numsol(U,F,Times2,'TVD',dt,dx,0.1);
% [U3,F3,Rho3,Ux3,Uy3,Uz3,Bx3,By3,Bz3,P3] = ShockTube1D_numsol(U,F,Times3,'TVD',dt,dx,0.1);

figure(1)
cor = {'#F3D9BE','#EF9163','#CD4432','#611E1E'};
set(gcf,'Position',[0,-200,800,1000])
subplot(4,2,1);
plot(x,Rho0,'color',cor{1});hold on;
plot(x,Rho1,'color',cor{2});hold on;
plot(x,Rho2,'color',cor{3});hold on;
plot(x,Rho3,'color',cor{4});hold on;
legend(['Time:',num2str(0)],['Time:',num2str(T1)],['Time:',num2str(T2)],['Time:',num2str(T3)],'Location','southwest')
xlabel('x')
ylabel('\rho');
ylim([0 1.1])
grid off

subplot(4,2,2);
plot(x,Ux0,'color',cor{1});hold on;
plot(x,Ux1,'color',cor{2});hold on;
plot(x,Ux2,'color',cor{3});hold on;
plot(x,Ux3,'color',cor{4});hold on;
xlabel('x')
ylabel('U_x');
ylim([-0.8 0.8])
grid off

subplot(4,2,3);
plot(x,Uy0,'color',cor{1});hold on;
plot(x,Uy1,'color',cor{2});hold on;
plot(x,Uy2,'color',cor{3});hold on;
plot(x,Uy3,'color',cor{4});hold on;
xlabel('x')
ylabel('U_y');
ylim([-2.2 0.2])
grid off

subplot(4,2,4);
plot(x,Uz0,'color',cor{1});hold on;
plot(x,Uz1,'color',cor{2});hold on;
plot(x,Uz2,'color',cor{3});hold on;
plot(x,Uz3,'color',cor{4});hold on;
xlabel('x')
ylabel('U_z');
ylim([-1e-17 1e-17])
grid off

subplot(4,2,5);
plot(x,Bx0,'color',cor{1});hold on;
plot(x,Bx1,'color',cor{2});hold on;
plot(x,Bx2,'color',cor{3});hold on;
plot(x,Bx3,'color',cor{4});hold on;
xlabel('x')
ylabel('B_x');
grid off

subplot(4,2,6);
plot(x,By0,'color',cor{1});hold on;
plot(x,By1,'color',cor{2});hold on;
plot(x,By2,'color',cor{3});hold on;
plot(x,By3,'color',cor{4});hold on;
xlabel('x')
ylabel('B_y');
grid off

subplot(4,2,7);
plot(x,Bz0,'color',cor{1});hold on;
plot(x,Bz1,'color',cor{2});hold on;
plot(x,Bz2,'color',cor{3});hold on;
plot(x,Bz3,'color',cor{4});hold on;
xlabel('x')
ylabel('B_z');
ylim([-1e-17 1e-17])
grid off

subplot(4,2,8);
plot(x,P0,'color',cor{1});hold on;
plot(x,P1,'color',cor{2});hold on;
plot(x,P2,'color',cor{3});hold on;
plot(x,P3,'color',cor{4});hold on;
xlabel('x')
ylabel('P');
ylim([0 1.1])
grid off