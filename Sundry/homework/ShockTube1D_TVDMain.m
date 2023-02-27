% Harten 2nd TVD?????¨®???????¡è?¡§??????
clear;clc;close all
%% ?¡§???????¨°?¡ã????????
% ?¡§???????¨°
dx=1; dt=0.15;  
x=[0:dx:1000];     % ??????????
U=zeros(7,1001);  % Euler¡¤???¡À???U
Up=zeros(7,1001); % Euler¡¤?????n+1?¡À??¡À???Up
F=zeros(7,1001);  % Euler¡¤???¡À???F
Fe=zeros(7,1001); % Euler¡¤???????¡À???Fe
g=zeros(7,1001);
ge=zeros(7,1001);
ae=zeros(7,1001);
Me=zeros(7,1001);
phe=zeros(7,1001);
Qk=zeros(7,1001);
Lambda=zeros(7,1001);
RR=zeros(7,7,1001);
rho=zeros(1,1001);    % ????
u=zeros(1,1001);      % x¡¤??¨°????¡¤???
v=zeros(1,1001);      % y¡¤??¨°????¡¤???
w=zeros(1,1001);      % z¡¤??¨°????¡¤???
Bx=zeros(1,1001);     % x¡¤??¨°????¡¤???
By=zeros(1,1001);     % y¡¤??¨°????¡¤???
Bz=zeros(1,1001);     % z¡¤??¨°????¡¤???
p=zeros(1,1001);      % ????
rho_=zeros(1,1001);   % ????Roe???¨´?¡§??????
u_=zeros(1,1001);      
v_=zeros(1,1001);      
w_=zeros(1,1001);      
Bx_=zeros(1,1001);     
By_=zeros(1,1001);     
Bz_=zeros(1,1001);    
p_=zeros(1,1001);
global R mu0
R=5/3;               % ¡À???¡À?
mu0=4*pi;              % ???????????¡§????10

for i=1:1001
    Bx(i)=0.75*sqrt(4*pi);
end
for i=1:500
    rho(i)=1;
    u(i)=0;
    v(i)=0;
    w(i)=0;
    p(i)=1;
    By(i)=sqrt(4*pi);
    Bz(i)=0;
end
for i=501:1001
    rho(i)=0.125;
    u(i)=0;
    v(i)=0;
    w(i)=0;
    p(i)=0.1;
    By(i)=-sqrt(4*pi);
    Bz(i)=0;
end
U(1,:)=rho;
U(2,:)=rho.*u;
U(3,:)=rho.*v;
U(4,:)=rho.*w;
U(5,:)=By;
U(6,:)=Bz;
U(7,:)=rho.*(p./((R-1)*rho)+0.5*(u.^2+v.^2+w.^2)+(Bx.^2+By.^2+Bz.^2)./(2*mu0*rho));
F(1,:)=rho.*u;
F(2,:)=rho.*u.^2+p+(-Bx.^2+By.^2+Bz.^2)/(2*mu0);
F(3,:)=rho.*u.*v-(Bx.*By)/mu0;
F(4,:)=rho.*u.*w-Bx.*Bz/mu0;
F(5,:)=u.*By-v.*Bx;
F(6,:)=u.*Bz-w.*Bx;
F(7,:)=(U(7,:)+p+(Bx.^2+By.^2+Bz.^2)/(2*mu0)).*u-Bx/mu0.*(u.*Bx+v.*By+w.*Bz);

T1 = 50;T2 = 150;T3 = 250;e = 0.1;
[U1,F1,Rho1,Ux1,Uy1,Uz1,Bx1,By1,Bz1,P1] = ShockTube1D_TVD(U,F,T1,dt,dx,Bx,e);
[U2,F2,Rho2,Ux2,Uy2,Uz2,Bx2,By2,Bz2,P2] = ShockTube1D_TVD(U,F,T2,dt,dx,Bx,e);
[U3,F3,Rho3,Ux3,Uy3,Uz3,Bx3,By3,Bz3,P3] = ShockTube1D_TVD(U,F,T3,dt,dx,Bx,e);

figure(1)
cor = {'#F3D9BE','#EF9163','#CD4432','#611E1E'};
set(gcf,'Position',[0,-200,800,1000])
subplot(4,2,1);
plot(x,rho,'color',cor{1});hold on;
plot(x,Rho1,'color',cor{2});hold on;
plot(x,Rho2,'color',cor{3});hold on;
plot(x,Rho3,'color',cor{4});hold on;
legend(['Time:',num2str(0)],['Time:',num2str(T1)],['Time:',num2str(T2)],['Time:',num2str(T3)],'Location','southwest')
xlabel('x')
ylabel('\rho');
ylim([0 1.1])
grid off

subplot(4,2,2);
plot(x,u,'color',cor{1});hold on;
plot(x,-Ux1,'color',cor{2});hold on;
plot(x,-Ux2,'color',cor{3});hold on;
plot(x,-Ux3,'color',cor{4});hold on;
xlabel('x')
ylabel('U_x');
ylim([-0.8 0.8])
grid off

subplot(4,2,3);
plot(x,v,'color',cor{1});hold on;
plot(x,-Uy1,'color',cor{2});hold on;
plot(x,-Uy2,'color',cor{3});hold on;
plot(x,-Uy3,'color',cor{4});hold on;
xlabel('x')
ylabel('U_y');
ylim([-2.2 0.2])
grid off

subplot(4,2,4);
plot(x,w,'color',cor{1});hold on;
plot(x,-Uz1,'color',cor{2});hold on;
plot(x,-Uz2,'color',cor{3});hold on;
plot(x,-Uz3,'color',cor{4});hold on;
xlabel('x')
ylabel('U_z');
% ylim([-1e-17 1e-17])
grid off

subplot(4,2,5);
plot(x,Bx,'color',cor{1});hold on;
plot(x,Bx1,'color',cor{2});hold on;
plot(x,Bx2,'color',cor{3});hold on;
plot(x,Bx3,'color',cor{4});hold on;
xlabel('x')
ylabel('B_x');
grid off

subplot(4,2,6);
plot(x,By,'color',cor{1});hold on;
plot(x,By1,'color',cor{2});hold on;
plot(x,By2,'color',cor{3});hold on;
plot(x,By3,'color',cor{4});hold on;
xlabel('x')
ylabel('B_y');
grid off

subplot(4,2,7);
plot(x,Bz,'color',cor{1});hold on;
plot(x,Bz1,'color',cor{2});hold on;
plot(x,Bz2,'color',cor{3});hold on;
plot(x,Bz3,'color',cor{4});hold on;
xlabel('x')
ylabel('B_z');
% ylim([-1e-17 1e-17])
grid off

subplot(4,2,8);
plot(x,p,'color',cor{1});hold on;
plot(x,P1,'color',cor{2});hold on;
plot(x,P2,'color',cor{3});hold on;
plot(x,P3,'color',cor{4});hold on;
xlabel('x')
ylabel('P');
ylim([0 1.1])
grid off
 