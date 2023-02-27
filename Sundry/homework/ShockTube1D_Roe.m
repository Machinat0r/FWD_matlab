function [U,F,Rho,Ux,Uy,Uz,Bx,By,Bz,P]=ShockTube1D_Roe(U,F,t,varargin)
narginchk(3,10)
try
dt = varargin{1};
dx = varargin{2};
Bx = varargin{3};
n = varargin{4};
catch
    error('no Enough Inputs')
end
global r mu
%%
rho=U(1,:);
u=U(2,:)./rho;
v=U(3,:)./rho;
w=U(4,:)./rho;
By=U(5,:);
Bz=U(6,:);
p=(r-1)*(U(7,:)-0.5.*rho.*(u.^2+v.^2+w.^2)-(Bx.^2+By.^2+Bz.^2)/(2*mu));
et = p./(r-1)./rho + 0.5*(u.^2 + v.^2 + w.^2)  + (Bx.^2 + By.^2 + Bz.^2)/(2*mu)./rho; 
H=r*p/(r-1)./rho+0.5*(u.^2+v.^2+w.^2)+(Bx.^2+By.^2+Bz.^2)/mu./rho;

for j = 0:dt :t
    %求A及数值通量f
    for i=1:n-1
        %求roe平均
        rho_=((sqrt(rho(i))+sqrt(rho(i+1)))/2).^2;
        u_=(sqrt(rho(i))*u(i)+sqrt(rho(i+1))*u(i+1))/(sqrt(rho(i))+sqrt(rho(i+1)));
        v_=(sqrt(rho(i))*v(i)+sqrt(rho(i+1))*v(i+1))/(sqrt(rho(i))+sqrt(rho(i+1)));
        w_=(sqrt(rho(i))*w(i)+sqrt(rho(i+1))*w(i+1))/(sqrt(rho(i))+sqrt(rho(i+1)));
        Bx_=0.25*(sqrt(rho(i))+sqrt(rho(i+1)))*(Bx(i+1)/sqrt(rho(i+1))+Bx(i)/sqrt(rho(i)));
        By_=0.25*(sqrt(rho(i))+sqrt(rho(i+1)))*(By(i+1)/sqrt(rho(i+1))+By(i)/sqrt(rho(i)));
        Bz_=0.25*(sqrt(rho(i))+sqrt(rho(i+1)))*(Bz(i+1)/sqrt(rho(i+1))+Bz(i)/sqrt(rho(i)));
        H_=(sqrt(rho(i))*H(i)+sqrt(rho(i+1))*H(i+1))/(sqrt(rho(i))+sqrt(rho(i+1)));
        
        %计算A,A_
         X=(Bx(i)-Bx(i+1))^2/2/(sqrt(rho(i))+sqrt(rho(i+1)))^2+(By(i)-By(i+1))^2/2/(sqrt(rho(i))+sqrt(rho(i+1)))^2+(Bz(i)-Bz(i+1))^2/2/(sqrt(rho(i))+sqrt(rho(i+1)))^2;
        delta21=-u_^2+(2-r)*X+(r-1)/2*(u_^2+v_^2+w_^2);
        delta22=2*u_-(r-1)*u_;
        delta23=-(r-1)*v_;
        delta24=-(r-1)*w_;
        delta25=(2-r)*By_/mu;
        delta26=(2-r)*Bz_/mu;
        delta27=r-1;
        delta71=-u_*H_+u_*(delta21+u_^2)+Bx_*(u_*Bx_+v_*By_+w_*Bz_)/rho_/mu;
        delta72=H_+u_*(delta22-2*u_);
        delta73=u_*delta23-Bx_*By_/rho_/mu;
        delta74=u_*delta24-Bx_*Bz_/rho_/mu;
        delta75=u_*delta25-Bx_*v_/mu;
        delta76=u_*delta26-Bx_*w_/mu;
        delta77=u_+u_*delta27;
        
        A=[0 1 0 0 0 0 0;
            delta21 delta22 delta23 delta24 delta25 delta26 delta27;
            -u_*v_ v_ u_ 0 -Bx_/mu 0 0;
            -u_*w_ w_ 0 u_ 0 -Bx_/mu 0;
            -By_*u_/rho_+Bx_*v_/rho_ By_/rho_ -Bx_/rho_ 0 u_ 0 0;
            -Bz_*u_/rho_+Bx_*w_/rho_ Bz_/rho_ 0 -Bx_/rho_ 0 u_ 0;
            delta71 delta72 delta73 delta74 delta75 delta76 delta77];
        [R,D]=eig(A);
        A_ = R*abs(D) * R^(-1);
        
        %计算数值通量
        f(:, i) = 0.5 * (F(:, i) + F(:, i+1)) - 0.5 * A_ * (U(:,i+1)-U(:,i));
    end
    
    %推进到下一个时间步的U
    for i = 2:n-1
        U(:, i) = U(:, i) - dt/dx * (f(:, i) - f(:, i-1));
    end
    U(:,1)=U(:,2);
    U(:,n)=U(:,n-1);
    
    %获得原始变量
    rho=U(1,:);
    u=U(2,:)./rho;
    v=U(3,:)./rho;
    w=U(4,:)./rho;
    By=U(5,:);
    Bz=U(6,:);
    p=(r-1)*(U(7,:)-0.5.*rho.*(u.^2+v.^2+w.^2)-(Bx.^2+By.^2+Bz.^2)/(2*mu));
    F(1, :) = rho.*u;
    F(2, :) = rho.*u.*u + p + (-Bx.^2 + By.^2 + Bz.^2)/(2*mu);
    F(3, :) = rho.*u.*v - Bx.*By/mu;
    F(4, :) = rho.*u.*w - Bx.*Bz/mu;
    F(5, :) = u.*By - v.*Bx;
    F(6, :) = u.*Bz - w.*Bx;
    F(7, :) = (U(7, :) + p + (Bx.^2 + By.^2 + Bz.^2)/(2*mu)).*u - Bx.*(u.*Bx + v.*By + w.*Bz)/mu;
    H=r*p/(r-1)./rho+0.5*(u.^2+v.^2+w.^2)+(Bx.^2+By.^2+Bz.^2)/mu./rho;
end
Rho = rho;
Ux = u;Uy = v;Uz = w;
Bx = Bx;By = By;Bz = Bz;
P = p;
end