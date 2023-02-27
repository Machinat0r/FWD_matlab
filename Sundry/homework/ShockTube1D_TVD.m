function [U,F,Rho,Ux,Uy,Uz,Bx,By,Bz,P]=ShockTube1D_TVD(U,F,t,varargin)
narginchk(3,10)
try
dt = varargin{1};
dx = varargin{2};
Bx = varargin{3};
e = varargin{4};
catch
    error('no Enough Inputs')
end
global R mu0
    rho=U(1,:);
    u=U(2,:)./rho;
    v=U(3,:)./rho;
    w=U(4,:)./rho;
    By=U(5,:);
    Bz=U(6,:);
    p=(R-1)*(U(7,:)-0.5.*rho.*(u.^2+v.^2+w.^2)-(Bx.^2+By.^2+Bz.^2)/(2*mu0));
    %%
for t=dt:dt:t
    % ?¡Á??????ge(7*1001),ae(7*1001)
    % ?¡Àj=1001?¡À
    j=1001;
    rho1002=2*rho(1001)-rho(1000);     % ?????¡Â????????Roe???¨´
    u1002=2*u(1001)-u(1000);
    v1002=2*v(1001)-v(1000);
    w1002=2*w(1001)-w(1000);
    Bx1002=2*Bx(1001)-Bx(1000);
    By1002=2*By(1001)-By(1000);
    Bz1002=2*Bz(1001)-Bz(1000);
    Ro=sqrt(rho1002);
    Lo=sqrt(rho(1001));
    rho_(1001)=(0.5*(Ro+Lo))^2;
    u_(1001)=(Ro*u1002+Lo*u(1001))/(Ro+Lo);
    v_(1001)=(Ro*v1002+Lo*v(1001))/(Ro+Lo);
    w_(1001)=(Ro*w1002+Lo*w(1001))/(Ro+Lo);
    Bx_(1001)=0.25*(Ro+Lo)*(Bx1002/Ro+Bx(1001)/Lo);
    By_(1001)=0.25*(Ro+Lo)*(By1002/Ro+By(1001)/Lo);
    Bz_(1001)=0.25*(Ro+Lo)*(Bz1002/Ro+Bz(1001)/Lo);
    VV=u_(j)^2+v_(j)^2+w_(j)^2;
    UW=[1       0       0       0       0       0       0;
        u_(j)  rho_(j)  0       0       0       0       0;
        v_(j)   0    rho_(j)    0       0       0       0;
        w_(j)   0       0    rho_(j)    0       0       0;
        0       0       0       0       1       0       0;
        0       0       0       0       0       1       0;
     0.5*VV rho_(j)*u_(j) rho_(j)*v_(j)  rho_(j)*w_(j) By_(j)/mu0 Bz_(j)/mu0 1/(R-1)];
 
    Aw=[u_(j) rho_(j)   0       0       0       0       0;
        0      u_(j)    0       0 By_(j)/mu0/rho_(j) Bz_(j)/mu0/rho_(j) 1/rho_(j);
        0      0      u_(j)     0 -Bx_(j)/mu0/rho_(j) 0 0;
        0      0        0     u_(j)     0 -Bx_(j)/mu0/rho_(j) 0;
        0   By_(j)   -Bx_(j)    0     u_(j)     0       0;
        0   Bz_(j)      0    -Bx_(j)    0     u_(j)     0;
        0  R*rho_(j)    0       0       0       0     u_(j)];
 
    WU=[1      0        0       0       0       0       0;
       -u_(j)/rho_(j) 1/rho_(j) 0   0   0       0       0;
       -v_(j)/rho_(j)  0  1/rho_(j) 0   0       0       0;
       -w_(j)/rho_(j) 0  0 1/rho_(j)    0       0       0;
        0      0        0       0       1       0       0;
        0      0        0       0       0       1       0;
        0.5*(R-1)*VV -(R-1)*u_(j) -(R-1)*v_(j) -(R-1)*w_(j) -(R-1)*By_(j)/mu0 -(R-1)*Bz_(j)/mu0 (R-1)];
    Aue=UW*Aw*WU;         % Roe???¨´???????¨®
    [Re,A]=eig(Aue);      % Re???????¡Â???¨®??A?????¡Â?????¨®
    Le=Re^(-1);           % Le??¡Á¨®???¡Â???¨®
    for L=1:7
        Lambda(L,j)=A(L,L);
    end
    RR(:,:,j)=Re;
    U1002=2*U(:,j)-U(:,j-1);
    ae(:,j)=Le*(U(:,j)-U1002);
    for L=1:7
        z=abs(dt/dx*A(L,L));
        if z>e
            ge(L,j)=0.5*(z-z^2)*ae(L,j);
        else
            ge(L,j)=0.5*((z^2+e^2)/(2*e)-z^2)*ae(L,j);
        end
    end
    % j??1??1000???¨¦??
    for j=1:1000
        Ro=sqrt(rho(j+1));     % ???¡Â????????Roe???¨´
        Lo=sqrt(rho(j));
        rho_(j)=(0.5*(Ro+Lo))^2;
        u_(j)=(Ro*u(j+1)+Lo*u(j))/(Ro+Lo);
        v_(j)=(Ro*v(j+1)+Lo*v(j))/(Ro+Lo);
        w_(j)=(Ro*w(j+1)+Lo*w(j))/(Ro+Lo);
        Bx_(j)=0.25*(Ro+Lo)*(Bx(j+1)/Ro+Bx(j)/Lo);
        By_(j)=0.25*(Ro+Lo)*(By(j+1)/Ro+By(j)/Lo);
        Bz_(j)=0.25*(Ro+Lo)*(Bz(j+1)/Ro+Bz(j)/Lo);
        VV=u_(j)^2+v_(j)^2+w_(j)^2;
        UW=[1       0       0       0       0       0       0;
            u_(j)  rho_(j)  0       0       0       0       0;
            v_(j)   0    rho_(j)    0       0       0       0;
            w_(j)   0       0    rho_(j)    0       0       0;
            0       0       0       0       1       0       0;
            0       0       0       0       0       1       0;
         0.5*VV rho_(j)*u_(j) rho_(j)*v_(j)  rho_(j)*w_(j) By_(j)/mu0 Bz_(j)/mu0 1/(R-1)];
 
        Aw=[u_(j) rho_(j)   0       0       0       0       0;
            0      u_(j)    0       0 By_(j)/mu0/rho_(j) Bz_(j)/mu0/rho_(j) 1/rho_(j);
            0      0      u_(j)     0 -Bx_(j)/mu0/rho_(j) 0 0;
            0      0        0     u_(j)     0 -Bx_(j)/mu0/rho_(j) 0;
            0   By_(j)   -Bx_(j)    0     u_(j)     0       0;
            0   Bz_(j)      0    -Bx_(j)    0     u_(j)     0;
            0  R*rho_(j)    0       0       0       0     u_(j)];
 
        WU=[1      0        0       0       0       0       0;
           -u_(j)/rho_(j) 1/rho_(j) 0   0   0       0       0;
           -v_(j)/rho_(j)  0  1/rho_(j) 0   0       0       0;
           -w_(j)/rho_(j) 0  0 1/rho_(j)    0       0       0;
            0      0        0       0       1       0       0;
            0      0        0       0       0       1       0;
            0.5*(R-1)*VV -(R-1)*u_(j) -(R-1)*v_(j) -(R-1)*w_(j) -(R-1)*By_(j)/mu0 -(R-1)*Bz_(j)/mu0 (R-1)];
        Aue=UW*Aw*WU;         % Roe???¨´???????¨®
        [Re,A]=eig(Aue);      % Re???????¡Â???¨®??A?????¡Â?????¨®
        Le=Re^(-1);           % Le??¡Á¨®???¡Â???¨®
        for L=1:7
            Lambda(L,j)=A(L,L);
        end
        RR(:,:,j)=Re;
        ae(:,j)=Le*(U(:,j)-U(:,j+1));
        for L=1:7
            z=abs(dt/dx*A(L,L));
            if z>e
                ge(L,j)=0.5*(z-z^2)*ae(L,j);
            else
                ge(L,j)=0.5*((z^2+e^2)/(2*e)-z^2)*ae(L,j);
            end
        end
    end
    % ????g(7*1001)
    ge0=2*ge(:,1)-ge(:,2);
    for L=1:7
        if ge(L,1)*ge0(L)>0
            g(L,1)=sign(ge(L,1))*min(abs(ge(L,1)),abs(ge0(L)));
        else
            g(L,1)=0;
        end
    end
    for j=2:1001
        for L=1:7
            if ge(L,j)*ge(L,j-1)>0
                g(L,j)=sign(ge(L,j))*min(abs(ge(L,j)),abs(ge(L,j-1)));
            else
                g(L,j)=0;
            end
        end
    end
    % ????Me(7*1001)
    for j=1:1000
        for L=1:7
            if ae(L,j)==0
                Me(L,j)=0;
            else
                Me(L,j)=(g(L,j+1)-g(L,j))/ae(L,j);
            end
        end
    end
    for L=1:7
        if ae(L,1001)==0
            Me(L,1001)=0;
        else
            Me(L,1001)=(g(L,1001)-g(L,1000))/ae(L,1001);
        end
    end
    % ????Qk(7*1001)
    for j=1:1001
        for L=1:7
            z=abs(dt/dx*Lambda(L,j)+Me(L,j));
            if z>e
                Qk(L,j)=z;
            else
                Qk(L,j)=(z^2+e^2)/(2*e);
            end
        end
    end
    % ????phe(7*1001)
    g1002=2*g(:,1001)-g(:,1000);
    phe(:,1001)=dx/dt*(g(:,1001)+g1002-Qk(:,1001).*ae(:,1001));
    for j=1:1000
        phe(:,j)=dx/dt*(g(:,j)+g(:,j+1)-Qk(:,j).*ae(:,j));
    end
    % ????????¡À???Fe(7*1001)
    F1002=2*F(:,1001)-F(:,1000);
    Fe(:,1001)=0.5*(F1002+F(:,1001)+RR(:,:,1001)*phe(:,1001));
    for j=1:1000
        Re=RR(:,:,j);
        Fe(:,j)=0.5*(F(:,j+1)+F(:,j)+Re*phe(:,j));
    end
    % ????Up(7*1001)
    Fe0=2*Fe(:,1)-Fe(:,2);
    Up(:,1)=U(:,1)-dt/dx*(Fe0-Fe(:,1));
    for j=2:1001
        Up(:,j)=U(:,j)-dt/dx*(Fe(:,j-1)-Fe(:,j));
    end
    U=Up;
    % ???????????¡Â???¨ª??
    rho=U(1,:);
    u=U(2,:)./rho;
    v=U(3,:)./rho;
    w=U(4,:)./rho;
    By=U(5,:);
    Bz=U(6,:);
    p=(R-1)*(U(7,:)-0.5.*rho.*(u.^2+v.^2+w.^2)-(Bx.^2+By.^2+Bz.^2)/(2*mu0));
    F(1,:)=rho.*u;
    F(2,:)=rho.*u.^2+p+(-Bx.^2+By.^2+Bz.^2)/(2*mu0);
    F(3,:)=rho.*u.*v-(Bx.*By)/mu0;
    F(4,:)=rho.*u.*w-Bx.*Bz/mu0;
    F(5,:)=u.*By-v.*Bx;
    F(6,:)=u.*Bz-w.*Bx;
    F(7,:)=(U(7,:)+p+(Bx.^2+By.^2+Bz.^2)/(2*mu0)).*u-Bx/mu0.*(u.*Bx+v.*By+w.*Bz);
end
Rho = rho;
Ux = u;Uy = v;Uz = w;P = p;
end