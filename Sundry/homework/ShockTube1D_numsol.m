function [U,F,Rho,Ux,Uy,Uz,Bx,By,Bz,P] = ShockTube1D_numsol(U,F,Times,varargin)
% This function calculate U & F with the specified iteration method for
% specified times.
% The default iteration method is "Lax"
%------written by Wending Fu, Dec.24.2022 in Beijing------------
%% Choose method
narginchk(3,10)
try
Method = varargin{1};
dt = varargin{2};
dx = varargin{3};
catch
    error('no Enough Inputs')
end
global gamma mu0
%% Lax-Friedrichs Method
switch Method
case 'Lax-Friedrichs'
% iteration
for t=1:Times 
    % transport boundary condition
    Uexpand = [U(:,1), U, U(:,end)];
    Fexpand = [F(:,1), F, F(:,end)];

    tempU = 0.5*(Uexpand(:,3:end)+Uexpand(:,1:end-2)) - 0.5*dt*(Fexpand(:,3:end)-Fexpand(:,1:end-2))/dx;
    U = tempU;
    
    % reverse solution
    [tempRho,tempUx,tempUy,tempUz,tempBx,tempBy,tempBz,tempP] = ReverseSolution(U);
    
    % solve F
    F = SolveF(tempRho,tempUx,tempUy,tempUz,tempBx,tempBy,tempBz,tempP,U);
end
%% MacCormack Method
case 'MacCormack'
alpha = varargin{4};
for t=1:Times 
    % transport boundary condition
    Uexpand = [U(:,1), U, U(:,end)];
    Fexpand = [F, F(:,end)];
    tempU = Uexpand(:,2:end-1)-dt/dx*(Fexpand(:,2:end)-Fexpand(:,1:end-1))+...
        alpha*dt/(dx^2)*(Uexpand(:,3:end)+Uexpand(:,1:end-2)-2*Uexpand(:,2:end-1));
    [tempRho,tempUx,tempUy,tempUz,tempBx,tempBy,tempBz,tempP] = ReverseSolution(tempU);
    tempF = SolveF(tempRho,tempUx,tempUy,tempUz,tempBx,tempBy,tempBz,tempP,tempU);
    
    Uexpand2 = [tempU(:,1), tempU, tempU(:,end)];
    Fexpand2 = [tempF(:,1), tempF];
    U = 0.5*(Uexpand(:,2:end-1)+Uexpand2(:,2:end-1))-0.5*dt/dx*(Fexpand2(:,2:end)-Fexpand2(:,1:end-1))...
        +alpha*dt/(dx^2)*(Uexpand2(:,3:end)+Uexpand2(:,1:end-2)-2*Uexpand2(:,2:end-1));
    [tempRho,tempUx,tempUy,tempUz,tempBx,tempBy,tempBz,tempP] = ReverseSolution(U);
    F = SolveF(tempRho,tempUx,tempUy,tempUz,tempBx,tempBy,tempBz,tempP,U);
end
%% Roe Method
case 'Roe'
for t=1:Times     
    [Rho,Ux,Uy,Uz,Bx,By,Bz,P] = ReverseSolution(U);
    F = SolveF(Rho,Ux,Uy,Uz,Bx,By,Bz,P,U);
    Uexpand = [U(:,:), U(:,end)];
    Fexpand = [F(:,:), F(:,end)];
    
    H=gamma/(gamma-1)*P./Rho+0.5*(Ux.^2+Uy.^2+Uz.^2)+(Bx.^2+By.^2+Bz.^2)/mu0./Rho;
    RhoEx = [Rho, Rho(end)];
    UxEx = [Ux, Ux(end)];
    UyEx = [Uy, Uy(end)];
    UzEx = [Uz, Uz(end)];
    BxEx = [Bx, Bx(end)];
    ByEx = [By, By(end)];
    BzEx = [Bz, Bz(end)];
    HEx = [H, H(end)];
    
    for i = 1:size(Rho,2)
    RhoR = sqrt(RhoEx(i+1)); RhoL = sqrt(RhoEx(i));
    Rhom = (0.5*(RhoR+RhoL)).^2;
    Uxm = (RhoR.*UxEx(i+1)+RhoL.*UxEx(i))./(RhoR+RhoL);
    Uym = (RhoR.*UyEx(i+1)+RhoL.*UyEx(i))./(RhoR+RhoL);
    Uzm = (RhoR.*UzEx(i+1)+RhoL.*UzEx(i))./(RhoR+RhoL);
    Bxm = 0.25*(RhoR+RhoL).*(BxEx(i+1)./RhoR+BxEx(i)./RhoL);
    Bym = 0.25*(RhoR+RhoL).*(ByEx(i+1)./RhoR+ByEx(i)./RhoL);
    Bzm = 0.25*(RhoR+RhoL).*(BzEx(i+1)./RhoR+BzEx(i)./RhoL);
    Hm = (RhoL.*HEx(i)+RhoR.*HEx(i+1))./(RhoR+RhoL);

    X = (BxEx(i)-BxEx(i+1))^2/2/(sqrt(RhoEx(i))+sqrt(RhoEx(i+1)))^2+(ByEx(i)-ByEx(i+1))^2/2/(sqrt(RhoEx(i))+sqrt(RhoEx(i+1)))^2+(BzEx(i)-BzEx(i+1))^2/2/(sqrt(RhoEx(i))+sqrt(RhoEx(i+1)))^2;
    a21 = -Uxm^2+(2-gamma)*X+(gamma-1)/2*(Uxm^2+Uym^2+Uzm^2);
    a22 = 2*Uxm-(gamma-1)*Uxm;
    a23 = -(gamma-1)*Uym;
    a24 = -(gamma-1)*Uzm;
    a25 = (2-gamma)*Bxm/mu0;
    a26 = (2-gamma)*Bym/mu0;
    a27 = (2-gamma)*Bzm/mu0;
    a28 = gamma-1;
    a81 = -Uxm*Hm+Uxm*(a21+Uxm^2)+Bxm*(Uxm*Bxm+Uym*Bym+Uzm*Bzm)/Rhom/mu0;
    a82 = Hm+Uxm*(a22-2*Uxm);
    a83 = Uxm*a23-Bxm*Bym/Rhom/mu0;
    a84 = Uxm*a24-Bxm*Bzm/Rhom/mu0;
    a85 = (-Uxm*Bxm-Uym*Bym-Uzm*Bzm)/mu0;
    a86 = Uxm*a26-Bxm*Uym/mu0;
    a87 = Uxm*a27-Bxm*Uzm/mu0;
    a88 = Uxm+Uxm*a28;

    D=[0 1 0 0 0 0 0 0;
       a21 a22 a23 a24 a25 a26 a27 a28;
       -Uxm*Uym Uym Uxm 0 -Bym/mu0 -Bxm/mu0 0 0;
       -Uxm*Uzm Uzm 0 Uxm -Bzm/mu0 0 -Bxm/mu0 0;
       0 0 0 0 0 0 0 0;
       -Bym*Uxm/Rhom+Bxm*Uym/Rhom Bym/Rhom -Bxm/Rhom 0 Uym Uxm 0 0;
       -Bzm*Uxm/Rhom+Bxm*Uzm/Rhom Bzm/Rhom 0 -Bxm/Rhom Uzm 0 Uxm 0;
       a81 a82 a83 a84 a85 a86 a87 a88];
    
    [Re,D]=eig(D);
    Am = Re*abs(D) * Re^(-1);
    c_eval('Ftemp(:, i) = 0.5*(Fexpand(:,i)+Fexpand(:,i+1)) - 0.5*Am*(Uexpand(:,i+1)-Uexpand(:,i));',1)
    end
    
    Fexpand = [Ftemp Ftemp(:,end)];
    U = Uexpand(:,1:end-1)-dt/dx*(Fexpand(:,2:end)-Fexpand(:,1:end-1));
end        
%% TVD Method
case 'TVD'
entro = varargin{4};
for t = 1:Times
    [Rho,Ux,Uy,Uz,Bx,By,Bz,P] = ReverseSolution(U);
    F = SolveF(Rho,Ux,Uy,Uz,Bx,By,Bz,P,U);
    Uexpand = [U, U(:,end)];
    Fexpand = [F, F(:,end)];
    
    H=gamma/(gamma-1)*P./Rho+0.5*(Ux.^2+Uy.^2+Uz.^2)+(Bx.^2+By.^2+Bz.^2)/mu0./Rho;
    RhoEx = [Rho, Rho(end)];
    UxEx = [Ux, Ux(end)];
    UyEx = [Uy, Uy(end)];
    UzEx = [Uz, Uz(end)];
    BxEx = [Bx, Bx(end)];
    ByEx = [By, By(end)];
    BzEx = [Bz, Bz(end)];
    HEx = [H, H(end)];
    
    RhoR = sqrt(RhoEx(i+1)); RhoL = sqrt(RhoEx(i));
    Rhom = (0.5*(RhoR+RhoL)).^2;
    Uxm = (RhoR.*UxEx(i+1)+RhoL.*UxEx(i))./(RhoR+RhoL);
    Uym = (RhoR.*UyEx(i+1)+RhoL.*UyEx(i))./(RhoR+RhoL);
    Uzm = (RhoR.*UzEx(i+1)+RhoL.*UzEx(i))./(RhoR+RhoL);
    Bxm = 0.25*(RhoR+RhoL).*(BxEx(i+1)./RhoR+BxEx(i)./RhoL);
    Bym = 0.25*(RhoR+RhoL).*(ByEx(i+1)./RhoR+ByEx(i)./RhoL);
    Bzm = 0.25*(RhoR+RhoL).*(BzEx(i+1)./RhoR+BzEx(i)./RhoL);
    Hm = (RhoL.*HEx(i)+RhoR.*HEx(i+1))./(RhoR+RhoL);
    
    X = (BxEx(i)-BxEx(i+1))^2/2/(sqrt(RhoEx(i))+sqrt(RhoEx(i+1)))^2+(ByEx(i)-ByEx(i+1))^2/2/(sqrt(RhoEx(i))+sqrt(RhoEx(i+1)))^2+(BzEx(i)-BzEx(i+1))^2/2/(sqrt(RhoEx(i))+sqrt(RhoEx(i+1)))^2;
    a21 = -Uxm^2+(2-gamma)*X+(gamma-1)/2*(Uxm^2+Uym^2+Uzm^2);
    a22 = 2*Uxm-(gamma-1)*Uxm;
    a23 = -(gamma-1)*Uym;
    a24 = -(gamma-1)*Uzm;
    a25 = (2-gamma)*Bxm/mu0;
    a26 = (2-gamma)*Bym/mu0;
    a27 = (2-gamma)*Bzm/mu0;
    a28 = gamma-1;
    a81 = -Uxm*Hm+Uxm*(a21+Uxm^2)+Bxm*(Uxm*Bxm+Uym*Bym+Uzm*Bzm)/Rhom/mu0;
    a82 = Hm+Uxm*(a22-2*Uxm);
    a83 = Uxm*a23-Bxm*Bym/Rhom/mu0;
    a84 = Uxm*a24-Bxm*Bzm/Rhom/mu0;
    a85 = (-Uxm*Bxm-Uym*Bym-Uzm*Bzm)/mu0;
    a86 = Uxm*a26-Bxm*Uym/mu0;
    a87 = Uxm*a27-Bxm*Uzm/mu0;
    a88 = Uxm+Uxm*a28;

    A=[0 1 0 0 0 0 0 0;
       a21 a22 a23 a24 a25 a26 a27 a28;
       -Uxm*Uym Uym Uxm 0 -Bym/mu0 -Bxm/mu0 0 0;
       -Uxm*Uzm Uzm 0 Uxm -Bzm/mu0 0 -Bxm/mu0 0;
       0 0 0 0 0 0 0 0;
       -Bym*Uxm/Rhom+Bxm*Uym/Rhom Bym/Rhom -Bxm/Rhom 0 Uym Uxm 0 0;
       -Bzm*Uxm/Rhom+Bxm*Uzm/Rhom Bzm/Rhom 0 -Bxm/Rhom Uzm 0 Uxm 0;
       a81 a82 a83 a84 a85 a86 a87 a88];

    [Re,D]=eig(A);
    Le=Re^(-1);
    
    c_eval('L(?,i) = D(?,?);',1:8);
    R{i}=Re;
    a(:,i)=Le*(Uexpand(:,i)-Uexpand(:,i+1));
    
    for ii=1:8
        z=abs(dt/dx*D(ii,ii));
        if z>entro
            ge(ii,i)=0.5*(z-z^2)*a(ii,i);
        else
            ge(ii,i)=0.5*((z^2+entro^2)/(2*entro)-z^2)*a(ii,i);
        end
    end


    ge0=2*ge(:,1)-ge(:,2);
    g = zeros(8,size(Rho,2)+1);
    c_eval('g(ge(?,1)*ge0(?)>0,1) = sign(ge(?,1)).*min(abs(ge(?,1)),abs(ge0(?)));',1:8);
    for ii=2:size(Rho,2)+1
        c_eval('g(ge(?,ii)*ge(?,ii)>0,ii) = sign(ge(?,ii)).*min(abs(ge(?,ii)),abs(ge(?,ii-1)));',1:8); 
    end
    
    g = [g, g(:,end)];
    M = zeros(8,size(Rho,2)+1);
    for ii=1:size(Rho,2)+1
        c_eval('M(?,a(?,ii)~=0) = (g(?,ii+1)-g(?,ii))/a(?,ii);',1:8)
        c_eval('z(?) = abs(dt/dx*L(?,ii)+M(?,ii));',1:8)
        for j = 1:8
            if z>entro
                Qk(ii,j)=z;
            else
                Qk(ii,j)=(z^2+entro^2)/(2*entro);
            end
        end
    end

    Phi = dx/dt*(g(:,1:end-1)+g(:,2:end)-Qk(:,1:end).*a(:,1:end));
    c_eval('Fe(:,?) = 0.5*(Fexpand(:,?+1)+Fexpand(:,?)+R{?}*Phi(:,?));',1:size(Rho,2));
    Feexpand = [Fe(:,1),Fe];
    U = U - dt/dx*(Feexpand(:,1:end-1)-Feexpand(:,2:end));
end
end
[Rho,Ux,Uy,Uz,Bx,By,Bz,P] = ReverseSolution(U);
end

function [Rho,Ux,Uy,Uz,Bx,By,Bz,P] = ReverseSolution(U)
global gamma mu0
Rho = U(1,:);
Ux = U(2,:)./Rho;
Uy = U(3,:)./Rho;
Uz = U(4,:)./Rho;
Bx = U(5,:);
By = U(6,:);
Bz = U(7,:);
P = (gamma-1)*(U(8,:)-0.5.*Rho.*(Ux.^2+Uy.^2+Uz.^2)-(Bx.^2+By.^2+Bz.^2)/(2*mu0));
end

function F = SolveF(Rho,Ux,Uy,Uz,Bx,By,Bz,P,U)
global mu0
F(1,:) = Rho.*Ux;
F(2,:) = Rho.*Ux.^2+P+(-Bx.^2+By.^2+Bz.^2)/(2*mu0);
F(3,:) = Rho.*Ux.*Uy-(Bx.*By)/mu0;
F(4,:) = Rho.*Ux.*Uz-Bx.*Bz/mu0;
F(5,:) = 0;
F(6,:) = Ux.*By-Uy.*Bx;
F(7,:) = Ux.*Bz-Uz.*Bx;
F(8,:) = (U(8,:)+P+(Bx.^2+By.^2+Bz.^2)/(2*mu0)).*Ux-Bx/mu0.*(Ux.*Bx+Uy.*By+Uz.*Bz);
end