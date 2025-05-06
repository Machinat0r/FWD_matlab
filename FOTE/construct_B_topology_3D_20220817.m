clc
clear 
Tcomputsta=clock; 
plottype = 'cline';
%--------------------------------------------
tintStr = '19-Sep-2015';
Tsta='2015-09-19T07:43:30.45Z'; 
Tend='2015-09-19T07:43:31.65Z'; 
tint=irf.tint('2015-09-19T07:43:30.12Z/2015-09-19T07:43:31.62Z');
Time='2015-09-19T07:43:31.115Z';

% Tsta='2018-08-27T12:15:30.000Z'; 
% Tend='2018-08-27T12:15:50.000Z'; 

mms.db_init('local_file_db','E:\MMS\')
                        
% tint=irf.tint('2018-08-27T12:15:30.000Z/2018-08-27T12:15:50.000Z');

% Tnull='2018-08-27T12:15:43.045Z';
% Tnull1='2018-08-27T12:15:43.045Z';
Tnull = '2015-09-19T07:43:31.115Z';
Tnull1 = '2015-09-19T07:43:31.115Z';
% Time1=Tnull1;
% time1 = irf_time(Time1,'utc>epochtt');
% Tnull2='2018-08-27T12:15:43.045Z';
Tnull2 = '2015-09-19T07:43:31.115Z';
Time1=Tnull1;
time1 = irf_time(Time1,'utc>epochtt');

ic=1:4;
    c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',tint,?);',ic);
    c_eval('B?=irf.ts2mat(Bxyz?);',ic);
    c_eval('Bt?=Bxyz?.abs;',ic);
%     c_eval('Rgse?=mms.get_data(''R_gse'',tint,?);',ic);
%     c_eval('R?=irf.ts2mat(Rgse?);',ic);
Pos = mms.get_data('R_gse',tint);
R_time = Pos.time.epoch;
c_eval('R? = Pos.gseR?;')
c_eval('R? = [Pos.time.epochUnix R?(:,1:3)];')


B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
tint=[iso2epoch(Tsta) iso2epoch(Tend)];

lth = 4;

for ic=1:4
    c_eval('R?=irf_resamp(R?,B?);',ic);
    c_eval('B?=irf_tlim(B?,tint);',ic);
    c_eval('R?=irf_tlim(R?,tint);',ic);
    c_eval('Ro?=irf_tlim(R?,tint);',ic);
    c_eval('Rts_sc?=ones(lth,3);',ic);
end

gradB=c_4_grad('R?','B?','grad');

idxnull_1=find(gradB(:,1,1)>=iso2epoch(Tnull1)); idxnull_1=idxnull_1(1);
idxnull_2=find(gradB(:,1,1)>=iso2epoch(Tnull2)); idxnull_2=idxnull_2(1);
% idxnull = 330;

%% axis_calculation_a
dB_null_1=reshape(gradB(idxnull_1,2:end),3,3);
[V,D] = eig(dB_null_1);
if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
    [m,~]=find(D==min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))])));
    [mm,nn]=find(D==-(min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]))));
    if isempty(m)
         m=mm;
    end
    if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
        [m_spine,n_spine]=find(D==-([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));
        [m_fan,n_fan]=find(D==([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));   
    else
        [m_spine,n_spine]=find(D==([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));
        [m_fan,n_fan]=find(D==-([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));
    end
else
    [m,~]=find(imag(D)==0&D~=0);  m_spine =m;
end  

    
e1=V(:,m)';
e1=irf_norm(e1);
% etmp=[0 1 0];
e_fan=cross(V(:,m_fan(1))',V(:,m_fan(2))');  
% e3=cross(e1,etmp);   
e3=irf_norm(e_fan);
e2=cross(e3,e1);  


for ic=1:4
    c_eval('B?_a=irf_newxyz(B?, e1,e2,e3);',ic);
    c_eval('R?_a=irf_newxyz(R?, e1,e2,e3);',ic);
end
V1 = V';
V = irf_newxyz(V1, e1,e2,e3);
% V([m,3],:) = V([3,m],:);
a_spine = V(m_spine,:);
a3 = V(m,:);
as=[0,a_spine(2),a_spine(3)];
e3=as/norm(as);
b3 = cross(e3,a3);  

e3_spine = irf_norm(V(m_spine,:));
e1_spine = cross(e3_spine,a3);
e2_spine = cross(e3_spine,e1_spine);
M =[e1_spine;e2_spine;e3_spine];

%change the coordinates to eigen vector corrdinate

%--------------------------------------------
% for ic=1:4
% c_eval('B?_a=irf_newxyz(B?_a, a3,b3,e3);',ic);
% c_eval('R?_a=irf_newxyz(R?_a, a3,b3,e3);',ic);
% end

gradB_a=c_4_grad('R?_a','B?_a','grad');

dB_null_a=reshape(gradB_a(idxnull_1,2:end),3,3);

dR1_a=dB_null_a\B1_a(idxnull_1,2:4)';   
R_null_a=R1_a(idxnull_1,2:4)-dR1_a';

Rsc1_a=R1_a(idxnull_1,2:end);
Rsc2_a=R2_a(idxnull_1,2:end);
Rsc3_a=R3_a(idxnull_1,2:end);
Rsc4_a=R4_a(idxnull_1,2:end);

Rsc1_a=Rsc1_a-R_null_a;
Rsc2_a=Rsc2_a-R_null_a;
Rsc3_a=Rsc3_a-R_null_a;
Rsc4_a=Rsc4_a-R_null_a;

%% axis_calculation_b
dB_null_2=reshape(gradB(idxnull_2,2:end),3,3);
[V,D] = eig(dB_null_2);
if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
    [m,~]=find(D==min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))])));
    [mm,nn]=find(D==-(min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]))));
    if isempty(m)
         m=mm;
    end
    if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
        [m_spine,n_spine]=find(D==-([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));
        [m_fan,n_fan]=find(D==([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));   
    else
        [m_spine,n_spine]=find(D==([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));
        [m_fan,n_fan]=find(D==-([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]));
    end
else
    [m,~]=find(imag(D)==0&D~=0);  m_spine =m;
end  

    
e1=V(:,m)';
e1=irf_norm(e1);
% etmp=[0 1 0];
e_fan=cross(V(:,m_fan(1))',V(:,m_fan(2))');  
% e3=cross(e1,etmp);   
e3=irf_norm(e_fan);
e2=cross(e3,e1);     
for ic=1:4
    c_eval('B?_b=irf_newxyz(B?, e1,e2,e3);',ic);
    c_eval('R?_b=irf_newxyz(R?, e1,e2,e3);',ic);
end
V1 = V';
V = irf_newxyz(V1, e1,e2,e3);
V([m,3],:) = V([3,m],:);
a_spine = V(m_spine,:);
a3 = V(m,:);
as=[0,a_spine(2),a_spine(3)];
e3=as/norm(as);
b3 = cross(e3,a3);      

% for ic=1:4
% c_eval(['B?_b=irf_newxyz(B?_b, a3,b3,e3);'],ic);
% c_eval(['R?_b=irf_newxyz(R?_b, a3,b3,e3);'],ic);
% end

gradB_b=c_4_grad('R?_b','B?_b','grad');

dB_null_b=reshape(gradB_b(idxnull_2,2:end),3,3);

dR1_b=dB_null_b \ B1_b(idxnull_2,2:4)';   
R_null_b=R1_b(idxnull_2,2:4)-dR1_b';

Rsc1_b=R1_b(idxnull_2,2:end);
Rsc2_b=R2_b(idxnull_2,2:end);
Rsc3_b=R3_b(idxnull_2,2:end);
Rsc4_b=R4_b(idxnull_2,2:end);


Rsc1_b=Rsc1_b-R_null_b;
Rsc2_b=Rsc2_b-R_null_b;
Rsc3_b=Rsc3_b-R_null_b;
Rsc4_b=Rsc4_b-R_null_b;

%% construct B around null
set(0,'defaultLineLineWidth', 0.5);
set(0,'defaultAxesFontSize', 12);
set(0,'defaultTextFontSize', 12);
set(0,'defaultAxesFontUnits', 'pixels');
fig1=figure( ...
'Name','Dataset coverage', ...
'Tag','XYXgsm');clf;
set(fig1,'PaperUnits','centimeters')
% xSize = 20; ySize = 20; coef=floor(min(800/xSize,800/ySize));
% xLeft = 0; yTop = -1;
% set(fig1,'PaperPosition',[xLeft yTop xSize ySize])
% set(fig1,'Position',[10 10 xSize*coef ySize*coef])

xSize = 800; ySize=800;
set(gcf,'position',[10,300,xSize,ySize])
set(gcf,'render','painters');
Paper_X = 20; Paper_Y = 20; 
coef=floor(max(xSize/Paper_X,ySize/Paper_Y));
FigSize_X = xSize/coef; FigSize_Y = ySize/coef;
xLeft2 = (Paper_X- FigSize_X)/2;  yTop2 = (Paper_Y- FigSize_Y)/2; 
set(gcf,'PaperSize', [Paper_X Paper_Y]); 
set(gcf,'PaperPosition',[xLeft2 yTop2 FigSize_X FigSize_Y])

h=[];

h(1)=axes('position',[0.1 0.1 0.8 0.8]); % [x y dx dy]
% h(1)=axes('position',[0.23 0.23 0.55 0.55]);
aaa=1;

BoxWid=150; 
Flim=100;   Flimx=100;   Flimy=100;
Flim_a=100; Flim_ax=100; Flim_ay=100;
Flim_b=50;  Flim_bx=50;  Flim_by=50;



cmin=0; cmax=75;
% BoxWid=200;
Xline_preva = []; Xline_nexta = [];
Yline_preva = []; Yline_nexta = [];
Zline_preva = []; Zline_nexta = [];
Bmline_preva = []; Bmline_nexta = [];

Xline_prevb = []; Xline_nextb = [];
Yline_prevb = []; Yline_nextb = [];
Zline_prevb = []; Zline_nextb = [];
Bmline_prevb = []; Bmline_nextb = [];
ns=1;

%% calculate spine region
% for Zgrid=[-25 -15 -1 1 15 25]  %31.32  31.3 55
for Zgrid=[-45 -15 15 45]
for theta=[0:5:355]*pi/180   %31.3  
    RRR = [13];
Xgrid=RRR*cos(theta); 
Ygrid=RRR*sin(theta);

XYZnew = [Xgrid,Ygrid,Zgrid]*M;
% plot3(gca, [XYZnew(1)], [XYZnew(2)],[XYZnew(3)], 'rs', 'Linewidth',1, ...
% 'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12); hold on;

%   for RRR = [0,5,10,20,30,40,50]
%       for theta = [0:30:360]*pi/180
%           Xgrid=RRR*cos(theta);
%           Ygrid=RRR*sin(theta); 
%           Zgrid=0;
%-----inverse trace-----

ii=1;
% Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;
Xprev=XYZnew(1);  Yprev=XYZnew(2);  Zprev=XYZnew(3);
while abs(Xprev)<BoxWid && abs(Yprev)<BoxWid && abs(Zprev)<BoxWid 
Xcurt=Xprev; 
Ycurt=Yprev;
Zcurt=Zprev;            
Br = dB_null_a * ([Xcurt Ycurt Zcurt])';            
%             Br = Br+([Bxbc Bybc Bzbc])';   %加上重心处磁场
Bxcurt=Br(1);
Bycurt=Br(2);
Bzcurt=Br(3); 
Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2); 

%             step=Bmcurt;
step=1;
stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;          
Xprev=Xcurt-stepvec(1); 
Yprev=Ycurt-stepvec(2);
Zprev=Zcurt-stepvec(3); 

Xline(ii)=Xcurt;  
Yline(ii)=Ycurt; 
Zline(ii)=Zcurt;
Bmline(ii)=Bmcurt; 
ii=ii+1;
end
    
Xline_preva(ns,1:length(Xline)) = Xline; 
Yline_preva(ns,1:length(Xline)) = Yline; 
Zline_preva(ns,1:length(Xline)) = Zline; 
Bmline_preva(ns,1:length(Xline)) = Bmline; 

clear Xline 
clear Yline 
clear Zline
clear Bmline 
clear Nlin

ii=1;
% Xprev2=Xgrid;  Yprev2=Ygrid;  Zprev2=Zgrid;
Xprev2=XYZnew(1);  Yprev2=XYZnew(2);  Zprev2=XYZnew(3);
while abs(Xprev2)<BoxWid && abs(Yprev2)<BoxWid && abs(Zprev2)<BoxWid 
    Xcurt2=Xprev2;
    Ycurt2=Yprev2;
    Zcurt2=Zprev2;


    Br2 = dB_null_b * ([Xcurt2 Ycurt2 Zcurt2])';
    %             Br = Br+([Bxbc Bybc Bzbc])';   %加上重心处磁场
    Bxcurt2=Br2(1);
    Bycurt2=Br2(2);
    Bzcurt2=Br2(3);         
    Bmcurt2=sqrt(Br2(1).^2+Br2(2).^2+Br2(3).^2);
    %             step=Bmcurt;
    step=1;
    stepvec2=[Bxcurt2 Bycurt2 Bzcurt2]/norm([Bxcurt2 Bycurt2 Bzcurt2])*step;
    Xprev2=Xcurt2-stepvec2(1);
    Yprev2=Ycurt2-stepvec2(2);
    Zprev2=Zcurt2-stepvec2(3);

    Xline2(ii)=Xcurt2; 
    Yline2(ii)=Ycurt2;
    Zline2(ii)=Zcurt2;
    Bmline2(ii)=Bmcurt2;
    ii=ii+1;
end

Xline_prevb(ns,1:length(Xline2)) = Xline2;
Yline_prevb(ns,1:length(Xline2)) = Yline2;
Zline_prevb(ns,1:length(Xline2)) = Zline2;
Bmline_prevb(ns,1:length(Xline2)) = Bmline2;

clear Xline2 
clear Yline2
clear Zline2
clear Bmline2 


%-----positive trace----one side
ii=1;
% Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
% Xnext2=Xgrid;  Ynext2=Ygrid;  Znext2=Zgrid;
Xnext=XYZnew(1);  Ynext=XYZnew(2);  Znext=XYZnew(3);

while abs(Xnext)<BoxWid && abs(Ynext)<BoxWid && abs(Znext)<BoxWid 
    Xcurt=Xnext;
    Ycurt=Ynext;
    Zcurt=Znext;

    Br = dB_null_a * ([Xcurt Ycurt Zcurt])';
    %             Br = Br+([Bxbc Bybc Bzbc])';  %加上重心处磁场
    Bxcurt=Br(1);
    Bycurt=Br(2);
    Bzcurt=Br(3);
    Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
    %             step=Bmcurt;
    step=1;
    stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
    Xnext=Xcurt+stepvec(1);
    Ynext=Ycurt+stepvec(2);
    Znext=Zcurt+stepvec(3);

    Xline(ii)=Xcurt;
    Yline(ii)=Ycurt;
    Zline(ii)=Zcurt;
    Bmline(ii)=Bmcurt;
    ii=ii+1;
end

Xline_nexta(ns,1:length(Xline)) = Xline;
Yline_nexta(ns,1:length(Xline)) = Yline;
Zline_nexta(ns,1:length(Xline)) = Zline;
Bmline_nexta(ns,1:length(Xline)) = Bmline;

clear Xline
clear Yline
clear Zline
clear Bmline
clear Nlin

%-----positive trace----the other side
ii=1;
% Xnext2=Xgrid;  Ynext2=Ygrid;  Znext2=Zgrid;
Xnext2=XYZnew(1);  Ynext2=XYZnew(2);  Znext2=XYZnew(3);
while abs(Xnext2)<BoxWid && abs(Ynext2)<BoxWid && abs(Znext2)<BoxWid 
    Xcurt2=Xnext2;
    Ycurt2=Ynext2;
    Zcurt2=Znext2;

    Br2 = dB_null_b * ([Xcurt2 Ycurt2 Zcurt2])';
    %             Br = Br+([Bxbc Bybc Bzbc])';  %加上重心处磁场
    Bxcurt2=Br2(1);
    Bycurt2=Br2(2);
    Bzcurt2=Br2(3);
    Bmcurt2=sqrt(Br2(1).^2+Br2(2).^2+Br2(3).^2);
    %             step=Bmcurt;
    step=1;
    stepvec2=[Bxcurt2 Bycurt2 Bzcurt2]/norm([Bxcurt2 Bycurt2 Bzcurt2])*step;
    Xnext2=Xcurt2+stepvec2(1);
    Ynext2=Ycurt2+stepvec2(2);
    Znext2=Zcurt2+stepvec2(3);

    Xline2(ii)=Xcurt2;
    Yline2(ii)=Ycurt2;
    Zline2(ii)=Zcurt2;
    Bmline2(ii)=Bmcurt2;
    ii=ii+1;
end

Xline_nextb(ns,1:length(Xline2)) = Xline2;
Yline_nextb(ns,1:length(Xline2)) = Yline2;
Zline_nextb(ns,1:length(Xline2)) = Zline2;
Bmline_nextb(ns,1:length(Xline2)) = Bmline2;
ns = ns +1;
clear Xline2
clear Yline2
clear Zline2
clear Bmline2
clear Nlin
end
end
%% calculate fan region
for RRR = [55]
    for theta = [-80:2.5:50,100:2.5:230]*pi/180
%       for theta = [0:2.5:362.5]*pi/180
        for Zgrid=[-0.3,0 0.3]
            Xgrid=RRR*cos(theta);
            Ygrid=RRR*sin(theta); 
%             Zgrid=RRR*sin(phi);
              
            %-----inverse trace-----

            ii=1;
            Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;
                while abs(Xprev)<BoxWid && abs(Yprev)<BoxWid && abs(Zprev)<BoxWid 
                Xcurt=Xprev; 
                Ycurt=Yprev;
                Zcurt=Zprev;            
                Br = dB_null_a * ([Xcurt Ycurt Zcurt])';            
                %             Br = Br+([Bxbc Bybc Bzbc])';   %加上重心处磁场
                Bxcurt=Br(1);
                Bycurt=Br(2);
                Bzcurt=Br(3); 
                Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2); 

                %             step=Bmcurt;
                step=1;
                stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;          
                Xprev=Xcurt-stepvec(1); 
                Yprev=Ycurt-stepvec(2);
                Zprev=Zcurt-stepvec(3); 

                Xline(ii)=Xcurt;  
                Yline(ii)=Ycurt; 
                Zline(ii)=Zcurt;
                Bmline(ii)=Bmcurt; 
                ii=ii+1;
                end

            Xline_preva(ns,1:length(Xline)) = Xline; 
            Yline_preva(ns,1:length(Xline)) = Yline; 
            Zline_preva(ns,1:length(Xline)) = Zline; 
            Bmline_preva(ns,1:length(Xline)) = Bmline; 

            clear Xline 
            clear Yline 
            clear Zline
            clear Bmline 
            clear Nlin

            ii=1;
            Xprev2=Xgrid;  Yprev2=Ygrid;  Zprev2=Zgrid;

            while abs(Xprev2)<BoxWid && abs(Yprev2)<BoxWid && abs(Zprev2)<BoxWid 
                Xcurt2=Xprev2;
                Ycurt2=Yprev2;
                Zcurt2=Zprev2;


                Br2 = dB_null_b * ([Xcurt2 Ycurt2 Zcurt2])';
    %             Br = Br+([Bxbc Bybc Bzbc])';   %加上重心处磁场
                Bxcurt2=Br2(1);
                Bycurt2=Br2(2);
                Bzcurt2=Br2(3);         
                Bmcurt2=sqrt(Br2(1).^2+Br2(2).^2+Br2(3).^2);
    %             step=Bmcurt;
                step=1;
                stepvec2=[Bxcurt2 Bycurt2 Bzcurt2]/norm([Bxcurt2 Bycurt2 Bzcurt2])*step;
                Xprev2=Xcurt2-stepvec2(1);
                Yprev2=Ycurt2-stepvec2(2);
                Zprev2=Zcurt2-stepvec2(3);

                Xline2(ii)=Xcurt2; 
                Yline2(ii)=Ycurt2;
                Zline2(ii)=Zcurt2;
                Bmline2(ii)=Bmcurt2;
                ii=ii+1;
            end

            Xline_prevb(ns,1:length(Xline2)) = Xline2;
            Yline_prevb(ns,1:length(Xline2)) = Yline2;
            Zline_prevb(ns,1:length(Xline2)) = Zline2;
            Bmline_prevb(ns,1:length(Xline2)) = Bmline2;

            clear Xline2 
            clear Yline2
            clear Zline2
            clear Bmline2 


            %-----positive trace----one side
            ii=1;
            Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
            Xnext2=Xgrid;  Ynext2=Ygrid;  Znext2=Zgrid;
            while abs(Xnext)<BoxWid & abs(Ynext)<BoxWid & abs(Znext)<BoxWid 
                Xcurt=Xnext;
                Ycurt=Ynext;
                Zcurt=Znext;

                Br = dB_null_a * ([Xcurt Ycurt Zcurt])';
                %             Br = Br+([Bxbc Bybc Bzbc])';  %加上重心处磁场
                Bxcurt=Br(1);
                Bycurt=Br(2);
                Bzcurt=Br(3);
                Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
                %             step=Bmcurt;
                step=1;
                stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
                Xnext=Xcurt+stepvec(1);
                Ynext=Ycurt+stepvec(2);
                Znext=Zcurt+stepvec(3);

                Xline(ii)=Xcurt;
                Yline(ii)=Ycurt;
                Zline(ii)=Zcurt;
                Bmline(ii)=Bmcurt;
                ii=ii+1;
            end

            Xline_nexta(ns,1:length(Xline)) = Xline;
            Yline_nexta(ns,1:length(Xline)) = Yline;
            Zline_nexta(ns,1:length(Xline)) = Zline;
            Bmline_nexta(ns,1:length(Xline)) = Bmline;

            clear Xline
            clear Yline
            clear Zline
            clear Bmline
            clear Nlin

            %-----positive trace----the other side
            ii=1;
            Xnext2=Xgrid;  Ynext2=Ygrid;  Znext2=Zgrid;
            while abs(Xnext2)<BoxWid && abs(Ynext2)<BoxWid && abs(Znext2)<BoxWid 
                Xcurt2=Xnext2;
                Ycurt2=Ynext2;
                Zcurt2=Znext2;

                Br2 = dB_null_b * ([Xcurt2 Ycurt2 Zcurt2])';
                %             Br = Br+([Bxbc Bybc Bzbc])';  %加上重心处磁场
                Bxcurt2=Br2(1);
                Bycurt2=Br2(2);
                Bzcurt2=Br2(3);
                Bmcurt2=sqrt(Br2(1).^2+Br2(2).^2+Br2(3).^2);
                %             step=Bmcurt;
                step=1;
                stepvec2=[Bxcurt2 Bycurt2 Bzcurt2]/norm([Bxcurt2 Bycurt2 Bzcurt2])*step;
                Xnext2=Xcurt2+stepvec2(1);
                Ynext2=Ycurt2+stepvec2(2);
                Znext2=Zcurt2+stepvec2(3);

                Xline2(ii)=Xcurt2;
                Yline2(ii)=Ycurt2;
                Zline2(ii)=Zcurt2;
                Bmline2(ii)=Bmcurt2;
                ii=ii+1;
            end

            Xline_nextb(ns,1:length(Xline2)) = Xline2;
            Yline_nextb(ns,1:length(Xline2)) = Yline2;
            Zline_nextb(ns,1:length(Xline2)) = Zline2;
            Bmline_nextb(ns,1:length(Xline2)) = Bmline2;
            ns = ns +1;
            clear Xline2
            clear Yline2
            clear Zline2
            clear Bmline2
            clear Nlin
        end
    end
end

%% plot magnetic field lines
pa=1;pb=1;na=1;nb=1;arrP=0.5;ap=15;
cmin1=0;
cmax1=75;


    for nk = 1:length(Xline_preva(:,1))
    ns0=find(Xline_preva(nk,:)~=0); a=size(ns0);a=a(2);
    if a~=0
    ns0 = ns0(length(ns0));  
    if Zline_preva(nk,ns0)<0
        switch plottype
            case('line')
                plot3(gca, Xline_preva(nk,1:ns0), Yline_preva(nk,1:ns0), Zline_preva(nk,1:ns0),'b'); hold on;
            otherwise
                cline(Xline_preva(nk,1:ns0), Yline_preva(nk,1:ns0), Zline_preva(nk,1:ns0), Bmline_preva(nk,1:ns0), cmin1, cmax1, jet);hold on;
        end
    % arrow3matrix_preva(pa,1:3)=[Xline_preva(nk,fix(arrP*ns0)),Yline_preva(nk,fix(arrP*ns0)),Zline_preva(nk,fix(arrP*ns0))];
    % arrow3matrix_preva(pa,4:6)=[Xline_preva(nk,fix(arrP*ns0)-ap),Yline_preva(nk,fix(arrP*ns0)-ap),Zline_preva(nk,fix(arrP*ns0)-ap)];
    pa = pa+1;
    end
    ns0=find(Xline_nexta(nk,:)~=0); ns0 = ns0(length(ns0));
    if Zline_nexta(nk,ns0)<0
        switch plottype
            case('line')
                plot3(gca, Xline_nexta(nk,1:ns0), Yline_nexta(nk,1:ns0), Zline_nexta(nk,1:ns0),'r'); hold on;
            otherwise
                cline(Xline_nexta(nk,1:ns0), Yline_nexta(nk,1:ns0), Zline_nexta(nk,1:ns0), Bmline_nexta(nk,1:ns0), cmin1, cmax1, jet);hold on;
        end
    % arrow3matrix_nexta(na,1:3)=[Xline_nexta(nk,fix(arrP*ns0)),Yline_nexta(nk,fix(arrP*ns0)),Zline_nexta(nk,fix(arrP*ns0))];
    % arrow3matrix_nexta(na,4:6)=[Xline_nexta(nk,fix(arrP*ns0)+ap),Yline_nexta(nk,fix(arrP*ns0)+ap),Zline_nexta(nk,fix(arrP*ns0)+ap)];
    na = na+1;
    end
    end  
    end

    for nk = 1:length(Xline_prevb(:,1))
    ns0=find(Xline_preva(nk,:)~=0); a=size(ns0);a=a(2);
    if a~=0
    ns0 = ns0(length(ns0));  
    if Zline_prevb(nk,ns0)>0
         switch plottype
             case('line')
                 plot3(gca, Xline_prevb(nk,1:ns0), Yline_prevb(nk,1:ns0), Zline_prevb(nk,1:ns0),'b'); hold on;
             otherwise
                 cline(Xline_prevb(nk,1:ns0), Yline_prevb(nk,1:ns0), Zline_prevb(nk,1:ns0), Bmline_prevb(nk,1:ns0), cmin1, cmax1, jet);hold on;
         end
    % arrow3matrix_prevb(pb,1:3)=[Xline_prevb(nk,fix(arrP*ns0)),Yline_prevb(nk,fix(arrP*ns0)),Zline_prevb(nk,fix(arrP*ns0))];
    % arrow3matrix_prevb(pb,4:6)=[Xline_prevb(nk,fix(arrP*ns0)-ap),Yline_prevb(nk,fix(arrP*ns0)-ap),Zline_prevb(nk,fix(arrP*ns0)-ap)];
    pb=pb+1;
    end
    ns0=find(Xline_nextb(nk,:)~=0); ns0 = ns0(length(ns0));
    if Zline_nextb(nk,ns0)>0
          switch plottype
              case('line')
                  plot3(gca, Xline_nextb(nk,1:ns0), Yline_nextb(nk,1:ns0), Zline_nextb(nk,1:ns0),'r'); hold on;
              otherwise
                  cline(Xline_nextb(nk,1:ns0), Yline_nextb(nk,1:ns0), Zline_nextb(nk,1:ns0), Bmline_nextb(nk,1:ns0), cmin1, cmax1, jet);hold on;
          end
    % arrow3matrix_nextb(nb,1:3)=[Xline_nextb(nk,fix(arrP*ns0)),Yline_nextb(nk,fix(arrP*ns0)),Zline_nextb(nk,fix(arrP*ns0))];
    % arrow3matrix_nextb(nb,4:6)=[Xline_nextb(nk,fix(arrP*ns0)+ap),Yline_nextb(nk,fix(arrP*ns0)+ap),Zline_nextb(nk,fix(arrP*ns0)+ap)];
    nb = nb+1;
    end
    end   
    end


% Jlinelngth=BoxWid;
% Jendp1=Jlinelngth*irf_norm(j_null);
% Jendp2=-Jendp1;
% arrow3(Jendp2, Jendp1,'s5',3,5); hold on;
% light('position',Jendp1); lighting gouraud;

%--------------------------------------------


%% PLOT
plot3(gca, [Rsc1_a(1) Rsc2_a(1) Rsc3_a(1) Rsc4_a(1)], [Rsc1_a(2) Rsc2_a(2) Rsc3_a(2) Rsc4_a(2)], ...
[Rsc1_a(3) Rsc2_a(3) Rsc3_a(3) Rsc4_a(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc2_a(1) Rsc4_a(1) Rsc1_a(1) Rsc3_a(1)], [Rsc2_a(2) Rsc4_a(2) Rsc1_a(2) Rsc3_a(2)], ...
[Rsc2_a(3) Rsc4_a(3) Rsc1_a(3) Rsc3_a(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc1_a(1)], [Rsc1_a(2)],[Rsc1_a(3)], 'ks', 'Linewidth',1, ...
'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12); hold on;
plot3(gca, [Rsc2_a(1)], [Rsc2_a(2)],[Rsc2_a(3)], 'rs', 'Linewidth',1, ...
'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12); hold on;
plot3(gca, [Rsc3_a(1)], [Rsc3_a(2)],[Rsc3_a(3)], 'gs', 'Linewidth',1, ...
'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12); hold on;
plot3(gca, [Rsc4_a(1)], [Rsc4_a(2)],[Rsc4_a(3)], 'bs', 'Linewidth',1, ...
'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12); hold on;

% text(0, -0.95*Flim,0.95*Flim,epoch2iso(gradB(idxnull_1,1)));      
% quiver3(0,0,0,100*V(1,1),100*V(1,2),100*V(1,3)); hold on
% quiver3(0,0,0,100*V(2,1),100*V(2,2),100*V(2,3)); hold on
% quiver3(0,0,0,100*V(3,1),100*V(3,2),100*V(3,3)); hold off
% k = irf_newxyz([0,0,1], e1_spine,e2_spine,e3_spine);
% quiver3(0,0,0,100*k(1),100*k(2),100*k(3),'b','linewidth',12)
% quiver3(0,0,0,100*e3_spine(1),100*e3_spine(2),100*e3_spine(3),'g','linewidth',12)
% quiver3(0,0,0,100*e2_spine(1),100*e2_spine(2),100*e2_spine(3),'g','linewidth',12)
% quiver3(0,0,0,100*e1_spine(1),100*e1_spine(2),100*e1_spine(3),'g','linewidth',12)

% maxB=max(Bmax);
% minB=min(Bmin);
box on       
% hcb=colorbar('peer',h(1),'North');
hcb=colorbar('peer',h(1),'East','Tickdirection','out');
posFig=get(h(1),'Position'); 
left=posFig(1)+posFig(3)*1; low=posFig(2)+posFig(4)*0.87; width=posFig(3)/40; height=0.1;
% left=posFig(1); low=posFig(4)*1; width=posFig(3); height=posFig(4)/35;
set(hcb,'Position',[left low width height]);
ylabel(hcb,'|B| (nT)');


caxis([cmin1, cmax1]);   %here is derived from minB & maxB. it should be consistent with cline range
set(gca,'Xlim',[-Flimx Flimx], 'Ylim',[-Flimy Flimy], 'Zlim',[-Flim Flim]);
set(gca,'xtick',[-Flim_ax:Flim_bx:Flim_ax], 'ytick',[-Flim_ay:Flim_by:Flim_ay], 'ztick',[-Flim_a:Flim_b:Flim_a],'fontsize',13);
set(gca,'DataAspectRatio',[1 1.0 1]);
xlabel(gca,'e_{1} [km]');
ylabel(gca,'e_{2} [km]');
zlabel(gca,'e_{3} [km]');
% grid off;
% box off
a=0.5;b=1.5;c=0;

% arrow3(arrow3matrix_preva(:,1:3),arrow3matrix_preva(:,4:6),'b-0',a)
% arrow3(arrow3matrix_prevb(:,1:3),arrow3matrix_prevb(:,4:6),'b-0',a)
% arrow3(arrow3matrix_nexta(:,1:3),arrow3matrix_nexta(:,4:6),'b-0',a)
% arrow3(arrow3matrix_nextb(:,1:3),arrow3matrix_nextb(:,4:6),'b-0',a)

% angles=get(gca,'view');
% angles=[-37.5,30];
% angles=[45,30];
% set(gca,'view',angles)
view(45,30)


%% save figure
 set(fig1,'render','painters');  %矢量
figname=['B_around_null'];
print(fig1, '-dpdf', [figname '.pdf']);


set(fig1,'renderer','opengl'); %位图

% axis off;
% figname=['Btopology_3D_e1e2e3'];
% print(fig1, '-dpng','-r400',[figname '.png']);
% print(fig1, '-dpdf','-r400',[figname '.pdf']);

% set(gca,'view',[0 90])
% figname=['Btopology_xy_plane'];
% print(fig1, '-dpng','-r400',[figname '.png']);

% set(gca,'view',[0 0])
% figname=['Btopology_xz_plane'];
% print(fig1, '-dpng','-r400',[figname '.png']);

% set(gca,'view',[135 26])
% zoom out
% figname=['Btopology_3D_plane_'];
% figname=[figname Tnull1(18:23)];
% print(fig1, '-dpng','-r400',[figname '.png']);

% set(gca,'view',[90 0])
% figname=['Btopology_2D_plane_'];
% figname=[figname Tnull1(18:23)];

% print(fig1, '-dpng','-r400',[figname '.png']);
% mark1=Blmn1(tInd,1);
% irf_pl_mark(h,mark1,'r','linewidth',0.8);


Tcomputend=clock;
Telaps=etime(Tcomputend, Tcomputsta)/60  %minute