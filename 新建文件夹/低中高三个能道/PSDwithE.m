clear 
clc

ic=1;

tint=irf.tint('2017-06-11T01:59:35Z/2017-06-11T01:59:55Z');

c_eval('fpiFile? = dataobj(''D:\Matlab\xy-matlab\MMS\mms_db\data\mms1\fpi\brst\l2\des-dist\2017\06\11\mms1_fpi_brst_l2_des-dist_20170611015853_v3.2.0.cdf'');',ic);
c_eval('diste_struct = get_variable(fpiFile?,''mms?_des_dist_brst'');',ic);
thetae=diste_struct.DEPEND_2.data;
energye0=diste_struct.DEPEND_3.data(1,1:end);
energye1=diste_struct.DEPEND_3.data(2,1:end);

c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint);',ic);
% c_eval('energy0=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy0_brst'',Tintr);',ic);
% c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_energy1_brst'',Tintr);',ic);
c_eval('phie=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_phi_brst'',tint);',ic);
% c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'',''mms?_des_theta_brst'',Tintr);',ic);
c_eval('stepTablee=mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_steptable_parity_brst'',tint);',ic);


diste.data=diste.data*1e30;

energyspec = ones(length(diste.time),1)*energye0;
for ii = 1:length(diste.time);
    if stepTablee.data(ii),
        energyspec(ii,:) = energye1;
    end
end


% define angles
dangle = pi/16;
lengthphi = 32;

z2 = ones(lengthphi,1)*sind(thetae);
solida = dangle*dangle*z2;
allsolide = zeros(size(diste.data));

for ii = 1:length(diste.time);
    for jj=1:length(energye0);
        allsolide(ii,jj,:,:) = solida;
    end
end

distes = diste.data.*allsolide;

% Electron analysis - OMNI
PSDomni = zeros(length(diste.time),length(energye0));
for ii = 1:length(diste.time);
    disttemp = squeeze(distes(ii,:,:,:));
    PSDomni(ii,:) = squeeze(irf.nanmean(irf.nanmean(disttemp,2),3))/(mean(mean(solida)));
end

timeomni=irf_time(phie.time,'epochTT>epoch');
time1='2017-06-11T01:59:35.00Z';
time1=iso2epoch(time1);

time2='2017-06-11T01:59:40.00Z';
time2=iso2epoch(time2);

time3='2017-06-11T01:59:45.00Z';
time3=iso2epoch(time3);

time4='2017-06-11T01:59:50.00Z';
time4=iso2epoch(time4);

idt1=find(timeomni>time1);
t=idt1(1);
chan=energyspec(t,:);
    psd=PSDomni(t,:);
    psd=transpose(psd);
    loglog(chan,psd,'b+');
    hold on;



for i=2:4
c_eval('idt?=find(timeomni>time?);',i);
c_eval('t?=idt?(1);',i);
c_eval('chan?=energyspec(t?,:);',i);

    c_eval('psd?=PSDomni(t?,:);',i);
    c_eval('psd?=transpose(psd?);',i);

end


loglog(chan2,psd2,'gs');
loglog(chan3,psd3,'ko');
loglog(chan4,psd4,'r^');

ylabel('PSD (s^3 km^{-6})');
xlabel('E (ev)');
% set(gca,'Xlim',[1700 6000])
irf_legend({'+ 01:59:35.00'},[0.2 0.55],'color','b')
irf_legend({'¡õ 01:59:40.00'},[0.2 0.45],'color','g')
irf_legend({'¡ð 01:59:45.00'},[0.2 0.35],'color','k')
irf_legend({'¦¤ 01:59:50.00'},[0.2 0.25],'color','r')






X=chan3;
Y=psd3
Y=log10(Y);
n=[-1e-30 -1e-26];
T=[11 100];

k=1.38e-23; 
m=9.1*10^-31;
e=1.6e-19;


Var=zeros(32,32);
n1=n(1);
n2=n(2);
T1=T(1);
T2=T(2);
N=linspace(n1,n2,32);
N=N';
T=linspace(T1,T2,32);
for i=1:32
    for ii=1:32
        Te=T(i);
        Ne=N(ii);
        F=Ne*1e24*power(m/(2*pi*k*Te*11605),3/2)*exp(-X*e/(k*Te*11605));
%        
        F=F'%ºó²¹
%         F=F-Y;
        Var(ii,i)=sum(abs(F));
    end
end

i_subplot=1;
plot(X,F)
set(gca, 'Xlim',[0 100000]);
set(gca, 'Ylim',[-1e-36 -1e-24]);
% end