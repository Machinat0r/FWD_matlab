function [Null_loc, Trans_mat, dB_null] = FOTE(R1,R2,R3,R4,B1,B2,B3,B4,Tnull)

% [ax,args,nargs] = axescheck(varargin{:});
% B1 = args{1}; B2 = args{2}; B3 = args{3}; B4 = args{4};
% R1 = args{5}; R2 = args{6}; R3 = args{7}; R4 = args{8};
% Tnull =  args{9};

if nargin~=9 && nargin~=3
	disp('Too few parameters. See usage:');
	help FOTE;     
	return
end

if nargin==3
	rs = R1;
	bs = R2;
    Tnull = R3;
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(rs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('R? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(bs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('B? =ttt;',cl_id)); clear ttt
    end    
	clear bs rs
end

%error of curolmeter
c_eval(['R?=irf_resamp(R?,B?);']);
[j,divB] = c_4_j('R?','B?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

%FOTE
gradB=c_4_grad('R?','B?','grad');
idxnull = find(gradB(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);

dB_null=reshape(gradB(idxnull,2:end),3,3);
[V,D] = eig(dB_null);

% eigVal_err=abs(real(D(1,1)+D(2,2)+D(3,3)))/max(abs([real(D(1,1)), real(D(2,2)), real(D(3,3))])) * 100;
% sumeigVal=D(1,1)+D(2,2)+D(3,3);
eigVal_err_v2=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;

if  err_4C(idxnull,2)>50 || eigVal_err_v2>50
    disp('error parameter is large!!');
    disp(err_4C(idxnull,2));
    disp(eigVal_err_v2);
end
    Axis_m = [1 0 0;0 1 0;0 0 1];
    Axis_m0 = [1 0 0;0 1 0;0 0 1];
    
if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
    [m,n]=find(D==min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))])));
    [mm,nn]=find(D==-(min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]))));
    if isempty(m)
        m=mm;
    end
    D_sel = D(m,m)*[D(1,1),D(2,2),D(3,3)];
    m0 = m;
    if sum(D_sel>0) == 1
        m =  find(logical((D_sel~=min(D_sel)).* (D_sel<0)));
    end
    
    eS=V(:,m)';
    e1=V(:,m)';
    e1=irf_norm(e1); etmp=[0 1 0]; e3=cross(e1,etmp); e3=irf_norm(e3);
    e2=cross(e3,e1);
    
    for ic=1:4
        c_eval(['Rp?=irf_newxyz(R?, e1,e2,e3);'],ic);
        c_eval(['Bp?=irf_newxyz(B?, e1,e2,e3);'],ic);
    end
    
    eS2 = irf_newxyz(eS, e1,e2,e3);
    Axis_m = irf_newxyz(Axis_m, e1,e2,e3);
    
    V1=V';
    V=irf_newxyz(V1, e1,e2,e3);
    V([m,3],:)=V([3,m],:);
    
    a1=V(1,:);
    a2=V(2,:);
    a3=V(3,:);
    as1=[0,a1(2),a1(3)];
    as2=[0,a2(2),a2(3)];
    
    b1=as1/norm(as1);
    b2=as2/norm(as2);
    b3=(b1+b2)/norm(b1+b2);
    c=dot(b1,b2/(norm(b1)*norm(b2)));
    c2=dot(a1,a2/(norm(a1)*norm(a2)));
    
    %     angle01=acosd(c);
    %     angle02=acosd(c2);
    %     angle11=acosd(abs(c));
    %     angle12=acosd(abs(c2));
    
else
    [m,n]=find(imag(D)==0 & D~=0);
    
    ii001 = 1;
    for ii00 = 1:3
        if ii00~=m
            MMn(ii001)=ii00;
            ii001 = ii001+1;
        end
    end
    
    e1=irf_norm((V(:,MMn(1))+V(:,MMn(2)))');
    e0= V(:,MMn(2))';
    eS= irf_norm(V(:,m)');
    
    for ii01 = 1:3
        if imag(e0(ii01))~=0
            e0(ii01)=imag(e0(ii01));
        else
            e0(ii01)=0; %e0(ii)
        end
    end
    e0 = irf_norm(e0);
    e3=cross(e1,e0);
    e2=cross(e3,e1);
    
    for ic=1:4
        c_eval(['Bp?=irf_newxyz(B?, e3,e1,e2);'],ic);
        c_eval(['Rp?=irf_newxyz(R?, e3,e1,e2);'],ic);
        
    end
    eS2=irf_newxyz(eS, e3,e1,e2);
    Axis_m = irf_newxyz(Axis_m, e1,e2,e3);
    c = abs(eS2(1));
    c2 = eS2(1);
end

Axis_m3 = irf_newxyz(Axis_m0, Axis_m(1,:),Axis_m(2,:),Axis_m(3,:));

gradB=c_4_grad('Rp?','Bp?','grad');
dB_null=reshape(gradB(idxnull,2:end),3,3);

dR1=inv(dB_null) * (Bp1(idxnull,2:4))';   
% R_null=R1(idxnull,2:4)-dR1';

Rsc1=Rp1(idxnull,2:end);
Rsc2=Rp2(idxnull,2:end);
Rsc3=Rp3(idxnull,2:end);
Rsc4=Rp4(idxnull,2:end);

Rbarycenter=(Rsc1+Rsc2+Rsc3+Rsc4)./[4.0 4.0 4.0];    %

R_null=Rsc1-dR1';
%  R_null(1)=0;
dRsc1=Rsc1-R_null;
dRsc2=Rsc2-R_null;
dRsc3=Rsc3-R_null;
dRsc4=Rsc4-R_null;
%save data
Trans_mat = struct('t',Tnull,'e1',Axis_m3(1,:),'e2',Axis_m3(2,:),'e3',Axis_m3(3,:),'es_xyz',eS,'es_e123',eS2);
Null_loc = struct('R_sc1', dRsc1, 'R_sc2', dRsc2,'R_sc3', dRsc3,'R_sc4', dRsc4,'Rbc', Rbarycenter);
end