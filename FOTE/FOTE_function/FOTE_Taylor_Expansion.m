function [Null_info_dis, Null_info_type, Null_info_err] = FOTE_Taylor_Expansion(R1,R2,R3,R4,B1,B2,B3,B4,smooth_span,smooth_method)

if nargin~=8 && nargin~=9 && nargin~=10 && nargin~=2 && nargin~=3 && nargin~=4
	disp('Too few parameters. See usage:');
	help FOTE;     
	return
end

if nargin==3 
    smooth_span = R3;
end
if nargin==4
    smooth_span = R3;
    smooth_method = R4;
end

if nargin == 2 || nargin == 3 || nargin == 4 
	rs = R1;
	bs = R2;
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(rs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('R? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(bs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('B? =ttt;',cl_id)); clear ttt
    end    
	clear bs rs
end


if nargin == 9|| nargin == 3
    for ic = 1:4
        c_eval('B?(:,2) = smooth(B?(:,2),smooth_span);',ic)
        c_eval('B?(:,3) = smooth(B?(:,3),smooth_span);',ic)
        c_eval('B?(:,4) = smooth(B?(:,4),smooth_span);',ic)
      
    end
end

if nargin == 10|| nargin == 4
    for ic = 1:4
        c_eval('B?(:,2) = smooth(B?(:,2),smooth_span,smooth_method);',ic)
        c_eval('B?(:,3) = smooth(B?(:,3),smooth_span,smooth_method);',ic)
        c_eval('B?(:,4) = smooth(B?(:,4),smooth_span,smooth_method);',ic)
        
    end
end
for ic = 1:4
    c_eval('Bo? = B?;',ic)
    c_eval('Ro? = R?;',ic)
end
gradB=c_4_grad('R?','B?','grad');
gradB2 = gradB;
d_dot_B=c_4_grad('R?','B?','div');
d_cros_B=c_4_grad('R?','B?','curl');
[j,divB] = c_4_j('R?','B?');

temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);
err_4C(:,2)=abs(err_4C(:,2))*100;

temp=irf_abs(d_cros_B);
ksd=ones(length(temp),2);
ksd(:,1)=temp(:,1);
angle=ones(length(temp),2);
angle(:,1)=temp(:,1);
angle2=ones(length(temp),2);
angle2(:,1)=temp(:,1);
d_cros_B_mag=temp(:,[1 5]);
err_curlmeter=irf_multiply(1,d_dot_B,1,d_cros_B_mag,-1);
err_curlmeter(:,2)=abs(err_curlmeter(:,2))*100;

%Null type identification
for ii=1:length(d_dot_B(:,1))
    mksizSim(ii)=4;  
    thresold=0.5; 
    deltB_null=reshape(gradB(ii,2:end),3,3);
    [V,D] = eig(deltB_null);
    
    
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        [m,n]=find(D==min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))])));
        [mm,nn]=find(D==-(min(([abs(D(1,1)) abs(D(2,2)) abs(D(3,3))]))));
        if isempty(m)
            m=mm;
        else m=m;
        end
        D_sel = D(m,m)*[D(1,1),D(2,2),D(3,3)];
        m0 = m;
        if sum(D_sel>0) == 1
            m =  find(logical((D_sel~=min(D_sel)).* (D_sel<0)));
        end
        
        e1=V(:,m)';
        e1=irf_norm(e1);
        etmp=[0 1 0];
        e3=cross(e1,etmp);          e3=irf_norm(e3);
        e2=cross(e3,e1);
        
        for ic=1:4
            c_eval(['Rp?=irf_newxyz(Ro?, e1,e2,e3);'],ic);
            c_eval(['Bp?=irf_newxyz(Bo?, e1,e2,e3);'],ic);
        end
        
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
        
%         angle01=acosd(c);
%         angle02=acosd(c2);
%         angle11=acosd(abs(c));
%         angle12=acosd(abs(c2));
        
        
       
        
    else [m,n]=find(imag(D)==0&D~=0);
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
            c_eval(['Bp?=irf_newxyz(Bo?, e3,e1,e2);'],ic);
            c_eval(['Rp?=irf_newxyz(Ro?, e3,e1,e2);'],ic);
            eS=irf_newxyz(eS, e3,e1,e2);
        end
        c = abs(eS(1));
        c2 = eS(1);
%         angle11=acosd(c);
    end
    
    gradB3=c_4_grad('Rp?','Bp?','grad');
    dB_null3=reshape(gradB3(ii,2:end),3,3);  
    dR1=inv(dB_null3) * (Bp1(ii,2:4))';
    R_null=Rp1(ii,2:4)-dR1';
    
    for nsk = 1:4
        c_eval(['Rts_sc?(ii,:) = Rp?(ii,2:end);'],nsk);
        c_eval(['dR2Dmag?(ii,1)=Rp1(ii,1);'],nsk);
        c_eval(['Rts_sc?(ii,:) = Rts_sc?(ii,:)-R_null;'],nsk);
        c_eval(['Rts_sc?(ii,1) = 0;'],nsk);
        c_eval(['dR2Dmag?=irf_abs(Rts_sc?);'],nsk);
        
    end
    
    dR2Dmin(ii,2)=min([dR2Dmag1(ii,4) dR2Dmag2(ii,4) dR2Dmag3(ii,4) dR2Dmag4(ii,4)], [], 2);
    dR2Dmin(ii,1)=dR2Dmag1(ii,1);

    angle(ii,2)=acosd(abs(c));
    angle2(ii,2)=acosd(abs(c2));
    
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
       ksd(ii,2) = 100*min(abs([D(1,1) D(2,2) D(3,3)]))/sum(abs([D(1,1) D(2,2) D(3,3)]));
    else ksd(ii,2) = 100*abs(D(m,m))/sum(abs([D(1,1) D(2,2) D(3,3)]));
    end
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='w';
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='w';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            type(ii)='X'; clr(ii)='k'; faceclr(ii)='w';
        end
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            type(ii)='>'; clr(ii)='b'; faceclr(ii)='b';
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                type(ii)='^'; clr(ii)='r'; faceclr(ii)='r';
            else
                type(ii)='s'; clr(ii)='k'; faceclr(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            type(ii)='o'; clr(ii)='k'; faceclr(ii)='w';
        end
    end
    %=========================================================
    
    %=========================================================
    if max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) == 0
        if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='w';mksizSim(ii)=8;
        else
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='w';mksizSim(ii)=8;
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if min(abs([D(1,1) D(2,2) D(3,3)]))/max(abs([D(1,1) D(2,2) D(3,3)]))<thresold
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 2
                typeSim(ii)='X'; clrSim(ii)='b'; faceclrSim(ii)='w'; mksizSim(ii)=13;
            end
            if length(find([D(1,1) D(2,2) D(3,3)]>0)) == 1
                typeSim(ii)='X'; clrSim(ii)='r'; faceclrSim(ii)='w'; mksizSim(ii)=13;
            end
        end
        if min(abs([D(1,1) D(2,2) D(3,3)]))==0
            typeSim(ii)='X'; clrSim(ii)='k'; faceclrSim(ii)='w'; mksizSim(ii)=13;
        end
        %------------------------------------------------------------------
    else
        if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
            typeSim(ii)='>'; clrSim(ii)='b'; faceclrSim(ii)='b';mksizSim(ii)=8;
        else
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='^'; clrSim(ii)='r'; faceclrSim(ii)='r';mksizSim(ii)=8;
            else
                typeSim(ii)='s'; clrSim(ii)='k'; faceclrSim(ii)='w';
            end
        end
        %------------Simplification Procedure------------------------------
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))/max(abs([imag(D(1,1)) imag(D(2,2)) imag(D(3,3))])) < (thresold*2)
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 2
                typeSim(ii)='o'; clrSim(ii)='b'; faceclrSim(ii)='w';
            end
            if length(find([real(D(1,1)) real(D(2,2)) real(D(3,3))]>0)) == 1
                typeSim(ii)='o'; clrSim(ii)='r'; faceclrSim(ii)='w';
            end
        end
        if max(abs([real(D(1,1)) real(D(2,2)) real(D(3,3))]))==0
            typeSim(ii)='o'; clrSim(ii)='k'; faceclrSim(ii)='w';
        end
        %------------------------------------------------------------------
    end
    %=========================================================
    
    eigVal_err(ii,2)=abs(real(D(1,1)+D(2,2)+D(3,3)))/max(abs([real(D(1,1)), real(D(2,2)), real(D(3,3))])) * 100;
    sumeigVal(ii,2)=D(1,1)+D(2,2)+D(3,3);
    eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;
end

eigVal_err(:,1)=d_dot_B(:,1);
sumeigVal(:,1)=d_dot_B(:,1);
eigVal_err_v2(:,1)=d_dot_B(:,1);


%find null position
for ii=1:length(B1(:,1))
    dBeach=reshape(gradB(ii,2:end),3,3);
    dR1(ii,2:4)=B1(ii,2:4)*inv(dBeach');
    dR2(ii,2:4)=B2(ii,2:4)*inv(dBeach');
    dR3(ii,2:4)=B3(ii,2:4)*inv(dBeach');
    dR4(ii,2:4)=B4(ii,2:4)*inv(dBeach');
end
dR1(:,1)=B1(:,1);
dR2(:,1)=B1(:,1);
dR3(:,1)=B1(:,1);
dR4(:,1)=B1(:,1);

dRmag1=irf_abs(dR1);
dRmag2=irf_abs(dR2);
dRmag3=irf_abs(dR3);
dRmag4=irf_abs(dR4);

dRmin(:,2)=min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1)=dRmag1(:,1);

Null_info_type = struct('t', B1(:,1) ,'type', typeSim, 'color', clrSim, 'f_color', faceclrSim, 'marksize', reshape(mksizSim,[length(mksizSim),1]), 'spine_angle', angle);
Null_info_dis = struct('dR1', dR1, 'dR2', dR2, 'dR3', dR3, 'dR4', dR4,'dRmin', dRmin);
Null_info_err = struct('err1', err_curlmeter,'err2', eigVal_err_v2, 'err3', eigVal_err);


end