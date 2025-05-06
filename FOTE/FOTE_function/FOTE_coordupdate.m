function [Null_loc, dB_null] = FOTE_coordupdate(R1,R2,R3,R4,B1,B2,B3,B4,Tnull,e1,e2,e3)
if nargin~= 12 && nargin~=6
	disp('Too few parameters. See usage:');
	help FOTE;     
	return
end

if nargin==6
	rs = R1;
	bs = R2;
    Tnull = R3;
    e1 = R4; e2 = R5; e3 = R6;
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(rs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('R? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(bs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('B? =ttt;',cl_id)); clear ttt
    end    
	clear bs rs
end

for ic=1:4
    c_eval(['Rp?=irf_newxyz(R?, e1,e2,e3);'],ic);
    c_eval(['Bp?=irf_newxyz(B?, e1,e2,e3);'],ic);
end
gradB=c_4_grad('Rp?','Bp?','grad');
idxnull = find(gradB(:,1,1)>=iso2epoch(Tnull)); idxnull=idxnull(1);
dB_null=reshape(gradB(idxnull,2:end),3,3);
dR1=inv(dB_null) * (Bp1(idxnull,2:4))';   

Rsc1=Rp1(idxnull,2:end);
Rsc2=Rp2(idxnull,2:end);
Rsc3=Rp3(idxnull,2:end);
Rsc4=Rp4(idxnull,2:end);

Rbarycenter=(Rsc1+Rsc2+Rsc3+Rsc4)./[4.0 4.0 4.0];    %

R_null=Rsc1-dR1';
dRsc1=Rsc1-R_null;
dRsc2=Rsc2-R_null;
dRsc3=Rsc3-R_null;
dRsc4=Rsc4-R_null;
Null_loc = struct('R_sc1', dRsc1, 'R_sc2', dRsc2,'R_sc3', dRsc3,'R_sc4', dRsc4,'Rbc', Rbarycenter);
end