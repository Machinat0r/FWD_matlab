function [Bx,By,Bz,Bt] = calculateB(Rx,Ry,Rz,Q,Loc)
d = [Rx,Ry,Rz]-Loc;
d2 = sum(d.^2,2);
Brec = Q/(4*pi*d2);
Brec = [Brec*Rx/sqrt(d2),Brec*Ry/sqrt(d2),Brec*Rz/sqrt(d2)];
% % %     gridBt(i) = sqrt(sum(gridB{i}.*gridB{i},2));
Bx = Brec(1); By = Brec(2); Bz = Brec(3);
Bt = sqrt(sum(Brec.^2,2));
end