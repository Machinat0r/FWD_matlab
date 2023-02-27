function curl = curl_Jac(Jac)
%input: Vector field Jacobian matrix(grad matrix)
%size: N*9 or N*10
%each line: (Time), dxBx, dxBy, dxBz, dyBx, dyBy, dyBz, dzBx, dzBy, dzBz
%output each line:(Time), curlBx,curlBy,curlBz
%Example: curlB=curl_Jac(gradB);
%------written by Wending Fu, Apr.14.2022 in Beijing------------
Time = [];
if  size(Jac,2) == 10
    Time = Jac(:,1);
    Jac = Jac(:,2:end);
elseif size(Jac,2) ~= 9
    disp('Jacobian matrix dimension error. See usage:')
    help curl_Jac
end

curl = zeros(size(Jac,1),3);
for rec = 1:size(Jac,1)
    curl(rec,1) = Jac(rec,6)-Jac(rec,8);
    curl(rec,2) = Jac(rec,3)-Jac(rec,7);
    curl(rec,3) = Jac(rec,2)-Jac(rec,4);
end

if ~isempty(Time)
    curl = [Time,curl];
end
end