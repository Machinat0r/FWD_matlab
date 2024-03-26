function [Br1,Br2,Br3,Br4,curlB] =  EliminateCurl(varargin)
% Input: 
% R:Location for 4 S/C, Format for each line: (Time) Rx Ry Rz
% B:Magnetic field strength for 4 S/C, Format for each line: (Time) Bx By Bz
% idx(optional parameters):Only a portion of the B is eliminated curl
% Output: Irrotational field B for 4 S/C
% Example: [Br1,Br2,Br3,Br4] =  EliminateCurl('R?_gsm','B?_gsm',1:length(B1));
%------written by Wending Fu, Apr.25.2022 in Beijing------------
if length(varargin) == 4 || length(varargin) == 10
    idx_R = varargin{end};
    varargin = varargin(1:end-1);
    idx_B = varargin{end};
    varargin = varargin(1:end-1);
if length(varargin) == 2
    Rs = varargin{1};
	Bs = varargin{2};
    for id=1:4
        temp = evalin('caller',irf_ssub(Rs,id));
        eval(irf_ssub('R? =temp;',id));
        temp = evalin('caller',irf_ssub(Bs,id));
        eval(irf_ssub('B? =temp;',id)); clear ttt
    end
elseif length(varargin) == 8
    c_eval('R? = varargin{?};',1:4);
    c_eval('B? = varargin{?+4};',1:4);
else
    disp('Incorrect number of input parameters. See usage:')
    help EliminateCurl
    return
end

TimeZero = datestr(datenum(1970,1,1,0,0,0)+B1(idx_B,1)/86400,'yyyy-mm-ddTHH:MM:SS.FFFZ');
[~, Trans_mat,~] = FOTE('R?','B?',TimeZero);
e = [Trans_mat.e1;Trans_mat.e2;Trans_mat.e3]';

c_eval('B? = irf_resamp(B?,B1);',2:4);
if size(R1,2) >= 4, c_eval('tR? = R?(:,1);');c_eval('R? = R?(:,2:4);');end
if size(B1,2) >= 4, c_eval('tB? = B?(:,1);');c_eval('B? = B?(:,2:4);');end
c_eval('B? = B?(idx_B,:);')
c_eval('R? = R?(idx_R,:);')

gradB=c_4_grad('R?','B?','grad');
deltB_null=reshape(gradB,3,3);
[V,D] = eig(deltB_null);
a = real(D(1,1));b = real(D(2,2));c = real(D(3,3));
Y = e*[a,0,0;0,b,0;0,0,c]*e';
% % % % x11 = deltB_null(1,1); x22 = deltB_null(2,2); x33 = deltB_null(3,3);
% % % syms a b c x11 x22 x33 y12 y13 y23
% % % Y = [x11,y12,y13;y12,x22,y23;y13,y23,x33];
% % % f1 = Y*e(:,1) == a*e(:,1);
% % % f2 = Y*e(:,2) == b*e(:,2);
% % % f3 = Y*e(:,3) == c*e(:,3);
% % % f4 = (a-x11)*(a-x22)*(a-x33)-2*y12*y13*y23-(a-x33)*y12^2-(a-x11)*y23^2-(a-x22)*y13^2;
% % % f5 = (b-x11)*(b-x22)*(b-x33)-2*y12*y13*y23-(b-x33)*y12^2-(b-x11)*y23^2-(b-x22)*y13^2;
% % % f6 = (c-x11)*(c-x22)*(c-x33)-2*y12*y13*y23-(c-x33)*y12^2-(c-x11)*y23^2-(c-x22)*y13^2;
% % % % % % f4 = norm(deltB_null,2) == norm(Y,2);
% % % [a,b,c,x11,x22,x33,y12,y23,y13] = vpasolve([f1,f2,f3],[a,b,c,x11,x22,x33,y12,y23,y13]);
% % % % [a,b,c,x11,x22,x33,y12,y23,y13] = solve([f1,f2,f3],[a,b,c,x11,x22,x33,y12,y23,y13]);
% % % y12 = double(y12);
% % % y13 = double(y13);
% % % y23 = double(y23);
% % % x11 = double(x11);
% % % x22 = double(x22);
% % % x33 = double(x33);
% % % for i = 1:length(y12)
% % %     tempy12 = y12(i);
% % %     tempy13 = y13(i);
% % %     tempy23 = y23(i);
% % % % % %     tempx11 = x11(i);
% % % % % %     tempx22 = x22(i);
% % % % % %     tempx33 = x33(i);
% % %     if isreal(tempy12) && isreal(tempy13) && isreal(tempy23)
% % %         Y = [x11,tempy12,tempy13;tempy12,x22,tempy23;tempy13,tempy23,x33];
% % %         [V1,D1] = eig(Y);
% % %         if norm(Y,2) == norm(deltB_null,2)
% % %             break
% % %         end
% % %     end
% % % end
% % % try 
% % %     Y;
% % % catch
% % %     id = 1;
% % %     y12 = y12(id);
% % %     y13 = y13(id);
% % %     y23 = y23(id);
% % %     Y = [x11,y12,y13;y12,x22,y23;y13,y23,x33];
% % %     Y = real(Y);
% % %     [V1,D1] = eig(Y);
% % % end

% Solve inhomogeneous linear equations
% Solve curl field strength
% [k1,k2,k3,k4]=c_4_k(R1,R2,R3,R4);
% c_eval('k? = irf_resamp([tR? k?],[tB? B?]);')
% c_eval('K? = zeros(length(k?),4);',1:3)
% % c_eval('K?(:,1:4) = [k1(:,?+1),k2(:,?+1),k3(:,?+1),k4(:,?+1)];',1:3);
% c_eval('resR? = irf_resamp([tR? R?],[tB? B?]);')
% c_eval('resR? = resR?(:,2:4);')
% % % gradB = c_4_grad('R?','B?','grad');
% % % gradB = gradB(:,1:9);
% % % deltB_null=reshape(gradB,3,3);
% % % curlB = curl_Jac(gradB);
% % % if size(curlB,2) >= 4,c_eval('curlB = curlB(:,2:4);');end
% % % anticurlB = zeros(size(curlB,1),9);
% % % for i = 1:size(curlB,1)
% % %     anticurlB(i,:) = -[0,curlB(i,3)/2,curlB(i,2)/2,-curlB(i,3)/2,0,curlB(i,1)/2,-curlB(i,2)/2,-curlB(i,1)/2,0];
% % % end
% % % JacBr = anticurlB +gradB;
% % % JacB_null = reshape(JacBr,3,3);

% if ~isempty(idx_B) && ~isempty(idx_R)
% c_eval('B? = B?(idx_B,:);')
% c_eval('R? = R?(idx_R,:);')
% end

c_eval('Br? = zeros(size(B?));')
for i = 1:size(B1,1)
    deltB_null = reshape(gradB(1,:),3,3);
    Rcenter = [R1(i,:)+R2(i,:)+R3(i,:)+R4(i,:)]/4;
    Bcenter = transpose(deltB_null * (Rcenter-R1(i,:))') + B1(i,:);
% % %     deltB_Br = reshape(JacBr(1,:),3,3);
    c_eval("Br?(i,:) = transpose(Y * (R?(i,:)-Rcenter)');")
% % %     c_eval("Br?(i,:) = transpose(deltB_Br * (R?(i,:)-Rcenter)');")
%     c_eval("Br?(i,:) = transpose(deltB_Br * (R?(i,:)-Rcenter)');")
end
curlB = deltB_null-Y;
c_eval('Br? = [tB?(idx_B),Br?];')
end