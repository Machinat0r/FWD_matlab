

%% FOTE 误差
gradB=c_4_grad('R?','B?','grad');
eigVal_err_v2=B1(:,1);

for ii=1:length(B1(:,1))  
deltB_null=reshape(gradB(ii,2:end),3,3);
[V,D] = eig(deltB_null);
eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;   %% Figure 1o    以最大特征值归一化
end
    

[j,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);          %% Figure 1n    受背景电流影响大
err_4C(:,2)=abs(err_4C(:,2))*100;    

%% 以上两个实际上都不是很合理，建议按照如下计算
gradB=c_4_grad('R?','B?','grad');
c_eval('B?=irf_abs(B?);',1:4);
B_mean=mean([B1(:,5) B2(:,5) B3(:,5) B4(:,5)],2);

divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)])];      %% 未归一化散度
divB(:,2)=divB(:,2)./B_mean;                                     %% 归一化散度



