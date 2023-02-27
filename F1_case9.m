

%% FOTE ���
gradB=c_4_grad('R?','B?','grad');
eigVal_err_v2=B1(:,1);

for ii=1:length(B1(:,1))  
deltB_null=reshape(gradB(ii,2:end),3,3);
[V,D] = eig(deltB_null);
eigVal_err_v2(ii,2)=abs(D(1,1)+D(2,2)+D(3,3))/max([abs(D(1,1)), abs(D(2,2)), abs(D(3,3))]) * 100;   %% Figure 1o    ���������ֵ��һ��
end
    

[j,divB,~,jxB,divTshear,divPb] = c_4_j('R?','B?');
temp=irf_abs(j);
jmag=temp(:,[1 5]);
err_4C=irf_multiply(1,divB,1,jmag,-1);          %% Figure 1n    �ܱ�������Ӱ���
err_4C(:,2)=abs(err_4C(:,2))*100;    

%% ��������ʵ���϶����Ǻܺ������鰴�����¼���
gradB=c_4_grad('R?','B?','grad');
c_eval('B?=irf_abs(B?);',1:4);
B_mean=mean([B1(:,5) B2(:,5) B3(:,5) B4(:,5)],2);

divB=[gradB(:,1) sum([gradB(:,2) gradB(:,6) gradB(:,10)])];      %% δ��һ��ɢ��
divB(:,2)=divB(:,2)./B_mean;                                     %% ��һ��ɢ��



