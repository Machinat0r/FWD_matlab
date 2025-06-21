function [R] = get_R_gse2gsm(epochUnix)
% 传入 ISDAT epoch 或普通 unix epoch 时间戳
% 返回 GSE->GSM 的 3x3 旋转矩阵 R

    % 三个基矢，在 GSE 下分别是 (1,0,0)、(0,1,0)、(0,0,1)
    ex_gse = [epochUnix 1 0 0];
    ey_gse = [epochUnix 0 1 0];
    ez_gse = [epochUnix 0 0 1];

    ex_gsm = irf.geocentric_coordinate_transformation(ex_gse,'gse>gsm');
    ey_gsm = irf.geocentric_coordinate_transformation(ey_gse,'gse>gsm');
    ez_gsm = irf.geocentric_coordinate_transformation(ez_gse,'gse>gsm');

    % 取 xyz 分量（去掉时间列），按列拼成 R
    R = [ex_gsm(2:4)' , ey_gsm(2:4)' , ez_gsm(2:4)'];
end
