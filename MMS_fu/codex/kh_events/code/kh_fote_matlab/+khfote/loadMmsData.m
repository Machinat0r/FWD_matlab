function raw = loadMmsData(tint,velocityField)
%KH FOTE 数据读取模块：读取MMS四星原生时间序列。
%
% 模块职责
%   1. 读取四颗卫星的FGM burst磁场B_GSE，单位 nT；
%   2. 读取FPI burst离子或电子速度V_GSE，单位 km/s；
%   3. 读取与速度种类一致的数密度N，单位 cm^-3；
%   4. 读取四星GSE位置R，单位 km。
%
% 输入
%   tint           IRFU时间区间
%   velocityField  'Vi'或'Ve'
%
% 输出
%   raw.B/V/N/R 均为4x1 cell。每个cell第一列为Unix秒，后续列为数据。
%   本模块保留每种产品的原生采样率，不执行插值或平滑。

raw = struct;
raw.B = cell(4,1);
raw.V = cell(4,1);
raw.N = cell(4,1);
raw.R = cell(4,1);

%% 1. 读取磁场、速度和密度
for sc = 1:4
    bts = mms.get_data('B_gse_brst_l2',tint,sc);
    if strcmpi(velocityField,'Ve')
        vts = mms.get_data('Ve_gse_fpi_brst_l2',tint,sc);
        nts = mms.get_data('Ne_fpi_brst_l2',tint,sc);
    else
        vts = mms.get_data('Vi_gse_fpi_brst_l2',tint,sc);
        nts = mms.get_data('Ni_fpi_brst_l2',tint,sc);
    end
    raw.B{sc} = tseriesToMatrix(bts,3,sprintf('MMS%d B',sc));
    raw.V{sc} = tseriesToMatrix(vts,3, ...
        sprintf('MMS%d %s',sc,velocityField));
    raw.N{sc} = tseriesToMatrix(nts,1,sprintf('MMS%d density',sc));
end

%% 2. 读取卫星位置
for sc = 1:4
    % burst FGM CDF内含与磁场配套的位置，可覆盖MEC/survey缺文件情况。
    dataset = sprintf('mms%d_fgm_brst_l2',sc);
    variable = sprintf('mms%d_fgm_r_gse_brst_l2',sc);
    rts = mms.db_get_ts(dataset,variable,tint);
    if isempty(rts)
        % 对仅保存MEC或survey星位的数据库使用兼容回退。
        rts = mms.get_data('R_gse',tint,sc);
    end
    raw.R{sc} = tseriesToMatrix(rts,3,sprintf('MMS%d R',sc));
end
end


function out = tseriesToMatrix(ts,nComponents,label)
% 将IRFU TSeries转换为普通double矩阵，便于后续模块独立处理。
if isempty(ts) || isempty(ts.time)
    error('%s returned no data.',label);
end
tEpoch = ts.time.epochUnix;
t = double(tEpoch(:));
values = squeeze(double(ts.data));
if isvector(values), values = values(:); end
if size(values,1)~=numel(t) && size(values,2)==numel(t)
    values = values.';
end
n = min(numel(t),size(values,1));
if n==0 || size(values,2)<nComponents
    error('%s has an unexpected data shape.',label);
end
out = [t(1:n),values(1:n,1:nComponents)];
end
