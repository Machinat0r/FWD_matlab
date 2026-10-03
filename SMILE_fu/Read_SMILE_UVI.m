clear;clc

%% input
% 文件地址：可以填写一个 CDF 的完整路径，也可以填写存放 CDF 的文件夹。
FilePath = 'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS';
% FilePath = 'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS\SMILE_UVI_L2_AURORA-GEO_20260720T230011(1).cdf';

FileType = 'AURORA-GEO'; % UVI 产品类型；依据文件内 Logical_source 核对
ic = 1;                % 从读取列表中选第几个文件，提取下方常用变量
IrfDir = 'C:\Users\Administrator\Documents\irfu-matlab-master';

% 输出：UVI{1}, UVI{2}, ... 保存各文件全部原始变量；CDFInfo 保存元数据。
% Data = UVI{ic}；Image、ImageGEO、ImageAACGM 等为所选文件的常用变量。
% 保留原数值类型、时间、填充值及数组排列，不平滑、插值、转置或重新定标。
% 当前已用两份 AURORA-GEO CDF 实测；其他产品仅按其原变量通用读取。

%% file list
addpath(fullfile(IrfDir,'contrib','nasa_cdf_patch'));
FilePath = char(FilePath);
FileType = upper(strtrim(char(FileType)));
assert(~isempty(FileType),'请填写 UVI 产品类型，例如 AURORA-GEO。');
if isfolder(FilePath)
    % 目录读取：按标准文件名筛选，不递归搜索；列表按文件名排序。
    Files = dir(fullfile(FilePath,['SMILE_UVI_*_' FileType '_*.cdf']));
    Files = Files(~[Files.isdir]);
    [~,Order] = sort({Files.name});
    Files = Files(Order);
elseif isfile(FilePath)
    [~,~,Ext] = fileparts(FilePath);
    assert(strcmpi(Ext,'.cdf'),'本程序读取 CDF 文件。');
    Files = dir(FilePath);
else
    error('SMILE_UVI:PathNotFound','文件或目录不存在：%s',FilePath);
end
assert(~isempty(Files),'指定目录中没有找到 %s 产品的 CDF 文件。',FileType);
assert(isscalar(ic) && isfinite(ic) && ic==fix(ic) && ic>=1 && ic<=numel(Files), ...
    'ic 必须为 1 到 %d 之间的整数。',numel(Files));
FileList = fullfile({Files.folder},{Files.name});
disp(string(FileList(:)));

%% load data
% 使用 MATLAB cdfinfo 读取元数据，避开本地 spdfcdfinfo.Variables
% 第10--12列已发现的摘要异常；FILLVAL 等属性见 CDFInfo.VariableAttributes。
% 使用 IRFU 自带 NASA spdfcdfread 逐变量读取，保留各变量自身的记录数。
UVI = cell(1,numel(Files));
CDFInfo = cell(1,numel(Files));
TimeUTC = cell(1,numel(Files));
TempDir = fullfile(tempdir,'SMILE_UVI_read');
if ~exist(TempDir,'dir'), mkdir(TempDir); end
OriginalDir = pwd;
cd(TempDir); % CDF 接口可能创建临时文件，统一留在 TEMP
try
    for iFile = 1:numel(Files)
        CDFInfo{iFile} = cdfinfo(FileList{iFile});
        Info = CDFInfo{iFile};

        % 核对真实产品标识；Logical_source 中的空项及 DOI 不参与匹配。
        assert(isfield(Info.GlobalAttributes,'Logical_source'), ...
            '文件缺少 Logical_source，无法核对产品类型：%s',Files(iFile).name);
        Source = Info.GlobalAttributes.Logical_source;
        Source = Source(cellfun(@ischar,Source));
        IsProduct = startsWith(Source,'SMILE_UVI_','IgnoreCase',true) & ...
            endsWith(Source,['_',FileType],'IgnoreCase',true);
        assert(any(IsProduct), ...
            '文件声明的产品与 FileType=%s 不符：%s',FileType,Files(iFile).name);

        for iVar = 1:size(Info.Variables,1)
            VarName = Info.Variables{iVar,1};
            if Info.Variables{iVar,3} == 0
                UVI{iFile}.(VarName) = []; % 原文件中该变量没有记录
                continue
            end
            Value = spdfcdfread(FileList{iFile},'Variables',{VarName}, ...
                'KeepEpochAsIs',true,'CombineRecords',true);
            if iscell(Value), Value=Value{1}; end
            UVI{iFile}.(VarName) = Value;
        end
        fprintf('已读取 %s：%d 个变量。\n',Files(iFile).name,size(Info.Variables,1));
    end
catch ME
    cd(OriginalDir);
    rethrow(ME);
end
cd(OriginalDir);

%% load time
% 原始时间仍保存在 UVI 中。这里只另外生成 UTC 字符串，不覆盖时间数值。
% 根据 CDF 数据类型识别 TT2000，沿用本机 CDF 闰秒表；不作星上时钟校正。
% 原时间的纳秒存储位数不代表已经验证了纳秒授时精度。
for iFile = 1:numel(Files)
    TimeUTC{iFile} = struct;
    Info = CDFInfo{iFile};
    for iVar = 1:size(Info.Variables,1)
        VarName = Info.Variables{iVar,1};
        if ~strcmpi(Info.Variables{iVar,4},'tt2000'), continue; end
        TT = UVI{iFile}.(VarName);
        if isempty(TT), continue; end
        assert(isa(TT,'int64'),'TT2000 必须保留 int64 类型。');
        TimeUTC{iFile}.(VarName) = string(spdfencodett2000(TT,'Format',4));
    end
end

%% select file
Data = UVI{ic};
Info = CDFInfo{ic};
UTC = TimeUTC{ic};
fprintf('\n当前文件：%s\n',FileList{ic});
disp(Info.Variables(:,1:4)); % 变量名、CDF维度、记录数、数据类型
% 属性示例：Info.VariableAttributes.UNITS / FILLVAL / CATDESC。
% UTC 中仅解码 TT2000；若产品使用其他时间格式，原时间仍完整保存在 Data。

%% load AURORA-GEO
% 三类图像及坐标直接引用文件变量；只对已核对结构的 AURORA-GEO 提取。
% 其他类型的数据直接查看 Data 或 UVI{iFile}，不套用此产品的变量名。
if strcmp(FileType,'AURORA-GEO')
    % load epoch
    Epoch = Data.EPOCH;                            % int64 TT2000
    EpochUTC = UTC.EPOCH;                          % UTC 字符串
    ExpoStart = Data.SMILE_UVI_EXPO_START;          % 原始曝光开始时间
    ExpoEnd = Data.SMILE_UVI_EXPO_END;              % 原始曝光结束时间
    Cadence = Data.SMILE_UVI_CADENCE;               % s，文件原值

    % load corrected image
    Image = Data.SMILE_UVI_IMAGE;                  % Rayleigh
    ProjGeoLat = Data.SMILE_UVI_PROJ_GEO_LAT;       % deg
    ProjGeoLon = Data.SMILE_UVI_PROJ_GEO_LON;       % deg
    ProjGeoAlt = Data.SMILE_UVI_PROJ_GEO_ALT;       % km，文件给定的投影高度
    ProjFlag = Data.SMILE_UVI_PROJ_FLAG;            % 文件注释：1-Earth, 0-Space

    % load GEO image
    ImageGEO = Data.SMILE_UVI_GRID_IMG_GEO;         % Rayleigh
    GeoLat = Data.SMILE_UVI_GRID_GEO_LAT;           % deg
    GeoLon = Data.SMILE_UVI_GRID_GEO_LON;           % deg

    % load AACGM image
    ImageAACGM = Data.SMILE_UVI_GRID_IMG_AACGM;     % Rayleigh
    AACGMLat = Data.SMILE_UVI_GRID_AACGM_LAT;        % degree
    AACGMMLT = Data.SMILE_UVI_GRID_AACGM_MLT;        % hour
    disp(EpochUTC);
end

%% notes
% 1. Image等保留读取器返回的行列与记录维度；不预设图像上方为北方或正午。
% 2. 本程序完全保留FILLVAL、0与-999。它们的具体含义须结合文件属性和标志。
% 3. 本批SC_POS的描述与坐标属性分别写GSE/J2000，SC_ATT定义也有冲突；
%    原值留在Data中，不执行位置坐标转换或姿态解算。
% 4. CDF已给出校正图像和投影网格；本程序不重复上游校正、投影或标定。
% 5. 不自动合并不同文件的时间/图像，不生成全量MAT/CSV，不修改原CDF。
% 6. NASA接口会将本批模式字符串末尾的NUL显示为空格；文本内容ITEN保持。