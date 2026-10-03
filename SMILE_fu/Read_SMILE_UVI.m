%% input
if ~exist('FilePath','var')
    FilePath = 'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS';
end
if ~exist('FileType','var')
    FileType = 'AURORA-GEO';
end
if ~exist('ic','var')
    ic = 1;
end
if ~exist('IrfDir','var')
    IrfDir = 'C:\Users\Administrator\Documents\irfu-matlab-master';
end

%% file list
addpath(fullfile(IrfDir,'contrib','nasa_cdf_patch'));
FilePath = char(FilePath);
FileType = upper(strtrim(char(FileType)));
assert(~isempty(FileType),'请填写 UVI 产品类型，例如 AURORA-GEO。');
if isfolder(FilePath)
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
UVI = cell(1,numel(Files));
CDFInfo = cell(1,numel(Files));
TimeUTC = cell(1,numel(Files));
TempDir = fullfile(tempdir,'SMILE_UVI_read');
if ~exist(TempDir,'dir'), mkdir(TempDir); end
OriginalDir = pwd;
cd(TempDir);
try
    for iFile = 1:numel(Files)
        CDFInfo{iFile} = cdfinfo(FileList{iFile});
        Info = CDFInfo{iFile};
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
                UVI{iFile}.(VarName) = [];
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
disp(Info.Variables(:,1:4));

%% load AURORA-GEO
clear Epoch EpochUTC ExpoStart ExpoEnd Cadence Image ProjGeoLat ProjGeoLon ...
    ProjGeoAlt ProjFlag ImageGEO GeoLat GeoLon ImageAACGM AACGMLat AACGMMLT
if strcmp(FileType,'AURORA-GEO')
    Epoch = Data.EPOCH;
    EpochUTC = UTC.EPOCH;
    ExpoStart = Data.SMILE_UVI_EXPO_START;
    ExpoEnd = Data.SMILE_UVI_EXPO_END;
    Cadence = Data.SMILE_UVI_CADENCE;

    Image = Data.SMILE_UVI_IMAGE;
    ProjGeoLat = Data.SMILE_UVI_PROJ_GEO_LAT;
    ProjGeoLon = Data.SMILE_UVI_PROJ_GEO_LON;
    ProjGeoAlt = Data.SMILE_UVI_PROJ_GEO_ALT;
    ProjFlag = Data.SMILE_UVI_PROJ_FLAG;

    ImageGEO = Data.SMILE_UVI_GRID_IMG_GEO;
    GeoLat = Data.SMILE_UVI_GRID_GEO_LAT;
    GeoLon = Data.SMILE_UVI_GRID_GEO_LON;

    ImageAACGM = Data.SMILE_UVI_GRID_IMG_AACGM;
    AACGMLat = Data.SMILE_UVI_GRID_AACGM_LAT;
    AACGMMLT = Data.SMILE_UVI_GRID_AACGM_MLT;
    disp(EpochUTC);
end