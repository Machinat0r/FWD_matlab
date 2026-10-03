%% SMILE UVI：读取两份 CDF，核对原始数组，预览文件内三种图像
% 输入：用户提供的两个原始 CDF。输出：工作区 UVI、UVIInfo、UVIAttr、TimeUTC。
% UVI 保留 spdfcdfread 返回的数值类型、数组及 int64 TT2000 时间。
% 使用本地 IRFU 自带 NASA CDF 接口；独立调用 MATLAB cdflib 逐项核对。
% 结果 PNG/TXT 保存到项目结果目录；临时文件在 TEMP。原 CDF 不改写。
% 本脚本针对已检查的每文件一帧结构；不平滑、不插值、不重新标定。

%% 路径
InputDir = 'C:\Users\Administrator\Documents\Recovery-Work_SMILE-MMS';
OutputDir = fullfile(InputDir,'SMILE_UVI_readout_20261003');
TempDir = fullfile(tempdir,'SMILE_UVI_CDF_inspect_20261003');
addpath('C:\Users\Administrator\Documents\irfu-matlab-master\contrib\nasa_cdf_patch');
FileNames = {'SMILE_UVI_L2_AURORA-GEO_20260720T230011(1).cdf', ...
             'SMILE_UVI_L2_AURORA-GEO_20260720T230041(1).cdf'};
if ~exist(OutputDir,'dir'), mkdir(OutputDir); end
if ~exist(TempDir,'dir'), mkdir(TempDir); end
OriginalDir = pwd;
cd(TempDir); % NASA CDF 读取器可能在工作目录创建临时文件
UVI = cell(1,2); UVIInfo = cell(1,2); UVIAttr = cell(1,2);
TimeUTC = cell(2,3); % 列依次为 EPOCH、曝光开始、曝光结束
Summary = cell(2,1);

%% 读取与独立核对；属性直接使用其实际 CDF 数据类型
% 本机 spdfcdfinfo.Variables 的第10--12列对部分 double 变量出现错误值。
% 因此 FILLVAL/VALIDMIN/MAX 等属性使用 cdflib.getAttrEntry 读取；
% 并逐项与 spdfcdfinfo.VariableAttributes 核对，不使用上述摘要列。
for iFile = 1:2
    FilePath = fullfile(InputDir,FileNames{iFile});
    Info = spdfcdfinfo(FilePath,'Validate',true);
    assert(all(cell2mat(Info.Variables(:,3))==1),'此脚本要求每变量一记录。');
    UVIInfo{iFile} = Info;
    CdfId = cdflib.open(FilePath);
    try
        cdflib.setReadOnlyMode(CdfId,'READONLYon');
        assert(strcmp(cdflib.getMajority(CdfId),'ROW_MAJOR'));
        Names = Info.Variables(:,1);
        AttrNames = fieldnames(Info.VariableAttributes);
        Lines = strings(numel(Names),1);
        for iVar = 1:numel(Names)
            Name = Names{iVar};
            A = spdfcdfread(FilePath,'Variables',{Name}, ...
                'KeepEpochAsIs',true,'CombineRecords',true);
            if iscell(A), A = A{1}; end
            UVI{iFile}.(Name) = A;
            VarId = cdflib.getVarNum(CdfId,Name);
            B = cdflib.getVarRecordData(CdfId,VarId,0);
            if isnumeric(A)
                % 两个现成接口对本文件 ROW_MAJOR 数组的排列相反。
                % 仅将独立验证数组转置；UVI 中保留 NASA 接口的原返回数组。
                assert(isequaln(A,B.'),'两接口数值不一致：%s',Name);
            else
                % CDF 字符串末尾 NUL 被 NASA 接口显示为空格；仅比较文本。
                assert(strcmp(strtrim(A),strtrim(strrep(B.',char(0),''))));
            end
            Attr = struct;
            for iAttr = 1:numel(AttrNames)
                AttrName = AttrNames{iAttr};
                Entries = Info.VariableAttributes.(AttrName);
                Row = find(strcmp(Entries(:,1),Name));
                if isempty(Row), continue; end
                AttrId = cdflib.getAttrNum(CdfId,AttrName);
                Value = cdflib.getAttrEntry(CdfId,AttrId,VarId);
                assert(isequaln(Value,Entries{Row,2}), ...
                    '属性核对失败：%s / %s',Name,AttrName);
                Attr.(AttrName) = Value;
            end
            UVIAttr{iFile}.(Name) = Attr;
            Units = ''; if isfield(Attr,'UNITS'), Units=Attr.UNITS; end
            if isnumeric(A)
                FillMask = false(size(A));
                if isfield(Attr,'FILLVAL') && ~strcmp(Name,'SMILE_UVI_PROJ_FLAG')
                    FillMask = (A == cast(Attr.FILLVAL,'like',A));
                end
                % FLAG 的0同时定义为空间与FILLVAL，保留0/1全部原值。
                Good = isfinite(A) & ~FillMask;
                if isa(A,'int64')
                    RangeText = "TT2000 int64 retained";
                else
                    RangeText = sprintf('min=%.12g max=%.12g',min(A(Good)),max(A(Good)));
                end
                Lines(iVar) = sprintf('%s | %s | %s | %s | fill=%d | %s', ...
                    Name,mat2str(size(A)),class(A),Units,nnz(FillMask),RangeText);
            else
                Lines(iVar) = sprintf('%s | %s | char | value=%s',Name,mat2str(size(A)),strtrim(A));
            end
        end
        cdflib.close(CdfId);
    catch ME
        cdflib.close(CdfId);
        cd(OriginalDir);
        rethrow(ME);
    end
    Summary{iFile} = Lines;
end

%% 时间：直接用现成 TT2000 接口，不先转换为 double 或 datenum
TimeNames = {'EPOCH','SMILE_UVI_EXPO_START','SMILE_UVI_EXPO_END'};
ExposureSeconds = zeros(2,1);
for iFile=1:2
    for iTime=1:3
        TT = UVI{iFile}.(TimeNames{iTime});
        TimeUTC{iFile,iTime} = char(string(spdfencodett2000(TT,'Format',4)));
        PartsNASA=spdfbreakdowntt2000(TT); PartsMATLAB=cdflib.breakdownTT2000(TT);
        assert(isequal(PartsNASA(:),PartsMATLAB(:))); 
    end
    S=UVI{iFile}.SMILE_UVI_EXPO_START; E=UVI{iFile}.SMILE_UVI_EXPO_END;
    assert(UVI{iFile}.EPOCH-S == E-UVI{iFile}.EPOCH);
    ExposureSeconds(iFile)=double(E-S)*1e-9;
    assert(isequal(UVI{iFile}.SMILE_UVI_PROJ_GEO_LAT==-999,UVI{iFile}.SMILE_UVI_PROJ_FLAG==0));
    assert(isequal(UVI{iFile}.SMILE_UVI_PROJ_GEO_LON==-999,UVI{iFile}.SMILE_UVI_PROJ_FLAG==0));
end
FrameSeparationSeconds = double(UVI{2}.EPOCH-UVI{1}.EPOCH)*1e-9;

%% 绘图副本：仅精确匹配各变量的 FILLVAL，原始 UVI 数组保持不变
ImageNames={'SMILE_UVI_IMAGE','SMILE_UVI_GRID_IMG_GEO','SMILE_UVI_GRID_IMG_AACGM'};
ImageTitles={'Corrected image','GEO grid image','AACGM grid image'};
ImagePlot=cell(2,3); Limits=zeros(3,2);
for iImage=1:3
    Both=[];
    for iFile=1:2
        Name=ImageNames{iImage}; Raw=UVI{iFile}.(Name);
        Missing=(Raw==cast(UVIAttr{iFile}.(Name).FILLVAL,'like',Raw));
        ImagePlot{iFile,iImage}=double(Raw);
        ImagePlot{iFile,iImage}(Missing)=NaN;
        Both=[Both; ImagePlot{iFile,iImage}(~Missing)]; %#ok<AGROW>
    end
    Limits(iImage,:)=[min(Both),max(Both)];
end

%% 原数组行列预览：行向下增加，不解释上下左右的地理含义
Fig=figure('Color','w','Position',[50 50 1550 1020],'Visible','off');
Layout=tiledlayout(Fig,2,3,'TileSpacing','compact','Padding','compact');
for iFile=1:2
    for iImage=1:3
        Ax=nexttile(Layout); C=ImagePlot{iFile,iImage};
        imagesc(Ax,C,'AlphaData',isfinite(C));
        set(Ax,'YDir','reverse','Color',[0.85 0.85 0.85],'FontSize',11);
        axis(Ax,'image'); colormap(Ax,parula(256)); clim(Ax,Limits(iImage,:));
        xlabel(Ax,'Array column'); ylabel(Ax,'Array row');
        title(Ax,{sprintf('%s (%d x %d)',ImageTitles{iImage},size(C,1),size(C,2)), ...
            ['EPOCH ' TimeUTC{iFile,1}(12:23) ' UTC']},'FontWeight','normal');
        CB=colorbar(Ax); CB.Label.String='Rayleigh';
    end
end
sgtitle(Layout,{'SMILE UVI | 20 July 2026', ...
    'Native array preview | rows increase downward | grey = declared FILLVAL | common linear scale within each column'}, ...
    'FontSize',15,'FontWeight','normal');
exportgraphics(Fig,fullfile(OutputDir,'SMILE_UVI_two_frames_array_preview.png'),'Resolution',150);
close(Fig);

%% 保存读取摘要与明确限制；不导出全量 MAT/CSV
Fid=fopen(fullfile(OutputDir,'SMILE_UVI_readout.txt'),'w','n','UTF-8');
fprintf(Fid,'SMILE UVI CDF MATLAB 读取结果\nMATLAB %s\n\n',version);
fprintf(Fid,'两个文件各有23个变量、每变量1记录。所有数值数组逐元素核对通过。\n');
fprintf(Fid,'NASA/IRFU spdfcdfread 与 MATLAB cdflib.getVarRecordData 对本文件存在转置关系；读取结果保留前者排列。\n');
fprintf(Fid,'全部已有变量属性经两个接口逐项核对。字符模式为 ITEN，原生尾部NUL在NASA接口中显示为空格。\n\n');
for iFile=1:2
    fprintf(Fid,'文件%d：%s\n',iFile,FileNames{iFile});
    fprintf(Fid,'EPOCH：%s\n曝光开始：%s\n曝光结束：%s\n曝光长度：%.9f s\n', ...
        TimeUTC{iFile,1},TimeUTC{iFile,2},TimeUTC{iFile,3},ExposureSeconds(iFile));
    fprintf(Fid,'原CADENCE：%.9f s；投影高度：%.9f km；CALIB_COEF：%.9f\n', ...
        UVI{iFile}.SMILE_UVI_CADENCE,UVI{iFile}.SMILE_UVI_PROJ_GEO_ALT,UVI{iFile}.SMILE_UVI_CALIB_COEF);
    fprintf(Fid,'原文件ROI：%.12g Rayleigh（ROI范围定义未给出）\n',UVI{iFile}.SMILE_UVI_ROI);
    fprintf(Fid,'%s\n',Summary{iFile}); fprintf(Fid,'\n');
end
fprintf(Fid,'两帧EPOCH间隔：%.9f s。每帧EPOCH经整数运算核对恰好为曝光中点。\n\n',FrameSeparationSeconds);
fprintf(Fid,'假设、显示约定和待确认事项\n');
fprintf(Fid,'1. 时间按CDF的TT2000类型及TIME_SCALE=Terrestrial Time解码；NASA与MATLAB解码结果一致。本机内置闰秒表与文件标记的最后更新均为2017-01-01。这里沿用文件时间基准和该表，未另作星上时钟校正，也未验证绝对授时精度。\n');
fprintf(Fid,'2. Rayleigh、GEO、AACGM和110 km投影高度沿用文件属性与变量值；它们是文件给出的信息。未重新计算投影、曝光归一化或辐射定标，也未独立验证既有L2处理。\n');
fprintf(Fid,'3. 预览仅按NASA读取后的数组行列绘制，行号向下增加；未赋予上方北向或正午含义，未旋转或重投影。每列两帧共用原数据全范围线性色标；灰色仅表示精确FILLVAL。零值全部保留。\n');
fprintf(Fid,'4. FILLVAL在本文件多为single类型的-1e31，实际值为-9.9999998482432073e30；比较时按真实属性类型转换到变量类型。仅绘图副本替换为NaN，未应用额外VALIDMIN/MAX筛选。\n');
fprintf(Fid,'5. PROJ_FLAG属性明示1-Earth,0-Space，同时FILLVAL也为0，存在含义重叠；保留原0/1值。PROJ_GEO_LAT/LON的-999与FLAG=0逐像素对应，未将这些-999当作真实经纬度。\n');
fprintf(Fid,'6. SC_POS的CATDESC/FIELDNAM写GSE，COORDINATE_SYSTEM写J2000；SC_ATT描述写GSE，坐标属性写Spacecraft Coordinate System，且数据是4分量、单位写deg。坐标与姿态含义待产品方确认，未据此做转换或解算。\n');
fprintf(Fid,'7. 全局Descriptor写MLAT-MLT-IMG，Logical_source/文件名写AURORA-GEO，而文件同时包含GEO和AACGM产品。保留各自变量级标签。\n');
fprintf(Fid,'8. KEOGRAM_NS标注沿23 MLT，EW标注沿66度，未提供各采样点的完整坐标轴；ROI区域定义和ITEN模式含义也未给出，未自行补定。\n');
fprintf(Fid,'9. 本机spdfcdfinfo.Variables的FILLVAL/VALIDMIN/MAX摘要列在部分double变量上返回异常值；使用经MATLAB cdflib核对的VariableAttributes真实属性，未修改IRFU安装文件。\n');
fprintf(Fid,'10. 两文件均无额外平滑、插值、背景扣除、质量门槛或人为填补。预览里的已有网格图像由原CDF直接读取；上游生成它们时采用的具体处理未由本次读取验证。\n');
fclose(Fid);
cd(OriginalDir);
fprintf('READOUT_COMPLETE\nEPOCH1=%s\nEPOCH2=%s\nFRAME_SEPARATION=%.9f s\nOUTPUT=%s\n', ...
    TimeUTC{1,1},TimeUTC{2,1},FrameSeparationSeconds,OutputDir);