function Tpos = KH_event_spatial_distribution_20261003
% MMS磁层顶KH事件的GSE空间分布：XY、XZ平面。
% 只保留此前CDF核查中四星B、Vi均有有效记录的74个目录事件。
% 每个点为目录事件中点时刻的四星几何中心，不代表已识别的单个涡旋。
% 原CDF只读取星位和对应Epoch；不另存CSV、MAT或中间数据。

%% 1. 目录、IRFU函数和模型参数
CodeDir = 'C:\Users\Administrator\Documents\FWD_matlab\MMS_fu';
DataDir = 'Z:\SPART-WORK\Data\MMS';
FigDir = 'C:\Users\Administrator\Documents\KH';
IrfDir = 'C:\Users\Administrator\Documents\irfu-matlab-master';
AuditDir = fullfile(DataDir,'derived','KH','catalog_audit_20260903');
addpath(CodeDir,IrfDir,fullfile(IrfDir,'irf'), ...
    fullfile(IrfDir,'contrib','nasa_cdf_patch'));
setenv('CDF_LEAPSECONDSTABLE', ...
    fullfile(IrfDir,'contrib','nasa_cdf_patch','CDFLeapSeconds.txt'));

Units = irf_units;
RE = Units.R_Earth/1000;                 % 星位km转为地球半径
Dp = 2;                                % nPa，仅用于示意磁层顶
Bz = 0;                                % nT，IMF Bz (GSM)
[Xmp,Ymp] = irf_magnetosphere('mp_shue1998',Dp,Bz);

%% 2. 按已有核查记录筛选四星B+Vi事件
T = readtable(fullfile(AuditDir,'catalog_90_audit.csv'), ...
    'TextType','string','VariableNamingRule','preserve');
A = jsondecode(fileread(fullfile(AuditDir,'cdf_presence.json')));
E = A.entries;
EventID = string({E.EventID});
Product = string({E.Product});
Status = string({E.Status});
Sc = [E.Spacecraft];
Keep = false(height(T),1);
for ie = 1:height(T)
    iB = EventID==T.EventID(ie) & Product=="B" & Status=="有有效记录";
    iV = EventID==T.EventID(ie) & Product=="Vi" & Status=="有有效记录";
    Keep(ie) = isequal(sort(unique(Sc(iB))),1:4) && ...
        isequal(sort(unique(Sc(iV))),1:4);
end
T = T(Keep,:);
assert(height(T)==74,'当前输入的四星B+Vi事件数发生变化，请核对目录。');
Tsta = datetime(T.StartUTC,'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','UTC');
Tend = datetime(T.EndUTC,'InputFormat','yyyy-MM-dd HH:mm:ss','TimeZone','UTC');
Tmid = Tsta+(Tend-Tsta)/2;

%% 3. 从原CDF读取四星位置，对齐到各事件中点
% 优先MEC；没有覆盖中点的MEC时，使用FGM survey CDF附带的星位。
% 只在相邻有效星位之间调用irf_resamp，不外推、不跨数据缺口。
Rsc = nan(height(T),3,4);
CDFfiles = strings(height(T),4);
Products = strings(height(T),4);
PositionCadence = nan(height(T),4);
for ie = 1:height(T)
    for ic = 1:4
        [R,CDFfiles(ie,ic),Products(ie,ic),PositionCadence(ie,ic)] = ...
            read_position(DataDir,ic,Tmid(ie));
        Rsc(ie,:,ic) = R;
    end
    fprintf('%s  %s UTC  four positions OK\n',T.EventID(ie), ...
        char(datetime(Tmid(ie),'Format','yyyy-MM-dd HH:mm:ss')));
end
assert(all(isfinite(Rsc),'all'),'存在无法读取的四星星位，不能静默丢弃事件。');
Rc = mean(Rsc,3)/RE;
Separation = max(sqrt(sum((Rsc-mean(Rsc,3)).^2,2)),[],3);
Tpos = table(T.EventID,Tsta,Tend,Tmid,Rc(:,1),Rc(:,2),Rc(:,3), ...
    Separation,CDFfiles,Products,PositionCadence, ...
    'VariableNames',{'EventID','StartUTC','EndUTC','MidUTC', ...
    'X_GSE_RE','Y_GSE_RE','Z_GSE_RE','MaxDistanceToCentre_km', ...
    'SourceCDF','SourceProduct','PositionCadence_s'});
fprintf('POSITION_RANGE_RE min=[%.3f %.3f %.3f] max=[%.3f %.3f %.3f]\n', ...
    min(Rc,[],1),max(Rc,[],1));
fprintf('SOURCE_COUNTS: MEC=%d FGM_srvy=%d\n', ...
    nnz(Products=="MEC"),nnz(Products=="FGM survey"));

%% 4. XY、XZ投影：同一坐标尺度，晨昏侧使用相同颜色
N = height(T);
iDusk = Rc(:,2)>=0;
iDawn = ~iDusk;
DuskColor = [0.85 0.325 0.098];
DawnColor = [0 0.447 0.741];
XYmax = max(20,5*ceil((max(abs(Rc(:,2:3)),[],'all')+0.5)/5));
Xmax = max(15,5*ceil((max(Rc(:,1))+2)/5));
Xmin = min(-25,5*floor((min(Rc(:,1))-2)/5));
Fig = figure('Color','w','Position',[80 80 1500 1000], ...
    'Visible','off','Name','MMS KH event spatial distribution');
Ax = gobjects(1,2);
for ip = 1:2
    Ax(ip) = axes('Parent',Fig,'Position',[0.075+(ip-1)*0.49 0.19 0.40 0.66]);
    hold(Ax(ip),'on');
    hMP = plot(Ax(ip),[fliplr(Xmp) Xmp],[fliplr(Ymp) -Ymp], ...
        '--','Color',[0.35 0.35 0.35],'LineWidth',1.6);
    hMP.HandleVisibility = 'off';
    add_earth_irfu(Ax(ip));
    if ip==1
        Coordinate = Rc(:,2);
        YLabel = 'Y_{GSE} (R_E)';
        PanelTitle = '(a) XY plane';
    else
        Coordinate = Rc(:,3);
        YLabel = 'Z_{GSE} (R_E)';
        PanelTitle = '(b) XZ plane';
    end
    hDusk = scatter(Ax(ip),Rc(iDusk,1),Coordinate(iDusk),58,DuskColor, ...
        'filled','MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',0.45, ...
        'MarkerFaceAlpha',0.85);
    hDawn = scatter(Ax(ip),Rc(iDawn,1),Coordinate(iDawn),58,DawnColor, ...
        'filled','MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',0.45, ...
        'MarkerFaceAlpha',0.85);
    add_event_datatips(hDusk,Tpos(iDusk,:));
    add_event_datatips(hDawn,Tpos(iDawn,:));
    axis(Ax(ip),'equal');
    xlim(Ax(ip),[Xmin Xmax]);
    ylim(Ax(ip),[-XYmax XYmax]);
    set(Ax(ip),'XDir','reverse','FontName','Arial','FontSize',14, ...
        'LineWidth',1,'Box','on','TickDir','out','Layer','top', ...
        'XGrid','on','YGrid','on','GridAlpha',0.12);
    xlabel(Ax(ip),'X_{GSE} (R_E)','FontSize',16);
    ylabel(Ax(ip),YLabel,'FontSize',16);
    title(Ax(ip),PanelTitle,'FontSize',17,'FontWeight','normal');
    text(Ax(ip),0.035,0.95,'Sunward (+X)','Units','normalized', ...
        'FontSize',11,'Color',[0.3 0.3 0.3],'VerticalAlignment','top');
    if ip==1
        Lg = legend(Ax(ip),[hDusk hDawn], ...
            {sprintf('Dusk: %d events',nnz(iDusk)), ...
            sprintf('Dawn: %d events',nnz(iDawn))}, ...
            'Orientation','horizontal','Box','off','FontSize',13);
        Lg.Position = [0.285 0.09 0.43 0.04];
    end
end
annotation(Fig,'textbox',[0.04 0.91 0.92 0.06], ...
    'String',sprintf('MMS magnetopause KH events (N = %d)',N), ...
    'HorizontalAlignment','center','FontName','Arial','FontSize',22, ...
    'FontWeight','bold','EdgeColor','none');
annotation(Fig,'textbox',[0.04 0.865 0.92 0.04], ...
    'String','Four-spacecraft B + Vi available | Event-midpoint barycentres | GSE', ...
    'HorizontalAlignment','center','FontName','Arial','FontSize',13, ...
    'Color',[0.3 0.3 0.3],'EdgeColor','none');
annotation(Fig,'textbox',[0.04 0.035 0.92 0.035], ...
    'String','Dashed: reference magnetopause, Shue et al. (1998); P_{dyn} = 2 nPa, IMF B_z = 0 nT.', ...
    'Interpreter','tex','HorizontalAlignment','center','FontName','Arial', ...
    'FontSize',12,'Color',[0.25 0.25 0.25],'EdgeColor','none');

%% 5. 只保存图件；来源路径和事件坐标嵌入可编辑FIG
Fig.UserData = struct('Positions',Tpos,'EarthRadius_km',RE, ...
    'MagnetopauseModel','IRFU irf_magnetosphere: mp_shue1998', ...
    'Pressure_nPa',Dp,'IMFBz_GSM_nT',Bz, ...
    'PositionDefinition','Four-spacecraft barycentre at catalogue interval midpoint', ...
    'Selection','Four-spacecraft B and Vi CDF records confirmed; not continuous coverage', ...
    'CatalogFile',fullfile(AuditDir,'catalog_90_audit.csv'), ...
    'PresenceFile',fullfile(AuditDir,'cdf_presence.json'));
FileBase = fullfile(FigDir,'MMS_KH_74_events_GSE_XY_XZ_20261003');
drawnow;
exportgraphics(Fig,[FileBase '.png'],'Resolution',200,'BackgroundColor','white');
exportgraphics(Fig,[FileBase '.pdf'],'ContentType','vector','BackgroundColor','white');
savefig(Fig,[FileBase '.fig']);
fprintf('FINISHED: %d events; dusk=%d dawn=%d; all 4-spacecraft positions finite.\n', ...
    N,nnz(iDusk),nnz(iDawn));
fprintf('OUTPUT: %s\n',FileBase);
close(Fig);
end


function [R,File,Product,Cadence] = read_position(DataDir,ic,Time,AllowDownload)
% 同日多版本按数值版本排序，只读取实际覆盖目标时刻的CDF。
if nargin<4, AllowDownload=true; end
Sc = sprintf('mms%d',ic);
YYYY = char(datetime(Time,'Format','yyyy'));
MM = char(datetime(Time,'Format','MM'));
Day = char(datetime(Time,'Format','yyyyMMdd'));
Target = posixtime(Time);
Folder = {fullfile(DataDir,Sc,'mec','srvy','l2','epht89d',YYYY,MM), ...
    fullfile(DataDir,Sc,'fgm','srvy','l2',YYYY,MM)};
Pattern = {[Sc '_mec_srvy_l2_epht89d_' Day '*_v*.cdf'], ...
    [Sc '_fgm_srvy_l2_' Day '*_v*.cdf']};
Var = {[Sc '_mec_r_gse'],[Sc '_fgm_r_gse_srvy_l2']};
Epoch = {'Epoch','Epoch_state'};
ProductNames = ["MEC","FGM survey"];
Problems = strings(0,1);
for ip = 1:2
    F = dir(fullfile(Folder{ip},Pattern{ip}));
    Version = zeros(numel(F),3);
    for jf = 1:numel(F)
        Token = regexp(F(jf).name,'_v(\d+)\.(\d+)\.(\d+)\.cdf$','tokens','once');
        if ~isempty(Token), Version(jf,:) = str2double(Token); end
    end
    [~,Order] = sortrows(Version,[-1 -2 -3]);
    for jf = Order(:)'
        File = string(fullfile(F(jf).folder,F(jf).name));
        try
            Raw = spdfcdfread(char(File),'Variables',{Epoch{ip},Var{ip}}, ...
                'CombineRecords',true,'ConvertEpochToDatenum',true);
            t = (double(Raw{1}(:))-datenum(1970,1,1))*86400;
            r = double(Raw{2});
            if size(r,1)~=numel(t), r=r'; end
            assert(size(r,1)==numel(t) && size(r,2)>=3,'CDF星位维度异常。');
            r = r(:,1:3);
            [t,I] = unique(t,'sorted');
            r = r(I,:);
            Dt = diff(t);
            Cadence = median(Dt(isfinite(Dt) & Dt>0));
            i0 = find(t<=Target,1,'last');
            i1 = find(t>=Target,1,'first');
            if isempty(i0) || isempty(i1), continue; end
            if t(i1)-t(i0)>3*Cadence, continue; end
            if ~all(isfinite(r(i0:i1,:)),'all') || ...
                    any(abs(r(i0:i1,:))>1e7,'all'), continue; end
            if i0==i1
                R = r(i0,:);
            else
                Rtmp = irf_resamp([t(i0:i1) r(i0:i1,:)],Target,'linear');
                R = Rtmp(2:4);
            end
            assert(norm(R)>6000 && norm(R)<5e5,'星位长度单位或数值异常。');
            Product = ProductNames(ip);
            return
        catch ME
            Problems(end+1) = File+": "+string(ME.message); %#ok<AGROW>
        end
    end
end
if AllowDownload
    download_mec(DataDir,ic,Time);
    [R,File,Product,Cadence] = read_position(DataDir,ic,Time,false);
    return
end
error('无法读取%s在%s的原CDF星位。%s',Sc,char(Time),strjoin(Problems,newline));
end


function download_mec(DataDir,ic,Time)
% 复用用户已有SDCFilenames查询；websave下载原CDF并核验后按标准目录归档。
% 仅补齐事件位置所需星历，下载失败时不把残缺文件留作正式CDF。
Day = char(datetime(Time,'Format','yyyyMMdd'));
Date1 = char(datetime(Time,'Format','yyyy-MM-dd'));
Date2 = char(datetime(Time+days(1),'Format','yyyy-MM-dd'));
Names = SDCFilenames([Date1 '/' Date2],ic,'inst','mec','drm','srvy','dpt','epht89d');
Pattern = sprintf('^mms%d_mec_srvy_l2_epht89d_%s_v(\\d+)\\.(\\d+)\\.(\\d+)\\.cdf$',ic,Day);
Match = regexp(Names,Pattern,'tokens','once');
Keep = ~cellfun(@isempty,Match);
Names = Names(Keep);
Match = Match(Keep);
assert(~isempty(Names),'官方SDC未返回指定日期的MEC星历。');
Version = cell2mat(cellfun(@str2double,Match,'UniformOutput',false)');
[~,Order] = sortrows(Version,[-1 -2 -3]);
Name = Names{Order(1)};
Folder = fullfile(DataDir,sprintf('mms%d',ic),'mec','srvy','l2','epht89d',Day(1:4),Day(5:6));
if ~isfolder(Folder), mkdir(Folder); end
File = fullfile(Folder,Name);
Part = [File '.partial.cdf'];
URL = ['https://lasp.colorado.edu/mms/sdc/public/files/api/v1/download/science?file=' Name];
fprintf('DOWNLOAD_MEC %s\n',Name);
try
    websave(Part,URL,weboptions('Timeout',60));
    Info = spdfcdfinfo(Part,'Variables',{'Epoch',sprintf('mms%d_mec_r_gse',ic)});
    assert(size(Info.Variables,1)==2 && all(cell2mat(Info.Variables(:,3))>1), ...
        '下载的CDF缺少星位或Epoch。');
    movefile(Part,File,'f');
catch ME
    if isfile(Part), delete(Part); end
    rethrow(ME);
end
end


function add_earth_irfu(h)
% 复用IRFU mms.mms4_pl_conf内部add_Earth的terminator画法。
% 该辅助函数不是独立公开接口，故将原有画法用于本图的两个axes。
theta = 0:pi/20:pi;
xEarth = sin(theta);
yEarth = cos(theta);
patch(-xEarth,yEarth,'k','EdgeColor','none','Parent',h,'HandleVisibility','off');
patch(xEarth,yEarth,'w','EdgeColor','k','Parent',h,'HandleVisibility','off');
end


function add_event_datatips(h,T)
% MATLAB打开FIG后，点击点即可查看事件编号、UTC和三维星位。
h.DataTipTemplate.DataTipRows = [ ...
    dataTipTextRow('Event',cellstr(T.EventID)); ...
    dataTipTextRow('Midpoint UTC',cellstr(string(T.MidUTC))); ...
    dataTipTextRow('X GSE (RE)',T.X_GSE_RE); ...
    dataTipTextRow('Y GSE (RE)',T.Y_GSE_RE); ...
    dataTipTextRow('Z GSE (RE)',T.Z_GSE_RE)];
end
