function res = hs_plot(varargin)
%PLOT([H], dobj, var_s, [comp], [options])  plot a variable
%
% OPTIONS - one of the following:
%	'AX'       - axis handles to use
%   'COMP'     - components to plot
%   'SUM_DIM1' - sum over first dimension (frequency, pitch angle)
%   'COMP_DIM1' - form subplots from that component
%   'nolabels' - only plot data, do not add any labels or text
%   'ColorbarLabel' - specify colorbar label
%   'FitColorbarLabel' - fit the text of colorbar label to colobar height
%   'ClusterColors' - use Cluster colors C1-black, C2-red, C3-green, C4-blue
%
%
% for common cluster variables see: http://spreadsheets.google.com/pub?key=0AjU2FPHMIluIdEcydzhHeXhRUzRiM3lsN1VzQWMySmc&hl=en&single=true&gid=0&output=html
% $Id: plot.m,v 1.45 2011/02/09 11:23:31 yuri Exp $

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,14,nargin))

[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax),
    ax=gca;
end
dobj=args{1};
var_s=args{2};
args=args(3:end);
nargs=nargs-2;

LCOMP = 3;

if ~ischar(var_s), error('VAR_S must be a stirng'), end

data = getv(dobj,var_s);
if isempty(data), error('VAR_S not found'), end
dim = length(data.variance(3:end));
dep = getdep(dobj,var_s);
units = corr_latex(getunits(dobj,var_s));
fieldnam = findva(dobj,'FIELDNAM',var_s);
ii = regexp(fieldnam,'_'); fieldnam(ii) = ' '; % Get rid of underscores
lablaxis = getlablaxis(dobj,var_s);
ii = regexp(lablaxis,'_'); lablaxis(ii) = ' '; % Get rid of underscores
cs = getcs(dobj,var_s);
if ~isempty(cs), cs = [' ' cs]; end
fillv = getfillval(dobj,var_s);
data.data(data.data==fillv) = NaN;


%% INPUT ARGUMENTS

% Default values that can be override by options
sum_dim = 0;  % along which dimension to sum
comp_dim = 0; % component dimension used to separate subplots 
use_comp = 0; % pick up separate component values in separate plots
comp = [];    % index of component vector values to pick up 
ydim = 0;     % default dimension of data used for y axis (0- data value itself);
plot_properties=cell(0);
create_axes = 1;
flag_lineplot = 0;
flag_spectrogram = 0;
flag_labels_is_on=1;
flag_colorbar_label_is_manually_specified=0;
flag_colorbar_label_fit_to_colorbar_height_is_on=0;
line_color=''; % default line color; can be changed with flags, e.g. clustercolors
flag_use_cluster_colors=0;

arg_pos = 0;
while ~isempty(args)
    arg_pos = arg_pos + 1;
    l = 1;
    if arg_pos==1 && isnumeric(args{1})
        use_comp = 1;
        comp = args{1};
    else
        switch(lower(args{1}))
            case 'ax'
                l = 2;
                if all(ishandle(args{2}))
                    ax = args{2};
                    create_axes = 0;
                else disp('invalid value for AX')
                end
            case 'comp'
                l = 2;
                if isnumeric(args{2})
                    use_comp = 1;
                    comp = args{2};
                else
                    disp('invalid value for COMP')
                end
            case 'clustercolors'
                flag_use_cluster_colors = 1;
                if findstr('C1',var_s), line_color='k';
                elseif findstr('C2',var_s), line_color='r';
                elseif findstr('C3',var_s), line_color='g';
                elseif findstr('C4',var_s), line_color='b';
                else flag_use_cluster_colors = 0;
                end
            case 'nolabels'
                flag_labels_is_on = 0;
            case 'colorbarlabel'
                l=2;
                if ischar(args{2})
                    flag_colorbar_label_is_manually_specified=1;
                    colorbar_label = args{2};
                else
                    disp('invalid value for ColorbarLabel in PLOT')
                end
            case 'fitcolorbarlabel'
                flag_colorbar_label_fit_to_colorbar_height_is_on=1;
            case 'sum_dim1'
                sum_dim = 1;
            case 'sum_dim2'
                sum_dim = 2;
            case 'sum_dim3'
                sum_dim = 3;
            case 'comp_dim1'
                comp_dim = 1;
            case 'comp_dim2'
                comp_dim = 2;
            case 'comp_dim3'
                comp_dim = 3;
            otherwise
                disp('unknown argument')
                disp('the rest or arguments are passed to plot routines');
                plot_properties=args; 
                args=cell(0);
                break
        end
    end
    args = args(l+1:end);
end

if comp_dim ~=0 && comp_dim == sum_dim
    error('SUM_DIM and COMP_DIM must be different')
end

%% DATA PROCESSING

    
% define summing dimension and component to plot when component not defined
if comp_dim == 0,
    if dim == 1, comp_dim = 1;
    elseif dim == 2 && sum_dim == 0,    comp_dim = 2; ydim = 1;
    elseif dim == 2 && sum_dim == 1,    comp_dim = 2; ydim = 0;
    elseif dim == 2 && sum_dim == 2,    comp_dim = 1;
    elseif dim == 3 && sum_dim == 0,    comp_dim = 3; sum_dim = 1;
    elseif dim == 3 && sum_dim == 1,    comp_dim = 3;
    elseif dim == 3 && sum_dim == 2,    comp_dim = 3;
    elseif dim == 3 && sum_dim == 3,    comp_dim = 2; ydim = 1;     
    end
end

if dim ==3,  
    switch comp_dim
        case 1
            if sum_dim == 2, ydim = 3; else ydim = 2; end
        case 2
            if sum_dim == 1, ydim = 3; else ydim = 1; end
        case 3
            if sum_dim == 1, ydim = 2; else ydim = 1; end
    end
end


if dim == 0 
    plot_data = {double(data.data)};
    flag_lineplot = 1;

elseif dim == 1   
    if use_comp,
        plot_data = cell(size(comp));
        for i=1:length(comp)
            plot_data{i} = double(data.data(:,comp(i)));
        end
        flag_lineplot = 1;
    else
        plot_data = {double(data.data)};
        if dim == 1
            if isfield(data,'DEPEND_1')
                flag_spectrogram = 1;
                ydim = 1;
            else
                flag_lineplot = 1;
            end
        end
    end
    
elseif dim == 2
    if sum_dim > 0
        % Do not count NaNs when summing
        tt = data.data; tt(~isnan(tt)) = 1; tt(isnan(tt)) = 0;
        data.data(isnan(data.data)) = 0;
        data.data = sum(data.data,sum_dim+1,'double')./...
            sum(tt,sum_dim+1,'double'); % The data is 2D from now on
        clear tt
        
        if sum_dim == 1, comp_dim = 2;
        else  comp_dim = 1;
        end
    end 
    
    if sum_dim > 0 && isfield(data,'DEPEND_1') && ~use_comp
        plot_data = {squeeze(data.data)};
        ydim = comp_dim;
    else
        ndim = data.dim(comp_dim);
        if ~use_comp, comp=1:ndim; end        
            plot_data = cell(size(comp));
            for i=1:length(comp)
                switch comp_dim
                    case 1
                        plot_data{i} = squeeze(data.data(:,comp(i),:));
                    case 2
                        plot_data{i} = squeeze(data.data(:,:,comp(i)));
                    otherwise
                        error('smth wrong')
                end
            end
            if comp_dim == 2, ydim =1; else ydim = 2; end
    end
    
    plot_f = 2;
    if ~isfield(data,'DEPEND_1'), plot_f = plot_f-1; end
    if sum_dim > 0, plot_f = plot_f-1; end
    if use_comp, plot_f = plot_f-1; end
    if plot_f > 0
        flag_spectrogram = 1;
    else
        flag_lineplot = 1;
    end
    
elseif dim == 3    
    if use_comp == 0 && sum_dim == 0
        if comp_dim==2
            sum_dim = 1;
        else
            sum_dim = 2;
        end
    end
    
    % Do not count NaNs when summing
    tt = data.data; tt(~isnan(tt)) = 1; tt(isnan(tt)) = 0;
    data.data(isnan(data.data)) = 0;
    data.data = sum(data.data,sum_dim+1,'double')./...
        sum(tt,sum_dim+1,'double'); % The data is 2D from now on
    clear tt
    
    ndim = data.dim(comp_dim);
    if ~use_comp, comp=1:ndim; end
    if ndim == 1
        plot_data = {data.data};
    else
        plot_data = cell(size(comp));
        for i=1:length(comp)
            switch comp_dim
                case 1
                    plot_data{i} = squeeze(data.data(:,comp(i),:,:));
                case 2
                    plot_data{i} = squeeze(data.data(:,:,comp(i),:));
                case 3
                    plot_data{i} = squeeze(data.data(:,:,:,comp(i)));
                otherwise
                    error('smth wrong')
            end
        end
    end    
    flag_spectrogram = 1;
else
    error('plotting not implememnted')
end


if flag_lineplot
    %% PLOTTING -- LINE PLOT
    if isfield(dep,'DEPEND_O')
        h = irf_plot(ax,[dep.DEPEND_O plot_data{:}],line_color,plot_properties{:});
    else
        h = plot(ax,data.data,line_color,plot_properties{:});
    end
    flab = getlablaxis(dobj,var_s);
    lab_1 = ''; 
    if ~isempty(dep.DEPEND_X) % 
        dep_x_s = dep.DEPEND_X{comp_dim,1};
        dep_x = getv(dobj,dep_x_s);
        if ~isempty(dep_x)
            if strcmp(dep_x.type,'char')
                if use_comp, % pick up components, data should be char
                    lab_1 = ['(' dep_x.data(1,:,comp) ')'];
                else % data are values, label under LABLAXIS
                    legend(num2cell(dep_x.data(1,:,:),2), 'Location','NorthWest')
                    legend('boxoff')
                end
            else
                lab_1 = ['(' num2str(dep_x.data(1,comp),'%6.2f') dep_x.UNITS ')'];
            end
        end
    end
    ylabel(sprintf('%s%s [%s]', flab, lab_1, units))
    
    text_s = [dobj.GlobalAttributes.OBSERVATORY{1} ' > ' ...
        dobj.GlobalAttributes.INSTRUMENT_NAME{1} ' > ' fieldnam];
    if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
    add_text(h,text_s);
    
elseif flag_spectrogram 
    %% PLOT -- SPECTRPGRAM
    
    dep_x=cell(size(dep.DEPEND_X,1));
    for d = 1:length(dep_x)
        dep_x{d} = getv(dobj,dep.DEPEND_X{d,1});
        dep_x{d}.s = dep.DEPEND_X{d,1};
        dep_x{d}.fillv = getfillval(dobj,dep_x{d}.s);
        if ~strcmp(dep_x{d}.type,'char')
            dep_x{d}.data(dep_x{d}.data==dep_x{d}.fillv) = NaN;
        end
        dep_x{d}.units = getunits(dobj,dep_x{d}.s);
        dep_x{d}.lab = getlablaxis(dobj,dep_x{d}.s);
        % check if DELTA_PLUS and  DELTA_MINUS are given
        if isfield(dep_x{d},'DELTA_PLUS') && isfield(dep_x{d},'DELTA_MINUS')
            dep_x{d}.df=struct('plus',dep_x{d}.DELTA_PLUS,'minus',dep_x{d}.DELTA_MINUS);
            if ischar(dep_x{d}.DELTA_PLUS)
                deltaplus= getv(dobj,dep_x{d}.DELTA_PLUS);
                deltaminus= getv(dobj,dep_x{d}.DELTA_MINUS);
                dep_x{d}.df.plus=deltaplus.data(1,:);
                dep_x{d}.df.minus=deltaminus.data(1,:);
            end
        else dep_x{d}.df=[];
        end
        
    end
    
    if sum_dim > 0
        fprintf('Summing over dimension %d (%s)\n', ...
            sum_dim, dep_x{sum_dim}.lab)
    end
    
    specrec = struct('t',dep.DEPEND_O,'f',dep_x{1}.data(1,:),'f_unit',dep_x{1}.units,'p',[],'df',dep_x{1}.df);
    
    lab_2 ='';
    if length(dep_x)>1 && ~isempty(dep_x{comp_dim})
        if strcmp(dep_x{comp_dim}.type,'char') && strcmp(dep_x{comp_dim}.variance,'F/T')...
                && findstr(dep_x{comp_dim}.s,'LABEL_2')
%            reclen = size(dep_x{comp_dim}.data,2)/length(dep.DEPEND_O);
            lab_2 = shiftdim(dep_x{comp_dim}.data(1,:,:),1)';
%            lab_2 = dep_x{comp_dim}.data(:,1:reclen);
        elseif strcmp(dep_x{comp_dim}.type,'single') && ...
                (strcmp(dep_x{comp_dim}.variance,'F/T') || ...
                strcmp(dep_x{comp_dim}.variance,'T/T'))
            lab_2 = num2str(dep_x{comp_dim}.data(1,comp)',['%.2f ' dep_x{comp_dim}.units '\n']);
        else
            error('BAD type for DEPEND_X')
        end
    end
    
    text_s = ['Cluster'];
    if ~isempty(cs), text_s = [text_s ' [' shorten_cs(cs) ']']; end
    
    if isempty(comp), comp = 1; end
    ncomp = length(comp);
    h = zeros(1,ncomp);
    if create_axes, ax = zeros(1, ncomp); end
    
    
    if ydim > 1
        specrec.f = dep_x{ydim}.data(1,:);
        specrec.f_unit = dep_x{ydim}.units;
        specrec.df = dep_x{ydim}.df;
    end
    
    % special case for degrees
    ytick = [];
    if strcmpi(dep_x{ydim}.units,'degrees') || strcmpi(dep_x{ydim}.units,'deg')
        frange = abs(max(specrec.f)-min(specrec.f));
        if frange > 80 && frange <=150, da = 15;
        elseif frange > 150 && frange <=200, da = 45;
        elseif frange > 200 && frange <=380, da = 90;
        else da = [];
        end
        if ~isempty(da)
            ytick = round(min(specrec.f)/da):round(max(specrec.f)/da);
            ytick = ytick*da;
        end
    end
    
    for i=1:ncomp
        specrec.p = plot_data(i);
        if create_axes, ax(i) = irf_subplot(length(comp),1,-i); end %#ok<AGROW>
        h(i) = caa_spectrogram(ax(i),specrec);
        if ~isempty(ytick), set(ax(i) ,'YTick',ytick), end
        if flag_labels_is_on,
            if ncomp<=LCOMP % Small number of components
                ylabel(sprintf('%s [%s]', dep_x{ydim}.lab, dep_x{ydim}.units))
                if ~isempty(lab_2), lab_2s = [text_s ' > ' lab_2(i,:)];
                else lab_2s = text_s;
                end
                add_text(h(i),lab_2s);
            else % Large number of components
                if i==1, title(text_s), end
                if i==fix(ncomp/2)+1, ylabel(sprintf('%s [%s]', dep_x{ydim}.lab, dep_x{ydim}.units))
                else ylabel('')
                end
                add_text(h(i),lab_2(i,:));
            end
        end
    end
    % Add colorbar
    i=fix(ncomp/2)+1;
    hcb = colorbar('peer',h(i));
    hcbl = get(hcb,'ylabel');
    dy = get(ax(i),'Position'); dy = dy(3);
    pcb = get(hcb,'Position');
    if ncomp>1, set(hcb,'Position',[pcb(1) pcb(2)-pcb(4)*(ncomp-fix(ncomp/2)-1) pcb(3) pcb(4)*ncomp]); end
    if flag_labels_is_on,
        if ~flag_colorbar_label_is_manually_specified
            colorbar_label=['Log ' lablaxis ' [' units ']' ];
        end
        ylabel(hcb,colorbar_label);
        if flag_colorbar_label_fit_to_colorbar_height_is_on
            fit_colorbarlabel_height(hcbl);
        end
    else
        ylabel(hcb,'');
    end
    
    %huishan add
    value=plot_data(i);
    Vlow=min(value{1}(:));  Vup=max(value{1}(:)); 
    caxis([log10(Vlow) log10(Vup)]); 
    
    % Resize all panels after addition of the colorbar
    if ~isempty(dy)
        for i=1:ncomp
            tt = get(ax(i),'Position');
            set(ax(i),'Position',[tt(1) tt(2) dy tt(4)])
        end
    end
    set(ax(1:ncomp-1),'XTickLabel',[]);
    for i=1:1:ncomp-1, xlabel(ax(i),'');end
    
end

if nargout > 0, res = h; end

%huishan add
res=plot_data;

function add_text(h,txt)
text(0.99, 0.97, [' ' txt],'HorizontalAlignment','right','VerticalAlignment','top',...
    'units','normalized','fontsize',6,'parent',h)

function cs = shorten_cs(cs)

if isempty(cs), return, end

% Remove leading spaces
while cs(1) == ' ', cs(1) = []; end

if strcmpi(cs(1:3),'GSE'), cs = 'GSE'; end

% Try to correct latex
function s = corr_latex(s)
expr = {'\^-[1-3]','\^[2-3]'};
exprl = [2 1];
for i=1:length(expr)
    while 1
        ii = regexp(s,expr{i});
        if isempty(ii), break, end
        ii = ii(1);
        l = length(s);
        s_tmp = [s(1:ii) '{' s(ii+1:ii+exprl(i)) '}'];
        if l > ii+2, s = [s_tmp s(ii+exprl(i)+1:end)];
        else s = s_tmp;
        end
    end
end

function fit_colorbarlabel_height(hcb)
hy=hcb;
colorbar_label_fontsize=get(hy,'fontsize');
units=get(hy,'units');
set(hy,'units','normalized');
temp=get(hy,'Extent');
colorbarlabelheight = temp(4);
while colorbarlabelheight>1.1,
    colorbar_label_fontsize=colorbar_label_fontsize*0.95;
    %set(hy,'fontsize',colorbar_label_fontsize,'position',labelposition);
    set(hy,'fontsize',colorbar_label_fontsize);
    temp=get(hy,'Extent');
    colorbarlabelheight=temp(4);
end

