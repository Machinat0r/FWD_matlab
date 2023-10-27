function data_out = Read2BData(data_path, label_path)
%------written by Wending Fu, Aug.2023 in Beijing------------
% input 1 nargin for data path, and the label path is defaulted
% or
% input 2 nargins for data path and label path
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       南无电子阿弥陀佛驱散仿生bug
%                                _ooOoo_
%                               o8888888o
%                               88" . "88
%                               (| -_- |)
%                               O\  =  /O
%                            ____/`---'\____
%                          .'  \\|     |//  `.
%                         /  \\|||  :  |||//  \
%                        /  _||||| -:- |||||-  \
%                        |   | \\\  -  /// |   |
%                        | \_|  ''\-/''  |   |
%                        \  .-\__  `-`  ___/-. /
%                      ___`. .'  /-.-\  `. . __
%                   ."" '<  `.___\_<|>_/___.'  >'"".
%                  | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%                  \  \ `-.   \_ __\ /__ _/   .-` /  /
%             ======`-.____`-.___\_____/___.-`____.-'======
% 	                   `=-='
%                 天地玄宗，万气本根。广修亿劫，证吾神通。
%                 三界内外，惟道独尊。体有金光，覆映吾身。
%                 视之不见，听之不闻。包罗天地，养育群生。
%                 受持万遍，身有光明。三界侍卫，五帝司迎。
%                 万神朝礼，役使雷霆。鬼妖丧胆，精怪忘形。
%                 内有霹雳，雷神隐名。洞慧交彻，五炁腾腾。
%                金光速现，覆护真人。急急如律令，bug全去除！
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% file path check
narginchk(1,2)
% data_path = '/Users/fwd/Documents/MATLAB/Code/fwd/TianWen/MEPA/HP_B/HX1-Or_GRAS_MEPA-HP_SCI_N_20201129000002_20201129235956_00130_B.2B';
if nargin == 1
label_path = '/Users/fwd/Documents/MATLAB/Code/fwd_matlab_patch/TianWen/MEPA/HP_B/HX1-Or_GRAS_MEPA-HP_SCI_N_20201129000002_20201129235956_00130_B.2BL';
warning('no label path input, using the default label path（ˉ﹃ˉ）')
end
data = importdata(data_path);
%% read label
label_cell = importdata(label_path); 
label_str = cat(2, label_cell{:});
label_groups= extractBetween(label_str,'<Group_Field_Character>','</Group_Field_Character>');
label_fields = extractBetween(label_str,'<Field_Character>','</Field_Character>');
label_struct = cell(length(label_fields), 1);
data_cell = cell(label_fields);
for i = 1:length(label_fields)
    field = label_fields{i};
    label_struct{i} = struct();
    repetitions = 1;
    if sum(contains(label_groups, field))
        label = find(contains(label_groups, label_fields{i}) == 1);
        label_struct{i}.group_location = str2double(extractBetween(label_groups{label},'<group_location unit="byte">', '</group_location>'));
        label_struct{i}.group_length = str2double(extractBetween(label_groups{label},'<group_length unit="byte">', '</group_length>'));
        label_struct{i}.repetitions = str2double(extractBetween(label_groups{label},'<repetitions>', '</repetitions>'));
        repetitions = label_struct{i}.repetitions;
    end 
    
    label_struct{i}.name = char(extractBetween(field, '<name>', '</name>'));
    label_struct{i}.name = replace(label_struct{i}.name,'-','_'); label_struct{i}.name = replace(label_struct{i}.name,'/','_');
    label_struct{i}.field_number = char(extractBetween(field, '<field_number>', '</field_number>'));
    label_struct{i}.field_location = str2double(extractBetween(field,'<field_location unit="byte">', '</field_location>'));
    label_struct{i}.data_type = char(extractBetween(field, '<data_type>', '</data_type>'));
    label_struct{i}.field_length = str2double(extractBetween(field, '<field_length unit="byte">', '</field_length>'));
    label_struct{i}.field_format = char(extractBetween(field, '<field_format>', '</field_format>'));
    label_struct{i}.unit = char(extractBetween(field, '<unit>', '</unit>'));
    label_struct{i}.description = char(extractBetween(field, '<description>', '</description>'));

    data_cell{i} =  zeros(length(data), repetitions);
end
%% read data
group_location = zeros(length(label_struct),1); repetitions = ones(length(label_struct),1); location_sep = zeros(length(label_struct), 1);
for label = 1:length(label_struct)
    if isfield(label_struct{label}, 'group_location')
        group_location(label) = label_struct{label}.group_location - 1;
        repetitions(label) = label_struct{label}.repetitions;
        location_sep(label) = label_struct{label}.group_length/label_struct{label}.repetitions;
    end
end

for i = 1:length(data)
    for label = 1:length(label_struct)
        for repetition = 1:repetitions(label)
            temp = data{i}(label_struct{label}.field_location + group_location(label) + location_sep(label)*(repetition-1):...
                label_struct{label}.field_location + group_location(label) + label_struct{label}.field_length-1 + location_sep(label)*(repetition-1));
            if label_struct{label}.name == "UTC", temp = num2str(EpochTT(temp).epoch);
            elseif contains(temp,'NUL'), temp = 'NaN';end

            data_cell{label}(i,repetition) = str2double(temp);
            % eval([label_struct{label}.name, '_data(i,repetition) = ', num2str(temp), ';']);
        end
    end
    clc;disp(['读取数据：',repmat('■',1,round(10*i/length(data))),repmat('□',1,10-round(10*i/length(data))),' (ง •̀_•́)ง'])
end
%% data output
data_out = cell(length(label_struct),1);
for label = 1:length(label_struct)
    data_out{label} = label_struct{label};
    data_out{label}.data = data_cell{label};
end