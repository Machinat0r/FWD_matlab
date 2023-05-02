clc
clear

filename = '/Users/fwd/Documents/MATLAB/新建文件夹/gzz/JED_090_HIERSESP_CDR_2019307_V04.TAB';
filename_var = '/Users/fwd/Documents/MATLAB/新建文件夹/gzz/JED_HIERSESP_CDR_V02.FMT';

data = readmatrix(filename, 'FileType', 'text','TrimNonNumeric',false,'OutputType','string');

% data_head = readmatrix(filename_var, delimitedTextImportOptions('DataLines',[1,Inf],'FileType', 'text', 'TrimNonNumeric',false,'OutputType','string'));
data_head = readmatrix(filename_var, 'FileType', 'text', 'TrimNonNumeric',false,'OutputType','string');
opts = detectImportOptions(filename_var, 'FileType', 'text');

data_head0 = readmatrix(filename, delimitedTextImportOptions('DataLines',[1,5]));

AA = opts.DataLines;
if AA(1)~=1
%     data_head_tmp = readmatrix(filename_var, delimitedTextImportOptions('DataLines',[1,Inf]));
    siz1 = size(data_head);
    data_head(AA(1):siz1(1)+AA(1)-1,:) = data_head;
    data_head(1:2,:) = data_head(10:11,:);
    data_head(2,3) = 'UTC';
end

idx1 = find(data_head(:,1) == 'OBJECT');
idx2 = find(data_head(:,1) == 'END_OBJECT');
for ii = 1:length(idx1)
    for ii0 = idx1(ii)+1:idx2(ii+1)-1
        v_tmp = data_head(ii0,1);
        switch v_tmp
            case 'NAME'
                Variable_names(ii,2) = data_head(ii0,3);
            case 'COLUMN_NUMBER'
                Variable_names(ii,1) = data_head(ii0,3);
            case 'UNIT'
                Variable_names(ii,3) = data_head(ii0,3);
            case 'DESCRIPTION'
                Variable_names(ii,4) = data_head(ii0,3);
        end
    end
end

Data_JED = cellstr(Variable_names);