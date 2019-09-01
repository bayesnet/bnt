function mdl_data = mdl_algorithm(data,class_index,data_need_convert_map,varargin)

% this function achieve Fayyad & Irani's MDL method and support data list
% from database.
%
% mdl_data = mdl_algorithm(data,class_index,data_need_convert_map,varargin)
%
% input:
%    data: a matrix data，just support num
%    class_index: class attribute index
%    data_need_convert_map: a map show data whether need convert int
%    varargin: support 'id_index' ,representative id index,'clear_index', .
%              representative clear attribute index ,The mdl algorithm will
%              delete unrelated attribute at the end.
% output:
%    mdl_data: the result of mdl algorithm process data
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个函数实现了Fayyad & Irani 的mdl算法，并且支持从数据库中获取的数据。
% 输入:
%    data: 一个矩阵数据,只支持数字类型
%    class_index: 类属性的索引位置
%    data_need_convert_map: 数据是否需要转换为int类型的map
%    varargin: 只支持'id_index'，代表id的索引位置,'clear_index',需要清楚变量的索引
%              位置。mdl算法会在最后删除无关的属性。
% 输出:
%    mdl_data: mdl 算法处理数据的结果
%
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)


id_index = [];
clear_index = [];
fix_data = 'no';
args = varargin;
nargs = length(args);
if ~isempty(args)
    if ischar(args{1})
        for i=1:2:nargs
            switch args{i}
                case 'id_index'
                    id_index = args{i+1};
                case 'clear_index'
                    clear_index = args{i+1};
                case 'fix_data'
                    fix_data = args{i+1};
            end
        end
    end
end

[row,col] = size(data);

if ~isempty(id_index)
    if length(id_index) ~= 1
        error('id_index must one length')
    end
    if id_index > col
        error('id_index must small than data col size')
    end
end
if ~isempty(clear_index)
    if clear_index > col
        error('clear_index must small than data col size')
    end
end
if class_index > col
    error('class_index must small than data col size')
end

if ~isempty(id_index)
    data_need_convert_map(id_index) = 1;
    clear_index = [id_index,clear_index];
end
if ~isempty(clear_index)
    for i = 1:length(clear_index)
        data_need_convert_map(clear_index(i)) = 1;
    end
end

data_need_convert_map(class_index) = 1;

mdl_data = zeros(row,col);

for i = 1: col
    % determine data whether need mdl process
    data_need_mdl_flag = data_need_convert_map(i) == 0;
    if ~data_need_mdl_flag
        mdl_data(:,i) = data(:,i);
    end
    
    if data_need_mdl_flag
        [single_mdl_data,~] = single_mdl_algorithm(data(:,i)',data(:,class_index)');
        mdl_data(:,i) = single_mdl_data';
    end
    
end

% if id index is not empty, delete unrelated attribute
if ~isempty(clear_index)
    mdl_data(:,clear_index) = [];
end

if isequal('yes',fix_data)
    
    need_clear_row = [];
    for i=1:row
        if isequal(any(isnan(mdl_data(i,:))),1)
            need_clear_row = [need_clear_row,i];
        end
    end
    if ~isempty(need_clear_row)
        mdl_data(need_clear_row,:) = [];
    end
    
    for i = 1:col
        flag = ~isempty(find(mdl_data(:,i) == 0, 1));
        if isequal(flag,1)
            mdl_data(:,i) = mdl_data(:,i) + 1;
        end
        
        unqiue_i_num = sort(unique(mdl_data(:,i)));
        if length(unqiue_i_num)<max(mdl_data(:,i))
            for j= 1:length(unqiue_i_num)
                [index,~] = find(mdl_data(:,i)==unqiue_i_num(j,:));
                mdl_data(index,i) = j;
            end
        end
    end
    
end

end
