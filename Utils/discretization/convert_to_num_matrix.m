function [data, data_need_convert_map] = convert_to_num_matrix(data)

% this function can make cell data convert to num matrix
%
% [data, data_need_convert_map] = convert_to_num_matrix(data)
%
% input:
%      data: cell data,from database.
% output:
%      data: the result, use this function processed.
%      data_need_convert_map: data need convert map，key is index, value is
%                             logical, 1 means need convert.
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个函数用来将cell数据转换为矩阵数据
%
% 输入:
%    data: cell数据，从数据库获取到的。
% 输出:
%    data: 用该函数处理得到的数据。
%    data_need_convert_map: 数据需要转换的map, key是位置索引，value是logical,
%                           1代表需要进行转换。
% 
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

% get data type map , the key is num , the value is string
[~, data_need_convert_map] = get_data_type_map(data);

% convert str to int in data

% this function convert data from char to int
for i = 1:length(keys(data_need_convert_map))
    
    if data_need_convert_map(i) == 1
        unique_data = unique(data(:,i));
        for j = 1: length(unique_data)
            for k = 1: size(data(:,i))
                if isequal(unique_data(j), data(k,i))
                    data{k,i} = j;
                end
            end
        end
    end
    
end

% all cell are double, so we use cell to mat
data = cell2mat(data);

end
