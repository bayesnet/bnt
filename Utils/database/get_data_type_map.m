function [data_type_map, data_need_convert_map] = get_data_type_map(data)

% this function can get data type map and data need convert map. 
%
% input:
%      data: cell data,from database.
% output:
%      data_type_map: data type map
%      data_need_convert_map: data need convert map，key is index, value is
%                             logical, 1 means need convert.
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个函数用来获取数据类型map和数据需要转换的map。
%
% 输入:
%    data: cell数据，从数据库获取到的。
% 输出:
%    data_type_map: 数据类型map
%    data_need_convert_map: 数据需要转换的map, key是位置索引，value是logical,
%                           1代表需要进行转换。
%
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

% get data row and col
[data_row,data_col] = size(data);

% make a empty map
key_set{data_col,1} = zeros(data_col,1);
data_type_map_value_set{data_col,1} = zeros(data_col,1);
data_need_convert_flag_set{data_col,1} = zeros(data_col,1);

% data from mysql and get data types
for i = 1:data_col
    flag_count = 0;
    key_set{i,1} = i;
    for j = 1: data_row
        
        % judge num
        num_flag = isnumeric(data{j,i});
        if num_flag
            
            integer_flag = is_integer(data{j,i});
            if integer_flag
                flag_count = flag_count + 1;
                if data_row == flag_count
                    data_type_map_value_set{i,1} = {'int'};
                end
            else
                data_type_map_value_set{i,1} = {'decimal'};
            end
            
        end
        
        % judge char
        char_flag = ischar(data{j,i});
        if char_flag
            flag_count = flag_count + 1;
            if data_row == flag_count
                data_type_map_value_set{i,1} = {'char'};
            end
        end
        
    end
    
    % check whether data needs to be converted
    value_need_convert_flag = isequal(data_type_map_value_set{i,1}(1,1),{'char'});
    if value_need_convert_flag
        data_need_convert_flag_set{i,1} = 1;
    else
        data_need_convert_flag_set{i,1} = 0;
    end
    
end

data_type_map = containers.Map(key_set,data_type_map_value_set);
data_need_convert_map = containers.Map(key_set,data_need_convert_flag_set);

end


function flag = is_integer(x)

% this function check whether data is integer

flag = 0;

if x == fix(x)
    flag = 1;
end

return

end
