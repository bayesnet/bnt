clear
clc

% mysql default configure
datasource = 'growlithe';
username ='growlithe';
password = 'growlithe';
url = 'jdbc:mysql://localhost:3306/growlithe';

% make sql
select_var = 'id, sepal_length, sepal_width, petal_length, petal_width, class_name, status';
sql = ['SELECT ', select_var ,' FROM uci_database.iris where status = 1'];

% set class index
id_index = 1;
% if id not null, set id index
class_index = 6;
% if clear var are not null, set clear index
clear_index = 7;

% get data from mysql
data = get_data_from_mysql(sql,datasource,username,password,url);

if isequal(0,data)
    disp('check you sql,the sql get null data')
    return
end

% convert data to num matrix
[data, data_need_convert_map] = convert_to_num_matrix(data);

% use mdl algorithm get data which bnt can used
mdl_data = mdl_algorithm(data,class_index,data_need_convert_map,'id_index',...
     id_index, 'clear_index',clear_index);

