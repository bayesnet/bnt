clear
clc

% mysql default configure
mdl_data = get_default_mysql_data('iris');
class_index = 1;
node_sizes = get_node_sizes(mdl_data);
% make a dag
[mdl_data_row,mdl_data_col] = size(mdl_data);

% make a naive struct
dag = mk_naive_struct(length(node_sizes),class_index);

% 
% order = catty_order_algorithm(mdl_data,'class_node',1,'subtract_type','sum');
% dag = learn_struct_K2(mdl_data,node_size,order);
% draw_graph(dag);
% 
% % make a bnet
bnet = get_bnet(mdl_data,dag,'tabular');
% % set classified times
times = 10;
class_test_algorithm = 'CV-5';
[correct_rate, confusion_matrix] = times_classified_test(bnet,mdl_data,...
    class_index,class_test_algorithm,times);
% 
% [avg_correct_rate,avg_confusion_matrix] ...
%     = computer_avg_classified(times,correct_rate,confusion_matrix);



function [avg_correct_rate,avg_confusion_matrix] ...
    = computer_avg_classified(times,correct_rate,confusion_matrix)

computer_length = times * 5;
confusion_matrix_data ...
    = reshape(confusion_matrix,computer_length,1);
[confusion_matrix_data_row,confusion_matrix_data_col] ...
    = size(cell2mat(confusion_matrix_data(1)));
sum_confusion_matrix_mat ...
    = zeros(confusion_matrix_data_row,confusion_matrix_data_col);
for i = 1 : computer_length
    sum_confusion_matrix_mat = sum_confusion_matrix_mat +...
        cell2mat(confusion_matrix_data(i));
end
avg_confusion_matrix = sum_confusion_matrix_mat/computer_length;
avg_correct_rate = sum(correct_rate)/computer_length;

end
