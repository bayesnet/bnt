function [avg_correct_rate,avg_confusion_matrix] ...
    = computer_avg_classified(correct_rate,confusion_matrix)

% this function computer avg correct rate from correct rate and computer 
% avg confusion matrix from confusion matrix.
%
% input:
%      correct_rate: the correct rate from classified test.
%      confusion_matrix: the confusion matrix from classified test.
% output:
%      avg_correct_rate: the avg correct rate from correct rate.
%      avg_confusion_matrix: the avg confusion matrix from confusion matrix.
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个函数用来计算平均正确率 和 平均含混矩阵。
%
% 输入:
%    correct_rate: 从分类检验中获取的正确率。
%    confusion_matrix: 从分类检验中获取的含混矩阵。
% 输出:
%    avg_correct_rate: 从正确率中得到的平均正确率。
%    avg_confusion_matrix: 从含混矩阵中得到的平均含混矩阵。
% 
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

[confusion_matrix_row,confusion_matrix_col] = size(confusion_matrix);
confusion_matrix_length = confusion_matrix_row * confusion_matrix_col;
confusion_matrix_data = reshape(confusion_matrix,confusion_matrix_length,1);
[confusion_matrix_data_row,confusion_matrix_data_col] ...
    = size(cell2mat(confusion_matrix_data(1)));
sum_confusion_matrix_mat = zeros(confusion_matrix_data_row,confusion_matrix_data_col);

for i = 1 : confusion_matrix_length
    sum_confusion_matrix_mat = sum_confusion_matrix_mat + cell2mat(confusion_matrix_data(i));
end

avg_correct_rate = sum(correct_rate)/length(correct_rate);
avg_confusion_matrix = sum_confusion_matrix_mat/confusion_matrix_length;

end
