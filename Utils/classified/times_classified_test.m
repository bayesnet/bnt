function [correct_rate_list,confusion_matrix_list] = times_classified_test(bnet,...
    data,class_index,cvpartition_algorithm_name,times)

% this function .
%
% input:
%      bnet: .
%      data: .
%      class_index: .
%      cvpartition_algorithm_name: .
%      times: .
% output:
%      correct_rate_list:  .
%      confusion_matrix_list:  .
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个函数。
%
% 输入:
%      bnet: .
%      data: .
%      class_index: .
%      cvpartition_algorithm_name: .
%      times: .
% 输出:
%    avg_correct_rate: 。
%    avg_confusion_matrix: 。
% 
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

correct_rate_list = zeros(1,times);
for i = 1: times
    [correct_rate,confusion_matrix] = classified_test(bnet,data,class_index,...
        cvpartition_algorithm_name);
    correct_rate_list(i) = correct_rate;
    if i == 1
        [~,col] = size(confusion_matrix);
        confusion_matrix_list = num2cell(zeros(times,col));
    end
    confusion_matrix_list(i,:) = confusion_matrix;
end

end

function [correct_rate,confusion_matrix_list] = classified_test(bnet,data,class_index,...
    cvpartition_algorithm_name)

class_data = data(class_index,:);
cv = cvpartition_package(class_data,cvpartition_algorithm_name);
unique_class_data_length = length(unique(class_data));
confusion_matrix_list = cell(1,cv.NumTestSets);
for i = 1:cv.NumTestSets
    train_data = data(:,cv.training(i));
    train_bnet = learn_params(bnet,train_data);
   
    test_data = data(:,cv.test(i));
    test_class_data = test_data(class_index,:);
    
    prob = bn_inference(train_bnet,unique_class_data_length,class_index,test_data);
    
    [~,max_prob_class] = max(prob);
    
    [correct_rate,confusion_matrix] = computer_correct(test_class_data,max_prob_class);
    confusion_matrix_list{i} = confusion_matrix;
end

end

function [correct_rate,confusion_matrix] = computer_correct(test_class,max_prob_class)

test_class_length = length(test_class);
class_unique_num = unique(test_class);
class_unique_num_length = length(class_unique_num);
confusion_matrix = zeros(class_unique_num_length,class_unique_num_length);
correct = zeros(1,length(test_class));

for i = 1:test_class_length
    correct(i) = isequal(max_prob_class(i),test_class(i));
    confusion_matrix_i = test_class(i);
    confusion_matrix_j = max_prob_class(i);
    confusion_matrix(confusion_matrix_i,confusion_matrix_j) ...
        = confusion_matrix(confusion_matrix_i,confusion_matrix_j) + 1;
end

correct_rate = sum(correct)/length(test_class);

end

function cv = cvpartition_package(class_data,name)
switch name
    case 'HOLD_OUT'
        cv = cvpartition(class_data, 'HoldOut',1/3);
    case 'CV-5'
        cv = cvpartition(class_data, 'kfold',5);
    case 'CV-10'
        cv = cvpartition(class_data, 'kfold',10);
    otherwise
        % error classify algorithm
        fprintf('%s classifier algorithm is not support',classify_algorithm);
        return;
end
end


function prob = bn_inference(bnet,unique_class_data_length,class_index,...
    test_data)

test_length = size(test_data,2);
prob = zeros(unique_class_data_length,test_length);

engine = jtree_inf_engine(bnet);

for i = 1:test_length
    evidence = num2cell(test_data(:,i));
    evidence(class_index,class_index) = {[]};
    [engine,~] = enter_evidence(engine, evidence);
    marg = marginal_nodes(engine, bnet.dnodes(class_index));
    prob(:,i) = marg.T;
end

end

