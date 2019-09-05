function [mdl_data,cut_point] = single_mdl_algorithm(row_data,class)

% this function achieve Fayyad & Irani's MDL method.
%
% [mdl_data,cut_point] = single_mdl_algorithm(row_data,class)
%
% input:
%     row_data : a row data, just support one row and num data
%     class: class data, of course one row and num data
% output:
%     mdl_data: the result of mdl algorithm process data
%     cut_point: row_data cut point
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个函数实现了Fayyad & Irani 的mdl算法
%
% 输入:
%    row_data : 一行数据，只能是一行数字形式的数据
%    class: 类别数据，也是一行数字形式的数据
% 输出:
%    mdl_data: mdl 算法处理数据的结果
%    cut_point: row_data 的切点
%
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

mdl_data = row_data;
cut_point = mdl_core(row_data,class);
if isempty(cut_point)
    return
end
cut_point = sort(cut_point);
for i = 1: length(cut_point)+1
    switch i
        case 1
            [~,col,~] = find(row_data <= cut_point(1));
        case length(cut_point)+1
            [~,col,~] = find(row_data > cut_point(length(cut_point)));
        otherwise
            [~,col,~] = find(row_data > cut_point(i-1) & row_data <= cut_point(i));
    end
    
    mdl_data(col) = i;
end

mdl_data = mdl_data';

end


function cut_point = mdl_core(data,class)

% mdl algorithm core
% cut_point = mdl_algorithm(data,class)
% Input
%     data
%     class
% Output:
%    cut_point

cut_point = [];
[~,min_cut_point,~] = get_min_cut_point(data,class);
if isempty(min_cut_point)
    return
end

[left_data,left_class,right_data,right_class] = split_data(data,class,min_cut_point);

cost_nt = computer_cost_nt(class);
cost_ht = computer_cost_ht(left_class,right_class);
if cost_nt < cost_ht
    return
end
if cost_nt >= cost_ht
    cut_point = [cut_point,min_cut_point];
   
    left_min_cut_point = mdl_core(left_data,left_class);
    if ~isempty(left_min_cut_point)
        cut_point = [cut_point,left_min_cut_point];
    end
    
    right_min_cut_point = mdl_core(right_data,right_class);
    if ~isempty(right_min_cut_point)
       cut_point = [cut_point,right_min_cut_point];
    end
   
    cut_point = unique(cut_point);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ent = computer_ent(class)

% computer ent
% input: class
% output: ent

% 1 data preprocess
[data_r,data_n] = size(class);
if (data_r ~= 1)
    error('data in Ent algorithm is not 1 row');
end

% 2 algorithm
class_num = unique(class);
every_entropy = zeros(length(class_num),1);
for i =1:length(class_num)
    [time,~] =find(class == class_num(i));
    times = sum(time);
    probability = (times/data_n);
    every_entropy(i) = probability * log2(probability);
end

ent = -sum(every_entropy);

% 3 return pass

end


% computer cost_nt
% input:
% output:
function cost_nt = computer_cost_nt(class)

% 1 data preprocess
[~,n] = size(class);

% 2 algorithm
k = length(unique(class));
cost_nt = (n+k) * computer_ent(class);

% 3 return pass

end


function cost_ht = computer_cost_ht(left_class,right_class)

% computer cost_ht
% input
% output

% 1 data preprocess pass

% 2 algorithm
class = [left_class right_class];

% get n k
[~,n] = size(class);
k = length(unique(class));

common_ht = log2(n-1) + log2(3^k-2);
left_ht = computer_cost_nt(left_class);
right_ht = computer_cost_nt(right_class);

cost_ht = common_ht + left_ht + right_ht;
% 3 return pass

end


function cut_entropy = get_cut_entropy(left_class,right_class)

%  algorithm
class = [left_class,right_class];
class_length = length(class);

left_cut_entropy = length(left_class)/class_length * computer_ent(left_class);
right_cut_entropy = length(right_class)/class_length * computer_ent(right_class);
cut_entropy = left_cut_entropy + right_cut_entropy;

end


function [min_cut_entropy,min_cut_point,min_cut_point_index] = get_min_cut_point(data,class)

point = unique(data);
point_length = length(point);
min_cut_point = [];
min_cut_entropy = [];
min_cut_point_index = [];
for i = 1:point_length-1
    cut_point = (point(i) + point(i+1))/2;
    [~,left_class,~,right_class] = split_data(data,class,cut_point);
    cut_entropy = get_cut_entropy(left_class,right_class);
    if isempty(min_cut_entropy) || min_cut_entropy > cut_entropy
        min_cut_entropy = cut_entropy;
        min_cut_point = cut_point;
        min_cut_point_index = i;
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [left_data,left_class,right_data,right_class] = split_data(data,class,cut_point)

% split two groups
% input:
%       data,class,cut_point
% output:
%       left_data,left_class,right_data,right_class


% 1 data preprocess
data_check(data,class);
if cut_point > max(data) || cut_point < min(data)
    error('cut point must in data');
end
not_cut_point = unique(data);
for i=1:length(not_cut_point)
    if cut_point == not_cut_point(i)
        error('cut point must not equal data value');
    end
end

% 2 algorithm
left_data_index = find(data<cut_point);
left_data = data(left_data_index);
left_class = class(left_data_index);
right_data_index = find(data>cut_point);
right_data = data(right_data_index);
right_class = class(right_data_index);

% 3.return

end

% just check data
function [data_r,data_c] = data_check(data,class)

% 1 data preprocess
[data_r,data_c] = size(data);
[index_r,index_c] = size(class);
if(data_r ~= 1 || index_r ~= 1)
    error('must a row data');
end
if (data_c ~= index_c)
    error('data and index are not equal size')
end

for i = 1:length(class)
    if(~is_integer(class(i)))
        error('class must integer array')
    end
end

end


function flag = is_integer(x)

% this function check whether data is integer

flag = 0;

if x == fix(x)
    flag = 1;
end

return

end