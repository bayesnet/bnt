function node_sizes = get_node_sizes(data)

% this function can get node_sizes from data.
%
% node_sizes = get_node_sizes(data)
%
% input:
%    data: a matrix data，just support num.
% output:
%    node_sizes: the node sizes
%
% if you find bug, you can send email to me.
%
% Written by WANGXin(growlithe1205@gmail.com)
%
%
% 这个可以从数据中获得节点大小。
% 输入:
%    data: 一个矩阵数据,只支持数字类型
% 输出:
%    node_sizes: 节点大小
%
% 如果你发现了bug，你可以给我发email
%
% Written by WANGXin(growlithe1205@gmail.com)

[col,~]= size(data);

node_sizes = zeros(1,col);
for i = 1:col
    col_data = data(i,:);
    different_num = unique(col_data);
    node_sizes(i) = length(different_num);
end

node_sizes = node_sizes';
end