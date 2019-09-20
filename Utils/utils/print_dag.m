function node2node = print_dag(dag,node_name)

% this function will print relationship of node in dag 
%
% node2node = print_dag(dag,node_name)
% Input：
%     dag  
%     node_name 
%
% Output：
%     printf 
% 
% example: dag = [0 1 1;0 0 1;0 0 0]; 
%          node_name = ["A1","A2","A3"];
%          print_dag(dag,node_name);
%
% Written by WANGXin(growlithe1205@gmail.com)
%
% 该函数用于打印输出一个图模型中各节点的关系，允许该图不为DAG。
% 如果发现BUG请联系作者本人(growlithe1205@gmail.com)。

[r,c] = size(dag);

if r~=c
    disp('The input dag is not a square matrix')
    return
end
if r ~= length(node_name)
    disp('The label length is not equal dag size')
    return
end

order_name = cell(2,length(node_name));

for i=1:r
    order_name{1,i} = i;
    order_name{2,i} = node_name(1,i);
end

[row,col] = find(dag==1);
edge = [row,col];

node2node = cell(size(edge,1),1);
for i=1:size(edge,1)
    a_node = order_name{2,edge(i,1)};
    b_node = order_name{2,edge(i,2)};
    node2node{i} = sprintf('dag(node.%s,node.%s)=true;',a_node,b_node);
    fprintf('dag(node.%s,node.%s)=true;',a_node,b_node);
    fprintf('\n');
end

end