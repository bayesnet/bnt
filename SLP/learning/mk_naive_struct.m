function S = mk_naive_struct(n,C)
%
% S = mk_naive_struct(Number_of_nodes, Class_node)
%
S = zeros(n);
S(C,setdiff(1:n,C)) = 1;