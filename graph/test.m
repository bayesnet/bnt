% A - B
%     |
% D - C - E

A=1;B=2;C=3;D=4;E=5;
dag = zeros(5,5);
dag(A,B)=1;
%dag(A,D)=1;
dag(B,C)=1;
dag(C,D)=1;
dag(E,C)=1;
[d, pre, post, cycle, f, pred] = dfs(dag, A, 0)

[T, pre, post, cycle] = mk_rooted_tree(dag, A)

%[T, pre, post, cycle] = mkRootedTree(dag, A)
