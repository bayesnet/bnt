function b = connected(adj_mat, directed)
%
% b = connected(adj_mat).
% Returns true iff the graph is connected.

n = length(adj_mat);
start = 1;
[d, pre] = dfs(adj_mat, start, directed);
b = (length(pre) == n);
