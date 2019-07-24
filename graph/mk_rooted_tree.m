function [T, pre, post, cycle] = mk_rooted_tree(G, root)
% MK_ROOTED_TREE Make a directed tree pointing away from root
% [T, pre, post, cycle] = mk_rooted_tree(G, root)

n = length(G);
T = sparse(n,n); % not the same as T = sparse(n) !
directed = 0;
[d, pre, post, cycle, f, pred] = dfs(G, root, directed);
[junk, pre2] = sort(d);
assert(isequal(pre, pre2))
[junk, post2] = sort(f);
assert(isequal(post, post2));
%[d, pre, post, cycle, f, pred] = dfs(G, [], directed);
for i=1:length(pred)
  if pred(i)>0
    T(pred(i),i)=1;
  end
end

