function [d, pre, post, height, cycle, pred] = dfs(adj_mat, start, directed)
% DFS Perform a depth-first search of the graph starting from 'start'.
% [d, pre, post, height, cycle, pred] = dfs(adj_mat, start, directed)
%
% d(i) is the time at which node i is first discovered.
% pre is a listing of the nodes in the order in which they are first encountered (opened).
% post is a listing of the nodes in the order in which they are last encountered (closed).
% A node is last encountered once we have explored all of its neighbors.
% If the graph is directed, i's neighbors are its children.
% If the graph is a tree, preorder is parents before children, and
% postorder is children before parents.
% For a DAG, topological order = reverse(postorder).
% height(i) is the height (distance) of node i from the start.
% 'cycle' is true iff a (directed) cycle is found.
% pred(i) is the parent of i in the dfs tree rooted at start.
% See Cormen, Leiserson and Rivest, "An intro. to algorithms" 1994, p478.

% We can detect undirected cycles by checking if we are about to visit a node n which we have
% already visited. To detect *directed* cycles, we need to know if n has been closed or is still open.
% For example (where arcs are directed down)
%   1    2
%   \   /
%     3
% Assume we visit 1, 3 and then 2 in order. The fact that a child of 2 (namely, 3) has
% already been visited is okay, because 3 has been closed.
% The algorithms in Aho, Hopcroft and Ullman, and Sedgewick, do not detect directed cycles.

n = length(adj_mat);

global white gray black
white = 0; gray = 1; black = 2;

color = white*ones(1,n);
d = zeros(1,n);
height = zeros(1,n);
pred = zeros(1,n);
pre = [];
post = [];
cycle = 0;
global count
count = 0;
h = 0;
[d, pre, post, height, cycle, color, pred] = ...
    dfs2(adj_mat, start, directed, h, d, pre, post, height, cycle, color, pred);



%%%%%%%%%%

function [d, pre, post, height, cycle, color, pred] = ...
    dfs2(adj_mat, i, directed, h, d, pre, post, height, cycle, color, pred)

global count
global white gray black

color(i) = gray;
count = count + 1;
d(i) = count;
pre = [pre i];
height(i) = h;
if directed
  ns = children(adj_mat, i);
else
  ns = neighbors(adj_mat, i);
end
for j=1:length(ns)
  n=ns(j);
  if ~directed & n==pred(i) % don't go back up the edge you just came down
    % continue
  else
    if color(n) == gray % going back to a non-closed vertex via a new edge
      %fprintf(1, 'cycle from %d to %d\n', i, n);
      cycle = 1;
    end
    if color(n) == white % not visited n before
      pred(n)=i;
      [d, pre, post, height, cycle, color, pred] = ...
	  dfs2(adj_mat, n, directed, h+1, d, pre, post, height, cycle, color, pred);
    end
  end
end
color(i) = black;
post = [post i];

