function [jtree, root, B, w] = cliques_to_jtree(cliques, ns)
% MK_JTREE Make an optimal junction tree.
% [jtree, root, B, w] = mk_jtree(cliques, ns)
% 
% A junction tree is a tree that satisfies the jtree property, which says:
% for each pair of cliques U,V with intersection S, all cliques on the path between U and V
% contain S. (This ensures that local propagation leads to global consistency.)
%
% We can create a junction tree by computing the maximal spanning tree of the junction graph.
% (The junction graph connects all cliques, and the weight of an edge (i,j) is
% |C(i) intersect C(j)|, where C(i) is the i'th clique.)
%
% The best jtree is the maximal spanning tree which minimizes the sum of the costs on each edge,
% where cost(i,j) = w(C(i)) + w(C(j)), and w(C) is the weight of clique C,
% which is the total number of values C can take on.
%
% For details, see
%  - Jensen and Jensen, "Optimal Junction Trees", UAI 94.
%
% Input:
%  cliques{i} = nodes in clique i
%  ns(i) = number of values node i can take on
% Output:
%  jtree(i,j) = 1 iff cliques i and j aer connected
%  root = the clique that should be used as root
%  B(i,j) = 1 iff node j occurs in clique i
%  w(i) = weight of clique i



num_cliques = length(cliques);
w = zeros(num_cliques, 1); 
B = sparse(num_cliques, 1);
for i=1:num_cliques
  B(i, cliques{i}) = 1;
  w(i) = prod(ns(cliques{i}));
end


% C1(i,j) = length(intersect(cliques{i}, cliques{j})); 
% The length of the intersection of two sets is the dot product of their bit vector representation.
C1 = B*B';
C1 = setdiag(C1, 0);

% C2(i,j) = w(i) + w(j)
num_cliques = length(w);
W = repmat(w, 1, num_cliques);
C2 = W + W';
C2 = setdiag(C2, 0);

jtree = sparse(minimum_spanning_tree(-C1, C2)); % Using -C1 gives *maximum* spanning tree

% The root is arbitrary, but since the first pass is towards the root,
% we would like this to correspond to going forward in time in a DBN.
root = num_cliques;


