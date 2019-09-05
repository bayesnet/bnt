function [jtree, root, cliques, B, w, elim_order] = graph_to_jtree(MG, ns, partial_order, stages, clusters)
% GRAPH_TO_JTREE Triangulate a graph and make a junction tree from its cliques.
% [jtree, root, cliques, B, w, elim_order] = ...
%    graph_to_jtree(graph, node_sizes, partial_order, stages, clusters)
%
% INPUT:
% graph(i,j) = 1 iff there is an edge between i,j
% node_weights(i) = num discrete values node i can take on [1 if observed]
% partial_order = {} if no constraints on elimination ordering
% stages{i} = nodes that must be eliminated at i'th stage (if porder is empty)
% clusters{i} = list of nodes that must get connected together in the moral graph
%
% OUTPUT:
% jtree(i,j) = 1 iff there is an arc between clique i and clique j 
% root = the root clique
% cliques{i} = the nodes in clique i
% B(i,j) = 1 iff node j occurs in clique i
% w(i) = weight of clique i

N = length(MG);

if nargin >= 5
  % Add extra arcs between nodes in each cluster to ensure they occur in the same clique
  for i=1:length(clusters)
    c = clusters{i};
    MG(c,c) = 1;
  end
end
MG = setdiag(MG, 0);

% Find an optimal elimination ordering (NP-hard problem!)
if nargin < 4
  stages = {1:N};
end
if nargin < 3
  partial_order = {};
end
if isempty(partial_order)
  strong = 0;
  elim_order = best_first_elim_order(MG, ns, stages);
else
  strong = 1;
  elim_order = strong_elim_order(MG, ns, partial_order);
end

[MTG, cliques, fill_in_edges]  = triangulate(MG, elim_order);

% Connect the cliques up into a jtree,
[jtree, root, B, w] = cliques_to_jtree(cliques, ns);

if 0
  disp('testing dag to jtree');
  % Find the cliques containing each node, and check they form a connected subtree
  clqs_con_node = cell(1,N);
  for i=1:N
    clqs_con_node{i} = find(B(:,i))';
  end
  check_jtree_property(clqs_con_node, jtree);
end
