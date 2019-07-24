function [jtree, root, cliques, B, w, elim_order, moral_edges, fill_in_edges, strong] = ...
    dag_to_jtree(dag, node_sizes, partial_order, stages, clusters)
% DAG_TO_JTREE Moralize and triangulate a DAG, and make a junction tree from its cliques.
% [jtree, root, cliques, B, w, elim_order, moral_edges, fill_in_edges, strong] = ...
%    dag_to_jtree(dag, node_sizes, partial_order, stages, clusters)
%
% Input:
% dag(i,j) 
% jtree(i,j) = 1 iff there is an arc between clique i and clique j 
% root = the root clique
% cliques{i} = the nodes in clique i
% B(i,j) = 1 iff node j occurs in clique i
% w(i) = weight of clique i

N = length(bnet.dag);
if nargin < 2, obs_nodes = []; end
if nargin < 3, stages = { 1:N }; end
if nargin < 4, clusters = {}; end

[MG, moral_edges]  = moralize(bnet.dag);

% Add extra arcs between nodes in each cluster to ensure they occur in the same clique
for i=1:length(clusters)
  c = clusters{i};
  MG(c,c) = 1;
end
MG = setdiag(MG, 0);

% Find an optimal elimination ordering (NP-hard problem!)
ns = bnet.node_sizes(:);
ns(obs_nodes) = 1; % observed nodes have only 1 possible value
partial_order = determine_elim_constraints(bnet, obs_nodes);

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
