function check_jtree_property(cliques_containing_node, jtree)
% CHECK_JTREE_PROPERTY Raise an error if the graph does not satisfy the join tree property.
% check_jtree_property(cliques_containing_node, jtree_adj_mat)
%
% The join tree property says:
% For each node n in the dag, compute the node-induced subgraph G by looking at all the cliques
% that contain n, and make sure G forms a connected graph. 
% This ensures that local propagation leads to global consistency.

num_bn_nodes = length(cliques_containing_node);
directed = 0;
for i=1:num_bn_nodes
  cs = cliques_containing_node{i};
  G = jtree(cs,cs);
  if ~connected_graph(G, directed)
    error(['node ' num2str(i) ' violates jtree property']);
  end
end
