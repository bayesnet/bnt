function gdl = bnet_to_gdl_graph(bnet)
% BNET_TO_GDL_GRAPH Convert a Bayesian network to a GDL graph
% gdl = bnet_to_gdl_graph(bnet)
%
% Each node in the BN gets converted to a single node in the GDL graph,
% representing its family; its kernel function is the corresponding CPD.

N = length(bnet.dag);
doms = cell(1,N);
for i=1:N
  doms{i} = family(bnet.dag, i);
end

U = mk_undirected(bnet.dag);
gdl = mk_gdl_graph(U, doms, bnet.node_sizes, bnet.CPD, 'equiv_class', bnet.equiv_class, ...
		   'discrete', bnet.dnodes, 'chance', bnet.chance_nodes, ...
		   'decision', bnet.decision_nodes, 'utility', bnet.utility_nodes);

