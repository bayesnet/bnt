function fg = bnet_to_fgraph(bnet)
% BNET_TO_FGRAPH Convert a Bayes net to a factor graph
% fg = bnet_to_fgraph(bnet)
%
% We create one factor per family, whose kernel is the CPD

nnodes = length(bnet.dag);
G = zeros(nnodes, nnodes);
for i=1:nnodes
  G(family(bnet.dag, i), i) = 1;
end

fg = mk_fgraph(G, bnet.node_sizes, bnet.CPD, 'equiv_class', bnet.equiv_class, 'discrete', bnet.dnodes);

	       

