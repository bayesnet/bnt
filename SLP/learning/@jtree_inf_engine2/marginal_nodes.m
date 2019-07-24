function marginal = marginal_nodes(engine, query, add_ev)
% MARGINAL_NODES Compute the marginal on the specified query nodes (jtree)
% marginal = marginal_nodes(engine, query, add_ev)
%
% 'query' must be a subset of some clique; an error will be raised if not.
% add_ev is an optional argument; if 1, we will "inflate" the marginal of observed nodes
% to their original size, adding 0s to the positions which contradict the evidence

if nargin < 3, add_ev = 0; end

c = clq_containing_nodes(engine, query);
if c == -1
  error(['no clique contains ' num2str(query)]);
end
marginal = pot_to_marginal(marginalize_pot(engine.clpot{c}, query, engine.maximize));

if add_ev
  bnet = bnet_from_engine(engine);
  %marginal = add_ev_to_dmarginal(marginal, engine.evidence, bnet.node_sizes);
  marginal = add_evidence_to_gmarginal(marginal, engine.evidence, bnet.node_sizes, bnet.cnodes);
end

