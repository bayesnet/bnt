function marginal = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified query nodes (belprop)
% marginal = marginal_nodes(engine, query)

assert(length(query)==1);
marginal = pot_to_marginal(engine.marginal_nodes{query});
