function marginal = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified query nodes (belprop)
% marginal = marginal_nodes(engine, query)
%
% query must be a single node

if length(query)>1
  error('can only handle single node marginals')
end
marginal = engine.bel{query};
