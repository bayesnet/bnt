function m = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified query nodes (quickscore)
% marginal = marginal_nodes(engine, query)
%
% 'query' must be a single disease (root) node.

assert(length(query)==1);
p = engine.post(query);
m.T = [1-p p]';
m.domain = query;
  
