function [marginal, pot] = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified query nodes (belprop)
% [marginal, pot] = marginal_nodes(engine, query)
%
% query must be a subset of a family

if isempty(query)
  big_pot = engine.marginal_domains{1}; % pick an arbitrary domain
else
  big_pot = engine.marginal_domains{query(end)};   
end
pot = marginalize_pot(big_pot, query);
marginal = pot_to_marginal(pot);
 
