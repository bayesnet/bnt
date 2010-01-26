function [m, pot] = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified nodes (jtree_limid)
% [m, pot] = marginal_nodes(engine, query)
%
% query should be a subset of a family of a decision node

if isempty(query)
  bnet = bnet_from_engine(engine);
  d = bnet.decision_nodes(1); % pick an arbitrary decision node
  [dummy, big_pot] = marginal_family(engine, d); 
else
  [dummy, big_pot] = marginal_family(engine, query);
end
pot = marginalize_pot(big_pot, query);
m = pot_to_marginal(pot);

  
