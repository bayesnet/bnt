function [m, pot] = marginal_nodes(engine, query)
% MARGINAL_NODES Compute the marginal on the specified set of nodes (global_joint)
% [m, pot] = marginal_nodes(engine, query)

pot = marginalize_pot(engine.jpot, query);
m = pot_to_marginal(pot);
%[m.T, lik] = normalize(m.T);
%loglik = log(lik);
