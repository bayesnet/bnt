function [marginal, pot] = marginal_family(engine, query)
% MARGINAL_NODES Compute the marginal on the family of the specified query node (belprop)
% [marginal, pot] = marginal_family(engine, query)

pot = engine.marginal_domains{query};
marginal = pot_to_marginal(pot);
