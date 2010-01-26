function pot = marginal_family_pot(engine, i)
% MARGINAL_FAMILY_POT Compute the marginal on i's family and return as a potentila (inf_engine)
% function pot = marginal_family_pot(engine,i)

% This function is only called by solve_limid.
% It requires that engine's marginal_family function return a potential.
% This is true for jtree_inf_engine, but not for, say, jtree_ndx_inf_engine.
% All limids must be solved using potentials,
% but this is not true for bnets.

%[m, pot] = marginal_family(engine, i);

bnet = bnet_from_engine(engine);
[m, pot] = marginal_nodes(engine, family(bnet.dag, i));
