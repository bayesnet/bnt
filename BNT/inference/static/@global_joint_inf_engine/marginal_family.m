function [m, pot] = marginal_family(engine, i)
% MARGINAL_FAMILY Compute the marginal on i's family (global_inf_engine)
% [m, pot] = marginal_family(engine, i)
%

bnet = bnet_from_engine(engine);
[m, pot] = marginal_nodes(engine, family(bnet.dag, i));
