function marginal = marginal_family(engine, i, t, add_ev)
% MARGINAL_FAMILY Compute the joint distribution on a set of family (smoother_engine)
% function marginal = marginal_family(engine, i, t, add_ev)

if nargin < 4, add_ev = 0; end
marginal = marginal_family(engine.tbn_engine, engine.b{t}, i, t, add_ev);
