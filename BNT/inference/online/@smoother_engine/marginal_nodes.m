function marginal = marginal_nodes(engine, nodes, t, add_ev)
% MARGINAL_NODES Compute the joint distribution on a set of nodes (smoother_engine)
% function marginal = marginal_nodes(engine, nodes, t, add_ev)

if nargin < 4, add_ev = 0; end

marginal = marginal_nodes(engine.tbn_engine, engine.b{t}, nodes, t, add_ev);
