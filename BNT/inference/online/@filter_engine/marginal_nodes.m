function marginal = marginal_nodes(engine, nodes, t, add_ev)
% MARGINAL_NODES Compute the joint distribution on a set of nodes (filter_engine)
% function marginal = marginal_nodes(engine, nodes, t, add_ev)

if nargin < 4, add_ev = 0; end

if t ~= engine.t
  error('mixed up time stamps')
end
marginal = marginal_nodes(engine.tbn_engine, engine.b, nodes, t, add_ev);
