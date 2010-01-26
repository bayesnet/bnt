function marginal = marginal_family(engine, i, t, add_ev)
% MARGINAL_FAMILY Compute the joint distribution on a set of family (filter_engine)
% function marginal = marginal_family(engine, i, t, add_ev)

if nargin < 4, add_ev = 0; end

if t ~= engine.t
  error('mixed up time stamps')
end

marginal = marginal_family(engine.tbn_engine, engine.b, i, t, add_ev);
