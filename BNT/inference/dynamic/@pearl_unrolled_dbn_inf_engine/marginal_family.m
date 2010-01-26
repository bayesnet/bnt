function marginal = marginal_family(engine, i, t, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family (jtree_unrolled_dbn)
% marginal = marginal_family(engine, i, t)

if nargin < 3, t = 1; end
if nargin < 4, add_ev = 0; end
assert(~add_ev);

%marginal = marginal_family(engine.unrolled_engine, i + (t-1)*engine.ss, add_ev);
marginal = marginal_family(engine.unrolled_engine, i + (t-1)*engine.ss);
              
