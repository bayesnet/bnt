function marginal = marginal_family(engine, i, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family (jtree)
% marginal = marginal_family(engine, i)

if nargin < 3, add_ev = 0; end
assert(~add_ev);

bnet = bnet_from_engine(engine);
fam = family(bnet.dag, i);
c = engine.clq_ass_to_node(i);
marginal = pot_to_marginal(marginalize_pot(engine.clpot{c}, fam));
