function marginal = marginal_family(engine, i, t)
% MARGINAL_FAMILY Compute the marginal on node i in slice t and its parents  (frontier)
% marginal = marginal_family(engine, i, t)

bnet = bnet_from_engine(engine);    
fam = family(bnet.dag, i, t);
marginal = pot_to_marginal(normalize_pot(marginalize_pot(engine.fwdback{i,t}, fam)));
