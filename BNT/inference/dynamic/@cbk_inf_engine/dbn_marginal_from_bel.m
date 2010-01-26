function marginal = dbn_marginal_from_bel(engine, i)
% DBN_MARGINAL_FROM_BEL Compute the marginal on a node given the current belief state (bk)
% marginal = dbn_marginal_from_bel(engine, i)
  
if engine.slice1
  j = i;
  c = clq_containing_nodes(engine.sub_engine1, j);
else
  bnet = bnet_from_engine(engine);
  ss = length(bnet.intra);
  j = i+ss;
  c = clq_containing_nodes(engine.sub_engine, j);
end
assert(c >= 1);
bigpot = engine.bel_clpot{c};

pot = marginalize_pot(bigpot, j);
marginal = pot_to_marginal(pot);
