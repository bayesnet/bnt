function m = marginal_family(engine, b, i, t, add_ev)
% MARGINAL_FAMILY Compute the marginal on the specified family (jtree_2TBN)
% marginal = marginal_family(engine, b, i, t, add_ev)

bnet = bnet_from_engine(engine);
if t==1
  m = marginal_nodes(engine, b, family(bnet.dag, i), t, add_ev, 1);
else
  ss = length(bnet.intra);
  fam = family(bnet.dag, i+ss);
  m = marginal_nodes(engine, b, fam, t, add_ev, 1);
end
