function marginal = marginal_nodes(engine, b, nodes, t, add_ev, is_fam)
% function marginal = marginal_nodes(engine, b, nodes, t, add_ev, is_fam) (jtree_2TBN)

if nargin < 6, is_fam = 0; end
ss = engine.slice_size;

if ~is_fam & (t > 1) & all(nodes<=ss)
  nodes = nodes + ss;
end

if t==1
  c = clq_containing_nodes(engine.jtree_engine1, nodes, is_fam);
else
  c = clq_containing_nodes(engine.jtree_engine, nodes, is_fam);
end
if c == -1
  error(['no clique contains ' nodes])
end
bigpot = b.clpot{c};
pot = marginalize_pot(bigpot, nodes, engine.maximize);
marginal = pot_to_marginal(pot);

% we convert the domain to the unrolled numbering system
% so that add_ev_to_dmarginal (maybe called in update_ess) extracts the right evidence.
if t > 1
  marginal.domain = nodes+(t-2)*engine.slice_size;
end
assert(~add_ev);




