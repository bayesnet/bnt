function [f, logscale] = fwd1(engine, ev, t)
% Forwards pass for slice 1.

bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

CPDpot = cell(1,ss);
for n=1:ss
  fam = family(bnet.dag, n, 1);
  e = bnet.equiv_class(n, 1);
  CPDpot{n} = convert_to_pot(bnet.CPD{e}, engine.pot_type, fam(:), ev);
end       
f.evidence = ev;
f.t = t;

pots = CPDpot;
slice1 = 1:ss;
CPDclqs = engine.clq_ass_to_node1(slice1);

[f.clpot, f.seppot] =  init_pot(engine.jtree_engine1, CPDclqs, CPDpot, engine.pot_type, engine.observed1);
[f.clpot, f.seppot] = collect_evidence(engine.jtree_engine1, f.clpot, f.seppot);
for c=1:length(f.clpot)
  if isa(f.clpot{c}, 'struct')
     domain = f.clpot{c}.domain;
     sizes = f.clpot{c}.sizes;
     T = f.clpot{c}.T;
     f.clpot{c} = dpot(domain, sizes, T);
  end
  [f.clpot{c}, ll(c)] = normalize_pot(f.clpot{c});
end
logscale = ll(engine.root1);

