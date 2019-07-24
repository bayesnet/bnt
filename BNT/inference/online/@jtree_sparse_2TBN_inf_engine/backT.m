function b = backT(engine, f, t)

if t==1
  [b.clpot, seppot] = distribute_evidence(engine.jtree_engine1, f.clpot, f.seppot);
else
  [b.clpot, seppot] = distribute_evidence(engine.jtree_engine, f.clpot, f.seppot);
end
for c=1:length(b.clpot)
  [b.clpot{c}, ll(c)] = normalize_pot(b.clpot{c});
end
b.t = t;
