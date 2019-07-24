function [b, mpe] = backT_mpe(engine, f, ev2, t)

bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

if t==1
  % ev2 is just the evidence on slice 1
  [mpe, b.clpot] = find_max_config(engine.jtree_engine1, f.clpot, f.seppot, ev2);
else
  [mpe, b.clpot] = find_max_config(engine.jtree_engine, f.clpot, f.seppot, ev2);
  mpe = mpe((1:ss)+ss); % extract values for slice 2
end
for c=1:length(b.clpot)
  [b.clpot{c}, ll(c)] = normalize_pot(b.clpot{c});
end
b.t = t;
