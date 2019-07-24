function b = back1(engine, bfuture, f, t)

if t ~= 1
  error('mixed up time stamps')
end
bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

int = engine.interface;
D = engine.in_clq; % from J2
C = engine.int_clq1; % from J1
phiD = marginalize_pot(bfuture.clpot{D}, int, engine.maximize);
phiC = marginalize_pot(f.clpot{C}, int, engine.maximize);
ratio = divide_by_pot(phiD, phiC);
f.clpot{C} = multiply_by_pot(f.clpot{C}, ratio);

[b.clpot, seppot] = distribute_evidence(engine.jtree_engine1, f.clpot, f.seppot);
for c=1:length(b.clpot)
  [b.clpot{c}, ll(c)] = normalize_pot(b.clpot{c});
end
b.t = t;
