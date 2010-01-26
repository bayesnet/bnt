function b = back(engine, bfuture, f, t)

if f.t ~= t
  error('mixed up time stamps')
end
if t==1
  b = back1(engine, bfuture, f, t);
  return;
end

bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

int = engine.interface;
D = engine.in_clq;
C = engine.out_clq;
phiD = marginalize_pot(bfuture.clpot{D}, int, engine.maximize);
phiD = set_domain_pot(phiD, int+ss); % shift to slice 2
phiC = marginalize_pot(f.clpot{C}, int+ss, engine.maximize);
ratio = divide_by_pot(phiD, phiC);
f.clpot{C} = multiply_by_pot(f.clpot{C}, ratio);

[b.clpot, seppot] = distribute_evidence(engine.jtree_engine, f.clpot, f.seppot);
for c=1:length(b.clpot)
  [b.clpot{c}, ll(c)] = normalize_pot(b.clpot{c});
end
b.t = t;
