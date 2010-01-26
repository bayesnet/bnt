function [b, mpe] = back_mpe(engine, bfuture, f, ev2, t)

if f.t ~= t
  error('mixed up time stamps')
end
if t==1
  error('should call back1_mpe')
end

maximize = 1;
bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

int = engine.interface;
D = engine.in_clq;
C = engine.out_clq;
phiD = marginalize_pot(bfuture.clpot{D}, int, maximize);
phiD = set_domain_pot(phiD, int+ss); % shift to slice 2
phiC = marginalize_pot(f.clpot{C}, int+ss, maximize);
ratio = divide_by_pot(phiD, phiC);
f.clpot{C} = multiply_by_pot(f.clpot{C}, ratio);

[mpe, b.clpot] = find_max_config(engine.jtree_engine, f.clpot, f.seppot, ev2);
mpe = mpe((1:ss)+ss); % extract values for slice 2
for c=1:length(b.clpot)
  [b.clpot{c}, ll(c)] = normalize_pot(b.clpot{c});
end
b.t = t;
