function b = backT(engine, f, t)

b.t = t;
b.obslik = f.obslik;
Q = length(f.alpha);
b.beta = ones(Q,1);
b.gamma = f.alpha;
if t > 1
  bb_t = b.obslik;
  b.xi = normalise((engine.transprob .* (f.past_alpha * bb_t'))); % T-1,T
end
