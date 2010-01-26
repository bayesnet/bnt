function b = back(engine, bfuture, f, t)

if f.t ~= t
  error('mixed up time stamps')
end

b.t = t;
b.obslik = f.obslik;
bb_future = bfuture.beta .* bfuture.obslik;
if engine.maximize
  B = repmat(bb_future(:)', length(bfuture.beta), 1);
  b.beta = normalise(max(engine.transprob .* B, [], 2));
else
  b.beta = normalise((engine.transprob * bb_future));
end
b.gamma = normalise(f.alpha .* b.beta);
if t > 1
  bb_t = b.beta .* b.obslik;
  b.xi = normalise((engine.transprob .* (f.past_alpha * bb_t'))); % t-1,t
end

