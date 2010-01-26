function [f, logscale] = fwd(engine, fpast, ev, t)
% Forwards pass.

f.obslik = mk_hmm_obs_lik_vec(engine, ev);
transmat = engine.transprob;
f.past_alpha = fpast.alpha;
if engine.maximize
  Q = length(fpast.alpha);
  A = repmat(fpast.alpha, [1 Q]);
  m = max(transmat .* A, [], 1);
  [f.alpha, scale] = normalise(m(:) .* f.obslik);
else
  [f.alpha, scale] = normalise((transmat' * fpast.alpha) .* f.obslik);
end
logscale = log(scale);
%f.xi = normalise((fpast.alpha * obslik') .* transmat); % t-1,t
f.t = t;
