function [f, logscale] = fwd1(engine, ev, t)
% Forwards pass for slice 1.

if t ~= 1
  error('mixed up time stamps')
end
prior = engine.startprob(:);
f.obslik = mk_hmm_obs_lik_vec(engine, ev);
[f.alpha, lik] = normalise(prior .* f.obslik);
logscale = log(lik);
f.t = t;
