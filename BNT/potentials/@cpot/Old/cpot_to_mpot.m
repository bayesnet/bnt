function mom = cpot_to_mpot(can)
% CPOT_TO_MPOT Convert a canonical potential to moment form.
% mom = cpot_to_mpot(can)

[logp, mu, Sigma] = canonical_to_moment(can.g, can.h, can.K);
mom = mpot(can.domain, can.sizes, logp, mu, Sigma);

%%%%%%%

function [logp, mu, Sigma] = canonical_to_moment(g, h, K)
% CANONICAL_TO_MOMENT Convert canonical characteristics to moment form.
% [logp, mu, Sigma] = canonical_to_moment(g, h, K)

if det(K)==0
  Sigma = inf*size(K);
else
  Sigma = inv(K);
end
mu = Sigma*h;
n = length(mu);
if isempty(mu)
  logp = g - 0.5*(log(det(K)) - n*log(2*pi));
else
  logp = g - 0.5*(log(det(K)) - n*log(2*pi) - mu'*K*mu);
end
