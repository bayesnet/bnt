function can = mpot_to_cpot(mom)
% MPOT_TO_CPOT Convert a moment potential to canonical form.
% mom = mpot_to_cpot(can)

[g, h, K] = moment_to_canonical(mom.logp, mom.mu, mom.Sigma);
can = cpot(mom.domain, mom.sizes, g, h, K);

%%%%%%%%%%%

function [g, h, K] = moment_to_canonical(logp, mu, Sigma)
% MOMENT_TO_CANONICAL Convert moment characteristics to canonical form.
% [g, h, K] = moment_to_canonical(logp, mu, Sigma)

K = inv(Sigma);
h = K*mu;
n = length(K);
if isempty(mu)
  g = logp + 0.5*(log(det(K)) - n*log(2*pi));
else
  g = logp + 0.5*(log(det(K)) - n*log(2*pi) - mu'*K*mu);
end        
