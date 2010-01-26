function mom = cpot_to_mpot(can)
% CPOT_TO_MPOT Convert a canonical potential to moment form.
% mom = cpot_to_mpot(can)

[logp, mu, Sigma] = canonical_to_moment(can.g, can.h, can.K);
mom = mpot(can.domain, can.sizes, logp, mu, Sigma);

%%%%%%%

function [logp, mu, Sigma] = canonical_to_moment(g, h, K)
% CANONICAL_TO_MOMENT Convert canonical characteristics to moment form.
% [logp, mu, Sigma] = canonical_to_moment(g, h, K)

n = length(K);
if isempty(K)
  logp = g - 0.5*(log(det(K)) - n*log(2*pi));
  Sigma = [];
  mu = [];
else
  if det(K)==0
    Sigma = inf*ones(n,n);
    mu = zeros(n,1); % if the precision is zero, the mean is arbitrary
    logp = g; % the scaling factor for the uniform distribution is 1
  else
    Sigma = inv(K);
    mu = Sigma*h;
    logp = g - 0.5*(log(det(K)) - n*log(2*pi) - mu'*K*mu);
  end
end
