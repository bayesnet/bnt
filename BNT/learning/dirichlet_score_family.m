function LL = dirichlet_score_family(counts, prior)
% DIRICHLET_SCORE Compute the log marginal likelihood of a single family
% LL = dirichlet_score(counts, prior)
%
% counts(a, b, ..., z) is the number of times parent 1 = a, parent 2 = b, ..., child = z
% prior is an optional multidimensional array of the same shape as counts.
% It defaults to a uniform prior.
% 
% We marginalize out the parameters:
% LL = log \int \prod_m P(x(i,m) | x(Pa_i,m), theta_i) P(theta_i) d(theta_i)


% LL = log[  prod_j gamma(alpha_ij)/gamma(alpha_ij + N_ij)  *
%            prod_k gamma(alpha_ijk + N_ijk)/gamma(alpha_ijk)  ]
% Call the prod_k term U and the prod_j term  V.
% We reshape all quantities into (j,k) matrices
% This formula was first derived by Cooper and Herskovits, 1992.
% See also "Learning Bayesian Networks", Heckerman, Geiger and Chickering, MLJ 95.

ns = mysize(counts);
ns_ps = ns(1:end-1);
ns_self = ns(end);

if nargin < 2, prior = normalise(myones(ns)); end


if 1
  prior = reshape(prior(:), [prod(ns_ps) ns_self]);
  counts = reshape(counts,  [prod(ns_ps) ns_self]);
  %U = prod(gamma(prior + counts) ./ gamma(prior), 2); % mult over k
  LU = sum(gammaln(prior + counts) - gammaln(prior), 2);
  alpha_ij = sum(prior, 2); % sum over k
  N_ij = sum(counts, 2);
  %V = gamma(alpha_ij) ./ gamma(alpha_ij + N_ij);
  LV = gammaln(alpha_ij) - gammaln(alpha_ij + N_ij);
  %L = prod(U .* V);
  LL = sum(LU + LV);
else
  CPT = mk_stochastic(prior + counts);
  LL = sum(log(CPT(:) .* counts(:)));
end

