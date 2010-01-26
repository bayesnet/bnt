function [S, LL] = bic_score(counts, CPT, ncases)
% BIC_SCORE Bayesian Information Criterion score for a single family
% [S, LL] = bic_score(counts, CPT, ncases)
%
% S is a large sample approximation to the log marginal likelihood,
% which can be computed using dirichlet_score.
%
% S  = \log [ prod_j _prod_k theta_ijk ^ N_ijk ]  - 0.5*d*log(ncases) 
% where counts encode N_ijk, theta_ijk is the MLE comptued from counts,
% and d is the num of free parameters.

%CPT = mk_stochastic(counts);
tiny = exp(-700);
LL = sum(log(CPT(:)  + tiny) .* counts(:));
% CPT(i) = 0 iff counts(i) = 0 so it is okay to add tiny

ns = mysize(counts);
ns_ps = ns(1:end-1);
ns_self = ns(end);
nparams = prod([ns_ps (ns_self-1)]);
% sum-to-1 constraint reduces the effective num. vals of the node by 1

S = LL - 0.5*nparams*log(ncases);
