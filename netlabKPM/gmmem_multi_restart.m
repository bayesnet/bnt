function [means, covs, weights, ll] = gmmem_multi_restart(K, data, varargin)
% GMMEM_MULTI_RESTART Multiple restart wrapper for gmmem_kpm
% function [means, covs, weights, ll] = gmmem_multi_restart(K, data, varargin)
%
% Input:
% K = number of mixture components
% data(i,:) is the i'th example (feature vector)
%
% Output:
% The parameters for the k'th mixture component, k=1:K, are
%   means(k,:), covs(:,:,k) and weights(k)
%
% [ ... ] = gmmem_multi_restart(..., 'param1',val1, 'param2',val2, ...)
% allows you to specify optional parameter name/value pairs.
% Parameters are below [default value in brackets]
%
% 'nrestarts' - number of EM restarts [2]
% 'cov_type' - 'full', 'diag' or 'spherical' ['full']
% 'init_cov' - the initial covariance matrix [0.1*cov(data) for each k]
% 'init_means' - [] means sample from randn(); otherwise, use
%               init_means(k,:,r) for the k'th comp. on the r'th restart [ [] ]
% 'restartfn'  - this function, if non-empty,  will be called before/after every restart
%               (e.g., to display the parameters as they evolve) [ [] ]
%               The fn is called as fn(mix{r}, data, restart_num, niter, outerfnargs)
%               where niter is the number of iterations performed (0 initially)
% 'restartfnargs' - additional arguments to be passed to restartfn [ {} ]
%
% Optional arguments for gmmem_kpm are passed through.
%
% Written by Kevin P Murphy, 30 Dec 2002

[ndata nfeatures] = size(data);

%Cinit = repmat(0.1*diag(diag(cov(data))), [1 1 K]);
Cinit = repmat(0.1*cov(data), [1 1 K]);

[nrestarts, init_cov, init_means, cov_type, ...
 restartfn, restartfnargs, unused_args] = ...
    process_options(varargin, ...
	'nrestarts', 2, 'init_cov',  Cinit,  'init_means', [], ...
	'cov_type', 'full', 'restartfn', [], 'restartfnargs', {});

mix = cell(1, nrestarts);
cost = inf*ones(1,nrestarts);

for r=1:nrestarts
  mix{r} = gmm(nfeatures, K, cov_type); % random centers
  if ~isempty(init_means), mix{r}.centres = init_means(:,:,r); end
  mix{r}.covars = init_cov;
  if ~isempty(restartfn)
    feval(restartfn, mix{r}, data, r, 0, restartfnargs{:});
  end
  [mix{r}, niter, ll] = gmmem_kpm(mix{r}, data, unused_args{:}); 
  cost(r) = -ll; %-sum(log(gmmprob(mix{r}, data)));
  if ~isempty(restartfn)
    feval(restartfn, mix{r}, data, r, niter, restartfnargs{:});
  end
end

[nll, bestr] = min(cost);
fprintf('best r = %d\n', bestr);
ll = -nll;
means = mix{bestr}.centres;
covs = mix{bestr}.covars;
weights = mix{bestr}.priors;
