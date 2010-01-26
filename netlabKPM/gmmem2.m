function [mix, num_iter, ll] = gmmem_kpm(mix, x, varargin)
%GMMEM_KPM Like GMMEM, but with additional optional arguments
% function [mix, num_iter, ll] = gmmem_kpm(mix, x, varargin)
%
% Input:
% mix - structure created by gmminit or gmmem_multi_restart
% data - each row is an example
%
% Output:
% mix - modified structure
% num_iter - number of iterations needed to reach convergence
% ll - final log likelihood
%
% [ ... ] = gmmem_kpm(..., 'param1',val1, 'param2',val2, ...) allows you to
% specify optional parameter name/value pairs.
% Parameters are below [default value in brackets]
%
% 'max_iter' - maximum number of EM iterations [10]
% 'll_thresh' - change in log-likelihood threshold for convergence [1e-2]
% 'verbose' - 1 means display output while running [0]
% 'prior_cov' - this will be added to each estimated covariance
%               to prevent singularities  [1e-3*eye(d)]
% 'fn'        - this function, if non-empty,  will be called at every iteration
%               (e.g., to display the parameters as they evolve) [ [] ]
%               The fn is called as fn(mix, x, iter_num, fnargs).
%               It is also called before the iteration starts as
%               fn(mix, x, -1, fnargs), which can be used to initialize things.
% 'fnargs'    - additional arguments to be passed to fn [ {} ]
%
% Modified by Kevin P Murphy, 29 Dec 2002


% Check that inputs are consistent
errstring = consist(mix, 'gmm', x);
if ~isempty(errstring)
  error(errstring);
end

[ndata, xdim] = size(x);

[max_iter, ll_thresh, verbose, prior_cov, fn, fnargs] = ...
    process_options(varargin, ...
	'max_iter', 10, 'll_thresh', 1e-2, 'verbose', 1, ...
	'prior_cov', 1e-3*eye(xdim), 'fn', [], 'fnargs', {});

options = foptions;
if verbose, options(1)=1; else options(1)=-1; end
options(14) = max_iter;
options(3) = ll_thresh;


% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
test = 0;
if options(3) > 0.0
  test = 1;	% Test log likelihood for termination
end

check_covars = 0;
if options(5) >= 1
  if display >= 0
    disp('check_covars is on');
  end
  check_covars = 1;	% Ensure that covariances don't collapse
  MIN_COVAR = eps;	% Minimum singular value of covariance matrix
  init_covars = mix.covars;
end

mix0 = mix; % save init values for debugging

if ~isempty(fn)
  feval(fn, mix, x, -1, fnargs{:});
end

% Main loop of algorithm
for n = 1:niters
  
  % Calculate posteriors based on old parameters
  [post, act] = gmmpost(mix, x);
  
  % Calculate error value if needed
  if (display |  test)
    prob = act*(mix.priors)';
    % Error value is negative log likelihood of data
    e = - sum(log(prob + eps));
    if display > 0
      fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
    end
    if test
      if (n > 1 & abs(e - eold) < options(3))
        options(8) = e;
	ll = -e;
	num_iter = n;
        return; %%%%%%%%%%%%%%%% Exit here if converged
      else
        eold = e;
      end
    end
  end

  if ~isempty(fn)
    feval(fn, mix, x, n, fnargs{:});
  end

  % Adjust the new estimates for the parameters
  new_pr = sum(post, 1);
  new_c = post' * x;
  
  % Now move new estimates to old parameter vectors
  mix.priors = new_pr ./ ndata;
  
  mix.centres = new_c ./ (new_pr' * ones(1, mix.nin));
  
  switch mix.covar_type
  case 'spherical'
    n2 = dist2(x, mix.centres);
    for j = 1:mix.ncentres
      v(j) = (post(:,j)'*n2(:,j));
    end
    mix.covars = ((v./new_pr) + sum(diag(prior_cov)))./mix.nin;
    if check_covars
      % Ensure that no covariance is too small
      for j = 1:mix.ncentres
        if mix.covars(j) < MIN_COVAR
          mix.covars(j) = init_covars(j);
        end
      end
    end
  case 'diag'
    for j = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(j,:));
      wts = (post(:,j)*ones(1, mix.nin));
      mix.covars(j,:) = sum((diffs.*diffs).*wts + prior_cov, 1)./new_pr(j);
    end
    if check_covars
      % Ensure that no covariance is too small
      for j = 1:mix.ncentres
        if min(mix.covars(j,:)) < MIN_COVAR
          mix.covars(j,:) = init_covars(j,:);
        end
      end
    end
  case 'full'
    for j = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(j,:));
      diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
      mix.covars(:,:,j) = (diffs'*diffs + prior_cov)/new_pr(j);
    end
    if check_covars
      % Ensure that no covariance is too small
      for j = 1:mix.ncentres
        if min(svd(mix.covars(:,:,j))) < MIN_COVAR
          mix.covars(:,:,j) = init_covars(:,:,j);
        end
      end
    end
  case 'ppca'
    for j = 1:mix.ncentres
      diffs = x - (ones(ndata, 1) * mix.centres(j,:));
      diffs = diffs.*(sqrt(post(:,j))*ones(1, mix.nin));
      [mix.covars(j), mix.U(:,:,j), mix.lambda(j,:)] = ...
        ppca((diffs'*diffs)/new_pr(j), mix.ppca_dim);
    end
    if check_covars
      if mix.covars(j) < MIN_COVAR
        mix.covars(j) = init_covars(j);
      end
    end
    otherwise
      error(['Unknown covariance type ', mix.covar_type]);               
  end
end

ll = sum(log(gmmprob(mix, x)));
num_iter = n;

%if (display >= 0)
%  disp('Warning: Maximum number of iterations has been exceeded');
%end
  
