function [mu, con, B] = clvm_Mstep(w, Y, YY, YTY, X, XX, XY, varargin)
% MSTEP_CLG Compute ML/MAP estimates for a conditional linear Gaussian
% [mu, Sigma, B] = Mstep_clg(w, Y, YY, YTY, X, XX, XY, varargin)
%
% We fit P(Y|X,Q=i) = N(Y; B_i X + mu_i, Sigma_i) 
% and w(i,t) = p(M(t)=i|y(t)) = posterior responsibility
% See www.ai.mit.edu/~murphyk/Papers/learncg.pdf.
%
% See process_options for how to specify the input arguments.
%
% INPUTS:
% w(i) = sum_t w(i,t) = responsibilities for each mixture component
%  If there is only one mixture component (i.e., Q does not exist),
%  then w(i) = N = nsamples,  and 
%  all references to i can be replaced by 1.
% Y(:,i) = sum_t w(i,t) y(:,t) = weighted observations
% YY(:,:,i) = sum_t w(i,t) y(:,t) y(:,t)' = weighted outer product
% YTY(i) = sum_t w(i,t) y(:,t)' y(:,t) = weighted inner product
%   You only need to pass in YTY if Sigma is to be estimated as spherical.
%
% In the regression context, we must also pass in the following
% X(:,i) = sum_t w(i,t) x(:,t) = weighted inputs
% XX(:,:,i) = sum_t w(i,t) x(:,t) x(:,t)' = weighted outer product
% XY(i) = sum_t w(i,t) x(:,t) y(:,t)' = weighted outer product
%
% Optional inputs (default values in [])
%
% 'cov_type' - 'full', 'diag' or 'spherical' ['full']
% 'tied_cov' - 1 (Sigma) or 0 (Sigma_i) [0]
% 'clamped_cov' - pass in clamped value, or [] if unclamped [ [] ]
% 'clamped_mean' - pass in clamped value, or [] if unclamped [ [] ]
% 'clamped_weights' - pass in clamped value, or [] if unclamped [ [] ]
% 'cov_prior' - added to Sigma(:,:,i) to ensure psd [0.01*eye(d,d,Q)]
%
% If cov is tied, Sigma has size d*d.
% But diagonal and spherical covariances are represented in full size.

[cov_type, tied_cov, ...
 clamped_cov, clamped_mean, clamped_weights,  cov_prior, ...
 xs, ys, post] = ...
    process_options(varargin, ...
		    'cov_type', 'full', 'tied_cov', 0,  'clamped_cov', [], 'clamped_mean', [], ...
		    'clamped_weights', [], 'cov_prior', [], ...
		    'xs', [], 'ys', [], 'post', []);

[Ysz Q] = size(Y);

if isempty(X) % no regression
    B2 = zeros(Ysz, 1, Q);
  for i=1:Q
    B(:,:,i) = B2(:,1:0,i); % make an empty array of size Ysz x 0 x Q
  end
  [mu, con] = mixvonMises_Mstep(w, Y, YY, YTY, varargin{:});
  return;
end


  
