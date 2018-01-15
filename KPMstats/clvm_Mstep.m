function [mu, con, B] = clvm_Mstep(w, Y, YY, YTY, X, XX, XY, varargin)
% MSTEP_CLG Compute ML/MAP estimates for a conditional linear Von Mises
% [mu, Con, B] = clvm_Mstep(w, Y, YY, YTY, X, XX, XY, varargin)

% This currently only accounts for conditions where there are no cts parents. 
% Dsc-->Cts is the structure that this can handle. 

% INPUTS:
% w(i) = sum_t w(i,t) = responsibilities for each mixture component
%  If there is only one mixture component (i.e., Q does not exist),
%  then w(i) = N = nsamples,  and 
%  all references to i can be replaced by 1.
% Y(:,i) = sum_t w(i,t) cos(y(:,t)) = weighted observations
% YY(:,:,i) = sum_t w(i,t) sin(y(:,t)) = weighted outer product

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


  
