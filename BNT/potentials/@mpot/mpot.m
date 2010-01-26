function pot = mpot(members, sizes, logp, mu, Sigma)
% MPOT Make a moment Gaussian potential.
% pot = mpot(members, sizes, logp, mu, Sigma)
%
% All params default to 0 if omitted.

n = sum(sizes);
if nargin < 3, logp = 0; end
if nargin < 4, mu = zeros(n,1); end
if nargin < 5, Sigma = zeros(n,n); end
  
pot.domain = members;
pot.sizes = sizes;
pot.logp = logp;
pot.mu = mu;
pot.Sigma = Sigma;zeros(n,n);      
pot = class(pot, 'mpot');
