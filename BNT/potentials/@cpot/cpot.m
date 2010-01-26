function pot = cpot(members, sizes, g, h, K)
% CPOT Make a canonical Gaussian potential.
% pot = cpot(members, sizes, g, h, K)
%
% All params default to 0 if omitted.

n = sum(sizes);
if nargin < 3, g = 0; end
if nargin < 4, h = zeros(n,1); end
if nargin < 5, K = zeros(n,n); end
  
pot.domain = members;
pot.sizes = sizes(:)';
pot.g = g;
pot.h = h;
pot.K = K;
pot = class(pot, 'cpot');
