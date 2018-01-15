function p = vonMises_prob( x, m, k, use_log )
%VONMIS_PROB Calculates the probability of x coming from a Von Mises
%distribution with mean mu and concentration parameter k. 
%   p = vonMises_prob( x, m, k, use_log )

if nargin < 4, use_log = 0; end
[d N] = size(x);

m = m(:);
M = m*ones(1,N);

denom = (2*pi)*besseli(0,k);
mahal = (k*cos(x-M)); 

if use_log
    p = mahal-log(denom);
else
    p = exp(mahal)/denom;
end

