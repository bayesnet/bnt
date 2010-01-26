function smallpot = marginalize_pot(bigpot, keep, maximize, useC)
% MARGINALIZE_POT Marginalize a mpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, keep, maximize, useC)
%
% The maximize argument is ignored - maxing out a Gaussian is the same as summing it out,
% since the mode and mean are equal.
% The useC argument is ignored.


node_sizes = sparse(1, max(bigpot.domain));
node_sizes(bigpot.domain) = bigpot.sizes;
sum_over = mysetdiff(bigpot.domain, keep);

[logp, mu, Sigma] = marginalize_gaussian(bigpot.logp, bigpot.mu, bigpot.Sigma, ...
					 keep, sum_over, node_sizes);
smallpot = mpot(keep, node_sizes(keep), logp, mu, Sigma);

%%%%%%

function [logpX, muX, SXX] = marginalize_gaussian(logp, mu, Sigma, X, Y, ns)
% MARGINALIZE_GAUSSIAN Compute Pr(X) from Pr(X,Y) where X and Y are jointly Gaussian.
% [logpX, muX, SXX] = marginalize_gaussian(logp, mu, Sigma, X, Y, ns)
%
% sizes(i) is the size of the i'th block in domain.
 
[muX, muY, SXX, SXY, SYX, SYY] = partition_matrix_vec(mu, Sigma, X, Y, ns);
logpX = logp; % Lauritzen (1996) p161          
