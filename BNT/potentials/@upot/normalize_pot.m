function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the probability part of a utility potential
% [pot, loglik] = normalize_pot(pot)

[pot.p, lik] = normalise(pot.p);
%pot.u = pot.u - sum(pot.u(:));
%pot.u = pot.u ./ sum(pot.u(:)); % same as normalise(pot.u)
%pot.u = normalise(pot.u);
%pot.u = pot.u / 726.8121;
pot.u = pot.u / 10;
loglik = log(lik + (lik==0)*eps);

      
