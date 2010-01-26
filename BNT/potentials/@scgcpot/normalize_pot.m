function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the element of stable conditional gaussian potential Pr(X,E) into Pr(X|E) and return log Pr(E).
% [pot, loglik] = normalize_pot(pot)

loglik = log(pot.p);
pot.p = 1;
