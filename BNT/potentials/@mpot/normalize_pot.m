function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the moment potential Pr(X,E) into Pr(X|E) and return log Pr(E).
% [pot, loglik] = normalize_pot(pot)

loglik = pot.logp;
pot.logp = 0;
