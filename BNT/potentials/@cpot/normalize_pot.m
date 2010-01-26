function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the canonical potential Pr(X,E) into Pr(X|E) and return log Pr(E).
% [pot, loglik] = normalize_pot(pot)

mom = cpot_to_mpot(pot);  % move the normalizing constant out of g, to reveal the coefficient          
%loglik = scaling_factor_pot(mom);
%loglik = mom.logp; 
[temp, loglik] = normalize_pot(mom);
pot.g = pot.g - loglik;

