function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the discrete potential Pr(X,E) into Pr(X|E) and return log Pr(E).
% [pot, loglik] = normalize_pot(pot)

if isempty(pot.T)  %add to process sparse
   loglik = 0;
   return;
end
[pot.T, lik] = normalise(pot.T);
loglik = log(lik + (lik==0)*eps);

      
