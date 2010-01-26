function [mom2, loglik] = normalize_pot(can)
% NORMALIZE_POT Convert the canonical potential Pr(X,E) into moment potential Pr(X|E) and return log Pr(E).
% [mom, loglik] = normalize_pot(can)

mom = cpot_to_mpot(can);
mom = struct(mom); % violate privacy of object
loglik = mom.logp;
%mom.logp = 0; % now represents Pr(X | E) instead of Pr(X, E). 
mom2 = mpot(mom.domain, mom.sizes, 0, mom.mu, mom.Sigma);
