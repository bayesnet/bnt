function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the SCG potential Pr(X,E) into Pr(X|E) and return log Pr(E).
% [pot, loglik] = normalize_pot(pot)

% Marginalize down to [], so that the normalizing constant becomes Pr(E)
temp = marginalize_pot(pot, []);
[temp2, loglik] = normalize_pot(temp.scgpotc{1});
  
% Adjust scale factor to reflect the fact that the pot now represents Pr(X | E) instead of Pr(X,E).

scale = -loglik;
if 1
    for i=1:pot.dsize
        pot.scgpotc{i} = rescale_pot( pot.scgpotc{i}, scale);
    end
end
