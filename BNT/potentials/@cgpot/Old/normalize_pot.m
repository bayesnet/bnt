function [pot, loglik] = normalize_pot(pot)
% NORMALIZE_POT Convert the CG potential Pr(X,E) into Pr(X|E) and return log Pr(E).
% [pot, loglik] = normalize_pot(pot)

% Marginalize down to [], so that the normalizing constant becomes Pr(E)
temp = marginalize_pot(cg_can_to_mom(pot), []); 
%loglik = temp.mom{1}.logp;
[temp2, loglik] = normalize_pot(temp.mom{1});
  
% Adjust scale factor to reflect the fact that the pot now represents Pr(X | E) instead of Pr(X,E).

scale = -loglik;
if 1
switch pot.subtype
  case 'm'
    for i=1:pot.dsize
      pot.mom{i} = rescale_pot(pot.mom{i}, scale);
    end
  case 'c'
    for i=1:pot.dsize
      pot.can{i} = rescale_pot(pot.can{i}, scale);
    end
end        
end
