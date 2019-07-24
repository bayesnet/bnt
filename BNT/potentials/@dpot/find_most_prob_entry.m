function [indices, pot] = find_most_prob_entry(pot)
% function [indices, pot] = find_most_prob_entry(pot)
% function [indices, pot] = find_most_prob_entry(pot)
% Find the indices of the argmax, and set all other enties to 0.

%indices = argmax(pot.T);
[m i] = max(pot.T(:));
indices = ind2subv(pot.sizes, i);
pot.T = 0*myones(pot.sizes);
pot.T(i) = m;
