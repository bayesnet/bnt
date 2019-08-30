function bigpot = divide_by_pot(bigpot, smallpot)
% DIVIDE_BY_POT bigpot /= smallpot for cpot
% bigpot = divide_by_pot(bigpot, smallpot)
%
% smallpot's domain must be a subset of bigpot's domain.
if bigpot.g ==-Inf && smallpot.g == -Inf %bug fix by Sacha Gunaratne 2017-10-23
    bigpot.g = 0;
else
    bigpot.g = bigpot.g - smallpot.g;
end
if sum(smallpot.sizes) > 0
  mask = find_equiv_posns(smallpot.domain, bigpot.domain);
  u = block(mask, bigpot.sizes);
  bigpot.h(u) = bigpot.h(u) - smallpot.h;
  bigpot.K(u, u) = bigpot.K(u, u) - smallpot.K;
end               
