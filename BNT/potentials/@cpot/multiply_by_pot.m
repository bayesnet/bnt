function bigpot = multiply_by_pot(bigpot, smallpot, varargin)
% MULTIPLY_BY_POT bigpot *= smallpot for cpot
% bigpot = multiply_by_pot(bigpot, smallpot)
%
% smallpot's domain must be a subset of bigpot's domain.

bigpot.g = bigpot.g + smallpot.g;
if sum(smallpot.sizes) > 0
  mask = find_equiv_posns(smallpot.domain, bigpot.domain);
  u = block(mask, bigpot.sizes);
  bigpot.h(u) = bigpot.h(u) + smallpot.h;
  bigpot.K(u, u) = bigpot.K(u, u) + smallpot.K;
end               
