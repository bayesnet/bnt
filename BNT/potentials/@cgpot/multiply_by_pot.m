function bigpot = multiply_by_pot(bigpot, smallpot, varargin)
% MULTIPLY_BY_POT bigpot *= smallpot for cgpot
% bigpot = multiply_by_pot(bigpot, smallpot)
%
% smallpot's domain must be a subset of bigpot's domain.

bigpot = cg_mom_to_can(bigpot);
smallpot = cg_mom_to_can(smallpot);

mask = find_equiv_posns(smallpot.ddom, bigpot.ddom);
for i=1:bigpot.dsize
  if isempty(smallpot.ddom)
    src = 1;
  else
    sub = ind2subv(bigpot.dsizes, i);
    src = subv2ind(smallpot.dsizes, sub(mask));
  end
  bigpot.can{i} = multiply_by_pot(bigpot.can{i}, smallpot.can{src});
end                   
