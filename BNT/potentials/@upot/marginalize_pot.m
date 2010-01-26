function smallpot = marginalize_pot(bigpot, onto, maximize)
% MARGINALIZE_POT Marginalize a upot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, onto, maximize)
%
% The maximize argument is ignored

numer = marg_table(bigpot.p .* bigpot.u, bigpot.domain, bigpot.sizes, onto);
denom = marg_table(bigpot.p, bigpot.domain, bigpot.sizes, onto);

p = denom;
% replace 0s by 1s before dividing. This is valid since demon(i) = 0 => numer(i) = 0
denom = denom + (denom == 0); 
u = numer ./ denom;

ns = zeros(1, max(bigpot.domain));
ns(bigpot.domain) = bigpot.sizes;

smallpot = upot(onto, ns(onto), p, u);
