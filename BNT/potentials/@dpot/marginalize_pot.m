function smallpot = marginalize_pot(bigpot, onto, maximize)
% MARGINALIZE_POT Marginalize a dpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, onto, maximize)
%
% 'onto' must be in ascending order.

if nargin < 3, maximize = 0; end

ns = zeros(1, max(bigpot.domain));
ns(bigpot.domain) = bigpot.sizes;
%assert(isequal(bigpot.sizes, mysize(bigpot.T))); % may fail if there are trailing dimensions of size 1
if issparse(bigpot.T)
   smallT = marg_sparse_table(bigpot.T, bigpot.domain, bigpot.sizes, onto, maximize);
else 
   smallT = marg_table(bigpot.T, bigpot.domain, bigpot.sizes, onto, maximize);
end
smallpot = dpot(onto, ns(onto), smallT);
