function Tbig = divide_by_pot(Tbig, Tsmall)
% DIVIDE_BY_POT Tbig /= Tsmall
% Tbig = divide_by_pot(Tbig, Tsmall)
%
% Tsmall's domain must be a subset of Tbig's domain.

%process sparse dpot, we do not concern only one of the two pots is sparse
if issparse(Tbig.T) & issparse(Tsmall.T)
   Tbig.T = divide_by_sparse_table(Tbig.T, Tbig.domain, Tbig.sizes, Tsmall.T, Tsmall.domain, Tsmall.sizes);
else
   Tbig.T = divide_by_table(Tbig.T, Tbig.domain, Tbig.sizes, Tsmall.T, Tsmall.domain, Tsmall.sizes);
end


