function Tbig = multiply_by_pot(Tbig, Tsmall)
% MULTIPLY_BY_POT Tbig *= Tsmall
% Tbig = multiply_by_pot(Tbig, Tsmall)
%
% Tsmall's domain must be a subset of Tbig's domain.

%process sparse dpot, we do not consider only one of the two pots is sparse
if issparse(Tbig.T) & issparse(Tsmall.T)
   Tbig.T = mult_by_sparse_table(Tbig.T, Tbig.domain, Tbig.sizes, Tsmall.T, Tsmall.domain, Tsmall.sizes);
else 
   Tbig.T = mult_by_table(Tbig.T, Tbig.domain, Tbig.sizes, Tsmall.T, Tsmall.domain, Tsmall.sizes);
end

