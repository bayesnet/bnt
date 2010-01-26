function Tbig = multiply_by_pot(Tbig, Tsmall)
% MULTIPLY_BY_POT Tbig *= Tsmall
% Tbig = multiply_by_pot(Tbig, Tsmall)
%
% Tsmall's domain must be a subset of Tbig's domain.

smallp = extend_domain_table(Tsmall.p, Tsmall.domain, Tsmall.sizes, Tbig.domain, Tbig.sizes);
Tbig.p = Tbig.p .* smallp;

smallu = extend_domain_table(Tsmall.u, Tsmall.domain, Tsmall.sizes, Tbig.domain, Tbig.sizes);
Tbig.u = Tbig.u + smallu;

