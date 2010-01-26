function Tbig = divide_by_pot(Tbig, Tsmall)
% DIVIDE_BY_POT Tbig /= Tsmall
% Tbig = divide_by_pot(Tbig, Tsmall)
%
% Tsmall's domain must be a subset of Tbig's domain.

smallp = extend_domain_table(Tsmall.p, Tsmall.domain, Tsmall.sizes, Tbig.domain, Tbig.sizes);
smallp = smallp + (smallp==0);
Tbig.p = Tbig.p ./ smallp;

smallu = extend_domain_table(Tsmall.u, Tsmall.domain, Tsmall.sizes, Tbig.domain, Tbig.sizes);
Tbig.u = Tbig.u - smallu;

