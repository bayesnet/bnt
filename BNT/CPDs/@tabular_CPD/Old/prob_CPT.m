function p = prob_CPT(CPD, x)
% PROB_CPT Lookup the prob. of a family value in a tabular CPD
% p = prob_CPT(CPD, x)
%
% This is a version of prob_CPD optimized for tables.

switch length(x)
 case 1, p = CPD.CPT(x);
 case 2, p = CPD.CPT(x(1), x(2));
 case 3, p = CPD.CPT(x(1), x(2), x(3));
 case 4, p = CPD.CPT(x(1), x(2), x(3), x(4));
 otherwise,
  ind = subv2ind(mysize(CPD.CPT), x);
  p = CPD.CPT(ind);
end

