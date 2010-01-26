function T = multiply_pots(T1, T2)
% MULTIPLY_POTS Multiply a pair of dpots together pointwise (cpot)
% T = multiply_pots(pots)

dom = myunion(T1.domain, T2.domain);
ns = sparse(1, max(dom));
ns(T1.domain) = T1.sizes;
ns(T2.domain) = T2.sizes;
T = cpot(dom, ns(dom));
T = multiply_by_pot(T, T1);
T = multiply_by_pot(T, T2);   
