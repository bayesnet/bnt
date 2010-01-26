function T = multiply_pots(T1, T2)
% MULTIPLY_POTS Multiply a pair of dpots together pointwise (cgpot)
% T = multiply_pots(pots)

ddom = myunion(T1.ddom, T2.ddom);
cdom = myunion(T1.cdom, T2.cdom);
dom = myunion(ddom, cdom);
ns = zeros(1, max(dom));
ns(T1.ddom) = T1.dsizes;
ns(T2.ddom) = T2.dsizes;
ns(T1.cdom) = T1.csizes;
ns(T2.cdom) = T2.csizes;

T = cgpot(ddom, cdom, ns);
T = multiply_by_pot(T, T1);
T = multiply_by_pot(T, T2);   
