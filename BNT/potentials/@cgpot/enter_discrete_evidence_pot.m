function pot = enter_discrete_evidence_pot(pot, Y, y)

%ns = sparse(1, max(pot.domain));
ns = zeros(1, max(pot.domain));
ns(pot.ddom) = pot.dsizes;
ns(pot.cdom) = pot.csizes;

ddom = pot.ddom;
S = prod(ns(ddom));
sub = ind2subv(ns(ddom), 1:S);
mask = find_equiv_posns(Y, ddom);
sub(mask) = y;
ndx = subv2ind(ns(ddom), sub);

pot.can = pot.can(ndx);
pot.mom = pot.mom(ndx);
