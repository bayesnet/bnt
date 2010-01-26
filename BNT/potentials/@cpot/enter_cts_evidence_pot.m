function pot = enter_cts_evidence_pot(pot, Y, y)
% function pot = enter_cts_evidence_pot(pot, Y, y) (cpot)

ns = sparse(1, max(pot.domain));
ns(pot.domain) = pot.sizes;

X = mysetdiff(pot.domain, Y);
[hx, hy, KXX, KXY, KYX, KYY] = partition_matrix_vec(pot.h, pot.K, X, Y, ns);
pot.g = pot.g + hy'*y - 0.5*y'*KYY*y;
if ~isempty(X)
  pot.h = hx - KXY*y;
  pot.K = KXX;
end

pot.sizes(find_equiv_posns(Y,pot.domain)) = 0;
