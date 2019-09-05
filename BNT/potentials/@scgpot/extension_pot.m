function pot = extension_pot(oldpot, ddom_u, dsizes, ctaildom_u, csizes)
% EXTENSION_POT Extense a stable CG potential.
% pot = extension_pot(oldpot, ddom_u, ctaildom_u, dsizes, csizes)
% ddom_u Added discrete nodes
% ctaildom_u Added continuous tail nodes
% csizes is the size of the tail nodes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A CG potential can be extended by adding discrete variables to its %
% domain of continuous variables to its tail                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddom = myunion(oldpot.ddom, ddom_u);
ctaildom = myunion(oldpot.ctaildom, ctaildom_u);
cheaddom = oldpot.cheaddom;
udom = myunion(ddom_u, ctaildom_u);
domain = myunion(oldpot.domain, udom);

ns = zeros(1,max(domain));
ns(ddom_u) = dsizes;
ns(ctaildom_u) = csizes;
ns(oldpot.ddom) = oldpot.dsizes;
ns(oldpot.cheaddom) = oldpot.cheadsizes;
ns(oldpot.ctaildom) = oldpot.ctailsizes;

dsizes = ns(ddom);
dsize = prod(ns(ddom));
cheadsizes = ns(cheaddom);
cheadsize = sum(ns(cheaddom));
ctailsizes = ns(ctaildom);
ctailsize = sum(ns(ctaildom));

BZ = zeros(cheadsize, ctailsize);
potarray = cell(1, dsize);
mask = find_equiv_posns(oldpot.ddom, ddom);

tmask = find_equiv_posns(oldpot.ctaildom, ctaildom);
tu = block(tmask, ctailsizes);

for i=1:dsize
    sub1 = ind2subv(dsizes, i);
    sub2 = sub1(mask);
    ind = subv2ind(oldpot.dsizes, sub2);
    if isempty(ind)
        ind = 1;
    end
    potc = struct(oldpot.scgpotc{ind});
    p = potc.p;
    B = BZ;
    if ~isempty(B)
        B(:, tu) = potc.B;
    end
    potarray{i} = scgcpot(cheadsize, ctailsize, p, potc.A, B, potc.C);
end

pot = scgpot(ddom, cheaddom, ctaildom, ns,potarray);
