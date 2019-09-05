function pot = direct_combine_pots(pot1, pot2)
% DIRECTED_COMBINE_POTS The combination operation corresponds to ordinary composition of conditional distributions. 
% In some sense is similar to that of forming disjoint union of set.
% pot = direct_combine_pots(pot1, pot2)

% directed combine can be performed under the conditon that the head node set of pot1 is disjoint from the domain of 
% pot2 or vice versa. if the last conditon was satisfied we exchange the pot1 and pot2 firstly then perform the operation.
% If neither of them was satified the directed combine is undifined.


if isempty( myintersect(pot1.domain, pot2.cheaddom) )
    pot1 = pot1;
    pot2 = pot2;
elseif  isempty( myintersect(pot2.domain, pot1.cheaddom))
    temppot = pot1;
    pot1 = pot2;
    pot2 = temppot;
else
    assert(0);
    return;
end

domain = myunion(pot1.domain, pot2.domain);
nodesizes = zeros(1,max(domain));
nodesizes(pot2.ctaildom) = pot2.ctailsizes;
nodesizes(pot2.cheaddom) = pot2.cheadsizes;
nodesizes(pot2.ddom) = pot2.dsizes;
nodesizes(pot1.ctaildom) = pot1.ctailsizes;
nodesizes(pot1.cheaddom) = pot1.cheadsizes;
nodesizes(pot1.ddom) = pot1.dsizes;

dom_u = mysetdiff(pot2.ctaildom, pot1.cheaddom);
if ~isempty(dom_u) & ~mysubset(dom_u, pot1.ctaildom)
    pot1 = extension_pot(pot1, [], [], dom_u, nodesizes(dom_u));
end

dom_u = myunion(pot1.cheaddom, pot1.ctaildom);
if ~isempty(dom_u) & ~mysubset(dom_u, pot2.ctaildom)
    pot2 = extension_pot(pot2, [], [], dom_u, nodesizes(dom_u));
end


cheaddom = myunion(pot1.cheaddom, pot2.cheaddom);
ctaildom = mysetdiff(myunion(pot1.ctaildom, pot2.ctaildom), cheaddom);
cdom = myunion(cheaddom, ctaildom);
ddom = mysetdiff(domain, cdom);
dsizes = nodesizes(ddom);
dsize = prod(nodesizes(ddom));
cheadsizes = nodesizes(cheaddom);
cheadsize = sum(nodesizes(cheaddom));
ctailsizes = nodesizes(ctaildom);
ctailsize = sum(nodesizes(ctaildom));

r1 = pot1.cheadsize;
s1 = pot1.ctailsize;
scpot = cell(1, dsize);
mask1 = [];
mask2 = [];
if ~isempty(pot1.ddom)
    mask1 = find_equiv_posns(pot1.ddom, ddom);
end
if ~isempty(pot2.ddom)
    mask2 = find_equiv_posns(pot2.ddom, ddom);
end
cmask1 = [];
cmask2 = [];
if ~isempty(pot1.cheaddom)
    cmask1 = find_equiv_posns(pot1.cheaddom, cheaddom);
end
if ~isempty(pot2.cheaddom)
    cmask2 = find_equiv_posns(pot2.cheaddom, cheaddom);
end

u1 = block(cmask1, cheadsizes);
u2 = block(cmask2, cheadsizes);

fmaskh = find_equiv_posns(pot1.cheaddom, pot2.ctaildom);
fmaskt = find_equiv_posns(pot1.ctaildom, pot2.ctaildom);

fh = block(fmaskh, pot2.ctailsizes);
ft = block(fmaskt, pot2.ctailsizes);

for i=1:dsize
    sub = ind2subv(dsizes, i);
    sub1 = sub(mask1);
    sub2 = sub(mask2);
    ind1 = subv2ind(pot1.dsizes, sub1);
    ind2 = subv2ind(pot2.dsizes, sub2);
    
    if isempty(ind1)
        ind1 = 1;
    end
    if isempty(ind2)
        ind2 = 1;
    end
    potc1 = struct(pot1.scgpotc{ind1});
    potc2 = struct(pot2.scgpotc{ind2});
    p = potc1.p;
    q = potc2.p;
    ro = p*q;
    
    A = potc1.A;
    B = potc1.B;
    C = potc1.C;
   
    E = potc2.A;
    F = potc2.B;
    G = potc2.C;
    
    F1 = F(:, fh);
    F2 = F(:, ft);
    
    if ~isempty(F1)
        K1 = F1*A;
        K2 = F1*B;
        FCF = F1*C*F1';
        FC = F1*C;
        CFT = C*F1';
    else
        K1 = zeros(size(E));
        K2 = zeros(size(F2));
        FCF = zeros(size(G));
        FC = zeros(size(C, 1), size(G, 2));
        CFT = zeros(size(G, 2), size(C, 1));
    end
    
    
    U = zeros(cheadsize,1); 
    W = zeros(cheadsize,cheadsize);
    V = zeros(cheadsize,ctailsize); 
    
    if cheadsize > 0
        U(u1) = A;
        U(u2) = E + K1;
        W(u1, u1) = C;
        W(u2, u2) = G + FCF;
        W(u1, u2) = CFT;
        W(u2, u1) = FC;
    else
        U = zeros(cheadsize,1); 
        W = zeros(cheadsize,cheadsize); 
    end
    if cheadsize > 0 | ctailsize > 0
        if ~isempty(u1)
            V(u1, :) = B;
        else
            V(u1, :) = zeros(potc1.cheadsize, ctailsize);
        end
        if ~isempty(u2)
            V(u2, :) = F2 + K2;
        else
            V(u2, :) = zeros(potc2.cheadsize, ctailsize);
        end
    else
        V = zeros(cheadsize,ctailsize); 
    end

    scpot{i} = scgcpot(cheadsize, ctailsize, ro, U, V, W);
end

pot = scgpot(ddom, cheaddom, ctaildom, nodesizes, scpot);
