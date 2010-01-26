function smallpot = marginalize_pot(bigpot, keep)
% MARGINALIZE_POT Marginalize a cgpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, keep)

sumover = mysetdiff(bigpot.domain, keep);
cdom = myunion(bigpot.cheaddom, bigpot.ctaildom);
csumover = myintersect(sumover, bigpot.cheaddom);
dsumover = myintersect(sumover, bigpot.ddom);

dkeep = myintersect(keep, bigpot.ddom);
ckeep = myintersect(keep, bigpot.cheaddom);
cheaddom = myintersect(keep, bigpot.cheaddom);

assert(isempty(myintersect(csumover,bigpot.ctaildom)));
ns = zeros(1, max(bigpot.domain));
ns(bigpot.ddom) = bigpot.dsizes;
ns(bigpot.cheaddom) = bigpot.cheadsizes;
ns(bigpot.ctaildom) = bigpot.ctailsizes;


if sum(ns(csumover)) > 0
    for i=1:bigpot.dsize
      bigpot.scgpotc{i} = marginalize_pot(bigpot.scgpotc{i}, ckeep, csumover, ns);
    end
end

if (isequal(csumover, cheaddom))
    bigpot.ctaildom = [];
end
% If we are not marginalizing over any discrete nodes, we are done.
if prod(ns(dsumover))==1
  smallpot = scgpot(dkeep, cheaddom, bigpot.ctaildom, ns, bigpot.scgpotc);
  return;
end

if (~isempty(bigpot.ctaildom))
    assert(0);
    return;
end

I = prod(ns(dkeep));
J = prod(ns(dsumover));
C = sum(ns(ckeep));   
sum_map = find_equiv_posns(dsumover, bigpot.ddom);
keep_map = find_equiv_posns(dkeep, bigpot.ddom);
iv = zeros(1, length(bigpot.ddom)); % index vector

p1 = zeros(I,J);
A1 = zeros(C,J,I);
C1 = zeros(C,C,J,I);
for i=1:I
  keep_iv = ind2subv(ns(dkeep), i);
  iv(keep_map) = keep_iv;
  for j=1:J
    sum_iv = ind2subv(ns(dsumover), j);
    iv(sum_map) = sum_iv;
    k = subv2ind(ns(bigpot.ddom), iv);
    pot = struct(bigpot.scgpotc{k}); % violate object privacy
    p1(i,j) = pot.p;
    if C > 0 % so mu1 and Sigma1 are non-empty
      A1(:,j,i) = pot.A;
      C1(:,:,j,i) = pot.C;
    end
  end
end

% Collapse the mixture of Gaussians
coef = mk_stochastic(p1); % coef must be convex combination
%keyboard
p2 = sum(p1,2);
if (all(p2 == 0))
    p2 = p2 + (p2==0)*eps;
end
A = [];
S = [];

pot = cell(1,I);
ctailsize = sum(ns(bigpot.ctaildom));
tB = zeros(C, ctailsize);
for i=1:I
  if C > 0
    [A, S] = collapse_mog(A1(:,:,i), C1(:,:,:,i), coef(i,:));
  end
  p = p2(i);
  pot{i} = scgcpot(C, ctailsize, p, A, tB, S);
end

smallpot = scgpot(dkeep, ckeep, bigpot.ctaildom, ns, pot);




