function smallpot = marginalize_pot(bigpot, keep)
% MARGINALIZE_POT Marginalize a cgpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, keep)

sumover = mysetdiff(bigpot.domain, keep);
csumover = myintersect(sumover, bigpot.cdom);
dsumover = myintersect(sumover, bigpot.ddom);
dkeep = myintersect(keep, bigpot.ddom);
ckeep = myintersect(keep, bigpot.cdom);
%ns = sparse(1, max(bigpot.domain)); % must be full, so I is an integer
ns = zeros(1, max(bigpot.domain));
ns(bigpot.ddom) = bigpot.dsizes;
ns(bigpot.cdom) = bigpot.csizes;

% sum(ns(csumover))==0 is like isempty(csumover) but handles observed nodes.
% Similarly, prod(ns(dsumover))==1 is like isempty(dsumover)

% Marginalize the cts parts.
% If we are in canonical form, we stay that way, since moment form might not exist.
% Besides, we would like to minimize the number of conversions.
if sum(ns(csumover)) > 0
  if bigpot.subtype == 'm'
    for i=1:bigpot.dsize
      bigpot.mom{i} = marginalize_pot(bigpot.mom{i}, ckeep);
    end
  else
    for i=1:bigpot.dsize
      bigpot.can{i} = marginalize_pot(bigpot.can{i}, ckeep);
    end
  end
end

% If we are not marginalizing over any discrete nodes, we are done.
if prod(ns(dsumover))==1
  smallpot = cgpot(dkeep, ckeep, ns, bigpot.can, bigpot.mom, bigpot.subtype);
  return;
end

% To marginalize the discrete parts, we must be in moment form.
bigpot = cg_can_to_mom(bigpot);

I = prod(ns(dkeep));
J = prod(ns(dsumover));
C = sum(ns(ckeep));

% Reshape bigpot into the form mu1(:,j,i), where i is in dkeep, j is in dsumover
T1 = zeros(I,J);
mu1 = zeros(C,J,I);
Sigma1 = zeros(C,C,J,I);
sum_map = find_equiv_posns(dsumover, bigpot.ddom);
keep_map = find_equiv_posns(dkeep, bigpot.ddom);
iv = zeros(1, length(bigpot.ddom)); % index vector
for i=1:I
  keep_iv = ind2subv(ns(dkeep), i);
  iv(keep_map) = keep_iv;
  for j=1:J
    sum_iv = ind2subv(ns(dsumover), j);
    iv(sum_map) = sum_iv;
    k = subv2ind(ns(bigpot.ddom), iv);
    mom = struct(bigpot.mom{k}); % violate object privacy
    T1(i,j) = exp(mom.logp);
    if C > 0 % so mu1 and Sigma1 are non-empty
      mu1(:,j,i) = mom.mu;
      Sigma1(:,:,j,i) = mom.Sigma;
    end
  end
end

% Collapse the mixture of Gaussians
coef = mk_stochastic(T1); % coef must be convex combination
T2 = sum(T1,2);
T2 = T2 + (T2==0)*eps;
%if C > 0, disp('collapsing onto '); disp(leep); end
mu = [];
Sigma = [];
mom = cell(1,I);
for i=1:I
  if C > 0
    [mu, Sigma] = collapse_mog(mu1(:,:,i), Sigma1(:,:,:,i), coef(i,:));
  end
  logp = log(T2(i));
  mom{i} = mpot(ckeep, ns(ckeep), logp, mu, Sigma);
end

smallpot = cgpot(dkeep, ckeep, ns, [], mom, 'm');

