function smallpot = marginalize_pot(bigpot, keep, maximize, useC)
% MARGINALIZE_POT Marginalize a cgpot onto a smaller domain.
% smallpot = marginalize_pot(bigpot, keep, maximize, useC)
%
% If maximize = 1, we raise an error.
% useC is ignored.

if nargin < 3, maximize = 0; end
assert(~maximize);


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

% To marginalize the discrete parts, we partition the cts parts into those that depend
% on dkeep (i) and those that depend on on dsumover (j).

I = prod(ns(dkeep));
J = prod(ns(dsumover));
C = sum(ns(ckeep));   
sum_map = find_equiv_posns(dsumover, bigpot.ddom);
keep_map = find_equiv_posns(dkeep, bigpot.ddom);
iv = zeros(1, length(bigpot.ddom)); % index vector

% If in canonical form, marginalize if possible, else convert to moment form.
if 0 & bigpot.subtype == 'c'
  p1 = zeros(I,J);
  h1 = zeros(C,J,I);
  K1 = zeros(C,C,J,I);
  for i=1:I
    keep_iv = ind2subv(ns(dkeep), i);
    iv(keep_map) = keep_iv;
    for j=1:J
      sum_iv = ind2subv(ns(dsumover), j);
      iv(sum_map) = sum_iv;
      k = subv2ind(ns(bigpot.ddom), iv);
      can = struct(bigpot.can{k}); % violate object privacy
      p1(i,j) = exp(can.g);
      if C > 0 % so mu1 and Sigma1 are non-empty
	h1(:,j,i) = can.h;
	K1(:,:,j,i) = can.K;
      end
    end
  end
  
  % If the cts parts do not depend on j, we can just marginalize the weighting coefficient g.
  jdepends = 0;
  for i=1:I
    for j=2:J
      if ~approxeq(h1(:,j,i), h1(:,1,i)) | ~approxeq(K1(:,:,j,i), K1(:,:,1,i))
	jdepends = 1;
	break
      end
    end
  end

  if ~jdepends
    %g2 = log(sum(p1, 2));
    g2 = zeros(I,1);
    for i=1:I
      s = sum(p1(i,:));
      if s > 0
	g2(i) = log(s);
      end
    end
    h2 = h1;
    K2 = K1;
    can = cell(1,I);
    j = 1; % arbitrary
    for i=1:I
      can{i} = cpot(ckeep, ns(ckeep), g2(i), h2(:,j,i), K2(:,:,j,i));
    end
    smallpot = cgpot(dkeep, ckeep, ns, can, [], 'c');  
    return;
  else
    % Since the cts parts depend on j, we must convert to moment form
    bigpot = cg_can_to_mom(bigpot);
  end
end


% Marginalize in moment form
bigpot = cg_can_to_mom(bigpot);

% Now partition the moment components.
T1 = zeros(I,J);
mu1 = zeros(C,J,I);
Sigma1 = zeros(C,C,J,I);
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

