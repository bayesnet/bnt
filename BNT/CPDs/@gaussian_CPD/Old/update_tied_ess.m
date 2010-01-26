function CPD = update_tied_ess(CPD, domain, engine, evidence, ns, cnodes)

if ~adjustable_CPD(CPD), return; end
nCPDs = size(domain, 2);
fmarginal = cell(1, nCPDs);
for l=1:nCPDs
  fmarginal{l} = marginal_family(engine, nodes(l));
end

[ss cpsz dpsz] = size(CPD.weights);
if const_evidence_pattern(engine)
  dom = domain(:,1);
  dnodes = mysetdiff(1:length(ns), cnodes);
  ddom = myintersect(dom, dnodes);
  cdom = myintersect(dom, cnodes);
  odom = dom(~isemptycell(evidence(dom)));
  hdom = dom(isemptycell(evidence(dom)));
  % If all hidden nodes are discrete and all cts nodes are observed 
  % (e.g., HMM with Gaussian output)
  % we can add the observed evidence in parallel
  if mysubset(ddom, hdom) & mysubset(cdom, odom)
    [mu, Sigma, T] = add_cts_ev_to_marginals(fmarginal, evidence, ns, cnodes);
  else
    mu = zeros(ss, dpsz, nCPDs);
    Sigma = zeros(ss, ss, dpsz, nCPDs);
    T = zeros(dpsz, nCPDs);
    for l=1:nCPDs
      [mu(:,:,l), Sigma(:,:,:,l), T(:,l)] = add_ev_to_marginals(fmarginal{l}, evidence, ns, cnodes);
    end
  end
end
CPD.nsamples = CPD.nsamples + nCPDs;            


if dpsz == 1 % no discrete parents
  w = 1;
else
  w = fullm.T(:);
end
CPD.Wsum = CPD.Wsum + w;
% Let X be the cts parent (if any), Y be the cts child (self).
xi = 1:cpsz;
yi = (cpsz+1):(cpsz+ss);
for i=1:dpsz
  muY = fullm.mu(yi, i);
  SYY = fullm.Sigma(yi, yi, i);
  CPD.WYsum(:,i) = CPD.WYsum(:,i) + w(i)*muY;
  CPD.WYYsum(:,:,i) = CPD.WYYsum(:,:,i) + w(i)*(SYY + muY*muY'); % E[X Y] = Cov[X,Y] + E[X] E[Y]
  if cpsz > 0
    muX = fullm.mu(xi, i);
    SXX = fullm.Sigma(xi, xi, i);
    SXY = fullm.Sigma(xi, yi, i);
    CPD.WXsum(:,i) = CPD.WXsum(:,i) + w(i)*muX;
    CPD.WXYsum(:,:,i) = CPD.WXYsum(:,:,i) + w(i)*(SXY + muX*muY');
    CPD.WXXsum(:,:,i) = CPD.WXXsum(:,:,i) + w(i)*(SXX + muX*muX');
  end
end                


%%%%%%%%%%%%%

function fullm = add_evidence_to_marginal(fmarginal, evidence, ns, cnodes)


dom = fmarginal.domain;

% Find out which values of the discrete parents (if any) are compatible with 
% the discrete evidence (if any).
dnodes = mysetdiff(1:length(ns), cnodes);
ddom = myintersect(dom, dnodes);
cdom = myintersect(dom, cnodes);
odom = dom(~isemptycell(evidence(dom)));
hdom = dom(isemptycell(evidence(dom)));

dobs = myintersect(ddom, odom);
dvals = cat(1, evidence{dobs});
ens = ns; % effective node sizes
ens(dobs) = 1;
S = prod(ens(ddom));
subs = ind2subv(ens(ddom), 1:S);
mask = find_equiv_posns(dobs, ddom);
subs(mask) = dvals;
supportedQs = subv2ind(ns(ddom), subs);

if isempty(ddom)
  Qarity = 1;
else
  Qarity = prod(ns(ddom));
end
fullm.T = zeros(Qarity, 1);
fullm.T(supportedQs) = fmarginal.T(:);

% Now put the hidden cts parts into their right blocks,
% leaving the observed cts parts as 0.
cobs = myintersect(cdom, odom);
chid = myintersect(cdom, hdom);
cvals = cat(1, evidence{cobs});
n = sum(ns(cdom));
fullm.mu = zeros(n,Qarity);
fullm.Sigma = zeros(n,n,Qarity);

if ~isempty(chid)
  chid_blocks = block(find_equiv_posns(chid, cdom), ns(cdom));
end
if ~isempty(cobs)
  cobs_blocks = block(find_equiv_posns(cobs, cdom), ns(cdom));
end

for i=1:length(supportedQs)
  Q = supportedQs(i);
  if ~isempty(chid)
    fullm.mu(chid_blocks, Q) = fmarginal.mu(:, i);
    fullm.Sigma(chid_blocks, chid_blocks, Q) = fmarginal.Sigma(:,:,i);
  end
  if ~isempty(cobs)
    fullm.mu(cobs_blocks, Q) = cvals(:);
  end
end
