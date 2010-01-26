function [mu, Sigma, W] = CPD_to_linear_gaussian(CPD, domain, ns, cnodes, evidence)

ps = domain(1:end-1);
dnodes = mysetdiff(1:length(ns), cnodes);
dps = myintersect(ps, dnodes); % discrete parents

if isempty(dps)
  Q = 1;
else
  assert(~any(isemptycell(evidence(dps))));
  dpvals = cat(1, evidence{dps});
  Q = subv2ind(ns(dps), dpvals(:)');
end

mu = CPD.mean(:,Q);
Sigma = CPD.cov(:,:,Q);
W = CPD.weights(:,:,Q);


