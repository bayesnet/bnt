function fullm = add_evidence_to_gmarginal(fmarginal, evidence, ns, cnodes)
% ADD_EVIDENCE_TO_GMARGINAL 'pump up' observed nodes back to their original size.
% function fullm = add_evidence_to_gmarginal(fmarginal, evidence, ns, cnodes)
%
% We introduce 0s into the array in positions which are incompatible with the evidence.
% for both discrete and continuous nodes.
%
% See also add_ev_to_dmarginal

dom = fmarginal.domain;
fullm.domain = fmarginal.domain;

% Find out which values of the discrete parents (if any) are compatible with 
% the discrete evidence (if any).
dnodes = mysetdiff(1:length(ns), cnodes);
ddom = myintersect(dom, dnodes);
cdom = myintersect(dom, cnodes);
odom = dom(~isemptycell(evidence(dom)));
hdom = dom(isemptycell(evidence(dom)));

% Find the entries in the big table that are compatible with the discrete evidence.
% (We will put the probabilities from the small inferred table into these positions.)
% We could use add_ev_to_dmarginal to do this.
dobs = myintersect(ddom, odom);
dvals = cat(1, evidence{dobs});
ens = ns; % effective node sizes
ens(dobs) = 1;
S = prod(ens(ddom));
subs = ind2subv(ens(ddom), 1:S);
mask = find_equiv_posns(dobs, ddom);
%subs(mask) = dvals; % bug fix by P. Brutti
for i=1:length(mask),
  subs(:,mask(i)) = dvals(i);
end       
supportedQs = subv2ind(ns(ddom), subs);

if isempty(ddom)
  Qarity = 1;
else
  Qarity = prod(ns(ddom));
end
fullm.T = zeros(Qarity, 1);
fullm.T(supportedQs) = fmarginal.T(:);
fullm.T = myreshape(fullm.T, ns(ddom));


if isempty(cdom)
  fullm.mu = [];
  fullm.sigma = [];
  return;
end

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
