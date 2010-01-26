function pot = CPD_to_scgpot(CPD, domain, ns, cnodes, evidence)
% CPD_TO_CGPOT Convert a Gaussian CPD to a CG potential, incorporating any evidence   
% pot = CPD_to_cgpot(CPD, domain, ns, cnodes, evidence)

self = CPD.self;
dnodes = mysetdiff(1:length(ns), cnodes);
odom = domain(~isemptycell(evidence(domain)));
cdom = myintersect(cnodes, domain);
cheaddom = myintersect(self, domain);
ctaildom = mysetdiff(cdom,cheaddom);
ddom = myintersect(dnodes, domain);
cobs = myintersect(cdom, odom);
dobs = myintersect(ddom, odom);
ens = ns; % effective node size
ens(cobs) = 0;
ens(dobs) = 1;

% Extract the params compatible with the observations (if any) on the discrete parents (if any)
% parents are all but the last domain element
ps = domain(1:end-1);
dps = myintersect(ps, ddom);
dops = myintersect(dps, odom);

map = find_equiv_posns(dops, dps);
dpvals = cat(1, evidence{dops});
index = mk_multi_index(length(dps), map, dpvals);

dpsize = prod(ens(dps));
cpsize = size(CPD.weights(:,:,1), 2); % cts parents size
ss = size(CPD.mean, 1); % self size
% the reshape acts like a squeeze
m = reshape(CPD.mean(:, index{:}), [ss dpsize]);
C = reshape(CPD.cov(:, :, index{:}), [ss ss dpsize]);
W = reshape(CPD.weights(:, :, index{:}), [ss cpsize dpsize]);


% Convert each conditional Gaussian to a canonical potential
pot = cell(1, dpsize);
for i=1:dpsize
  %pot{i} = linear_gaussian_to_scgcpot(m(:,i), C(:,:,i), W(:,:,i), cdom, ns, cnodes, evidence);
  pot{i} = scgcpot(ss, cpsize, 1, m(:,i), W(:,:,i), C(:,:,i));
end

pot = scgpot(ddom, cheaddom, ctaildom, ens, pot);


function pot = linear_gaussian_to_scgcpot(mu, Sigma, W, domain, ns, cnodes, evidence)
% LINEAR_GAUSSIAN_TO_CPOT Convert a linear Gaussian CPD  to a stable conditional potential element.
% pot = linear_gaussian_to_cpot(mu, Sigma, W, domain, ns, cnodes, evidence)

p = 1;
A = mu;
B = W;
C = Sigma;
ns(odom) = 0;
%pot = scgcpot(, ns(domain), p, A, B, C);


