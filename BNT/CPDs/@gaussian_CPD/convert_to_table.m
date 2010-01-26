function T = convert_to_table(CPD, domain, evidence)
% CONVERT_TO_TABLE Convert a Gaussian CPD to a table
% T = convert_to_table(CPD, domain, evidence)


sz = CPD.sizes;
ns = zeros(1, max(domain));
ns(domain) = sz;

odom = domain(~isemptycell(evidence(domain)));
ps = domain(1:end-1);
cps = ps(CPD.cps);
dps = ps(CPD.dps);
self = domain(end);
cdom = [cps(:)' self];
ddom = dps;
cnodes = cdom;

[m, C, W] = gaussian_CPD_params_given_dps(CPD, domain, evidence);


ns(odom) = 1;
dpsize = prod(ns(dps));
self = domain(end);
assert(myismember(self, odom));
self_val = evidence{self};
T = zeros(dpsize, 1);
if length(cps) > 0 
  assert(~any(isemptycell(evidence(cps))));
  cps_vals = cat(1, evidence{cps});
  for i=1:dpsize
    T(i) = gaussian_prob(self_val, m(:,i) + W(:,:,i)*cps_vals, C(:,:,i));
  end
else
  for i=1:dpsize
    T(i) = gaussian_prob(self_val, m(:,i), C(:,:,i));
  end
end
