function T = convert_to_table(CPD, domain, evidence)
%CONVERT_TO_TABLE Convert a Von Mises CPD to a table
%    T = convert_to_table(CPD, domain, evidence)
sz = CPD.sizes;
ns = zeros(1, max(domain));
ns(domain) = sz;

odom = domain(~isemptycell(evidence(domain)));
ps = domain(1:end-1);
dps = ps(CPD.dps);

[m, k, ~] = vonMises_CPD_params_given_dps(CPD, domain, evidence);

ns(odom) = 1;
dpsize = prod(ns(dps));
self = domain(end);
assert(myismember(self, odom));
self_val = evidence{self};
T = zeros(dpsize, 1);
for i=1:dpsize
    T(i) = vonMises_prob(self_val, m(:,i), k(:,:,i));
end

end

