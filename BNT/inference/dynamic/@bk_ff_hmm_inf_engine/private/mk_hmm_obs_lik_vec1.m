function obslik = mk_hmm_obs_lik_vec1(bnet, evidence)
% MK_HMM_OBS_LIK_VEC1 Make the observation likelihood vector for the first slice
% obslik = mk_hmm_obs_lik_vec1(engine, evidence)
%
% obslik(i) = Pr(y(1) | X(1)=i)
% evidence{i} contains the evidence on node i in slice 1

ns = bnet.node_sizes;
ss = length(ns);
onodes = find(~isemptycell(evidence(:)));
hnodes = find(isemptycell(evidence(:)));
ens = ns;
ens(onodes) = 1;
Q = prod(ens(hnodes));
obslik = zeros(1,Q);
dom = (1:ss);
bigpot = dpot(dom, ens(dom));
for i=onodes(:)'
  e = bnet.equiv_class(i,1);
  fam = family(bnet.dag, i);
  pot = convert_to_pot(bnet.CPD{e}, 'd', fam(:), evidence);
  bigpot = multiply_by_pot(bigpot, pot);
end
m = pot_to_marginal(bigpot);
obslik = m.T(:);
