function obslik = mk_hmm_obs_lik_vec(bnet, evidence)
% MK_HMM_OBS_LIK_VEC Make the observation likelihood vector for one slice
% obslik = mk_obs_lik(bnet, evidence)
%
% obslik(i) = Pr(y(t) | X(t)=i)
% evidence{i,1} contains the evidence on node i in slice t-1
% evidence{i,2} contains the evidence on node i in slice t

ns = bnet.node_sizes;
ss = length(bnet.intra);
onodes = find(~isemptycell(evidence(:)));
hnodes = find(isemptycell(evidence(:)));
ens = ns;
ens(onodes) = 1;
Q = prod(ens(hnodes));
obslik = zeros(1,Q);
dom = (1:ss)+ss;
bigpot = dpot(dom, ens(dom));
onodes1 = find(~isemptycell(evidence(:,1)));
for i=onodes1(:)'
  e = bnet.equiv_class(i,2);
  fam = family(bnet.dag, i, 2);
  pot = convert_to_pot(bnet.CPD{e}, 'd', fam, evidence);
  bigpot = multiply_by_pot(bigpot, pot);
end
m = pot_to_marginal(bigpot);
obslik = m.T(:);
