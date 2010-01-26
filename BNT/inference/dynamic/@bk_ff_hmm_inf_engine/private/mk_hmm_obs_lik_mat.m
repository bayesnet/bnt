function obslik = mk_hmm_obs_lik_mat(bnet, onodes, evidence)
% MK_HMM_OBS_LIK_MAT Make the observation likelihood matrix for all slices
% obslik = mk_hmm_obs_lik_mat(bnet, onodes, evidence)
%
% obslik(i,t) = Pr(Y(t) | X(t)=i)

[ss T] = size(evidence);

hnodes = mysetdiff(1:ss, onodes);
ns = bnet.node_sizes_slice;
ns(onodes) = 1;
Q = prod(ns(hnodes));
obslik = zeros(Q,T);

dom = 1:ss;
for t=1:T
  bigpot = dpot(dom, ns(dom));
  for i=onodes(:)'
    if t==1
      e = bnet.equiv_class(i,1);
      fam = family(bnet.dag, i);
    else
      e = bnet.equiv_class(i,2);
      fam = family(bnet.dag, i, 2) + ss*(t-2);
    end
    pot = convert_to_pot(bnet.CPD{e}, 'd', fam(:), evidence);
    pot = set_domain_pot(pot, family(bnet.dag, i));
    bigpot = multiply_by_pot(bigpot, pot);
  end
  m = pot_to_marginal(bigpot);
  obslik(:,t) = m.T(:);
end


