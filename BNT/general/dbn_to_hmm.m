function [startprob, transprob, obsprob] = dbn_to_hmm(bnet)
% DBN_TO_HMM % Convert DBN params to HMM params
% [startprob, transprob, obsprob] = dbn_to_hmm(bnet, onodes)
% startprob(i)
% transprob(i,j)
% obsprob{k}.big_CPT(i,o) if k'th observed node is discrete
% obsprob{k}.big_mu(:,i), .big_Sigma(:,:,i) if k'th observed node is Gaussian
% Big means the domain contains all the hidden discrete nodes, not just the parents.

% Called by constructor and by update_engine

ss = length(bnet.intra);
onodes = bnet.observed;
hnodes = mysetdiff(1:ss, onodes);
evidence = cell(ss, 2);
ns = bnet.node_sizes(:);
Qh = prod(ns(hnodes));
tmp = dpot_to_table(compute_joint_pot(bnet, hnodes, evidence));
startprob = reshape(tmp, Qh, 1);

tmp = dpot_to_table(compute_joint_pot(bnet, hnodes+ss, evidence, [hnodes hnodes+ss]));
transprob = mk_stochastic(reshape(tmp, Qh, Qh));

% P(o|ps) is used by mk_hmm_obs_lik_vec for a single time slice
% P(o|h) (the big version), where h = all hidden nodes, is used by enter_evidence

obsprob = cell(1, length(onodes));
for i=1:length(onodes)
  o = onodes(i);
  if bnet.auto_regressive(o)
    % We assume the parents of this node are all the hidden nodes in the slice,
    % so the params already are "big". Also, we assume we regress only on our old selves.
    % slice 1
    e = bnet.equiv_class(o);
    CPD = struct(bnet.CPD{e});
    O = ns(o);
    ps = bnet.parents{o};
    Qps = prod(ns(ps));
    obsprob{i}.big_mu0 = reshape(CPD.mean, [O Qps]);
    obsprob{i}.big_Sigma0 = reshape(CPD.cov, [O O Qps]);

    % slice t>1
    e = bnet.equiv_class(o+ss);
    CPD = struct(bnet.CPD{e});
    O = ns(o);
    dps = mysetdiff(bnet.parents{o+ss}, o);
    Qdps = prod(ns(dps));
    obsprob{i}.big_mu = reshape(CPD.mean, [O Qdps]);
    obsprob{i}.big_Sigma = reshape(CPD.cov, [O O Qdps]);
    obsprob{i}.big_W = reshape(CPD.weights, [O O Qdps]);
  else
    e = bnet.equiv_class(o+ss);
    CPD = struct(bnet.CPD{e});
    O = ns(o);
    ps = bnet.parents{o};
    Qps = prod(ns(ps));
    % We make a big potential, replicating the params if necessary
    % e.g., for a 2 chain coupled HMM, mu(:,Q1) becomes mu(:,Q1,Q2)
    bigpot = pot_to_marginal(compute_joint_pot(bnet, onodes(i), evidence, [hnodes onodes(i)]));

    if myismember(o, bnet.dnodes)
      obsprob{i}.CPT = reshape(CPD.CPT, [Qps O]);
      obsprob{i}.big_CPT = reshape(bigpot.T, Qh, O); 
    else
      obsprob{i}.big_mu = bigpot.mu;
      obsprob{i}.big_Sigma = bigpot.Sigma;
      
      if 1
      obsprob{i}.mu = reshape(CPD.mean, [O Qps]);
      C = reshape(CPD.cov, [O O Qps]);
      obsprob{i}.Sigma = C;
      d = size(obsprob{i}.mu, 1);
      for j=1:Qps
	obsprob{i}.inv_Sigma(:,:,j) = inv(C(:,:,j));
	obsprob{i}.denom(j) = (2*pi)^(d/2)*sqrt(abs(det(C(:,:,j))));
      end
      end
      
    end % if discrete
  end % if ar
end % for
