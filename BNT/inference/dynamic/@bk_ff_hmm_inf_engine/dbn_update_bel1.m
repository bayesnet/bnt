function [engine, loglik] = dbn_update_bel1(engine, evidence)
% DBN_UPDATE_BEL Update the initial belief state (bk_ff_hmm)
% [engine, loglik] = dbn_update_bel(engine, evidence)
%
%  evidence{i} contains the evidence on node i in slice 1

oldbel = engine.bel;
bnet = bnet_from_engine(engine);
obslik = mk_hmm_obs_lik_vec1(bnet, evidence);
[newbel, lik] = normalise(oldbel .* obslik);
loglik = log(lik);

hnodes = engine.hnodes;
ns = bnet.node_sizes;
[marginals, marginalsT] = project_joint_onto_marginals(newbel, hnodes, ns);
newbel = combine_marginals_into_joint(marginalsT, hnodes, ns);          
engine.bel_marginals = marginalsT;
engine.bel = newbel;
