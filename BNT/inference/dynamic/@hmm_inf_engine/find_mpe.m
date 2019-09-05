function mpe = find_mpe(engine, evidence)
% FIND_MPE Find the most probable explanation (Viterbi)
% mpe = enter_evidence(engine, evidence, ...)
%
% evidence{i,t} = [] if if X(i,t) is hidden, and otherwise contains its observed value (scalar or column vector)
%

obslik = mk_hmm_obs_lik_matrix(engine, evidence);
path = viterbi_path(engine.startprob, engine.transprob, obslik);
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes_slice;
ns(bnet.observed) = 1;
ass = ind2subv(ns, path);
mpe = num2cell(ass');
mpe(bnet.observed,:) = evidence(bnet.observed,:);

