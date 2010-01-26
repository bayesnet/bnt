function engine = bk_ff_hmm_inf_engine(bnet)
% BK_FF_HMM_INF_ENGINE Naive (HMM-based) implementation of fully factored form of Boyen-Koller 
% engine = bk_ff_hmm_inf_engine(bnet)
%
% This is implemented on top of the forwards-backwards algo for HMMs,
% so it is *less* efficient than exact inference! However, it is good for educational purposes,
% because it illustrates the BK algorithm very clearly.

[persistent_nodes, transient_nodes] = partition_dbn_nodes(bnet.intra, bnet.inter);
assert(isequal(sort(bnet.observed), transient_nodes));
[engine.prior, engine.transmat] = dbn_to_hmm(bnet);

ss = length(bnet.intra);

engine.bel = [];
engine.bel_marginals = [];
engine.marginals = [];


engine = class(engine, 'bk_ff_hmm_inf_engine', inf_engine(bnet));

