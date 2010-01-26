function engine = kalman_inf_engine(bnet)
% KALMAN_INF_ENGINE Inference engine for Linear-Gaussian state-space models.
% engine = kalman_inf_engine(bnet)
%
% 'onodes' specifies which nodes are observed; these must be leaves.
% The remaining nodes are all hidden. All nodes must have linear-Gaussian CPDs.
% The hidden nodes must be persistent, i.e., they must have children in
% the next time slice. In addition, they may not have any children within the current slice,
% except to the observed leaves. In other words, the topology must be isomorphic to a standard LDS.
%
% There are many derivations of the filtering and smoothing equations for Linear Dynamical
% Systems in the literature. I particularly like the following
% - "From HMMs to LDSs", T. Minka, MIT Tech Report, (no date), available from
%    ftp://vismod.www.media.mit.edu/pub/tpminka/papers/minka-lds-tut.ps.gz

[engine.trans_mat, engine.trans_cov, engine.obs_mat, engine.obs_cov, engine.init_state, engine.init_cov] = ...
    dbn_to_lds(bnet);

% This is where we will store the results between enter_evidence and marginal_nodes
engine.one_slice_marginal = [];
engine.two_slice_marginal = [];

engine = class(engine, 'kalman_inf_engine', inf_engine(bnet));
