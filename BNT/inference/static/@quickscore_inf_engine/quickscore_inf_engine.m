function engine = quickscore_inf_engine(inhibit, leak, prior)
% QUICKSCORE_INF_ENGINE Exact inference for the QMR network
% engine = quickscore_inf_engine(inhibit, leak, prior)
%
% We create an inference engine for QMR-like networks.
% QMR is a bipartite graph, where the top layer contains hidden disease nodes,
% and the bottom later contains observed finding nodes.
% The diseases have Bernoulli CPDs, the findings noisy-or CPDs.
% The original QMR (Quick Medical Reference) network has specific parameter values which we are not
% allowed to release, for commercial reasons.
%
% inhibit(f,d) = inhibition probability on f->d arc for disease d, finding f
% If inhibit(f,d) = 1, there is effectively no arc from d->f
% leak(j) = inhibition prob. on leak node -> finding j arc
% prior(i) = prob. disease i is on
%
% We use exact inference, which takes O(2^P) time, where P is the number of positive findings.
% For details, see
% - Heckerman, "A tractable inference algorithm for diagnosing multiple diseases", UAI 89.
% - Rish and Dechter, "On the impact of causal independence", UCI tech report, 1998.
% Note that this algorithm is numerically unstable, since it adds a large number of positive and
% negative terms and hopes that some of them exactly cancel.
%
% For an interesting variational approximation, see
% - Jaakkola and Jordan, "Variational probabilistic inference and the QMR-DT network", JAIR 10, 1999.
%
% See also 
% - "Loopy belief propagation for approximate inference: an empirical study",
%      K. Murphy, Y. Weiss and M. Jordan, UAI 99.

engine.inhibit = inhibit;
engine.leak = leak;
engine.prior = prior;

% store results here between enter_evidence and marginal_nodes
engine.post = [];

engine = class(engine, 'quickscore_inf_engine'); % not a child of the inf_engine class!
