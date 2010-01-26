function engine = dhmm_inf_engine(bnet, onodes)
% DHMM_INF_ENGINE Inference engine for discrete DBNs which uses the forwards-backwards algorithm.
% engine = dhmm_inf_engine(bnet, onodes)
%
% 'onodes' specifies which nodes are observed; these must be leaves, and can be discrete or continuous.
% The remaining nodes are all hidden, and must be discrete.
% The DBN is converted to an HMM, with a single meganode, but which may have factored obs.

ss = length(bnet.intra);
hnodes = mysetdiff(1:ss, onodes);
evidence = cell(ss, 2);
ns = bnet.node_sizes;
Q = prod(ns(hnodes));
tmp = dpot_to_table(compute_joint_pot(bnet, hnodes, evidence));
engine.startprob = reshape(tmp, Q, 1);
tmp = dpot_to_table(compute_joint_pot(bnet, [hnodes hnodes+ss], evidence));
engine.transprob = mk_stochastic(reshape(tmp, Q, Q));
engine.obsprob = cell(1, length(onodes));
for i=1:length(onodes)
  tmp = dpot_to_table(compute_joint_pot(bnet, [hnodes onodes(i)], evidence));
  O = ns(onodes(i));
  engine.obsprob{i} = mk_stochastic(reshape(tmp, Q, O));
end

% This is where we will store the results between enter_evidence and marginal_nodes
engine.gamma = [];
engine.xi = [];

engine.onodes = onodes;
engine.hnodes = hnodes;
engine.maximize = [];

engine = class(engine, 'dhmm_inf_engine', inf_engine(bnet));

