function engine = cond_gauss_inf_engine(bnet)
% COND_GAUSS_INF_ENGINE Conditional Gaussian inference engine
% engine = cond_gauss_inf_engine(bnet)
%
% Enumerates all the discrete roots, and runs jtree on the remaining Gaussian nodes.

dnodes = mysetdiff(1:length(bnet.dag), bnet.cnodes);

%onodes = dnodes; % all the discrete ndoes will be observed
%engine.sub_engine = jtree_inf_engine(bnet, onodes);
bnet2 = bnet;
bnet2.observed = dnodes;
engine.sub_engine = jtree_inf_engine(bnet2);

% This is where we will store the results between enter_evidence and marginal_nodes
engine.T = [];
engine.mu = [];
engine.Sigma = [];
engine.joint_dmarginal = [];
engine.onodes = []; % needed for marginal_nodes
engine.evidence = []; % needed for marginal_nodes add_ev

engine = class(engine, 'cond_gauss_inf_engine', inf_engine(bnet));
