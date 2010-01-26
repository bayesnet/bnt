% make a factor graph corresponding to an  HMM with Gaussian outputs, where we absorb the
% evidence up front 

seed = 1;
rand('state', seed);
randn('state', seed);

T = 3;
Q = 3;
O = 2;
cts_obs = 1;
param_tying = 1;
bnet = mk_hmm_bnet(T, Q, O, cts_obs, param_tying);
N = 2*T;
onodes = bnet.observed;
hnodes = mysetdiff(1:N, onodes);

data = sample_bnet(bnet);

init_factor = bnet.CPD{1};
obs_factor = bnet.CPD{3};
edge_factor = bnet.CPD{2}; % trans matrix

nfactors = T;
nvars = T; % hidden only
G = zeros(nvars, nfactors);
G(1,1) = 1;
for t=1:T-1
  G(t:t+1, t+1)=1;
end

node_sizes = Q*ones(1,T);

% We tie params as follows:
% the first hidden node use init_factor (number 1)
% all hidden nodes on the backbone use edge_factor (number 2)
% all observed nodes use the same factor, namely obs_factor

small_fg = mk_fgraph_given_ev(G, node_sizes, {init_factor, edge_factor}, {obs_factor}, data(onodes), ...
			 'equiv_class', [1 2*ones(1,T-1)], 'ev_equiv_class', ones(1,T));

small_bnet = fgraph_to_bnet(small_fg);

% don't pre-process evidence
% big_fg = bnet_to_fgraph(bnet); % can't handle Gaussian node


engine = {};
engine{1} = jtree_inf_engine(bnet);
engine{2} = belprop_fg_inf_engine(small_fg, 'max_iter', 2*T);
engine{3} = jtree_inf_engine(small_bnet);
nengines = length(engine);


% on BN, use the original evidence
evidence = cell(1, 2*T);
evidence(onodes) = data(onodes);
tic; [engine{1}, ll(1)] = enter_evidence(engine{1}, evidence); toc


% on small_fg, we have already included the evidence
evidence = cell(1,T);
tic; [engine{2}, ll(2)] = enter_evidence(engine{2}, evidence); toc


% on small_bnet, we must add evidence to the dummy nodes 
V = small_fg.nvars;
dummy = V+1:V+small_fg.nfactors;
N = max(dummy);
evidence = cell(1, N);
evidence(dummy) = {1};
tic; [engine{3}, ll(3)] = enter_evidence(engine{3}, evidence); toc



marg = zeros(T, nengines, Q); % marg(t,e,:)
for t=1:T
  for e=1:nengines
    m = marginal_nodes(engine{e}, t);
    marg(t,e,:) = m.T;
  end
end
marg(:,:,1)
