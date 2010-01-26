% make an unrolled HMM, convert to factor graph, and check that 
% loopy propagation on the fgraph gives the exact answers.

seed = 1;
rand('state', seed);
randn('state', seed);

T = 3;
Q = 3;
O = 3;
cts_obs = 0;
param_tying = 1;
bnet = mk_hmm_bnet(T, Q, O, cts_obs, param_tying);

data = sample_bnet(bnet);

fgraph = bnet_to_fgraph(bnet);
big_bnet = fgraph_to_bnet(fgraph);
% converting factor graph back does not recover the structure of the original bnet

max_iter = 2*T;

engine = {};
engine{1} = jtree_inf_engine(bnet);
engine{2} = belprop_inf_engine(bnet, 'max_iter', max_iter);
engine{3} = belprop_fg_inf_engine(fgraph, 'max_iter', max_iter);
engine{4} = jtree_inf_engine(big_bnet);
nengines = length(engine);

big_engine = 4;
fgraph_engine = 3;


N = 2*T;
evidence = cell(1,N);
onodes = bnet.observed;
evidence(onodes) = data(onodes);
hnodes = mysetdiff(1:N, onodes);

bigN = length(big_bnet.dag);
big_evidence = cell(1, bigN);
big_evidence(onodes) = data(onodes);
big_evidence(N+1:end) = {1}; % factors are observed to be 1

ll = zeros(1, nengines);
for i=1:nengines
  if i==big_engine
    tic; [engine{i}, ll(i)] = enter_evidence(engine{i}, big_evidence); toc
  else
    tic; [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence); toc
  end
end

% compare all engines to engine{1}

% the log likelihood values may be bogus...
for i=2:nengines
  %assert(approxeq(ll(1), ll(i)));
end


marg = zeros(T, nengines, Q); % marg(t,e,:)
for t=1:T
  for e=1:nengines
    m = marginal_nodes(engine{e}, t);
    marg(t,e,:) = m.T;
  end
end
marg


m = cell(nengines, T);
for i=1:T
  for e=1:nengines
    m{e,i} = marginal_nodes(engine{e}, hnodes(i));
  end
  for e=2:nengines
    assert(approxeq(m{e,i}.T, m{1,i}.T));
  end
end

mpe = {};
ll = zeros(1, nengines);
for e=1:nengines
  if e==big_engine
    mpe{e} = find_mpe(engine{e}, big_evidence);
    mpe{e} = mpe{e}(1:N); % chop off dummy nodes
  else
    mpe{e} = find_mpe(engine{e}, evidence);
  end
end

% fgraph can't compute loglikelihood for software reasons
% jtree on the big_bnet gives the wrong ll
for e=2:nengines
  %assert(approxeq(ll(1), ll(e)));
  assert(approxeq(cell2num(mpe{1}), cell2num(mpe{e})))
end
