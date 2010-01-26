N = 1; % regular HMM
Q = 2;
ss = 2;
hnodes = 1;
onodes = 2;

rand('state', 0);
randn('state', 0);
O = 2;
discrete_obs = 1;
bnet = mk_chmm(N, Q, O, discrete_obs);
ns = bnet.node_sizes_slice;

engine = hmm_inf_engine(bnet, onodes);

T = 4;
ev = cell(ss,T);
ev(onodes,:) = num2cell(sample_discrete([0.5 0.5], N, T));


engine = dbn_init_bel(engine);
for t=1:T
  if t==1
    [engine, ll(t)] = dbn_update_bel1(engine, ev(:,t));
  else
    [engine, ll(t)] = dbn_update_bel(engine, ev(:,t-1:t));
  end
  % one-step ahead prediction
  lag = 1;
  engine2 = dbn_predict_bel(engine, lag);  
  marg = dbn_marginal_from_bel(engine2, 1)
  marg = dbn_marginal_from_bel(engine2, 2)
end
