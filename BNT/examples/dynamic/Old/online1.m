% Check that online inference gives same results as filtering for various algorithms

N = 3;
Q = 2;
ss = N*2;

rand('state', 0);
randn('state', 0);


obs_size = 1;
discrete_obs = 0;
bnet = mk_chmm(N, Q, obs_size, discrete_obs);
ns = bnet.node_sizes_slice;

engine = {};
engine{end+1} = hmm_inf_engine(bnet);
E = length(engine);

onodes = (1:N)+N;

T = 4;
ev = cell(ss,T);
ev(onodes,:) = num2cell(randn(N, T));


filter = 1;
loglik2 = zeros(1,E);
for e=1:E
  [engine2{e}, loglik2(e)] = enter_evidence(engine{e}, ev, 'filter', filter);
end

loglik = zeros(1,E);
marg1 = cell(E,N,T);
for e=1:E
  ll = zeros(1,T);
  engine{e} = dbn_init_bel(engine{e});
  for t=1:T
    [engine{e}, ll(t)] = dbn_update_bel(engine{e}, ev(:,t), t);
    for i=1:N
      marg1{e,i,t} = dbn_marginal_from_bel(engine{e}, i);
    end
  end
  loglik1(e) = sum(ll);
end

assert(approxeq(loglik1, loglik2))

a = zeros(E,N,T);
for e=1:E
  for t=1:T
    for i=1:N
      marg2{e,i,t} = marginal_nodes(engine2{e}, i, t);
      a(e,i,t) = (approxeq(marg2{e,i,t}.T(:), marg1{e,i,t}.T(:)));
    end
  end
end

assert(all(a(:)==1))
