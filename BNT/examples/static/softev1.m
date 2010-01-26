% Check that adding soft evidence to a hidden node is equivalent to evaluating its leaf CPD.

% Make an HMM
T = 3; Q = 2; O = 2; cts_obs = 0; param_tying = 0;
bnet = mk_hmm_bnet(T, Q, O, cts_obs, param_tying);
N = 2*T;
onodes = bnet.observed;
hnodes = mysetdiff(1:N, onodes);
for i=1:N
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

ev = sample_bnet(bnet);
evidence = cell(1,N);
evidence(onodes) = ev(onodes);

engine = jtree_inf_engine(bnet);

[engine, ll] = enter_evidence(engine, evidence);
query = 1;
m = marginal_nodes(engine, query);


% Make a Markov chain with the same backbone
bnet2 = mk_markov_chain_bnet(T, Q);
for i=1:T
  S = struct(bnet.CPD{hnodes(i)}); % violate object privacy
  bnet2.CPD{i} = tabular_CPD(bnet2, i, S.CPT);
end

% Evaluate the observed leaves of the HMM
soft_ev = cell(1,T);
for i=1:T
  S = struct(bnet.CPD{onodes(i)}); % violate object privacy
  dist = S.CPT(:, evidence{onodes(i)});
  soft_ev{i} = dist;
end

% Use the leaf potentials as soft evidence
engine2 = jtree_inf_engine(bnet2);
[engine2, ll2] = enter_evidence(engine2, cell(1,T), 'soft', soft_ev);
m2 = marginal_nodes(engine2, query);

assert(approxeq(m2.T, m.T))
assert(approxeq(ll2, ll))



% marginal on node 1 without evidence
[engine2, ll2] = enter_evidence(engine2, cell(1,T));
m2 = marginal_nodes(engine2, 1);

% add soft evidence
soft_ev=cell(1,T);
soft_ev{1}=[0.7 0.3]; 
[engine2, ll2] = enter_evidence(engine2, cell(1,T), 'soft', soft_ev);
m3 = marginal_nodes(engine2, 1);

assert(approxeq(normalise(m2.T .* [0.7 0.3]'), m3.T))

