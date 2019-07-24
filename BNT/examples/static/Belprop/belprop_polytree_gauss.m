% Do the example from Satnam Alag's PhD thesis, UCB ME dept 1996 p46

% Make the following polytree, where all arcs point down

% 1   2
%  \ /
%   3
%  / \
% 4   5

N = 5;
dag = zeros(N,N);
dag(1,3) = 1;
dag(2,3) = 1;
dag(3, [4 5]) = 1;

ns = [2 1 2 1 2];

bnet = mk_bnet(dag, ns, 'discrete', []);

bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', [1 0]', 'cov', [4 1; 1 4]);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', 1, 'cov', 1);
B1 = [1 2; 1 0]; B2 = [2 1]';
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', [0 0]', 'cov', [2 1; 1 1], ...
			   'weights', [B1 B2]);
H1 = [1 1];
bnet.CPD{4} = gaussian_CPD(bnet, 4, 'mean', 0, 'cov', 1, 'weights', H1);
H2 = [1 0; 1 1];
bnet.CPD{5} = gaussian_CPD(bnet, 5, 'mean', [0 0]', 'cov', eye(2), 'weights', H2);

engine = {};
engine{end+1} = jtree_inf_engine(bnet);
engine{end+1} = pearl_inf_engine(bnet, 'protocol', 'tree');
engine{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel');
E = length(engine);

if 1
% no evidence
evidence = cell(1,N);
ll = zeros(1,E);
for e=1:E
  [engine{e}, ll(e)] = enter_evidence(engine{e}, evidence);
  add_ev = 1;
  m = marginal_nodes(engine{e}, 3, add_ev);
  assert(approxeq(m.mu, [3 2]'))
  assert(approxeq(m.Sigma, [30 9; 9 6]))

  m = marginal_nodes(engine{e}, 4, add_ev);
  assert(approxeq(m.mu, 5))
  assert(approxeq(m.Sigma, 55))

  m = marginal_nodes(engine{e}, 5, add_ev);
  assert(approxeq(m.mu, [3 5]'))
  assert(approxeq(m.Sigma, [31 39; 39 55]))
end
end

if 1
% evidence on leaf 5
evidence = cell(1,N);
evidence{5} = [5 5]';
for e=1:E
  [engine{e}, ll(e)] = enter_evidence(engine{e}, evidence);
  add_ev = 1;
  m = marginal_nodes(engine{e}, 3, add_ev);
  assert(approxeq(m.mu, [4.4022 1.0217]'))
  assert(approxeq(m.Sigma, [0.7011 -0.4891; -0.4891 1.1087]))

  m = marginal_nodes(engine{e}, 4, add_ev);
  assert(approxeq(m.mu, 5.4239))
  assert(approxeq(m.Sigma, 1.8315))

  m = marginal_nodes(engine{e}, 1, add_ev);
  assert(approxeq(m.mu, [0.3478 1.1413]'))
  assert(approxeq(m.Sigma, [1.8261 -0.1957; -0.1957 1.0924]))

  m = marginal_nodes(engine{e}, 2, add_ev);
  assert(approxeq(m.mu, 0.9239))
  assert(approxeq(m.Sigma, 0.8315))

  m = marginal_nodes(engine{e}, 5, add_ev);
  assert(approxeq(m.mu, evidence{5}))
  assert(approxeq(m.Sigma, zeros(2)))
end
end

if 1
% evidence on leaf 4 (non-info-state version is uninvertible)
evidence = cell(1,N);
evidence{4} = 10;
for e=1:E
  [engine{e}, ll(e)] = enter_evidence(engine{e}, evidence);
  add_ev = 1;
  m = marginal_nodes(engine{e}, 3, add_ev);
  assert(approxeq(m.mu, [6.5455 3.3636]'))
  assert(approxeq(m.Sigma, [2.3455 -1.6364; -1.6364 1.9091]))

  m = marginal_nodes(engine{e}, 5, add_ev);
  assert(approxeq(m.mu, [6.5455 9.9091]'))
  assert(approxeq(m.Sigma, [3.3455 0.7091; 0.7091 1.9818]))

  m = marginal_nodes(engine{e}, 1, add_ev);
  assert(approxeq(m.mu, [1.9091 0.9091]'))
  assert(approxeq(m.Sigma, [2.1818 -0.8182; -0.8182 2.1818]))

  m = marginal_nodes(engine{e}, 2, add_ev);
  assert(approxeq(m.mu, 1.2727))
  assert(approxeq(m.Sigma, 0.8364))
end
end


if 1
% evidence on leaves 4,5 and root 2
evidence = cell(1,N);
evidence{2} = 0;
evidence{4} = 10;
evidence{5} = [5 5]';
for e=1:E
  [engine{e}, ll(e)] = enter_evidence(engine{e}, evidence);
  add_ev = 1;
  m = marginal_nodes(engine{e}, 3, add_ev);
  assert(approxeq(m.mu, [4.9964 2.4444]'));
  assert(approxeq(m.Sigma, [0.6738 -0.5556; -0.5556 0.8889]));

  m = marginal_nodes(engine{e}, 1, add_ev);
  assert(approxeq(m.mu, [2.2043 1.2151]'));
  assert(approxeq(m.Sigma, [1.2903 -0.4839; -0.4839 0.8065]));
end
end

if 1
  [time, engine] = cmp_inference_static(bnet, engine, 'maximize', 0, 'check_ll', 0, ...
				     'singletons_only', 0, 'observed', [1 3 5]);
end
