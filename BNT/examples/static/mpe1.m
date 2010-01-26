seed = 1;
rand('state', seed);
randn('state', seed);

N = 4; 
dag = zeros(N,N);
C = 1; S = 2; R = 3; W = 4;
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;

false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

bnet = mk_bnet(dag, ns);
if 0
  bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
  bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
  bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
  bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);
else
  for i=1:N, bnet.CPD{i} = tabular_CPD(bnet, i); end
end



evidence = cell(1,N);
onodes = [1 3];
data = sample_bnet(bnet);
evidence(onodes) = data(onodes);

clear engine;
engine{1} = belprop_inf_engine(bnet);
engine{2} = jtree_inf_engine(bnet);
engine{3} = global_joint_inf_engine(bnet);
engine{4} = var_elim_inf_engine(bnet);
E = length(engine);

clear mpe;
for e=1:E
  mpe{e} = find_mpe(engine{e}, evidence);
end
for e=2:E
  assert(isequal(mpe{1}, mpe{e}))
end
