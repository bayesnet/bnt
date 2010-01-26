% Evaluate effectiveness of likelihood weighting on the lawn sprinkler example

N = 4; 
dag = zeros(N,N);
C = 1; R = 2; S = 3; W = 4;
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;

false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

bnet = mk_bnet(dag, ns);
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);


clear engine;
engine{1} = jtree_inf_engine(bnet);
engine{2} = likelihood_weighting_inf_engine(bnet);

nengines = length(engine);
m = cell(1, nengines);
ll = zeros(1, nengines);

evidence = cell(1,N);
%evidence{C} = true; % evidence at the top is the easiest
evidence{W} = true; % evidence at the bottom is the hardets

query = [R];

i=1;
engine{i}  = enter_evidence(engine{i}, evidence);
exact_m = marginal_nodes(engine{i}, query);

i=2;
samples = 100:100:500;
err = zeros(1, length(samples));
for j=1:length(samples)
  nsamples = samples(j);
  engine{i}  = enter_evidence(engine{i}, evidence, nsamples);
  approx_m = marginal_nodes(engine{i}, query);
  a1=approxeq(approx_m.T,exact_m.T,1e-1);
  a2=approxeq(approx_m.T,exact_m.T,1e-2);
  a3=approxeq(approx_m.T,exact_m.T,1e-3);
  e = sum(abs(approx_m.T(:) - exact_m.T(:)));
  fprintf('%d samples, 1dp %d, 2dp %d, 3dp %d,  err %f\n', nsamples, a1, a2, a3, e);
  err(j) = e;
end
