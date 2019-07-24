% Lawn sprinker example from Russell and Norvig p454
% See www.cs.berkeley.edu/~murphyk/Bayes/usage.html for details.

rand('state', 0);
randn('state', 0);

N = 4;
dag = zeros(N,N);
C = 1; S = 2; R = 3; W = 4;
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

nsamples = 500;
samplesM = cell(N, nsamples);
for i=1:nsamples
  samplesM(:,i) = sample_bnet(bnet);
end

hide = rand(N, nsamples) > 0.9;
[I,J]=find(hide);
for k=1:length(I)
  samplesM{I(k), J(k)} = [];
end

% Make a initial chain like dag
G0 = zeros(N,N);
for i=1:N-1
   G0(i, i+1) = 1;
end

figure;
draw_graph(G0);

B0 = mk_bnet(G0, ns);
% use random params
for i=1:N
  B0.CPD{i} = tabular_CPD(B0, i, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
end

max_loop = 30;
%profile on -detail mmex
[B0, order, best_score] = learn_struct_EM(B0, samplesM, max_loop);
%profile report

dag1 = B0.dag;
dag1 = dag1(order,order);
