% learn the structure of "discrete1" network.

rand('state', 0);
randn('state', 0);
N = 9;
dag = zeros(N,N);
dag(1,2)=1; dag(1,3)=1; dag(1,4)=1;
dag(2,5)=1; dag(3,6)=1; dag(4,7)=1;
dag(5,8)=1; dag(6,8)=1; dag(6,9)=1; dag(7,9) = 1;


dnodes = 1:N;
false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

bnet = mk_bnet(dag, ns);
% use random params
for i=1:N
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

nsamples = 500;
samplesM = cell(N, nsamples);
for i=1:nsamples
  samplesM(:,i) = sample_bnet(bnet);
end

hide = rand(N, nsamples) > 0.8;
[I,J]=find(hide);
for k=1:length(I)
  samplesM{I(k), J(k)} = [];
end

engine = jtree_inf_engine(bnet);
[bnet, LL, engine] = learn_params_em(engine, samplesM, 1);
LL

G0 = zeros(N,N);
for i=1: N-1
   G0(i, i+1) = 1;
end

figure;
draw_graph(G0);

B0 = mk_bnet(G0, ns);
% use random params
for i=1:N
  B0.CPD{i} = tabular_CPD(B0, i, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
end

max_loop = 10;
profile on -detail mmex
[B0, order, best_score] = learn_struct_EM(B0, samplesM, max_loop);
profile report
G1 = B0.dag;

G0 = G1(order, order);

