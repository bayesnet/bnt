clear all
close all

bnet = mk_asia_bnet();

N = length(bnet.dag);

base_proba = 0.3;
bnet_miss = gener_MCAR_net(bnet, base_proba);

m=500;
[data, comp_data, bnet_miss, taux, bnet_orig, notok] = gener_data_from_bnet_miss(bnet_miss, m, base_proba);

dag0 = zeros(N);
bnet0 = mk_bnet(dag0, bnet.node_sizes);
for node = 1:N
  bnet0.CPD{node} = tabular_CPD(bnet0, node, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
end

max_loop = 6;
[bnet0, cpdag, BIC_score, nloop] = learn_struct_ges_EM(bnet0, data, max_loop);

[XX, YY] = draw_graph(bnet.dag);
figure;
draw_graph(bnet0.dag, {'1','2','3','4','5','6','7','8'}, zeros(1,N), XX, YY);
