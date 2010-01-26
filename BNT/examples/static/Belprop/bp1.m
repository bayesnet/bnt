% Compare different loopy belief propagation algorithms on a graph with a single loop.
% LBP should give exact results if it converges.

seed = 0;
rand('state', seed);
randn('state', seed);

N = 2;
dag = zeros(N,N);
dag(1,2)=1;
ns = ones(1,N); 
bnet = mk_bnet(dag, ns, 'discrete', []);
for i=1:N
  %bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', 0);
  bnet.CPD{i} = gaussian_CPD(bnet, i);
end

engines = {};
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = pearl_inf_engine(bnet, 'protocol', 'tree');

[time, engines] = cmp_inference_static(bnet, engines, 'maximize', 0, 'exact', 1:2, 'observed', [2], ...
				   'check_ll', 0, 'singletons_only', 1, 'check_converged', []);

