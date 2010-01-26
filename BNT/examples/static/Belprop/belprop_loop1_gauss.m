% Compare different loopy belief propagation algorithms on a graph with a single loop.
% LBP should give exact results if it converges.

N = 4;
dag = zeros(N,N);
C = 1; S = 2; R = 3; W = 4;
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;
ns = 2*ones(1,N); 
bnet = mk_bnet(dag, ns, 'discrete', []);
for i=1:N
  bnet.CPD{i} = gaussian_CPD(bnet, i);
end

engines = {};
engines{end+1} = jtree_inf_engine(bnet);
engines{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel', 'max_iter', 20);
%engines{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel', 'max_iter', 20, 'filename', ...
%				  '/home/eecs/murphyk/matlab/gausspearl.txt', 'tol', 1e-5);

% pearl gaussian does not compute loglik
[time, engines] = cmp_inference_static(bnet, engines, 'maximize', 0, 'exact', 1, 'observed', [2], ...
				   'check_ll', 0, 'singletons_only', 0, 'check_converged', [2]);

