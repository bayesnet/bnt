% Inference on a conditional Gaussian model

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

dnodes = 1;
%onodes = [1 5];
bnet = mk_bnet(dag, ns, 'discrete', dnodes, 'observed', dnodes);

bnet.CPD{1} = tabular_CPD(bnet, 1);
for i=2:N
  bnet.CPD{i} = gaussian_CPD(bnet, i);
end

engine = {};
engine{end+1} = jtree_inf_engine(bnet);
engine{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel');

[time, engine] = cmp_inference_static(bnet, engine, 'maximize', 0, 'check_ll', 0, ...
				      'singletons_only', 0, 'observed', [1 3]);
