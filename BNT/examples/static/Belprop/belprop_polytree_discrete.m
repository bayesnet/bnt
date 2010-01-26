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

ns = 2*ones(1,N); % binary nodes

onodes = [1 5];

bnet = mk_bnet(dag, ns, 'observed', onodes);

if 0
seed = 0;
rand('state', seed);
randn('state', seed);
end

for i=1:N
  %bnet.CPD{i} = tabular_CPD(bnet, i);
  bnet.CPD{i} = noisyor_CPD(bnet, i);
end

engine = {};
engine{end+1} = jtree_inf_engine(bnet);
engine{end+1} = pearl_inf_engine(bnet, 'protocol', 'tree');
engine{end+1} = pearl_inf_engine(bnet, 'protocol', 'parallel');

[err, time] = cmp_inference_static(bnet, engine, 'maximize', 0, 'check_ll', 1);

