% Computing most probable explanation.

% If you don't break ties consistently, loopy can give wrong mpe
% even though the graph has no cycles, and even though the max-marginals are the same.
% This example was contributed by Wentau Yih <wtyih@yahoo.com> 29 Jan 02.

% define loop-free graph structure (all edges point down)
%
% Xe1   Xe2
%  |    |
%  E1   E2
%    \ /
%     R1
%     |
%    Xr1

N = 6;
dag = zeros(N,N);
Xe1 = 1; Xe2 = 2; E1 = 3; E2 = 4; R1 = 5; Xr1 = 6;
dag(Xe1, E1) = 1;
dag(Xe2, E2) = 1;
dag([E1 E2], R1) = 1;
dag(R1, Xr1) = 1;

node_sizes = [ 1 1 2 2 2 1 ];

% create BN

bnet = mk_bnet(dag, node_sizes, 'observed', [Xe1 Xe2 Xr1]);

% fill in CPT

bnet.CPD{Xe1} = tabular_CPD(bnet, Xe1, [1]);
bnet.CPD{Xe2} = tabular_CPD(bnet, Xe2, [1]);
bnet.CPD{E1} = tabular_CPD(bnet, E1, [0.2 0.8]);
bnet.CPD{E2} = tabular_CPD(bnet, E2, [0.3 0.7]);
bnet.CPD{R1} = tabular_CPD(bnet, R1, [1 1 1 0.8 0 0 0 0.2]);
bnet.CPD{Xr1} = tabular_CPD(bnet, Xr1, [0.15 0.85]);

clear engine;
engine{1} = belprop_inf_engine(bnet);
engine{2} = jtree_inf_engine(bnet);
engine{3} = global_joint_inf_engine(bnet);
engine{4} = var_elim_inf_engine(bnet);

evidence = cell(1,N);
evidence{Xe1} = 1;  evidence{Xe2} = 1;  evidence{Xr1} = 1;

mpe = find_mpe(engine{1}, evidence, 'break_ties', 0) % gives wrong results
mpe = find_mpe(engine{1}, evidence)
for i=2:4
  mpe = find_mpe(engine{i}, evidence)
end
