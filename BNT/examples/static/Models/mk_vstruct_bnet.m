function oracle = mk_vstruct_bnet()
% MK_VSTRUCT_BNET Make a simple V-structured 3-node noisy-AND Bayes net
% oracle = mk_vstruct_bnet()

N = 3;
dag = zeros(N);
A = 1; B = 2; C = 3;
dag(A,C)=1;
dag(B,C)=1;
ns = 2*ones(1,N);

oracle = mk_bnet(dag, ns);
oracle.CPD{1} = tabular_CPD(oracle, 1, [0.5 0.5]);
oracle.CPD{2} = tabular_CPD(oracle, 2, [0.5 0.5]);
pnoise = 0.1; % degree of noise
oracle.CPD{3} = boolean_CPD(oracle, 3, 'named', 'all', pnoise);
