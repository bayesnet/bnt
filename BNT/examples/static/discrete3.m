% Compare various inference engines on the following network (from Jensen (1996) p84 fig 4.17)
%    1
%  / | \
% 2  3  4
% |  |  |
% 5  6  7
%  \/ \/
%  8   9
% where all arcs point downwards

N = 9;
dag = zeros(N,N);
dag(1,2)=1; dag(1,3)=1; dag(1,4)=1;
dag(2,5)=1; dag(3,6)=1; dag(4,7)=1;
dag(5,8)=1; dag(6,8)=1; dag(6,9)=1; dag(7,9) = 1;

dnodes = 1:N;
false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

onodes = [1];
evidence = cell(1,N);
evidence(onodes) = num2cell(1);
bnet = mk_bnet(dag, ns, 'observed', onodes);
% use random params
%for i=1:N
%  bnet.CPD{i} = tabular_CPD(bnet, i);
%end
bnet.CPD{1} = tabular_CPD(bnet, 1, 'sparse', 1, 'CPT', [0.8, 0.2]);
bnet.CPD{2} = tabular_CPD(bnet, 2, 'sparse', 1, 'CPT', [1 0 0 1]);
bnet.CPD{3} = tabular_CPD(bnet, 3, 'sparse', 1, 'CPT', [0 1 1 0]);
bnet.CPD{4} = tabular_CPD(bnet, 4, 'sparse', 1, 'CPT', [1 1 0 0]);
bnet.CPD{5} = tabular_CPD(bnet, 5, 'sparse', 1, 'CPT', [0 0 1 1]);
bnet.CPD{6} = tabular_CPD(bnet, 6, 'sparse', 1, 'CPT', [1 0 0 1]);
bnet.CPD{7} = tabular_CPD(bnet, 7, 'sparse', 1, 'CPT', [0 1 1 0]);
bnet.CPD{8} = tabular_CPD(bnet, 8, 'sparse', 1, 'CPT', [1 1 0 0 0 0 1 1]);
bnet.CPD{9} = tabular_CPD(bnet, 9, 'sparse', 1, 'CPT', [0 1 0 1 1 0 1 0]);

engine = jtree_sparse_inf_engine(bnet);
tic
[engine, ll] = enter_evidence(engine, evidence);
toc

