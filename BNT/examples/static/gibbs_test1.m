function gibbs_test1()

disp('gibbs test 1')

rand('state', 0);
randn('state', 0);

%[bnet onodes hnodes qnodes] = gibbs_ex_1;
[bnet onodes hnodes qnodes] = gibbs_ex_2;

je = jtree_inf_engine(bnet);
ge = gibbs_sampling_inf_engine (bnet, 'T', 50, 'burnin', 0, ...
				'order', [2 2 1 2 1]);

ev = sample_bnet(bnet);

evidence = cell(length(bnet.dag), 1);
evidence(onodes) = ev(onodes);
[je lj] = enter_evidence(je, evidence);
[ge lg] = enter_evidence(ge, evidence);


mj = marginal_nodes(je, qnodes);

[mg ge] = marginal_nodes (ge, qnodes);
for t = 1:100
  [mg ge] = marginal_nodes (ge, qnodes, 'reset_counts', 0);
  diff = mj.T - mg.T;
  err(t) = norm (diff(:), 1);
end
clf
plot(err);
%title('error vs num. Gibbs samples')


%%%%%%%

function [bnet, onodes, hnodes, qnodes] = gibbs_ex_1
% bnet = gibbs_ex_1
% a simple network to test the gibbs sampling engine
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

onodes = 8:9;
hnodes = 1:7;
qnodes = [1 2 6];
ns = [2 3 4 3 5 2 4 3 2];

eclass = [1 2 3 2 4 5 6 7 8];

bnet = mk_bnet (dag, ns, 'equiv_class', eclass);

for i = 1:3
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

for i = 4:8
  bnet.CPD{i} = tabular_CPD(bnet, i+1);
end



%%%%%%%

function [bnet, onodes, hnodes, qnodes] = gibbs_ex_2
% bnet = gibbs_ex_2
% a very simple network
%
% 1   2
%  \ /
%   3

N = 3;
dag = zeros(N,N);
dag(1,3)=1; dag(2,3)=1;

onodes = 3;
hnodes = 1:2;
qnodes = 1:2;
ns = [2 4 3];

eclass = [1 2 3];

bnet = mk_bnet (dag, ns, 'equiv_class', eclass);

for i = 1:3
  bnet.CPD{i} = tabular_CPD(bnet, i);
end







