% Compare the speeds of various inference engines on the DBN in Kjaerulff
% "dHugin: A computational system for dynamic time-sliced {B}ayesian networks",
% Intl. J. Forecasting 11:89-111, 1995.
%
% The intra structure is (all arcs point downwards)
%
%  1 -> 2
%   \  /
%     3
%     |
%     4
%    / \
%   5   6
%   \  /
%     7
%     |
%     8
%
% The inter structure is 1->1, 4->4, 8->8

seed = 0;
rand('state', seed);
randn('state', seed);

ss = 8;
intra = zeros(ss);
intra(1,[2 3])=1;
intra(2,3)=1;
intra(3,4)=1;
intra(4,[5 6])=1;
intra([5 6], 7)=1;
intra(7,8)=1;

inter = zeros(ss);
inter(1,1)=1;
inter(4,4)=1;
inter(8,8)=1;

ns = 2*ones(1,ss);
onodes = 2;
bnet = mk_dbn(intra, inter, ns, 'observed', onodes, 'eclass2', (1:ss)+ss);
for i=1:2*ss
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

T = 4;

engine = {};
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
engine{end+1} = jtree_dbn_inf_engine(bnet);
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
%engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet)); % observed nodes have children

inf_time = cmp_inference_dbn(bnet, engine, T)
learning_time = cmp_learning_dbn(bnet, engine, T)
