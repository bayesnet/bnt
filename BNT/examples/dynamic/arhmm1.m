% Make an HMM with autoregressive Gaussian observations (switching AR model)
%   X1 -> X2
%   |     | 
%   v     v
%   Y1 -> Y2 

seed = 0;
rand('state', seed);
randn('state', seed);

intra = zeros(2);
intra(1,2) = 1;
inter = zeros(2);
inter(1,1) = 1;
inter(2,2) = 1;
n = 2;

Q = 2; % num hidden states
O = 2; % size of observed vector

ns = [Q O];
dnodes = 1;
onodes = [2];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes);

bnet.CPD{1} = tabular_CPD(bnet, 1);
bnet.CPD{2} = gaussian_CPD(bnet, 2);
bnet.CPD{3} = tabular_CPD(bnet, 3);
bnet.CPD{4} = gaussian_CPD(bnet, 4);


T = 10; % fixed length sequences

engine = {};
%engine{end+1} = hmm_inf_engine(bnet);
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
%engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet));
%engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));

inf_time = cmp_inference_dbn(bnet, engine, T, 'check_ll',1);
learning_time = cmp_learning_dbn(bnet, engine, T, 'check_ll', 1);

