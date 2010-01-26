% Compare the speeds of various inference engines on a coupled HMM

N = 3;
Q = 2;
rand('state', 0);
randn('state', 0);
discrete = 0;
if discrete
  Y = 2; % size of output alphabet
else
  Y = 3; % size of observed vectors
end
coupled = 1;
bnet = mk_chmm(N, Q, Y, discrete, coupled); 
%bnet = mk_fhmm(N, Q, Y, discrete);  % factorial HMM
ss = length(bnet.node_sizes_slice);

T = 3;

USEC = exist('@jtree_C_inf_engine/collect_evidence','file');

engine = {};
engine{end+1} = jtree_dbn_inf_engine(bnet);
%engine{end+1} = jtree_ndx_dbn_inf_engine(bnet, 'ndx_type', 'SD');
%engine{end+1} = jtree_ndx_dbn_inf_engine(bnet, 'ndx_type', 'D');
%engine{end+1} = jtree_ndx_dbn_inf_engine(bnet, 'ndx_type', 'B');
if USEC, engine{end+1} = jtree_C_dbn_inf_engine(bnet); end
engine{end+1} = hmm_inf_engine(bnet);
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T); 

% times in matlab N=4 Q=4 T=5 (* = winner)
%     jtree    SD       B          hmm       dhmm      unrolled
%    0.6266    1.1563    8.3815    0.3069    0.1948*    0.8654  inf
%    0.9057*   2.1522   12.6314    2.6847    2.3107    3.1905  learn

%engine{end+1} = bk_inf_engine(bnet, 'ff', onodes);
%engine{end+1} = pearl_unrolled_dbn_inf_engine(bnet, T);

inf_time = cmp_inference_dbn(bnet, engine, T)
learning_time = cmp_learning_dbn(bnet, engine, T)

