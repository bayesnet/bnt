N = 1; % single chain = HMM - should give exact answers
Q = 2;
rand('state', 0);
randn('state', 0);
discrete = 1;
if discrete
  Y = 2; % size of output alphabet
else
  Y = 1;
end
coupled = 1;
bnet  = mk_chmm(N, Q, Y, discrete, coupled);
ss = N*2;

T = 3;

engine = {};
engine{end+1} = jtree_dbn_inf_engine(bnet); 
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T); 
engine{end+1} = pearl_unrolled_dbn_inf_engine(bnet, 'protocol', 'tree');  

inf_time = cmp_inference_dbn(bnet, engine, T)
learning_time = cmp_learning_dbn(bnet, engine, T)

