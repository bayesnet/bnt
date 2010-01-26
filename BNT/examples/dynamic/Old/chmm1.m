% Compare the speeds of various inference engines on a coupled HMM

N = 2;
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
[bnet, onodes] = mk_chmm(N, Q, Y, discrete, coupled);
ss = N*2;

T = 3;


engine = {};
tic; engine{end+1} = jtree_dbn_inf_engine(bnet, 'observed', onodes);  toc
%tic; engine{end+1} = jtree_ndxSD_dbn_inf_engine(bnet, onodes);  toc
%tic; engine{end+1} = jtree_ndxB_dbn_inf_engine(bnet, onodes);  toc
engine{end+1} = hmm_inf_engine(bnet, onodes);
%engine{end+1} = dhmm_inf_engine(bnet, onodes);
tic; engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T, onodes);  toc

%engine{end+1} = bk_inf_engine(bnet, 'ff', onodes);
%engine{end+1} = loopy_dbn_inf_engine(bnet, onodes);

exact = [1 2 3];

filter = 0;
single = 0;
maximize = 0;

[err, time, engine] = cmp_inference(bnet, onodes, engine, exact, T, filter, single, maximize);
%err = cmp_learning(bnet, onodes, engine, exact, T);


