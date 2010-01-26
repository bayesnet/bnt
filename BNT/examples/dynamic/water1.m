% Compare the speeds of various inference engines on the water DBN
seed = 0;
rand('state', seed);
randn('state', seed);

bnet = mk_water_dbn;

T = 3;
engine = {};
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet));
engine{end+1} = jtree_dbn_inf_engine(bnet);
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T); 

%engine{end+1} = bk_inf_engine(bnet, 'ff', onodes);
%engine{end+1} = bk_inf_engine(bnet, { [1 2], [3 4 5 6], [7 8] }, onodes);

inf_time = cmp_inference_dbn(bnet, engine, T)
learning_time = cmp_learning_dbn(bnet, engine, T)

