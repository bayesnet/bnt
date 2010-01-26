bnet = mk_square_hhmm(1, 1);

engine = {};
engine{end+1} = hmm_inf_engine(bnet);
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));

exact = 1:length(engine);
filter = 0;
single = 0;
maximize = 0;
T = 4;

[err, inf_time, engine] = cmp_inference(bnet, engine, exact, T, filter, single, maximize);
