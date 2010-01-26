% Compare the speeds of various inference engines on the water DBN
seed = 0;
rand('state', seed);
randn('state', seed);

%bnet = mk_water_dbn;
bnet = mk_orig_water_dbn;

T = 3;
engine = {};
%engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
%engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet));
%engine{end+1} = jtree_dbn_inf_engine(bnet);
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
engine{end+1} = cbk_inf_engine(bnet, 'clusters', {[1],[2],[3],[4],[5],[6],[7],[8]}); %ff
engine{end+1} = cbk_inf_engine(bnet, 'clusters', {[1 2],[3 4 5 6],[7 8]}); %manually designed marginally independent by BK
engine{end+1} = cbk_inf_engine(bnet, 'clusters', {[1:5], [3:7], [7:8]}); %manually designed conditionally independent by BK
engine{end+1} = cbk_inf_engine(bnet, 'clusters', {[1 3], [2 3 7], [3 5], [3 4 7], [6 7 8]}); %automatically found using TJTs offline 
engine{end+1} = cbk_inf_engine(bnet, 'clusters', {[1 3 5], [2 3 5 7], [3 4 7], [4 6 7], [6 7 8]}); %automatically found using TJTs offline 
engine{end+1} = cbk_inf_engine(bnet, 'clusters', {[1 3 4 5], [2 3 4 7 8], [4 6 7 8]}); %automatically found using TJTs offline 

% bk_inf_engine yields exactly the same results for the marginally independent cases. 
%engine{end+1} = bk_inf_engine(bnet, 'clusters', 'ff');
%engine{end+1} = bk_inf_engine(bnet, 'clusters', { [1 2], [3 4 5 6], [7 8] });


inf_time = cmp_inference_dbn(bnet, engine, T, 'exact', 1)
learning_time = cmp_learning_dbn(bnet, engine, T, 'exact', 1)