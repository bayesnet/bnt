% Compare the speeds of various inference engines on the BAT DBN
[bnet, names] = mk_bat_dbn;

T = 3; % fixed length sequence - we make it short just for speed

USEC = exist('@jtree_C_inf_engine/collect_evidence','file');

disp('constructing engines for BAT');
engine = {}; % time in seconds for inference
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T, 'useC', USEC);  % 0.39
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet)); % 4.89
engine{end+1} = jtree_dbn_inf_engine(bnet); % 4.45
if 0
engine{end+1} = jtree_ndx_dbn_inf_engine(bnet, 'ndx_type', 'SD'); % 2.98
engine{end+1} = jtree_ndx_dbn_inf_engine(bnet, 'ndx_type', 'D'); % 3.52
engine{end+1} = jtree_ndx_dbn_inf_engine(bnet, 'ndx_type', 'B'); % 2.40
if USEC, engine{end+1} = jtree_C_dbn_inf_engine(bnet); end % 3.54
%engine{end+1} = hmm_inf_engine(bnet, onodes); % too big
end

%tic; engine{end+1} = frontier_inf_engine(bnet); toc % very slow
% The frontier engine thrashes badly on the BAT network
%tic; engine{end+1} = bk_inf_engine(bnet, 'exact', onodes); toc  % SLOW!

%tic; engine{end+1} = bk_inf_engine(bnet, 'ff', onodes); toc

%clusters{1} = [stringmatch({'LeftClr', 'RightClr', 'LatAct', 'Xdot', 'InLane'}, names)];
%clusters{2} = [stringmatch({'FwdAct', 'Ydot', 'Stopped', 'EngStatus', 'FBStatus'}, names)];      

%tic; engine{end+1} = bk_inf_engine(bnet, clusters, onodes); toc

disp('inference')
time = cmp_inference_dbn(bnet, engine, T)

disp('learning')
time = cmp_learning_dbn(bnet, engine, T)








