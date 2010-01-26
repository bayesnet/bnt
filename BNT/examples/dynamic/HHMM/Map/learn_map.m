seed = 1;
rand('state', seed);
randn('state', seed);

obs_model = 'unique';  % each cell has a unique label (essentially fully observable)
%obs_model = 'four'; % each cell generates 4 observations, NESW

% Generate the true network, and a randomization of it
realnet = mk_map_hhmm('p', 0.9, 'obs_model', obs_model);
rndnet = mk_rnd_map_hhmm('obs_model', obs_model);
eclass = realnet.equiv_class;
U = 1; A = 2; C = 3; F = 4; onodes = 5;

ss = realnet.nnodes_per_slice;
T = 100;
evidence = sample_dbn(realnet, 'length', T);
ev = cell(ss,T);
ev(onodes,:) = evidence(onodes,:);

infeng = jtree_dbn_inf_engine(rndnet);

if 0
% suppose we do not observe the final finish node, but only know 
% it is more likely to be on that off
ev2 = ev;
infeng = enter_evidence(infeng, ev2, 'soft_evidence_nodes', [F T], 'soft_evidence',  {[0.3 0.7]'});
end


learnednet = learn_params_dbn_em(infeng, {evidence}, 'max_iter', 5);

disp('real model')
disp_map_hhmm(realnet)

disp('learned model')
disp_map_hhmm(learnednet)

disp('rnd model')
disp_map_hhmm(rndnet)

