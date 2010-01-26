BNT_HOME = '/home/ai2/murphyk/matlab/FullBNT'; % edit this

d = fullfile(BNT_HOME, 'BNT');
%PC = (strncmp(computer,'PC',2));

cd(sprintf('%s/potentials/Tables', d))
mex marg_table.c % used by @dpot/marginalize_pot.m
mex marg_sparse_table.c %used by sparse jtree
mex mult_by_table.c
mex mult_by_sparse_table.c
mex divide_by_table.c
mex divide_by_sparse_table.c
mex rep_mult.c

% Written by Wei Hu
cd(sprintf('%s/CPDs/@discrete_CPD', d))
mex convert_to_sparse_table.c

% Written by Wei Hu
cd(sprintf('%s/inference/static/@jtree_sparse_inf_engine', d))
mex init_pot.c
mex collect_evidence.c
mex distribute_evidence.c

% written by Bhaskara Marthi 
cd(sprintf('%s/inference/static/@gibbs_sampling_inf_engine/private', d))
mex compute_posterior.c
mex get_slice_dbn.c
mex sample_single_discrete.c



