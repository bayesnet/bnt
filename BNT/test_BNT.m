% If all is well, all of these scripts should run without errors.


% bnets
cg1
cg2
discrete1
fa1
gaussian1
gaussian2
if exist('@gibbs_sampling_inf_engine/private/compute_posterior','file')
  % only exists if installC has been run
  gibbs_test1
end
learn1
lw1
mfa1
mixexp1
mixexp2
mixexp3
mog1
mpe1
mpe2
qmr1
qmr2
sample1
softev1
softmax1
sprinkler1


% belief propagation
belprop_polytree_discrete
belprop_polytree_gauss % alag
belprop_polytree_cg
belprop_loop1_discrete
belprop_loop1_gauss
belprop_loopy_discrete
belprop_loopy_gauss
belprop_loopy_cg % like cg1


% factor graphs
%fg1   failed since marginals were not exact

fg2
fg3
fg_mrf1
fg_mrf2


% Structure learning
bic1
cooper_yoo
k2demo1
mcmc1
model_select1
pc1
%pc2   failed due to numerical problems in KPMstats/cond_indep_fisher_z




% limids
asia_dt1
id1
oil1
pigs1


% dbns
arhmm1
bat1
bkff1
chmm1
dhmm1
filter_test1
ghmm1
kalman1
kjaerulff1
loopy_dbn1
mhmm1
mildew1
reveal1
viterbi1
water1


% HHMMs
abcd_hhmm
sample_square_hhmm_discrete
%learn_square_hhmm_cts
sample_motif_hhmm

%sparse jtree engine & ndx 2TBN engine
if exist('@jtree_sparse_inf_engine/init_pot','file')
  % only exists if installC has been run
  discrete2
  discrete3 
  filter_test1
  water2
end

%find . -path '*.m' -exec wc -l {} \; | ~/count.pl

% we cannot use tic;toc to time test_BNT, since functions within this script
% reset the tic;toc timer. Hence we use the following:
%clock0=clock; cpu0 = cputime; test_BNT; cpu=cputime-cpu0; elapsed=etime(clock, clock0)
