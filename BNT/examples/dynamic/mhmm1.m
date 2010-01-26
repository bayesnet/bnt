% Make an HMM with mixture of Gaussian observations
%    Q1 ---> Q2
%  /  |   /  |
% M1  |  M2  | 
%  \  v   \  v
%    Y1     Y2 
% where Pr(m=j|q=i) is a multinomial and Pr(y|m,q) is a Gaussian     

%seed = 3;
%rand('state', seed);
%randn('state', seed);

intra = zeros(3);
intra(1,[2 3]) = 1;
intra(2,3) = 1;
inter = zeros(3);
inter(1,1) = 1;
n = 3;

Q = 2; % num hidden states
O = 2; % size of observed vector
M = 2; % num mixture components per state

ns = [Q M O];
dnodes = [1 2];
onodes = [3];
eclass1 = [1 2 3];
eclass2 = [4 2 3];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
	      'observed', onodes);

prior0 = normalise(rand(Q,1));
transmat0 = mk_stochastic(rand(Q,Q));
mixmat0 = mk_stochastic(rand(Q,M));
mu0 = rand(O,Q,M);
Sigma0 = repmat(eye(O), [1 1 Q M]);
bnet.CPD{1} = tabular_CPD(bnet, 1, prior0);
bnet.CPD{2} = tabular_CPD(bnet, 2, mixmat0);
%% we set the cov prior to 0 to give same results as HMM toolbox
%bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', mu0, 'cov', Sigma0, 'cov_prior_weight', 0);
% new version of HMM toolbox uses the same default prior on Gaussians as BNT
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean', mu0, 'cov', Sigma0);
bnet.CPD{4} = tabular_CPD(bnet, 4, transmat0);



T = 5; % fixed length sequences

engine = {};
engine{end+1} = hmm_inf_engine(bnet);
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet));
if 0
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
%engine{end+1} = frontier_inf_engine(bnet);
engine{end+1} = bk_inf_engine(bnet, 'clusters', 'exact');
engine{end+1} = jtree_dbn_inf_engine(bnet);
end

inf_time = cmp_inference_dbn(bnet, engine, T);

ncases = 2;
max_iter = 2;
[learning_time, CPD, LL, cases] = cmp_learning_dbn(bnet, engine, T, 'ncases', ncases, 'max_iter', max_iter);

% Compare to HMM toolbox

data = zeros(O, T, ncases);
for i=1:ncases
  data(:,:,i) = reshape(cell2num(cases{i}(onodes,:)), [O T]);
end
tic;
[LL2, prior2, transmat2, mu2, Sigma2, mixmat2] = ...
    mhmm_em(data, prior0, transmat0,  mu0, Sigma0, mixmat0, 'max_iter', max_iter);
t=toc;
disp(['HMM toolbox took ' num2str(t) ' seconds '])

for e = 1:length(engine)
  assert(approxeq(prior2, CPD{e,1}.CPT))
  assert(approxeq(mixmat2, CPD{e,2}.CPT))
  assert(approxeq(mu2, CPD{e,3}.mean))
  assert(approxeq(Sigma2, CPD{e,3}.cov))
  assert(approxeq(transmat2, CPD{e,4}.CPT))
  assert(approxeq(LL2, LL{e}))
end
