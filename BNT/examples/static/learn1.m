% Lawn sprinker example from Russell and Norvig p454
% See www.cs.berkeley.edu/~murphyk/Bayes/usage.html for details.

N = 4; 
dag = zeros(N,N);
C = 1; S = 2; R = 3; W = 4;
dag(C,[R S]) = 1;
dag(R,W) = 1;
dag(S,W)=1;

false = 1; true = 2;
ns = 2*ones(1,N); % binary nodes

bnet = mk_bnet(dag, ns);
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);

CPT = cell(1,N);
for i=1:N
  s=struct(bnet.CPD{i});  % violate object privacy
  CPT{i}=s.CPT;
end

% Generate training data
nsamples = 50;
samples = cell(N, nsamples);
for i=1:nsamples
  samples(:,i) = sample_bnet(bnet);
end
data = cell2num(samples);

% Make a tabula rasa
bnet2 = mk_bnet(dag, ns);
seed = 0;
rand('state', seed);
bnet2.CPD{C} = tabular_CPD(bnet2, C, 'clamped', 1, 'CPT', [0.5 0.5], ...
			   'prior_type', 'dirichlet', 'dirichlet_weight', 0);
bnet2.CPD{R} = tabular_CPD(bnet2, R, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
bnet2.CPD{S} = tabular_CPD(bnet2, S, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);
bnet2.CPD{W} = tabular_CPD(bnet2, W, 'prior_type', 'dirichlet', 'dirichlet_weight', 0);


% Find MLEs from fully observed data
bnet4 = learn_params(bnet2, samples);

% Bayesian updating with 0 prior is equivalent to ML estimation
bnet5 = bayes_update_params(bnet2, samples);

CPT4 = cell(1,N);
for i=1:N
  s=struct(bnet4.CPD{i});  % violate object privacy
  CPT4{i}=s.CPT;
end

CPT5 = cell(1,N);
for i=1:N
  s=struct(bnet5.CPD{i});  % violate object privacy
  CPT5{i}=s.CPT;
  assert(approxeq(CPT5{i}, CPT4{i}))
end


if 1
% Find MLEs from partially observed data

% hide 50% of the nodes
samplesH = samples;
hide = rand(N, nsamples) > 0.5;
[I,J]=find(hide);
for k=1:length(I)
  samplesH{I(k), J(k)} = [];
end

engine = jtree_inf_engine(bnet2);
max_iter = 5;
[bnet6, LL] = learn_params_em(engine, samplesH, max_iter);

CPT6 = cell(1,N);
for i=1:N
  s=struct(bnet6.CPD{i});  % violate object privacy
  CPT6{i}=s.CPT;
end

end
