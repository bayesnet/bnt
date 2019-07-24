% Compute Viterbi path discrete HMM by different methods

intra = zeros(2);
intra(1,2) = 1;
inter = zeros(2);
inter(1,1) = 1;
n = 2;

Q = 2; % num hidden states
O = 2; % num observable symbols

ns = [Q O];
dnodes = 1:2;
onodes = [2];
eclass1 = [1 2];
eclass2 = [3 2];
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'eclass1', eclass1, 'eclass2', eclass2, ...
	      'observed', onodes);

for seed=1:10
rand('state', seed);
prior = normalise(rand(Q,1));
transmat = mk_stochastic(rand(Q,Q));
obsmat = mk_stochastic(rand(Q,O));
bnet.CPD{1} = tabular_CPD(bnet, 1, prior);
bnet.CPD{2} = tabular_CPD(bnet, 2, obsmat);
bnet.CPD{3} = tabular_CPD(bnet, 3, transmat);


% Create a sequence
T = 5; 
ev = sample_dbn(bnet, T);
evidence = cell(2,T);
evidence(2,:) = ev(2,:); % extract observed component
data = cell2num(ev(2,:));

%obslik = mk_dhmm_obs_lik(data, obsmat);
obslik = multinomial_prob(data, obsmat);
path = viterbi_path(prior, transmat, obslik);

engine = {};
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));

mpe = find_mpe(engine{1}, evidence);

assert(isequal(cell2num(mpe(1,:)), path)) % extract values of hidden nodes
end
