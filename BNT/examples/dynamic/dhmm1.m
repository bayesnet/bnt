% Make an HMM with discrete observations
%   X1 -> X2
%   |     | 
%   v     v
%   Y1    Y2 

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

rand('state', 0);
prior1 = normalise(rand(Q,1));
transmat1 = mk_stochastic(rand(Q,Q));
obsmat1 = mk_stochastic(rand(Q,O));
bnet.CPD{1} = tabular_CPD(bnet, 1, prior1);
bnet.CPD{2} = tabular_CPD(bnet, 2, obsmat1);
bnet.CPD{3} = tabular_CPD(bnet, 3, transmat1);


T = 5; % fixed length sequences

engine = {};
engine{end+1} = jtree_unrolled_dbn_inf_engine(bnet, T);
engine{end+1} = hmm_inf_engine(bnet);
engine{end+1} = smoother_engine(hmm_2TBN_inf_engine(bnet));
engine{end+1} = smoother_engine(jtree_2TBN_inf_engine(bnet));
if 1
%engine{end+1} = frontier_inf_engine(bnet); % broken
engine{end+1} = bk_inf_engine(bnet, 'clusters', {[1]});
engine{end+1} = jtree_dbn_inf_engine(bnet);
end

inf_time = cmp_inference_dbn(bnet, engine, T);

ncases = 2;
max_iter = 2;
[learning_time, CPD, LL, cases] = cmp_learning_dbn(bnet, engine, T, 'ncases', ncases, 'max_iter', max_iter);

% Compare to HMM toolbox

data = zeros(ncases, T);
for i=1:ncases
  %data(i,:) = cat(2, cases{i}{onodes,:});
  data(i,:) = cell2num(cases{i}(onodes,:));
end
[LL2, prior2, transmat2, obsmat2] = dhmm_em(data, prior1, transmat1, obsmat1, 'max_iter', max_iter);

e = 1;
assert(approxeq(prior2, CPD{e,1}.CPT))
assert(approxeq(obsmat2, CPD{e,2}.CPT))
assert(approxeq(transmat2, CPD{e,3}.CPT))
assert(approxeq(LL2, LL{e}))        

