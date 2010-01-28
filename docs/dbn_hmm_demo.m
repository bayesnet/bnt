 Example due to Wang Hee Lin" <engp1622@nus.edu.sg


intra = zeros(2);
intra(1,2) = 1; 
inter = zeros(2);
inter(1,1) = 1; 

Q = 2; % num hidden states
O = 2; % num observable symbols
ns = [Q O];%number of states
dnodes = 1:2;
%onodes = [1:2]; % only possible with jtree, not hmm
onodes = [2]; 
bnet = mk_dbn(intra, inter, ns, 'discrete', dnodes, 'observed', onodes);
for i=1:4
  bnet.CPD{i} = tabular_CPD(bnet, i);
end

prior0 = normalise(rand(Q,1));
transmat0 = mk_stochastic(rand(Q,Q));
obsmat0 = mk_stochastic(rand(Q,O));

%engine = smoother_engine(hmm_2TBN_inf_engine(bnet));
engine = smoother_engine(jtree_2TBN_inf_engine(bnet));

ss = 2;%slice size(ss)
ncases = 10;%number of examples
T=10;
max_iter=2;%iterations for EM
cases = cell(1, ncases);
for i=1:ncases
  ev = sample_dbn(bnet, T);
  cases{i} = cell(ss,T);
  cases{i}(onodes,:) = ev(onodes, :);
end
[bnet2, LLtrace] = learn_params_dbn_em(engine, cases, 'max_iter', 4);
