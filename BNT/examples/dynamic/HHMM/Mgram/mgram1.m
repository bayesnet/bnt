% a multigram is a degenerate 2HHMM where the bottom level HMMs emit deterministic strings
% and the the top level abstract states are independent of each other
% cf. HSMM/test_mgram2 

words = {'the', 't', 'h', 'e'};
data = 'the';
nwords = length(words);
word_len = zeros(1, nwords);
word_prob = normalise(ones(1,nwords));
word_logprob = log(word_prob);
for wi=1:nwords
  word_len(wi)=length(words{wi});
end
D = max(word_len);

alphasize = 26;
data = letter2num(data);
T = length(data);

% node numbers
W = 1; % top level state = word id
L = 2; % bottom level state = letter position within word
F = 3;
O = 4;

ss = 4;
intra = zeros(ss,ss);
intra(W,[F L O])=1;
intra(L,[O F])=1;

inter = zeros(ss,ss);
inter(W,W)=1;
inter(L,L)=1;
inter(F,[W L])=1;

% node sizes
ns = zeros(1,ss);
ns(W) = nwords;
ns(L) = D;
ns(F) = 2;
ns(O) = alphasize;


% Make the DBN
bnet = mk_dbn(intra, inter, ns, 'observed', O);
eclass = bnet.equiv_class;



% uniform start distrib over words, uniform trans mat
Wstart = normalise(ones(1,nwords));
Wtrans = mk_stochastic(ones(nwords,nwords));

% always start in state 1 for each bottom level HMM
delta1_start = zeros(1, D);
delta1_start(1) = 1;
Lstart = repmat(delta1_start, nwords, 1);
LRtrans = mk_leftright_transmat(D, 0); % 0 self loop prob
Ltrans = repmat(LRtrans, [1 1 nwords]);

% Finish in the last letter of each word
Fprob = zeros(nwords, D, 2);
Fprob(:,:,1)=1;
for i=1:nwords
  Fprob(i,length(words{i}),2)=1;
  Fprob(i,length(words{i}),1)=0;
end

% Each state uniquely emits a letter
Oprob = zeros(nwords, D, alphasize);
for i=1:nwords
  for l=1:length(words{i})
    a = double(words{i}(l))-96;
    Oprob(i,l,a)=1;
  end
end


% Define CPDs for slice 
bnet.CPD{eclass(W,1)} = tabular_CPD(bnet, W, 'CPT', Wstart);
bnet.CPD{eclass(L,1)} = tabular_CPD(bnet, L, 'CPT', Lstart);
bnet.CPD{eclass(F,1)} = tabular_CPD(bnet, F, 'CPT', Fprob);
bnet.CPD{eclass(O,1)} = tabular_CPD(bnet, O, 'CPT', Oprob);

% Define CPDs for slice 2
bnet.CPD{eclass(W,2)} = hhmmQ_CPD(bnet, W+ss, 'Fbelow', F, 'startprob', Wstart,  'transprob', Wtrans);
bnet.CPD{eclass(L,2)} = hhmmQ_CPD(bnet, L+ss, 'Fself', F, 'Qps', W+ss, 'startprob', Lstart, 'transprob', Ltrans);

evidence = cell(ss,T);
evidence{W,1}=1;
sample = cell2num(sample_dbn(bnet, 'length', T, 'evidence', evidence));
str = lower(sample(4,:))

engine = jtree_dbn_inf_engine(bnet);
evidence = cell(ss,T);
evidence(O,:) = num2cell(data);
[engine, ll_dbn] = enter_evidence(engine, evidence);

gamma = zeros(nwords, T);
for t=1:T
  m = marginal_nodes(engine, [W F], t);
  gamma(:,t) = m.T(:,2);
end
gamma

xidbn = zeros(nwords, nwords);
for t=1:T-1
  m = marginal_nodes(engine, [W F W+ss], t);
  xidbn = xidbn + squeeze(m.T(:,2,:));
end

% thee
% xidbn(1,4)  = 0.9412  the->e
% (2,3)=0.0588 t->h
% (3,4)=0.0588 h-e
% (4,4)=0.0588 e-e
