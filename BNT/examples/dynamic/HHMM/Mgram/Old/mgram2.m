% like mgram1, except we use a durational HMM instead of an HHMM2

past = 0;

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


alphasize = 26*2;
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
inter(F,[W L O])=1;

% node sizes
ns = zeros(1,ss);
ns(W) = nwords;
ns(L) = D;
ns(F) = 2;
ns(O) = alphasize;
ns2 = [ns ns];

% Make the DBN
bnet = mk_dbn(intra, inter, ns, 'observed', O);
eclass = bnet.equiv_class;

% uniform start distrib over words, uniform trans mat
Wstart = normalise(ones(1,nwords));
Wtrans = mk_stochastic(ones(nwords,nwords));

% always start in state d = length(word) for each bottom level HMM
Lstart = zeros(nwords, D);
for i=1:nwords
  l = length(words{i});
  Lstart(i,l)=1;
end

% make downcounters
RLtrans = mk_rightleft_transmat(D, 0); % 0 self loop prob
Ltrans = repmat(RLtrans, [1 1 nwords]);

% Finish when downcoutner = 1
Fprob = zeros(nwords, D, 2);
Fprob(:,1,2)=1;
Fprob(:,2:end,1)=1;


% Define CPDs for slice 
bnet.CPD{eclass(W,1)} = tabular_CPD(bnet, W, 'CPT', Wstart);
bnet.CPD{eclass(L,1)} = tabular_CPD(bnet, L, 'CPT', Lstart);
bnet.CPD{eclass(F,1)} = tabular_CPD(bnet, F, 'CPT', Fprob);


% Define CPDs for slice 2
bnet.CPD{eclass(W,2)} = hhmmQ_CPD(bnet, W+ss, 'Fbelow', F, 'startprob', Wstart,  'transprob', Wtrans);
bnet.CPD{eclass(L,2)} = hhmmQ_CPD(bnet, L+ss, 'Fself', F, 'Qps', W+ss, 'startprob', Lstart, 'transprob', Ltrans);


if 0
% To test it is generating correctly, we create an artificial
% observation process that capitalizes at the start of a new segment
% Oprob(Ft-1,Qt,Dt,Yt)
Oprob = zeros(2,nwords,D,alphasize);
Oprob(1,1,3,letter2num('t'),1)=1;
Oprob(1,1,2,letter2num('h'),1)=1;
Oprob(1,1,1,letter2num('e'),1)=1;
Oprob(2,1,3,letter2num('T'),1)=1;
Oprob(2,1,2,letter2num('H'),1)=1;
Oprob(2,1,1,letter2num('E'),1)=1;
Oprob(1,2,1,letter2num('a'),1)=1;
Oprob(2,2,1,letter2num('A'),1)=1;
Oprob(1,3,1,letter2num('b'),1)=1;
Oprob(2,3,1,letter2num('B'),1)=1;
Oprob(1,4,1,letter2num('c'),1)=1;
Oprob(2,4,1,letter2num('C'),1)=1;

% Oprob1(Qt,Dt,Yt)
Oprob1 = zeros(nwords,D,alphasize);
Oprob1(1,3,letter2num('t'),1)=1;
Oprob1(1,2,letter2num('h'),1)=1;
Oprob1(1,1,letter2num('e'),1)=1;
Oprob1(2,1,letter2num('a'),1)=1;
Oprob1(3,1,letter2num('b'),1)=1;
Oprob1(4,1,letter2num('c'),1)=1;

bnet.CPD{eclass(O,2)} = tabular_CPD(bnet, O+ss, 'CPT', Oprob);
bnet.CPD{eclass(O,1)} = tabular_CPD(bnet, O, 'CPT', Oprob1);

evidence = cell(ss,T);
%evidence{W,1}=1;
sample = cell2num(sample_dbn(bnet, 'length', T, 'evidence', evidence));
str = num2letter(sample(4,:))
end




[log_obslik, obslik, match] = mk_mgram_obslik(lower(data), words, word_len, word_prob);
% obslik(j,t,d)
softCPDpot = cell(ss,T);
ens = ns;
ens(O)=1;
ens2 = [ens ens];
for t=2:T
  dom = [F W+ss L+ss O+ss];
  % tab(Ft-1, Q2, Dt)
  tab = ones(2, nwords, D);
  if past
    tab(1,:,:)=1; % if haven't finished previous word, likelihood is 1
    tab(2,:,:) = squeeze(obslik(:,t,:)); % otherwise likelihood of this segment
  else
    for d=1:max(1,min(D,T+1-t))
      tab(2,:,d) = squeeze(obslik(:,t+d-1,d));
    end
  end
  softCPDpot{O,t} = dpot(dom, ens2(dom), tab);
end
t = 1;
dom = [W L O];
% tab(Q2, Dt)
tab = ones(nwords, D);
if past
  tab = squeeze(obslik(:,t,:));
else
  for d=1:min(D,T-t)
    tab(:,d) = squeeze(obslik(:,t+d-1,d));
  end
end
softCPDpot{O,t} = dpot(dom, ens(dom), tab);


%bnet.observed = [];
% uniformative observations
%bnet.CPD{eclass(O,2)} = tabular_CPD(bnet, O+ss, 'CPT', mk_stochastic(ones(2,nwords,D,alphasize)));
%bnet.CPD{eclass(O,1)} = tabular_CPD(bnet, O, 'CPT', mk_stochastic(ones(nwords,D,alphasize)));

engine = jtree_dbn_inf_engine(bnet);
evidence = cell(ss,T);
% we add dummy data to O to force its effective size to be 1.
% The actual values have already been incorporated into softCPDpot 
evidence(O,:) = num2cell(ones(1,T));
[engine, ll_dbn] = enter_evidence(engine, evidence, 'softCPDpot', softCPDpot);


%evidence(F,:) = num2cell(2*ones(1,T));
%[engine, ll_dbn] = enter_evidence(engine, evidence);


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


