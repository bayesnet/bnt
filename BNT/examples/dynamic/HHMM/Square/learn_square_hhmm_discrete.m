% Try to learn a 3 level HHMM similar to mk_square_hhmm
% from synthetic discrete sequences


discrete_obs = 1;
supervised = 0;
obs_finalF2 = 0;

seed = 1;
rand('state', seed);
randn('state', seed);

bnet_init = mk_square_hhmm(discrete_obs, 0);

ss = 6;
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];

if supervised
  bnet_init.observed = [Q1 Q2 Onode];
else
  bnet_init.observed = [Onode];
end

if obs_finalF2
  engine_init = jtree_dbn_inf_engine(bnet_init);
  % can't use ndx version because sometimes F2 is hidden, sometimes observed
  error('can''t observe F when learning')
  % It is not possible to observe F2 if we learn
  % because the update_ess method for hhmmF_CPD and hhmmQ_CPD assume
  % the F nodes are always hidden (for speed).
  % However, for generating, we might want to set the final F2=true
  % to force all subroutines to finish.
else
  if supervised
    engine_init = jtree_ndx_dbn_inf_engine(bnet_init);
  else
    engine_init = hmm_inf_engine(bnet_init);
  end
end
  
% generate some synthetic data (easier to debug)
chars = ['L', 'l', 'U', 'u', 'R', 'r', 'D', 'd'];
L=find(chars=='L'); l=find(chars=='l');
U=find(chars=='U'); u=find(chars=='u');
R=find(chars=='R'); r=find(chars=='r');
D=find(chars=='D'); d=find(chars=='d');

cases = {};

T = 8;
ev = cell(ss, T);
ev(Onode,:) = num2cell([L l U u R r D d]);
if supervised
  ev(Q1,:) = num2cell(1*ones(1,T));
  ev(Q2,:) = num2cell( [1 1 2 2 3 3 4 4]);
end
cases{1} = ev;
cases{3} = ev;

T  = 8;
ev = cell(ss, T);
%we start with R then r, even though we are running the model 'backwards'!
ev(Onode,:) = num2cell([R r U u L l D d]);

if supervised
  ev(Q1,:) = num2cell(2*ones(1,T));
  ev(Q2,:) = num2cell( [3 3 2 2 1 1 4 4]);
end

cases{2} = ev;
cases{4} = ev;

if obs_finalF2
  for i=1:length(cases)
    T = size(cases{i},2);
    cases{i}(F2,T)={2}; % force F2 to be finished at end of seq
  end
end


% startprob should be shared for t=1:T,
% but in the DBN it is shared for t=2:T,
% so we train using a single long sequence.
long_seq = cat(2, cases{:});
[bnet_learned, LL, engine_learned] = ...
    learn_params_dbn_em(engine_init, {long_seq}, 'max_iter', 200);

% figure out which subsequence each model is responsible for
mpe = calc_mpe_dbn(engine_learned, long_seq);
pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, chars);


% The "true" segmentation of the training sequence  is
% Q1: 1                 2
% O:  L l U u R r D d | R r U u L l D d | etc.
% 
% When we learn in a supervised fashion, we recover the "truth".

% When we learn in an unsupervised fashion with seed=1, we get
% Q1: 2                       1
% O:  L l U u R r D d  R r | U u L l D d | etc.
%
% This means for model 1:
% starts in state 2
% transitions 2->1, 1->4, 4->e, 3->2
%
% For model 2,
% starts in state 1
% transitions 1->2, 2->3, 3->4 or e, 4->3

% examine the params
eclass = bnet_learned.equiv_class;
CPDQ1=struct(bnet_learned.CPD{eclass(Q1,2)});
CPDQ2=struct(bnet_learned.CPD{eclass(Q2,2)});
CPDQ3=struct(bnet_learned.CPD{eclass(Q3,2)});
CPDF2=struct(bnet_learned.CPD{eclass(F2,1)});
CPDF3=struct(bnet_learned.CPD{eclass(F3,1)});
CPDO=struct(bnet_learned.CPD{eclass(Onode,1)});

A_learned =add_hhmm_end_state(CPDQ2.transprob, CPDF2.termprob(:,:,2));
squeeze(A_learned(:,1,:))
squeeze(A_learned(:,2,:))


% Does the "true" model have higher likelihood than the learned one?
% i.e., Does the unsupervised method learn the wrong model because
% we have the wrong cost fn, or because of local minima?

bnet_true = mk_square_hhmm(discrete_obs,1);

% examine the params
eclass = bnet_learned.equiv_class;
CPDQ1_true=struct(bnet_true.CPD{eclass(Q1,2)});
CPDQ2_true=struct(bnet_true.CPD{eclass(Q2,2)});
CPDQ3_true=struct(bnet_true.CPD{eclass(Q3,2)});
CPDF2_true=struct(bnet_true.CPD{eclass(F2,1)});
CPDF3_true=struct(bnet_true.CPD{eclass(F3,1)});

A_true =add_hhmm_end_state(CPDQ2_true.transprob, CPDF2_true.termprob(:,:,2));
squeeze(A_true(:,1,:))


if supervised
  engine_true = jtree_ndx_dbn_inf_engine(bnet_true);
else
  engine_true = hmm_inf_engine(bnet_true);
end

%[engine_learned, ll_learned] = enter_evidence(engine_learned, long_seq);
%[engine_true, ll_true] = enter_evidence(engine_true, long_seq);
[engine_learned, ll_learned] = enter_evidence(engine_learned, cases{2});
[engine_true, ll_true] = enter_evidence(engine_true, cases{2});
ll_learned
ll_true


% remove concatentation artefacts
ll_learned = 0;
ll_true = 0;
for m=1:length(cases)
  [engine_learned, ll_learned_tmp] = enter_evidence(engine_learned, cases{m});
  [engine_true, ll_true_tmp] = enter_evidence(engine_true, cases{m});
  ll_learned = ll_learned + ll_learned_tmp;
  ll_true = ll_true + ll_true_tmp;
end
ll_learned
ll_true

% In both cases, ll_learned >> ll_true
% which shows we are using the wrong cost function!
