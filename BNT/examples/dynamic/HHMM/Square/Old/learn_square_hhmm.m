% Learn a 3 level HHMM similar to mk_square_hhmm

% Because startprob should be shared for t=1:T,
% but in the DBN is shared for t=2:T, we train using a single long sequence.

discrete_obs = 0;
supervised = 1;
obs_finalF2 = 0;
% It is not possible to observe F2 if we learn
% because the update_ess method for hhmmF_CPD and hhmmQ_CPD assume
% the F nodes are always hidden (for speed).
% However, for generating, we might want to set the final F2=true
% to force all subroutines to finish.

ss = 6;
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];

seed = 1;
rand('state', seed);
randn('state', seed);

if discrete_obs
  Qsizes = [2 4 2];
else
  Qsizes = [2 4 1];
end

D = 3;
Qnodes = 1:D;
startprob = cell(1,D);
transprob = cell(1,D);
termprob = cell(1,D);

startprob{1} = 'unif';
transprob{1} = 'unif';

% In the unsupervised case, it is essential that we break symmetry
% in the initial param estimates.
%startprob{2} = 'unif';
%transprob{2} = 'unif';
%termprob{2} = 'unif';
startprob{2} = 'rnd';
transprob{2} = 'rnd';
termprob{2} = 'rnd';

leftright = 0;
if leftright
  % Initialise base-level models as left-right.
  % If we initialise with delta functions,
  % they will remain delat funcitons after learning
  startprob{3} = 'leftstart';
  transprob{3}  = 'leftright';
  termprob{3} = 'rightstop';
else
  % If we want to be able to run a base-level model backwards...
  startprob{3} = 'rnd';
  transprob{3}  = 'rnd';
  termprob{3} = 'rnd';
end

if discrete_obs
  % Initialise observations of lowest level primitives in a way which we can interpret
  chars = ['L', 'l', 'U', 'u', 'R', 'r', 'D', 'd'];
  L=find(chars=='L'); l=find(chars=='l');
  U=find(chars=='U'); u=find(chars=='u');
  R=find(chars=='R'); r=find(chars=='r');
  D=find(chars=='D'); d=find(chars=='d');
  Osize = length(chars);
  
  p = 0.9;
  obsprob = (1-p)*ones([4 2 Osize]);
  %       Q2 Q3 O
  obsprob(1, 1, L) =  p;
  obsprob(1, 2, l) =  p;
  obsprob(2, 1, U) =  p;
  obsprob(2, 2, u) =  p;
  obsprob(3, 1, R) =  p;
  obsprob(3, 2, r) =  p;
  obsprob(4, 1, D) =  p;
  obsprob(4, 2, d) =  p;
  obsprob = mk_stochastic(obsprob);
  Oargs = {'CPT', obsprob};

else
  % Initialise means of lowest level primitives in a way which we can interpret
  % These means are little vectors in the east, south, west, north directions.
  % (left-right=east, up-down=south, right-left=west, down-up=north)
  Osize = 2;
  mu = zeros(2, Qsizes(2), Qsizes(3));
  noise = 0;
  scale = 3;
  for q3=1:Qsizes(3)
    mu(:, 1, q3) = scale*[1;0] + noise*rand(2,1);
  end
  for q3=1:Qsizes(3)
    mu(:, 2, q3) = scale*[0;-1] + noise*rand(2,1);
  end
  for q3=1:Qsizes(3)
    mu(:, 3, q3) = scale*[-1;0] + noise*rand(2,1);
  end
  for q3=1:Qsizes(3)
    mu(:, 4, q3) = scale*[0;1] + noise*rand(2,1);
  end
  Sigma = repmat(reshape(scale*eye(2), [2 2 1 1 ]), [1 1 Qsizes(2) Qsizes(3)]);
  Oargs = {'mean', mu, 'cov', Sigma, 'cov_type', 'diag'};
end

bnet = mk_hhmm('Qsizes', Qsizes, 'Osize', Osize', 'discrete_obs', discrete_obs,...
	       'Oargs', Oargs, 'Ops', Qnodes(2:3), ...
	       'startprob', startprob, 'transprob', transprob, 'termprob', termprob);

if supervised
  bnet.observed = [Q1 Q2 Onode];
else
  bnet.observed = [Onode];
end

if obs_finalF2
  engine = jtree_dbn_inf_engine(bnet);
  % can't use ndx version because sometimes F2 is hidden, sometimes observed
  error('can''t observe F when learning')
else
  if supervised
    engine = jtree_ndx_dbn_inf_engine(bnet);
  else
    engine = jtree_hmm_inf_engine(bnet);
  end
end
  
if discrete_obs
  % generate some synthetic data (easier to debug)
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
  if leftright % base model is left-right
    ev(Onode,:) = num2cell([R r U u L l D d]);
  else
    ev(Onode,:) = num2cell([r R u U l L d D]);
  end
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

  if 0
    ev = cases{4};
    engine2 = enter_evidence(engine2, ev);
    T = size(ev,2);
    for t=1:T
      m=marginal_family(engine2, F2, t);
      fprintf('t=%d\n', t);
      reshape(m.T, [2 2])
    end
  end
  
  %  [bnet2, LL] = learn_params_dbn_em(engine, cases, 'max_iter', 10);
  long_seq = cat(2, cases{:});
  [bnet2, LL, engine2] = learn_params_dbn_em(engine, {long_seq}, 'max_iter', 200);
  
  % figure out which subsequence each model is responsible for
  mpe = calc_mpe_dbn(engine2, long_seq);
  pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, chars);

else
  load 'square4_cases' % cases{seq}{i,t} for i=1:ss 
  %plot_square_hhmm(cases{1})
  %long_seq = cat(2, cases{:});
  train_cases = cases(1:2);
  long_seq = cat(2, train_cases{:});
  if ~supervised
    T = size(long_seq,2);
    for t=1:T
      long_seq{Q1,t} = [];
      long_seq{Q2,t} = [];
    end
  end
  [bnet2, LL, engine2] = learn_params_dbn_em(engine, {long_seq}, 'max_iter', 100);

  CPDO=struct(bnet2.CPD{eclass(Onode,1)});
  mu = CPDO.mean;
  Sigma = CPDO.cov;
  CPDO_full = CPDO;
  
  % force diagonal covs after training
  for k=1:size(Sigma,3)
    Sigma(:,:,k) = diag(diag(Sigma(:,:,k)));
  end
  bnet2.CPD{6} = set_fields(bnet.CPD{6}, 'cov', Sigma);
  
  if 0
  % visualize each model by concatenating means for each model for nsteps in a row
  nsteps = 5;
  ev = cell(ss, nsteps*prod(Qsizes(2:3)));
  t = 1;
  for q2=1:Qsizes(2)
    for q3=1:Qsizes(3)
      for i=1:nsteps
	ev{Onode,t} = mu(:,q2,q3);
	ev{Q2,t} = q2;
	t = t + 1;
      end
    end
  end
  plot_square_hhmm(ev)      
  end

  % bnet3 is the same as the learned model, except we will use it in testing mode
  if supervised
    bnet3 = bnet2;
    bnet3.observed = [Onode];
    engine3 = hmm_inf_engine(bnet3);
    %engine3 = jtree_ndx_dbn_inf_engine(bnet3);
  else
    bnet3 = bnet2;
    engine3 = engine2;
  end
  
  if 0
  % segment whole sequence
  mpe = calc_mpe_dbn(engine3, long_seq);
  pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, []);
  end
  
  % segment each sequence
  test_cases = cases(3:4);
  for i=1:2
    ev = test_cases{i};
    T = size(ev, 2);
    for t=1:T
      ev{Q1,t} = [];
      ev{Q2,t} = [];
    end
    mpe = calc_mpe_dbn(engine3, ev);
    subplot(1,2,i)
    plot_square_hhmm(mpe)      
    %pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, []);
    q1s = cell2num(mpe(Q1,:));
    h = hist(q1s, 1:Qsizes(1));
    map_q1 = argmax(h);
    str = sprintf('test seq %d is of type %d\n', i, map_q1);
    title(str)
  end

end

if 0
% Estimate gotten by couting transitions in the labelled data
% Note that a self transition shouldnt count if F2=off.
Q2ev = cell2num(ev(Q2,:));
Q2a = Q2ev(1:end-1);
Q2b = Q2ev(2:end);
counts = compute_counts([Q2a; Q2b], [4 4]);
end

eclass = bnet2.equiv_class;
CPDQ1=struct(bnet2.CPD{eclass(Q1,2)});
CPDQ2=struct(bnet2.CPD{eclass(Q2,2)});
CPDQ3=struct(bnet2.CPD{eclass(Q3,2)});
CPDF2=struct(bnet2.CPD{eclass(F2,1)});
CPDF3=struct(bnet2.CPD{eclass(F3,1)});


A=add_hhmm_end_state(CPDQ2.transprob, CPDF2.termprob(:,:,2));
squeeze(A(:,1,:))
squeeze(A(:,2,:))
CPDQ2.startprob
 
if 0
S=struct(CPDF2.sub_CPD_term);
S.nsamples
reshape(S.counts, [2 4 2])
end
