% Try to learn a 3 level HHMM similar to mk_square_hhmm
% from hand-drawn squares.

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

seed = 1;
rand('state', seed);
randn('state', seed);

bnet = mk_square_hhmm(discrete_obs, 0);
 
ss = 6;
Q1 = 1; Q2 = 2; Q3 = 3; F3 = 4; F2 = 5; Onode = 6;
Qnodes = [Q1 Q2 Q3]; Fnodes = [F2 F3];
Qsizes = [2 4 1];
  
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
[bnet2, LL, engine2] = learn_params_dbn_em(engine, {long_seq}, 'max_iter', 2);

eclass = bnet2.equiv_class;
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
  %mpe = calc_mpe_dbn(engine3, ev);
  mpe = find_mpe(engine3, ev)
  subplot(1,2,i)
  plot_square_hhmm(mpe)      
  %pretty_print_hhmm_parse(mpe, Qnodes, Fnodes, Onode, []);
  q1s = cell2num(mpe(Q1,:));
  h = hist(q1s, 1:Qsizes(1));
  map_q1 = argmax(h);
  str = sprintf('test seq %d is of type %d\n', i, map_q1);
  title(str)
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
squeeze(A(:,1,:));
CPDQ2.startprob;
 
if 0
S=struct(CPDF2.sub_CPD_term);
S.nsamples
reshape(S.counts, [2 4 2])
end
