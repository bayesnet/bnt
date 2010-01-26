
seed = 0;
rand('state', seed);
randn('state', seed);

chars = ['a', 'c', 'g', 't'];
motif = 'accca';
motif_length = length(motif);
motif_code = zeros(1, motif_length);
for i=1:motif_length
  motif_code(i) = find(chars == motif(i));
end

[bnet_init, Qnodes, Fnodes, Onode] = mk_motif_hhmm('motif_length', length(motif));
%[bnet_init, Qnodes, Fnodes, Onode] = mk_motif_hhmm('motif_pattern', motif);
ss = bnet_init.nnodes_per_slice;



% We generate a training set by creating uniform sequences,
% and inserting a single motif at a random location.
ntrain = 100;
T = 20;
cases = cell(1, ntrain);

if 1
  % uniform background 
  background_dist = normalise(ones(1, length(chars)));
end
if 0
  % use a constant background
  background_dist = zeros(1, length(chars));
  m = find(chars=='t');
  background_dist(m) = 1.0;
end
if 0
  % use a background skewed away from the motif
  p = 0.01; q = (1-(2*p))/2;
  background_dist = [p p q q];
end

unif_pos = normalise(ones(1, T-length(motif)));
cases = cell(1, ntrain);
data = zeros(1,T);
for i=1:ntrain
  data = sample_discrete(background_dist, 1, T);
  L = sample_discrete(unif_pos, 1, 1);
  data(L:L+length(motif)-1) = motif_code;
  cases{i} = cell(ss, T);
  cases{i}(Onode,:) = num2cell(data);
end
disp('sample training cases')
for i=1:5
  chars(cell2num(cases{i}(Onode,:)))
end

engine_init = hmm_inf_engine(bnet_init);

[bnet_learned, LL, engine_learned] = ...
    learn_params_dbn_em(engine_init, cases, 'max_iter', 100, 'thresh', 1e-2);
%			'anneal', 1, 'anneal_rate', 0.7);

% extract the learned motif profile
eclass = bnet_learned.equiv_class;
CPDO=struct(bnet_learned.CPD{eclass(Onode,1)});
fprintf('columns = chars, rows = states\n');
profile_learned = squeeze(CPDO.CPT(2,:,:))
[m,ndx] = max(profile_learned, [], 2);
map_motif_learned = chars(ndx)
back_learned = squeeze(CPDO.CPT(1,1,:))'
%map_back_learned = chars(argmax(back_learned))

CPDO_init = struct(bnet_init.CPD{eclass(Onode,1)});
profile_init = squeeze(CPDO_init.CPT(2,:,:));
back_init = squeeze(CPDO_init.CPT(1,1,:))';
