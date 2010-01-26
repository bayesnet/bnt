function dag = learn_struct_K2(data, ns, order, varargin)
% LEARN_STRUCT_K2 Greedily learn the best structure compatible with a fixed node ordering
% best_dag = learn_struct_K2(data, node_sizes, order, ...)
%
% data(i,m) = value of node i in case m (can be a cell array).
% node_sizes(i) is the size of node i.
% order(i) is the i'th node in the topological ordering.
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% max_fan_in - this the largest number of parents we allow per node [N]
% scoring_fn - 'bayesian' or 'bic' [ 'bayesian' ]
%              Currently, only networks with all tabular nodes support Bayesian scoring.
% type       - type{i} is the type of CPD to use for node i, where the type is a string
%              of the form 'tabular', 'noisy_or', 'gaussian', etc. [ all cells contain 'tabular' ]
% params     - params{i} contains optional arguments passed to the CPD constructor for node i,
%              or [] if none.  [ all cells contain {'prior', 1}, meaning use uniform Dirichlet priors ]
% discrete   - the list of discrete nodes [ 1:N ]
% clamped    - clamped(i,m) = 1 if node i is clamped in case m [ zeros(N, ncases) ]
% verbose    - 'yes' means display output while running [ 'no' ]
%
% e.g., dag = learn_struct_K2(data, ns, order, 'scoring_fn', 'bic', 'params', [])
%
% To be backwards compatible with BNT2, you can also specify arguments as follows
%   dag = learn_struct_K2(data, node_sizes, order, max_fan_in)    
%
% This algorithm is described in
% - Cooper and Herskovits,  "A Bayesian method for the induction of probabilistic
%      networks from data", Machine Learning Journal 9:308--347, 1992

[n ncases] = size(data);

% set default params
type = cell(1,n);
params = cell(1,n);
for i=1:n
  type{i} = 'tabular';
  %params{i} = { 'prior', 1 };
  params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 };
end
scoring_fn = 'bayesian';
discrete = 1:n;
clamped = zeros(n, ncases);

max_fan_in = n;
verbose = 0;

args = varargin;
nargs = length(args);
if length(args) > 0 
  if isstr(args{1})
    for i=1:2:nargs
      switch args{i},
       case 'verbose',    verbose = strcmp(args{i+1}, 'yes');
       case 'max_fan_in', max_fan_in = args{i+1}; 
       case 'scoring_fn', scoring_fn = args{i+1};
       case 'type',       type = args{i+1}; 
       case 'discrete',   discrete = args{i+1}; 
       case 'clamped',    clamped = args{i+1}; 
       case 'params',     if isempty(args{i+1}), params = cell(1,n); else params = args{i+1};  end
      end
    end
  else
    max_fan_in = args{1};
  end
end

dag = zeros(n,n);

for i=1:n
  ps = [];
  j = order(i);
  u = find(clamped(j,:)==0);    
  score = score_family(j, ps, type{j}, scoring_fn, ns, discrete, data(:,u), params{j});
  if verbose, fprintf('\nnode %d, empty score %6.4f\n', j, score); end
  done = 0;
  while ~done & (length(ps) <= max_fan_in)
    pps = mysetdiff(order(1:i-1), ps); % potential parents
    nps = length(pps);
    pscore = zeros(1, nps);
    for pi=1:nps
      p = pps(pi);
      pscore(pi) = score_family(j, [ps p], type{j}, scoring_fn, ns, discrete, data(:,u), params{j});
      if verbose, fprintf('considering adding %d to %d, score %6.4f\n', p, j, pscore(pi)); end
    end
    [best_pscore, best_p] = max(pscore);
    best_p = pps(best_p);
    if best_pscore > score
      score = best_pscore;
      ps = [ps best_p];
      if verbose, fprintf('* adding %d to %d, score %6.4f\n', best_p, j, best_pscore); end
    else
      done = 1;
    end
  end
  if ~isempty(ps) % need this check for matlab 5.2
    dag(ps, j) = 1;
  end
end




