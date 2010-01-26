function score = score_dags(data, ns, dags, varargin)
% SCORE_DAGS Compute the score of one or more DAGs
% score = score_dags(data, ns, dags, varargin)
%
% data{i,m} = value of node i in case m (can be a cell array).
% node_sizes(i) is the number of size of node i.
% dags{g} is the g'th dag
% score(g) is the score of the i'th dag
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% scoring_fn - 'bayesian' or 'bic' [ 'bayesian' ]
%              Currently, only networks with all tabular nodes support Bayesian scoring.
% type       - type{i} is the type of CPD to use for node i, where the type is a string
%              of the form 'tabular', 'noisy_or', 'gaussian', etc. [ all cells contain 'tabular' ]
% params     - params{i} contains optional arguments passed to the CPD constructor for node i,
%              or [] if none.  [ all cells contain {'prior', 1}, meaning use uniform Dirichlet priors ]
% discrete   - the list of discrete nodes [ 1:N ]
% clamped    - clamped(i,m) = 1 if node i is clamped in case m [ zeros(N, ncases) ]
%
% e.g., score = score_dags(data, ns, mk_all_dags(n), 'scoring_fn', 'bic', 'params', []);
%
% If the DAGs have a lot of families in common, we can cache the sufficient statistics,
% making this potentially more efficient than scoring the DAGs one at a time.
% (Caching is not currently implemented, however.)

[n ncases] = size(data);

% set default params
type = cell(1,n);
params = cell(1,n);
for i=1:n
  type{i} = 'tabular';
  params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 };
end
scoring_fn = 'bayesian';
discrete = 1:n;

u = [1:ncases]'; % DWH
isclamped = 0; %DWH
clamped = zeros(n, ncases);

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'scoring_fn', scoring_fn = args{i+1};
   case 'type',       type = args{i+1}; 
   case 'discrete',   discrete = args{i+1}; 
   case 'clamped',    clamped = args{i+1}, isclamped = 1; %DWH
   case 'params',     if isempty(args{i+1}), params = cell(1,n); else params = args{i+1};  end
  end
end

NG = length(dags);
score = zeros(1, NG);
for g=1:NG
  dag = dags{g};
  for j=1:n
    if isclamped %DWH
        u = find(clamped(j,:)==0);    
    end
    ps = parents(dag, j);
    score(g) = score(g) + score_family(j, ps, type{j}, scoring_fn, ns, discrete, data(:,u), params{j});
  end
end
