function [score, cache] = score_dags(data, ns, dags, varargin)
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
% cache      - data structure used to memorize local score computations (cf. SCORE_INIT_CACHE function) [ [] ]
%
% e.g., score = score_dags(data, ns, mk_all_dags(n), 'scoring_fn', 'bic', 'params', [],'cache',cache);
%
% (Caching implementation : francois.olivier.c.h@gmail.com, philippe.leray@univ-nantes.fr)
% ("Clamped" optimisation : Derek Hoiem <dhoiem@cs.cmu.edu>)

[n ncases] = size(data);

% set default params
type = cell(1,n);
params = cell(1,n);
cache=[];
for i=1:n
  type{i} = 'tabular';
  params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 };
end
scoring_fn = 'bayesian';
discrete = 1:n;

isclamped = 0; % DWH
clamped = zeros(n, ncases);
u = [1:ncases]'; % DWH

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'scoring_fn', scoring_fn = args{i+1};
   case 'type',       type = args{i+1};
   case 'discrete',   discrete = args{i+1};
   case 'clamped',    clamped = args{i+1}, isclamped = 1; %DWH
   case 'params',     if isempty(args{i+1}), params = cell(1,n); else params = args{i+1};  end
   case 'cache',      cache=args{i+1} ;
  end
end

NG = length(dags);
score = zeros(1, NG);

for j=1:n
    if isclamped %DWH
        u = find(clamped(j,:)==0);
    end
    for g=1:NG
        if isempty(dags{g})
            score(g)=-Inf;
        else
            ps = parents(dags{g}, j);
            [scor cache] = score_family(j, ps, type{j}, scoring_fn, ns, discrete, data(:,u), params{j}, cache);
            score(g) = score(g) + scor;
        end
    end
end