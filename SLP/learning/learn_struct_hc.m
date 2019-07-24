function [dag,best_score] = learn_struct_hc(data, nodesizes, seeddag, varargin)
%
% LEARN_STRUCT_HC(data,seeddag) learns a structure of Bayesian net by Hill Climbing.
% dag = learn_struct_hc(data, nodesizes, seeddag)
%
% dag: the final structurre matrix
% Data : training data, data(i,m) is the m obsevation of node i
% Nodesizes: the size array of different nodes
% seeddag: given seed Dag for hill climbing, optional
%
% by Gang Li @ Deakin University (gli73@hotmail.com)

[N ncases] = size(data);
if (nargin < 3 ) 
    seeddag = zeros(N,N); % mk_rnd_dag(N); %call BNT function
elseif ~acyclic(seeddag)
    seeddag = mk_rnd_dag(N); %zeros(N,N);
end;

% set default params
scoring_fn = 'bic';
verbose  = 'yes';

% get params
args = varargin;
nargs = length(args);
if length(args) > 0
    if isstr(args{1})
    	for i = 1:2:nargs
    		switch args{i}
    		case 'scoring_fn', scoring_fn = args{i+1};
    		case 'verbose',  verbose  = strcmp(args{i+1},'yes');
    		end;
    	end;
    end;
end;

done = 0;
best_score = score_dags(data,nodesizes, {seeddag},'scoring_fn',scoring_fn);
while ~done
    [dags,op,nodes] = mk_nbrs_of_dag(seeddag);
    nbrs = length(dags);
    scores = score_dags(data, nodesizes, dags,'scoring_fn',scoring_fn);
    max_score = max(scores);
    new = find(scores == max_score );
    if ~isempty(new) & (max_score > best_score)
        p = sample_discrete(normalise(ones(1, length(new))));
        best_score = max_score;
        seeddag = dags{new(p)};
    else
        done = 1;
    end;
end;

dag = seeddag;

