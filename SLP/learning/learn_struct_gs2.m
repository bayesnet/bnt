function [dag, best_score, cache] = learn_struct_gs2(data, nodesizes, seeddag, varargin)
%
% LEARN_STRUCT_GS2(data,seeddag) learns a structure of Bayesian net by Greedy Search.
% dag = learn_struct_gs(data, nodesizes, seeddag)
%
% dag: the final structure matrix
% Data : training data, data(i,m) is the m obsevation of node i
% Nodesizes: the size array of different nodes
% seeddag: given seed Dag for hill climbing, optional
% cache : data structure used to memorize local score computations 
%   (cf. SCORE_INIT_CACHE function)
%
% by Gang Li @ Deakin University (gli73@hotmail.com)
% (use mk_nbrs_of_dag_topo, developped by Wei Hu, instead of mk_nbrs_of_dag)
% (Caching implementation : ofrancois.olivier.c.h@gmail.com, philippe.leray@univ-nantes.fr)
%

[N ncases] = size(data);
if (nargin < 3 ) 
    seeddag = zeros(N,N); % mk_rnd_dag(N); %call BNT function
elseif ~acyclic(seeddag)
    seeddag = mk_rnd_dag(N); %zeros(N,N);
end;

% set default params
scoring_fn = 'bic';
verbose  = 'yes';
cache=[];

% get params
args = varargin;
nargs = length(args);
if length(args) > 0
    if isstr(args{1})
    	for i = 1:2:nargs
    		switch args{i}
    		case 'scoring_fn', scoring_fn = args{i+1};
    		case 'verbose',  verbose  = strcmp(args{i+1},'yes');
    		case 'cache',  cache=args{i+1} ;
    		end;
    	end;
    end;
end;

done = 0;
[best_score cache] = score_dags(data,nodesizes, {seeddag},'scoring_fn',scoring_fn,'cache',cache);
while ~done
    [dags,op,nodes] = mk_nbrs_of_dag_topo(seeddag);
    nbrs = length(dags);
    [scores cache] = score_dags(data, nodesizes, dags,'scoring_fn',scoring_fn,'cache',cache);
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

outcount = 0; 
[best_score cache] = score_dags(data,nodesizes, {seeddag},'scoring_fn',scoring_fn,'cache',cache);
while outcount < 2
    innercount = 0;
    for i=1:N
        for j=1:N
           if i==j, continue;    end;
           if seeddag(i,j) == 0  % No edge i-->j, then try to add it
               tempdag = seeddag;
               tempdag(i,j) = 1;
               if acyclic(tempdag)
                    [temp_score cache] = score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn,'cache',cache);
                    if temp_score > best_score
                        seeddag = tempdag
                        best_score= temp_score;
                        innercount = innercount +1;
                    end;
               end
           else  % exists edge i--j, then try reverse it or remove it
               tempdag = seeddag;
               tempdag(i,j) = 0; tempdag(j,i) = 1; 
               if acyclic(tempdag)
                   [temp_score cache] = score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn,'cache',cache);
                   if temp_score > best_score
                       seeddag = tempdag;
                       best_score = temp_score;
                       innercount = innercount +1;
                   else
                       tempdag = seeddag;
                       tempdag(i,j) = 0;
                       [temp_score cache] = score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn,'cache',cache);
                       if temp_score > best_score
                           seeddag = tempdag;
                           best_score= temp_score;
                           innercount = innercount +1;
                       end;
                   end;
               else
                   tempdag = seeddag;
                   tempdag(i,j)=0;
                   [temp_score cache] = score_dags(data,nodesizes, {tempdag},'scoring_fn',scoring_fn,'cache',cache);
                   if temp_score > best_score
                       seeddag = tempdag
                       best_score= temp_score;
                       innercount = innercount +1;
                   end;
               end;
           end;
        end; % end for j
    end; % end for i
    if innercount == 0
        outcount = outcount +1;
    end;
end;  % end while

dag = seeddag;
