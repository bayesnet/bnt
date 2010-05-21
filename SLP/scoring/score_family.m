function [score, cache] = score_family(j, ps, node_type, scoring_fn, ns, discrete, data, args, cache)
% SCORE_FAMILY_COMPLETE Compute the score of a node and its parents given completely observed data
% score = score_family(j, ps, node_type, scoring_fn, ns, discrete, data, args, cache)
%
% data(i,m) is the value of node i in case m (can be a cell array, if contain missing value, it uses only available complete cases for an entry)
% args is a cell array containing optional arguments passed to the constructor,
% or is [] if none
% cache is a data structure used to memorize local score computations
% (cf. SCORE_INIT_CACHE function)
%
% We create a whole Bayes net which only connects parents to node,
% where node has a CPD of the specified type (with default parameters).
% We then evaluate its score ('bic' or 'bayesian')
% We should use a cache to avoid unnecessary computation.
% In particular, log_marginal_prob_node for tabular CPDs calls gammaln
% and compute_counts, both of which are slow.
%
% (Caching implementation : ofrancois.olivier.c.h@gmail.com, philippe.leray@univ-nantes.fr)
%

if (nargin<9 | isempty(cache)) , c=0; else c=1; end
%tic
if c==1
    [b,score]=score_find_in_cache(cache,j,ps,scoring_fn);
else
    b=0;
end
%Tfind=toc
ccc = iscell(data);
ps = unique(ps);

if b==0
    misv = -9999;
    if ccc, data = bnt_to_mat(data,misv); end
    [n ncases] = size(data);
    dag = zeros(n,n);
    % SML added to sort ps b/c mk_bnet, learn_params use sorted ps to make
    % CPTs    % Kevin had: if ~isempty(ps), dag(ps, j) = 1; end
    if ~isempty(ps), dag(ps, j) = 1;, ps = sort(ps);, end

    bnet = mk_bnet(dag, ns, 'discrete', discrete);
    fname = sprintf('%s_CPD', node_type);
    if isempty(args)
        bnet.CPD{j} = feval(fname, bnet, j);
    else
        bnet.CPD{j} = feval(fname, bnet, j, args{:});
    end
    %tic
    switch scoring_fn
    case 'bic',
        fam = [ps j];
        if ccc,
	    [tmp, available_case] = find(data(fam,:)==misv);
	    available_case = mysetdiff(1:ncases, available_case);
	else available_case = 1:ncases;
	end
        bnet.CPD{j} = learn_params(bnet.CPD{j}, fam, data(:,available_case), ns, bnet.cnodes);
        %L = log_prob_node(bnet.CPD{j}, data(j,:), data(ps,:));
	L = log_prob_node(bnet.CPD{j}, data(j,available_case), data(ps,available_case));
        S = struct(bnet.CPD{j}); % violate object privacy
        score = L - 0.5*S.nparams*log(length(available_case));
    case 'bicmod',
        fam = [ps j];
        if ccc,
	    [tmp, available_case] = find(data(fam,:)==misv);
	    available_case = mysetdiff(1:ncases, available_case);
	else available_case = 1:ncases;
	end
        bnet.CPD{j} = learn_params(bnet.CPD{j}, fam, data(:,available_case), ns, bnet.cnodes);
	L = log_prob_node(bnet.CPD{j}, data(j,available_case), data(ps,available_case));
        S = struct(bnet.CPD{j}); % violate object privacy
        score = L - S.nparams*log(length(available_case));
    case 'bayesian',
        fam = [ps j];
        if ccc,
	    [tmp, available_case] = find(data(fam,:)==misv);
	    available_case = mysetdiff(1:ncases, available_case);
	else available_case = 1:ncases;
	end
        score = log_marg_prob_node(bnet.CPD{j}, data(j,available_case), data(ps,available_case));
    otherwise,
        error(['unrecognized scoring fn ' scoring_fn]);
    end
    %Tcalc=toc
    %tic
    if c==1
%        fprintf('a\n')
        cache=score_add_to_cache(cache,j,ps,score,scoring_fn);
    end
    %Tecr=toc
% else
%     fprintf('*\n')
end

%===========================Inner functions

function [cache, place] = score_add_to_cache(cache,j,ps,score,scoring_fn)
% [cache place] = score_add_to_cache(cache,j,ps,score,scoring_fn)
%
% j is the son node,
% ps is the list of parents of j, for example [12 5 7],
% score is the score to add for this familly.
% scoring_fn is 'bic' or 'bayesian'.
%
% place = where the entry was add.
%
% example for 2 nodes with cache of size 5 :
%
% cache =
%   Nw  b        0      0      0 --> Nw=number of writings in cache (+1) and b==1 iff the cache is full
%   0   0        1   -239.12   1  --> 1st familly in the cache (node 1 without parents) calculate with bic
%   0   0        2   -318.98   1
%   1   0        2   -189.23   2  --> 3rd familly in the cache (node 2 with 1 as parent) calculate with bayésian
%   0   1        1   -251.09   1
% .ps2bool.      j    score  1or2 --> new entry
%   |   |        |      |      |
%   |   |        |      |      |___> scoring function : 1 for 'bic' or 2 for 'bayesian'
%   |   |        |      |__________> score of the familly
%   |   |        |_________________> son node of the familly
%   |   |__________________________> ==1 iff node 2 is parent of son node
%   |______________________________> ==1 iff node 1 is parent of son node
%
% If the cache is FULL then the new place is RanDomly choose.
%
% V1.1 : 5 may 2003 (O. Francois, Ph. Leray)

N=size(cache,2)-3;
L=size(cache,1)-1;
cache_full=cache(1,2) ;

place=0;

if ismember(j,ps)
    disp('This is a cyclic entry, nothing was done.');
elseif j>N | j<=0
    disp('This entry is not valid, nothing was done.');
else
    switch scoring_fn
    case 'bic',
        fn=1;
    case 'bayesian',
        fn=2;
    otherwise,
        fn=3;
        %error(['unrecognized scoring fn ' scoring_fn]);
    end

    if ~cache_full
        place=cache(1,1);
    else
    [ignore place]=max(rand(1,L)); place=place+1;
    end

    cache(place,:)=0;
    cache(place,ps)=1;
    cache(place,N+1)=j;
    cache(place,N+2)=score;
    cache(place,N+3)=fn;

    cache(1,1)=place+1;
    if place>L | cache(1,2)~=0
        cache(1,2)=1;
    end
end


%=========================================================================================
function [bool, score] = score_find_in_cache(cache,j,ps,scoring_fn)
% cache = score_find_in_cache(cache,j,ps,scoring_fn)
%
% V1.1 : 5 may 2003 (O. Francois, Ph. Leray)


%tic
L=size(cache,1)-1;
N=size(cache,2)-3;

if N<1
    bool=0;
    score=0;
    return
end

parents=zeros(1,N);
parents(ps)=1;
%parents(N+1)=j;

switch scoring_fn
case 'bic',
    fn=1;
case 'bayesian',
    fn=2;
otherwise,
    fn=3;
    %error(['unrecognized scoring fn ' scoring_fn]);
end

tmp=find(cache(2:L+1,N+3)==fn);
tmp=tmp+1;
tmp2=find(cache(tmp,N+1)==j);
candidats=tmp(tmp2);

i=1;
while i<=N & ~isempty(candidats)
    tmp=find(cache(candidats,i)==parents(i));
    candidats=candidats(tmp);
    i=i+1;
end

%Tpre=toc

bool=~isempty(candidats);

if bool
    score=cache(candidats(1),N+2);
else
    score=0;
end
