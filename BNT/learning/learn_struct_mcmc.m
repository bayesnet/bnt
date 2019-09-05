function [sampled_graphs, accept_ratio, num_edges] = learn_struct_mcmc(data, ns, varargin)
% MY_LEARN_STRUCT_MCMC  Monte Carlo Markov Chain search over DAGs assuming fully observed data
% [sampled_graphs, accept_ratio, num_edges] = learn_struct_mcmc(data, ns, ...)
% 
% data(i,m) is the value of node i in case m.
% ns(i) is the number of discrete values node i can take on.
%
% sampled_graphs{m} is the m'th sampled graph.
% accept_ratio(t) = acceptance ratio at iteration t
% num_edges(t) = number of edges in model at iteration t
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
% nsamples   - number of samples to draw from the chain after burn-in [ 100*N ]
% burnin     - number of steps to take before drawing samples [ 5*N ]
% init_dag   - starting point for the search [ zeros(N,N) ]
%
% e.g., samples = my_learn_struct_mcmc(data, ns, 'nsamples', 1000);
%
% Modified by Sonia Leach (SML) 2/4/02, 9/5/03



[n ncases] = size(data);


% set default params
type = cell(1,n);
params = cell(1,n);
for i=1:n
 type{i} = 'tabular';
 %params{i} = { 'prior', 1};
 params{i} = { 'prior_type', 'dirichlet', 'dirichlet_weight', 1 };
end
scoring_fn = 'bayesian';
discrete = 1:n;
clamped = zeros(n, ncases);
nsamples = 100*n;
burnin = 5*n;
dag = zeros(n);

args = varargin;
nargs = length(args);
for i=1:2:nargs
 switch args{i},
  case 'nsamples',   nsamples = args{i+1};
  case 'burnin',     burnin = args{i+1};
  case 'init_dag',   dag = args{i+1};
  case 'scoring_fn', scoring_fn = args{i+1};
  case 'type',       type = args{i+1}; 
  case 'discrete',   discrete = args{i+1}; 
  case 'clamped',    clamped = args{i+1}; 
  case 'params',     if isempty(args{i+1}), params = cell(1,n); else params = args{i+1};  end
 end
end

% We implement the fast acyclicity check described by P. Giudici and R. Castelo,
% "Improving MCMC model search for data mining", submitted to J. Machine Learning, 2001.

% SML: also keep descendant matrix C
use_giudici = 1;
if use_giudici
 [nbrs, ops, nodes, A] = mk_nbrs_of_digraph(dag);
else
 [nbrs, ops, nodes] = mk_nbrs_of_dag(dag);
 A = [];
end

num_accepts = 1;
num_rejects = 1;
T = burnin + nsamples;
accept_ratio = zeros(1, T);
num_edges = zeros(1, T);
sampled_graphs = cell(1, nsamples);
%sampled_bitv = zeros(nsamples, n^2);

for t=1:T
 [dag, nbrs, ops, nodes, A, accept] = take_step(dag, nbrs, ops, ...
                    nodes, ns, data, clamped, A, ...
                      scoring_fn, discrete, type, params);
 num_edges(t) = sum(dag(:));
 num_accepts = num_accepts + accept;
 num_rejects = num_rejects + (1-accept);
 accept_ratio(t) =  num_accepts/num_rejects;
 if t > burnin
   sampled_graphs{t-burnin} = dag;
   %sampled_bitv(t-burnin, :) = dag(:)';
 end
end


%%%%%%%%%


function [new_dag, new_nbrs, new_ops, new_nodes, A,  accept] = ...
   take_step(dag, nbrs, ops, nodes, ns, data, clamped, A,  ...
     scoring_fn, discrete, type, params, prior_w)


use_giudici = ~isempty(A);
if use_giudici
 [new_dag, op, i, j, new_A] =  pick_digraph_nbr(dag, nbrs, ops, nodes,A); % updates A
 [new_nbrs, new_ops, new_nodes] =  mk_nbrs_of_digraph(new_dag, new_A);
else
 d = sample_discrete(normalise(ones(1, length(nbrs))));
 new_dag = nbrs{d};
 op = ops{d};
 i = nodes(d, 1); j = nodes(d, 2);
 [new_nbrs, new_ops, new_nodes] = mk_nbrs_of_dag(new_dag);
end

bf =  bayes_factor(dag, new_dag, op, i, j, ns, data, clamped, scoring_fn, discrete, type, params);

%R = bf * (new_prior / prior) * (length(nbrs) / length(new_nbrs)); 
R = bf * (length(nbrs) / length(new_nbrs)); 
u = rand(1,1);
if u > min(1,R) % reject the move
 accept = 0;
 new_dag = dag;
 new_nbrs = nbrs;
 new_ops = ops;
 new_nodes = nodes;
else
 accept = 1;
 if use_giudici
A = new_A; % new_A already updated in pick_digraph_nbr
 end
end


%%%%%%%%%

function bfactor = bayes_factor(old_dag, new_dag, op, i, j, ns, data, clamped, scoring_fn, discrete, type, params)

u = find(clamped(j,:)==0);
LLnew = score_family(j, parents(new_dag, j), type{j}, scoring_fn, ns, discrete, data(:,u), params{j});
LLold = score_family(j, parents(old_dag, j), type{j}, scoring_fn, ns, discrete, data(:,u), params{j});
bf1 = exp(LLnew - LLold);

if strcmp(op, 'rev')  % must also multiply in the changes to i's family
 u = find(clamped(i,:)==0);
 LLnew = score_family(i, parents(new_dag, i), type{i}, scoring_fn, ns, discrete, data(:,u), params{i});
 LLold = score_family(i, parents(old_dag, i), type{i}, scoring_fn, ns, discrete, data(:,u), params{i});
 bf2 = exp(LLnew - LLold);
else
 bf2 = 1;
end
bfactor = bf1 * bf2;


%%%%%%%% Giudici stuff follows %%%%%%%%%%


% SML: This now updates A as it goes from digraph it choses
function [new_dag, op, i, j, new_A] = pick_digraph_nbr(dag, digraph_nbrs, ops, nodes, A)

d = sample_discrete(normalise(ones(1, length(digraph_nbrs))));
%d = myunidrnd(length(digraph_nbrs),1,1);
i = nodes(d, 1); j = nodes(d, 2);
new_dag = digraph_nbrs(:,:,d);
op = ops{d};
new_A = update_ancestor_matrix(A, op, i, j, new_dag); 


%%%%%%%%%%%%%%


function A = update_ancestor_matrix(A,  op, i, j, dag)

switch op
case 'add',
 A = do_addition(A,  op, i, j, dag);
case 'del', 
 A = do_removal(A,  op, i, j, dag);
case 'rev', 
 A = do_removal(A,  op, i, j, dag);
 A = do_addition(A,  op, j, i, dag);
end

 
%%%%%%%%%%%%

function A = do_addition(A, op, i, j, dag)

A(j,i) = 1; % i is an ancestor of j
anci = find(A(i,:));
if ~isempty(anci)
 A(j,anci) = 1; % all of i's ancestors are added to Anc(j)
end
ancj = find(A(j,:));
descj = find(A(:,j)); 
if ~isempty(ancj)
 for k=descj(:)'
   A(k,ancj) = 1; % all of j's ancestors are added to each descendant of j
 end
end

%%%%%%%%%%%
function A = do_removal(A, op, i, j, dag)

% find all the descendants of j, and put them in topological order

% SML: originally Kevin had the next line commented and the %* lines
% being used but I think this is equivalent and much less expensive
% I assume he put it there for debugging and never changed it back...?
descj = find(A(:,j));
%*  R = reachability_graph(dag);
%*  descj = find(R(j,:));

order = topological_sort(dag);

% SML: originally Kevin used the %* line but this was extracting the
% wrong things to sort
%* descj_topnum = order(descj);
[junk, perm] = sort(order); %SML:node i is perm(i)-TH in order
descj_topnum = perm(descj); %SML:descj(i) is descj_topnum(i)-th in order

% SML: now re-sort descj by rank in descj_topnum
[junk, perm] = sort(descj_topnum);
descj = descj(perm);

% Update j and all its descendants
A = update_row(A, j, dag);
for k = descj(:)'
   A = update_row(A, k, dag);
end

%%%%%%%%%%%

function A = old_do_removal(A, op, i, j, dag)

% find all the descendants of j, and put them in topological order
% SML: originally Kevin had the next line commented and the %* lines
% being used but I think this is equivalent and much less expensive
% I assume he put it there for debugging and never changed it back...?
descj = find(A(:,j)); 
%*  R = reachability_graph(dag);
%*  descj = find(R(j,:)); 

order = topological_sort(dag);
descj_topnum = order(descj);
[junk, perm] = sort(descj_topnum);
descj = descj(perm);
% Update j and all its descendants
A = update_row(A, j, dag);
for k = descj(:)'
 A = update_row(A, k, dag);
end

%%%%%%%%%

function A = update_row(A, j, dag)

% We compute row j of A
A(j, :) = 0;
ps = parents(dag, j);
if ~isempty(ps)
 A(j, ps) = 1;
end
for k=ps(:)'
 anck = find(A(k,:));
 if ~isempty(anck)
   A(j, anck) = 1;
 end
end
 
%%%%%%%%

function A = init_ancestor_matrix(dag)

order = topological_sort(dag);
A = zeros(length(dag));
for j=order(:)'
 A = update_row(A, j, dag);
end
