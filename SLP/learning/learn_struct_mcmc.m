function [sampled_graphs, accept_ratio, num_edges] = learn_struct_mcmc(data, ns, varargin)
% LEARN_STRUCT_MCMC  Monte Carlo Markov Chain search over DAGs assuming fully observed data
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
% 
% Modified by Mingyi Wang  (mingyiwang@hotmail.com) Sep 18, 2006 (based on Sonia Leach (SML)'s version ( 2/4/02, 9/5/03))
%
% Some bugs in update_ancestor_matrix() were fixed. This function can call mk_nbrs_of_digraph properly
% 

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
  case 'gconstraint', gconstraint=args{i+1};  %Added by mingyi
  case 'params',     if isempty(args{i+1}), params = cell(1,n); else params = args{i+1};  end
    
 end
end

% We implement the fast acyclicity check described by P. Giudici and R. Castelo,
% "Improving MCMC model search for data mining", submitted to J. Machine Learning, 2001.

% SML: also keep descendant matrix C
use_giudici = 1;
%use_giudici = 0; %Revised by MIngyi
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
 fprintf('MCMC: %d/%d\n',t,T);
end


%%%%%%%%%


function [new_dag, new_nbrs, new_ops, new_nodes, A,  accept] = ...
   take_step(dag, nbrs, ops, nodes, ns, data, clamped, A,  ...
     scoring_fn, discrete, type, params, prior_w)

global gconstraint;   %Added by Mingyi
use_giudici = ~isempty(A);
if use_giudici
 [new_dag, op, i, j, new_A] =  pick_digraph_nbr(dag, nbrs, ops, nodes,A); % updates A
 [new_nbrs, new_ops, new_nodes] =  mk_nbrs_of_digraph(new_dag,new_A);  
else
 d = sample_discrete(normalise(ones(1, length(nbrs))));
 new_dag = nbrs{d};
 op = ops{d};
 i = nodes(d, 1); j = nodes(d, 2);
 [new_nbrs, new_ops, new_nodes] = mk_nbrs_of_dag1(new_dag);   
end
%For debug
% fprintf('op:%s,i:%d,j:%d\n',op,i,j);
% if ~acyclic(new_dag)
%     error('new dag must be acyclic!')
% end
% if size(find(diag(new_A)),1)>0
%   A=A
%   new_A=new_A
%   error('new A must be acyclic!')
%  end
%debug ends

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
new_A = update_ancestor_matrix(A, op, i, j, dag); 
%for debug
% if op=='add'
%     if ~(dag(i,j)==0 & new_dag(i,j)==1)
%         fprintf('error add\n');
%     end
% end
% if op=='del'
%     if ~(dag(i,j)==1 & new_dag(i,j)==0)
%       fprintf('new dag del calculation is error!\n')
%     end
% end
% if op=='rev'
%     if ~(dag(i,j)==1 & dag(j,i)==0 & new_dag(i,j)==0 & new_dag(j,i)==1)
%       fprintf('new dag rev calculation is error!\n')
%     end
% end
% new_AA = reachability_graph(new_dag');
% if find(diag(new_AA)==1)
%     fprintf('cyclic\n');
% end
% if ~isequal(new_A,new_AA)
%    fprintf('new A calculation is error!\n')
% end
%debug ends

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

A(j,i) = 1;     % i is an ancestor of j
anci = find(A(i,:));
if ~isempty(anci)
 A(j,anci) = 1;   % all of i's ancestors are added to Anc(j)
end

descj = find(A(:,j));  %all the descendants of j are selected 
if ~isempty(descj)
 for k=descj(:)'
   A(k,i) = 1;        % i is the ancestor of descj
   if ~isempty(anci)  % all of i's ancestors are also the ancestor of each descendant of j
       A(k,anci)=1;
   end   
 end
end


%%%%%%%%%%%

function A = do_removal(A, op, i, j, dag)
descj = find(A(:,j)); 
A = update_row(A,i, j, dag);   % compute the A(j,:) row for dag i->j removal

if ~isempty(descj) 
  order = topological_sort(dag);  %all the parent nodes are before to the children nodes
  [junk, perm] = sort(order);     %node i is perm(i)-TH in order
  descj_topnum = perm(descj);     %descj(i) is descj_topnum(i)-th in order

% SML: now re-sort descj by rank in descj_topnum
  [junk, perm] = sort(descj_topnum);
  descj = descj(perm); 
  for k = descj(:)'
    A = old_update_row(A, k, dag);
  end
end

%%%%%%%%%

function A = update_row(A, i,j, dag)
% We compute row j of A
A(j, :) = 0;
ps = parents(dag, j);
ps=setdiff(ps,i);  % All the parents except i
if ~isempty(ps)
 A(j, ps) = 1;
end
for k=ps(:)'
 anck = find(A(k,:));
 if ~isempty(anck)
   A(j, anck) = 1;
 end
end

%%%%%%%%%

function A = old_update_row(A, j, dag)

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
