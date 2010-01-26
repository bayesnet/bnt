function engine = pearl_inf_engine(bnet, varargin)
% PEARL_INF_ENGINE Pearl's algorithm (belief propagation)
% engine = pearl_inf_engine(bnet, ...)
%
% If the graph has no loops (undirected cycles), you should use the tree protocol,
% and the results will be exact.
% Otherwise, you should use the parallel protocol, and the results may be approximate.
%
% Optional arguments [default in brackets]
% 'protocol' - tree or parallel ['parallel']
%
% Optional arguments for the loopy case
% 'max_iter' - specifies the max num. iterations to perform [2*num nodes]
% 'tol' - convergence criterion on messages  [1e-3]
% 'momentum' - msg = (m*old + (1-m)*new). [m=0]
% 'filename' -  msgs will be printed to this file, so you can assess convergence while it runs [[]]
% 'storebel' - 1 means save engine.bel{n,t} for every iteration t and hidden node n [0]
%
% If there are discrete and cts nodes, we assume all the discretes are observed. In this
% case, you must use the parallel protocol, and the evidence pattern must be fixed.


N = length(bnet.dag);
protocol = [];
max_iter = 2*N;
% We use N+2 for the following reason:
% In N iterations, we get the exact answer for a tree.
% In the N+1st iteration, we notice that the results are the same as before, and terminate.
% In loopy_converged, we see that N+1 < max = N+2, and declare convergence.
tol = 1e-3;
momentum = 0;
filename = [];
storebel = 0;

args = varargin;
for i=1:2:length(args)
  switch args{i},
   case 'protocol', protocol = args{i+1};
   case 'max_iter', max_iter = args{i+1};
   case 'tol', tol = args{i+1};
   case 'momentum', momentum = args{i+1};
   case 'filename', filename = args{i+1};
   case 'storebel', storebel = args{i+1};
  end
end

engine.filename = filename;
engine.storebel = storebel;
engine.bel = [];

if strcmp(protocol, 'tree')
  % We first send messages up to the root (pivot node), and then back towards the leaves.
  % If the bnet is a singly connected graph (no loops), choosing a root induces a directed tree.
  % Peot and Shachter discuss ways to pick the root so as to minimize the work,
  % taking into account which nodes have changed.
  % For simplicity, we always pick the root to be the last node in the graph.
  % This means the first pass is equivalent to going forward in time in a DBN.

  engine.root = N;
  [engine.adj_mat, engine.preorder, engine.postorder, loopy] = ...
    mk_rooted_tree(bnet.dag, engine.root);
  % engine.adj_mat might have different edge orientations from bnet.dag
  if loopy
    error('can only apply tree protocol to loop-less graphs')
  end
else
  engine.root = [];
  engine.adj_mat = [];
  engine.preorder = [];
  engine.postorder = [];
end

engine.niter = [];
engine.protocol = protocol;
engine.max_iter = max_iter;
engine.tol = tol;
engine.momentum = momentum;
engine.maximize = [];

%onodes = find(~isemptycell(evidence));
onodes = bnet.observed;
engine.msg_type = determine_pot_type(bnet, onodes, 1:N); % needed also by marginal_nodes
if strcmp(engine.msg_type, 'cg')
  error('messages must be discrete or Gaussian')
end
[engine.msg_dag, disconnected_nodes] = mk_msg_dag(bnet, engine.msg_type, onodes);
engine.disconnected_nodes_bitv = zeros(1,N);
engine.disconnected_nodes_bitv(disconnected_nodes) = 1;


% this is where we store stuff between enter_evidence and marginal_nodes
engine.marginal = cell(1,N);
engine.evidence = []; 
engine.msg = [];

[engine.parent_index, engine.child_index] = mk_loopy_msg_indices(engine.msg_dag);

engine = class(engine, 'pearl_inf_engine', inf_engine(bnet));
 

%%%%%%%%%

function [dag, disconnected_nodes] = mk_msg_dag(bnet, msg_type, onodes)

% If we are using Gaussian msgs, all discrete nodes must be observed;
% they are then disconnected from the graph, so we don't try to send
% msgs to/from them: their observed value simply serves to index into
% the right set of parameters for the Gaussian nodes (which use CPD.ps
% instead of parents(dag), and hence are unaffected by this "surgery").

disconnected_nodes = [];
switch msg_type
 case 'd', dag = bnet.dag;
 case 'g',
  disconnected_nodes = bnet.dnodes;
  dag = bnet.dag;
  for i=disconnected_nodes(:)'
    ps = parents(bnet.dag, i);
    cs = children(bnet.dag, i);
    if ~isempty(ps), dag(ps, i) = 0; end
    if ~isempty(cs), dag(i, cs) = 0; end
  end
end


%%%%%%%%%%
function [parent_index, child_index] = mk_loopy_msg_indices(dag)
% MK_LOOPY_MSG_INDICES Compute "port numbers" for message passing
% [parent_index, child_index] = mk_loopy_msg_indices(bnet)
%
% child_index{n}(c) = i means c is n's i'th child, i.e., i = find_equiv_posns(c, children(n))
% child_index{n}(c) = 0 means c is not a child of n.
% parent_index{n}{p} is defined similarly.
% We need to use these indices since the pi_from_parent/ lambda_from_child cell arrays
% cannot be sparse, and hence cannot be indexed by the actual number of the node.
% Instead, we use the number of the "port" on which the message arrived.

N = length(dag);
child_index = cell(1,N);
parent_index = cell(1,N);
for n=1:N
  cs = children(dag, n);
  child_index{n} = sparse(1,N);
  for i=1:length(cs)
    c = cs(i);
    child_index{n}(c) = i;
  end
  ps = parents(dag, n);
  parent_index{n} = sparse(1,N);
  for i=1:length(ps)
    p = ps(i);
    parent_index{n}(p) = i;
  end
end




