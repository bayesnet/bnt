function engine = jtree_sparse_inf_engine(bnet, varargin)
% JTREE_SPARSE_INF_ENGINE Junction tree inference engine when CPTs and Potentials are sparse
% engine = jtree_sparse_inf_engine(bnet, ...)
% It differs from jtree_inf_engine with all CPTs and potentials are 1D sparse arrays.
%
% The following optional arguments can be specified in the form of name/value pairs:
% [default value in brackets]
%
% clusters  - a cell array of sets of nodes we want to ensure are in the same clique (in addition to families) [ {} ]
% root      - the root of the junction tree will be a clique that contains this set of nodes [N]
% stages    - stages{t} is a set of nodes we want to eliminate before stages{t+1}, ... [ {1:N} ]
%
% e.g., engine = jtree_inf_engine(bnet, 'maximize', 1);
%
% For more details on the junction tree algorithm, see
% - "Probabilistic networks and expert systems", Cowell, Dawid, Lauritzen and Spiegelhalter, Springer, 1999
% - "Inference in Belief Networks: A procedural guide", C. Huang and A. Darwiche, 
%      Intl. J. Approximate Reasoning, 15(3):225-263, 1996.


% set default params
N = length(bnet.dag);
clusters = {};
root = N;
stages = { 1:N };

if nargin >= 2
  args = varargin;
  nargs = length(args);
  if ~isstr(args{1})
    error('the interface to jtree has changed; now, onodes is not allowed and all optional params must be passed by name')
  end
  for i=1:2:nargs
    switch args{i},
     case 'clusters', clusters = args{i+1}; 
     case 'root',     root = args{i+1}; 
     case 'stages',   stages = args{i+1}; 
     otherwise,  
      error(['invalid argument name ' args{i}]);       
    end
  end
end

engine = init_fields;
engine = class(engine, 'jtree_sparse_inf_engine', inf_engine(bnet));

onodes = bnet.observed;
%[engine.jtree, dummy, engine.cliques, B, w] = dag_to_jtree(bnet, onodes, stages, clusters);

porder = determine_elim_constraints(bnet, onodes);
strong = ~isempty(porder);
ns = bnet.node_sizes(:);
ns(onodes) = 1; % observed nodes have only 1 possible value
[engine.jtree, root2, engine.cliques, B, w] = ...
    graph_to_jtree(moralize(bnet.dag), ns, porder, stages, clusters);

engine.cliques_bitv = B;
engine.clique_weight = w;
C = length(engine.cliques);
engine.clpot = cell(1,C);

% Compute the separators between connected cliques.
[is,js] = find(engine.jtree > 0);
engine.separator = cell(C,C);
for k=1:length(is)
  i = is(k); j = js(k);
  engine.separator{i,j} = find(B(i,:) & B(j,:)); % intersect(cliques{i}, cliques{j});
end

% A node can be a member of many cliques, but is assigned to exactly one, to avoid
% double-counting its CPD. We assign node i to clique c if c is the "lightest" clique that
% contains i's family, so it can accomodate its CPD.

engine.clq_ass_to_node = zeros(1, N);
for i=1:N
  %c = clq_containing_nodes(engine, family(bnet.dag, i));
  clqs_containing_family = find(all(B(:,family(bnet.dag, i)), 2)); % all selected columns must be 1
  c = clqs_containing_family(argmin(w(clqs_containing_family)));  
  engine.clq_ass_to_node(i) = c; 
end

% Make the jtree rooted, so there is a fixed message passing order.
engine.root_clq = clq_containing_nodes(engine, root);
if engine.root_clq <= 0
  error(['no clique contains ' num2str(root)]);
end

[engine.jtree, engine.preorder, engine.postorder] = mk_rooted_tree(engine.jtree, engine.root_clq);

% collect 
engine.postorder_parents = cell(1,length(engine.postorder));
for n=engine.postorder(:)'
  engine.postorder_parents{n} = parents(engine.jtree, n);
end
% distribute
engine.preorder_children = cell(1,length(engine.preorder));
for n=engine.preorder(:)'
  engine.preorder_children{n} = children(engine.jtree, n);
end

ns = bnet.node_sizes;
engine.actual_node_sizes = ns;
 

%%%%%%%%

function engine = init_fields()

engine.jtree = [];
engine.cliques = [];
engine.separator = [];
engine.cliques_bitv = [];
engine.clique_weight = [];
engine.clpot = [];
engine.clq_ass_to_node = [];
engine.root_clq = [];
engine.preorder = [];
engine.postorder = [];
engine.preorder_children = [];
engine.postorder_parents = [];
engine.maximize = [];
engine.evidence = [];
engine.actual_node_sizes = [];



